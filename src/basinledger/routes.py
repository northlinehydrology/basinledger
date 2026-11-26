"""Routes: where a day's water goes in the models, and what each model actually took or gave.

A route is a string naming a source or sink:

    well:<name>                   a MODFLOW 6 well (withdrawal)
    swat:place_of_use:<name>      a SWAT+ water allocation place of use, taking from its channel (withdrawal)
    swat:point_source:<name>      a SWAT+ daily point source into its channel (discharge)
    ground:<name>                 water into an aquifer cell, through a WEL entry injecting it (discharge): a septic
                                  drainfield, a system's leaks ([[ground]])
    land:<name>                   water put on SWAT+'s fields, HRUs by weight, some of it running off ([[land]];
                                  swatplus-bmi 0.15.0's irrigation applied from outside)
    outfall:<name>                recorded in the ledger only (discharge)
    offsite:<name>                water from or to outside the basin, recorded in the ledger only (either way)

(``utility:<name>`` routes go to the utility chain, which then requests these.)

Each day the driver adds requests (m3), applies them to the models (SWAT+ before its day, MODFLOW 6 after preparing
its day), and reads back what each model did. For every route the ledger gets ``route:<r>:requested`` and
``route:<r>:delivered``; withdrawals from SWAT+ also get ``route:<r>:shortfall`` (what the channel couldn't supply).
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from .coupling import swatplus as S
from .ledger import FLOAT32
from .units import DAY

KINDS = ("well", "swat:place_of_use", "swat:point_source", "ground", "land", "outfall", "offsite")
LAND_SOIL = "land_irrigation_water__infiltration_volume_flow_rate"  # m3 s-1 per HRU, set for the next day, read back
LAND_RUNOFF = "land_irrigation_water__runoff_volume_flow_rate"     # m3 s-1 per HRU, the same
LAND_MASS = "land_irrigation_water_{x}__mass_flow_rate"            # kg s-1 per HRU into the top layer (no3, solp)


def kind(route: str) -> tuple[str, str]:
    for k in KINDS:
        if route.startswith(k + ":"):
            return k, route[len(k) + 1:]
    if route.startswith("utility:"):
        return "utility", route[len("utility:"):]
    raise ValueError(f"unknown route {route!r} (known: {', '.join(KINDS)}, utility)")


class Router:
    def __init__(self, swat: S.Swat, mf6, wells: dict[str, tuple[str, int]], ledger,
                 grounds: dict[str, tuple[str, int]] | None = None, lands: dict | None = None):
        self.swat, self.mf6, self.wells, self.ledger = swat, mf6, wells, ledger
        self.grounds = grounds or {}  # returns to the ground: name -> (WEL package, node)
        # water onto the fields: name -> (LandConfig, [(HRU index in SWAT+'s hru grid, weight)])
        self.lands = lands or {}
        self.hrus = swat.names(LAND_SOIL) if self.lands else []
        # a well or ground return over several cells: its route -> [(each cell's own route, weight)]; each cell is
        # an entry of its own above (<name>#1, #2, ...), and the route's request is split among them by weight
        self.groups: dict[str, list[tuple[str, float]]] = {}
        self.places = swat.names(S.DEMAND) if swat.has(S.DEMAND) else []
        self.points = swat.names(S.POINT) if swat.has(S.POINT) else []
        self.req: dict[str, float] = defaultdict(float)
        self.used: set[str] = set()

    def check(self, route: str):
        if route in self.groups:  # booked as one route: requested, and delivered over its cells
            if route not in self.used:
                self.used.add(route)
                self.ledger.pair(f"route:{route}:requested", f"route:{route}:delivered")
            return
        k, name = kind(route)
        known = {"well": self.wells, "swat:place_of_use": self.places, "swat:point_source": self.points,
                 "ground": self.grounds, "land": self.lands}.get(k)
        if known is not None and name not in known:
            raise KeyError(f"route {route!r}: no {k} named {name!r}")
        if k in KINDS and route not in self.used:
            self.used.add(route)
            r = f"route:{route}:"
            if k == "swat:place_of_use":
                self.ledger.pair(r + "requested", [r + "delivered", r + "shortfall"], rtol=FLOAT32)
                self.ledger.zero(r + "overdraw")
            elif k == "swat:point_source":
                self.ledger.pair(r + "requested", r + "delivered", rtol=FLOAT32)
            elif k == "land":  # what SWAT+ reports went in: into the soil, and running off
                self.ledger.pair(r + "requested", r + "delivered", rtol=FLOAT32)
                self.ledger.pair(r + "delivered", [r + "into_soil", r + "runoff"])
            else:
                self.ledger.pair(r + "requested", r + "delivered")

    def request(self, route: str, m3: float):
        self.check(route)
        self.req[route] += m3
        for sub, w in self.groups.get(route, ()):
            self.req[sub] += m3 * w

    def land_split(self) -> dict[str, np.ndarray]:
        """Each land route's water today on each of SWAT+'s HRUs (m3), by its weights."""
        out = {}
        for name, (_, rows) in self.lands.items():
            v = np.zeros(len(self.hrus))
            for i, w in rows:
                v[i] += self.req.get(f"land:{name}", 0.0) * w
            out[name] = v
        return out

    def apply_swat(self):
        if self.lands:  # the water onto the fields: into the soil, and running off by each route's share
            soil, off = np.zeros(len(self.hrus)), np.zeros(len(self.hrus))
            for name, v in self.land_split().items():
                share = self.lands[name][0].runoff
                soil += (1.0 - share) * v
                off += share * v
            self.swat.set(LAND_SOIL, soil / DAY)
            self.swat.set(LAND_RUNOFF, off / DAY)
        if self.places:
            d = np.zeros(len(self.places))
            for i, n in enumerate(self.places):
                d[i] = self.req.get(f"swat:place_of_use:{n}", 0.0) / DAY
            self.swat.set(S.DEMAND, d)
        if self.points:
            p = np.zeros(len(self.points))
            for i, n in enumerate(self.points):
                p[i] = self.req.get(f"swat:point_source:{n}", 0.0) / DAY
            self.swat.set(S.POINT, p)

    def read_swat(self):
        L = self.ledger
        if self.lands:  # what went in, each HRU's shared among the routes on it by what each asked there
            soil, off = self.swat.get(LAND_SOIL) * DAY, self.swat.get(LAND_RUNOFF) * DAY
            split = self.land_split()
            asked_soil = sum((1.0 - self.lands[n][0].runoff) * v for n, v in split.items())
            asked_off = sum(self.lands[n][0].runoff * v for n, v in split.items())
            f_soil = np.divide(soil, asked_soil, out=np.zeros_like(soil), where=asked_soil > 0)
            f_off = np.divide(off, asked_off, out=np.zeros_like(off), where=asked_off > 0)
            for name, v in split.items():
                r = f"land:{name}"
                if r in self.used:
                    share = self.lands[name][0].runoff
                    into, run = float(((1.0 - share) * v * f_soil).sum()), float((share * v * f_off).sum())
                    L.enter(f"route:{r}:requested", self.req.get(r, 0.0))
                    L.enter(f"route:{r}:delivered", into + run)
                    L.enter(f"route:{r}:into_soil", into)
                    L.enter(f"route:{r}:runoff", run)
        if self.places:
            got = self.swat.get(S.WITHDRAWN) * DAY
            for i, n in enumerate(self.places):
                r = f"swat:place_of_use:{n}"
                if r in self.used:
                    want = self.req.get(r, 0.0)
                    L.enter(f"route:{r}:requested", want)
                    L.enter(f"route:{r}:delivered", got[i])
                    L.enter(f"route:{r}:shortfall", max(want - got[i], 0.0))
                    L.enter(f"route:{r}:overdraw", max(got[i] - want * (1 + FLOAT32) - 1e-6, 0.0))
        if self.points:
            got = self.swat.get(S.POINT) * DAY
            for i, n in enumerate(self.points):
                r = f"swat:point_source:{n}"
                if r in self.used:
                    L.enter(f"route:{r}:requested", self.req.get(r, 0.0))
                    L.enter(f"route:{r}:delivered", got[i])

    def _maw(self, node):
        """A multi-aquifer well's (model view, 0-based number), or None for a WEL well's node."""
        if isinstance(node, tuple) and node[0] == "maw":
            view = self.mf6.view_of(node[2]) if hasattr(self.mf6, "view_of") else self.mf6
            return view, node[1]
        return None

    def apply_wells(self):
        if self.mf6 is None:
            self._entries, self._ground_entries = {}, {}
            return
        self._entries = {name: (pkg, node if self._maw(node) else self.mf6.entry(pkg, node))
                         for name, (pkg, node) in self.wells.items()}
        self._ground_entries = {name: (pkg, self.mf6.entry(pkg, node)) for name, (pkg, node) in self.grounds.items()}
        for name, (pkg, i) in self._ground_entries.items():  # injected: positive into the aquifer
            self.mf6.set_well(pkg, i, self.req.get(f"ground:{name}", 0.0))
        for name, (pkg, i) in self._entries.items():
            q = -self.req.get(f"well:{name}", 0.0)
            maw = self._maw(i)
            if maw:
                maw[0].set_maw_rate(pkg, maw[1], q)
            else:
                self.mf6.set_well(pkg, i, q)

    def read_wells(self):
        got = {}
        for name, (pkg, i) in self._ground_entries.items():
            got[f"ground:{name}"] = self.mf6.simvals(pkg)[i]
        for name, (pkg, i) in self._entries.items():
            if not self._maw(i):
                got[f"well:{name}"] = -self.mf6.simvals(pkg)[i]
        for route, subs in self.groups.items():
            if route in self.used:
                self.ledger.enter(f"route:{route}:requested", self.req.get(route, 0.0))
                self.ledger.enter(f"route:{route}:delivered", float(sum(got[s] for s, _ in subs)))
        for name, (pkg, i) in self._ground_entries.items():
            r = f"ground:{name}"
            if r in self.used:
                self.ledger.enter(f"route:{r}:requested", self.req.get(r, 0.0))
                self.ledger.enter(f"route:{r}:delivered", self.mf6.simvals(pkg)[i])
        for name, (pkg, i) in self._entries.items():
            r = f"well:{name}"
            if r in self.used:
                self.ledger.enter(f"route:{r}:requested", self.req.get(r, 0.0))
                maw = self._maw(i)
                if maw:  # what it pumped, and where that came from: the aquifer by its connections, and its wellbore
                    v, k = maw
                    self.ledger.enter(f"route:{r}:delivered", -v.maw_pumped(pkg, k))
                    self.ledger.enter(f"route:{r}:from_aquifer", -float(v.simvals(pkg)[v.maw_connections(pkg, k)].sum()))
                    self.ledger.enter(f"route:{r}:from_wellbore", v.maw_wellbore(pkg, k))
                else:
                    self.ledger.enter(f"route:{r}:delivered", -self.mf6.simvals(pkg)[i])

    def close_day(self):
        for r in self.used:
            if r.startswith(("outfall:", "offsite:")):
                self.ledger.enter(f"route:{r}:requested", self.req.get(r, 0.0))
                self.ledger.enter(f"route:{r}:delivered", self.req.get(r, 0.0))
        self.req = defaultdict(float)
