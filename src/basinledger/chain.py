"""Water systems: sources -> distribution -> customers -> sewer, septic and consumption -> plant -> outfall.

A system is accounting, not a model (docs/design.md, "Water systems with customers"). Each day it takes its
customers' demand and its users' withdrawals, asks its sources for the production that needs, and places the water in
the hydrologic models through routes:

    demand      = the customer classes' demand + the users' withdrawals from the tap
    asked       = demand / (1 - losses)               from the sources, by share or in priority up to capacity
    produced    = what the sources delivered          known after the models' day
    share       = produced / asked                    every customer, user, loss and return scales by it
    delivered   = share x demand                      the rest unmet
    losses      = produced - delivered                leaks (to the ground) + unaccounted (consumed)
    each class  = sewered + septic + consumed          the class's delivered water, by its shares
    sewer       = the classes' sewered + the users' returns (+ inflow and infiltration, if given)

Returns leave the next day: the sewer becomes tomorrow's plant outflow, and the septic water and leaks go to the
ground tomorrow (the delivery is only known after the day, and a return enters the models before it). The plant's
outflow goes to its outfall, less what other systems take from it as a source (``effluent:<system>``, reuse).

Config (one ``[[utilities]]`` table per system):

    name = "town"
    sources = [{route = "swat:place_of_use:intake", share = 0.7}, {route = "well:pw2", share = 0.3}]
                               # share: a constant, or 12 monthly values; or priority = 1, 2, ... and capacity
                               # (m3 s-1, the last may leave it out); a route, or "effluent:<system>"
    # shares = "town_shares.csv"   # or dated shares: date, then a column per source's route (its sources then
                               # {route} alone), each row summing to 1 and held until the next date
    losses = 0.10              # share of production
    leaks = {share = 0.5, to = "ground:town_leaks"}   # of the losses, to the ground; the rest unaccounted
    customers = [{name = "homes", demand = "homes.csv", sewered = 0.8, consumed = 0.2},
                 {name = "edge", demand = 0.001, septic = 0.9, septic_to = "ground:edge", consumed = 0.1}]
                               # demand: m3 s-1, 12 monthly values (m3 s-1), or a CSV of date and m3 by day
    outfall = "swat:point_source:wwtp"

The first form still reads: ``demand`` (m3 s-1) is a class of existing customers, ``sewered`` its sewered share
(the rest consumed), ``infiltration`` (m3 s-1) water from outside into the sewer.

Ledger accounts, all m3 a day, under ``utility:<name>:``; the pairs are in ``declare``.
"""

from __future__ import annotations

import bisect
import csv
import datetime as dt
from dataclasses import dataclass, field
from pathlib import Path

from .units import DAY

EFFLUENT = "effluent:"


def _monthly(where: str, v) -> list[float]:
    """A constant or 12 monthly values, as 12 values."""
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return [float(v)] * 12
    if isinstance(v, (list, tuple)) and len(v) == 12:
        return [float(x) for x in v]
    raise ValueError(f"{where}: a constant or 12 monthly values; got {v!r}")


@dataclass
class Source:
    route: str
    share: float | list[float] | None = None   # a constant or 12 monthly values
    priority: int | None = None
    capacity: float | None = None              # m3 s-1, with priority

    @property
    def effluent(self) -> str | None:
        return self.route[len(EFFLUENT):] if self.route.startswith(EFFLUENT) else None


@dataclass
class Customers:
    """A class of customers: its demand and where its water goes."""
    name: str
    demand: list[float] | dict[dt.date, float]  # 12 monthly m3 s-1, or m3 by date
    sewered: float = 1.0
    septic: float = 0.0
    consumed: float = 0.0
    septic_to: str | None = None

    def on(self, day: dt.date) -> float:
        """The class's demand on the day, m3."""
        if isinstance(self.demand, dict):
            return self.demand.get(day, 0.0)
        return self.demand[day.month - 1] * DAY


def _customers(system: str, c: dict, base: Path | None) -> Customers:
    where = f"utility {system}: customers {c.get('name')!r}"
    if set(c) - {"name", "demand", "sewered", "septic", "septic_to", "consumed"} or "name" not in c or "demand" not in c:
        raise ValueError(f"{where}: takes name, demand, sewered, septic, septic_to and consumed; got {sorted(c)}")
    d = c["demand"]
    if isinstance(d, str):
        path = Path(d) if Path(d).is_absolute() or base is None else Path(base) / d
        with open(path, newline="") as f:
            demand = {dt.date.fromisoformat(r["date"]): float(r["m3"]) for r in csv.DictReader(f)}
    else:
        demand = _monthly(where + " demand", d)
    sewered, septic = float(c.get("sewered", 0.0)), float(c.get("septic", 0.0))
    consumed = float(c.get("consumed", 1.0 - sewered - septic))
    if min(sewered, septic, consumed) < 0 or abs(sewered + septic + consumed - 1.0) > 1e-9:
        raise ValueError(f"{where}: sewered, septic and consumed must be shares that add up to 1")
    if septic > 0 and not str(c.get("septic_to", "")).startswith("ground:"):
        raise ValueError(f"{where}: septic water needs septic_to, a ground: route")
    return Customers(c["name"], demand, sewered, septic, consumed, c.get("septic_to"))


@dataclass
class UtilityChain:
    name: str
    sources: list[Source]
    outfall: str
    losses: float = 0.0
    demand: float = 0.0          # the first form: m3 s-1 of existing customers
    sewered: float = 1.0
    infiltration: float = 0.0    # m3 s-1 into the sewer, water from outside
    customers: list[Customers] = field(default_factory=list)
    leaks: float = 0.0           # share of the losses to the ground
    leaks_to: str | None = None
    shares: list = field(default_factory=list)  # [(date, {route: share}), ...] sorted, each held to the next date
    # the day, between plan and settle
    _withdrawals: float = field(default=0.0, repr=False)
    _returns: float = field(default=0.0, repr=False)          # users' returns to this sewer, scaled (settle)
    _day: dict = field(default_factory=dict, repr=False)
    # held overnight: tomorrow's plant outflow, septic water and leaks; reuse asked of this plant today
    held: dict = field(default_factory=lambda: {"sewer": 0.0, "ground": {}}, repr=False)
    _reuse_asked: dict = field(default_factory=dict, repr=False)
    _reuse_given: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        if self.demand and not any(c.name == "existing" for c in self.customers):
            self.customers.insert(0, Customers("existing", [self.demand] * 12, self.sewered, 0.0,
                                               1.0 - self.sewered))
        if not self.sources:
            raise ValueError(f"utility {self.name}: no sources")
        prio = [s.priority is not None for s in self.sources]
        if any(prio) and not all(prio):
            raise ValueError(f"utility {self.name}: sources take shares or priorities, not both")
        if self.shares:
            if any(prio) or any(s.share is not None for s in self.sources):
                raise ValueError(f"utility {self.name}: with shares from a table, sources take neither a share nor "
                                 "a priority")
            routes = {s.route for s in self.sources}
            for date, row in self.shares:
                if set(row) != routes:
                    raise ValueError(f"utility {self.name}: the shares' columns {sorted(row)} aren't its sources "
                                     f"{sorted(routes)}")
                if abs(sum(row.values()) - 1.0) > 1e-9 or min(row.values()) < 0:
                    raise ValueError(f"utility {self.name}: shares on {date} sum to {sum(row.values())}, not 1")
        elif not any(prio):
            for m in range(12):
                total = sum(_monthly(f"utility {self.name}: {s.route} share", s.share)[m] for s in self.sources)
                if abs(total - 1.0) > 1e-9:
                    raise ValueError(f"utility {self.name}: source shares sum to {total}, not 1"
                                     + (f" (month {m + 1})" if m else ""))
        if not 0 <= self.losses < 1 or not 0 <= self.sewered <= 1 or not 0 <= self.leaks <= 1:
            raise ValueError(f"utility {self.name}: losses must be in [0, 1), sewered and leaks in [0, 1]")
        if self.leaks > 0 and not str(self.leaks_to or "").startswith("ground:"):
            raise ValueError(f"utility {self.name}: leaks need to, a ground: route")
        for s in self.sources:
            if s.route.startswith("utility:"):
                raise ValueError(f"utility {self.name}: a source can't be another utility; its plant's effluent "
                                 f"can (effluent:<system>)")
            if s.effluent == self.name:
                raise ValueError(f"utility {self.name}: a system can't reuse its own effluent")

    @classmethod
    def from_config(cls, c: dict, base: Path | None = None) -> "UtilityChain":
        name = c["name"]
        sources = []
        for s in c["sources"]:
            if "priority" in s:
                sources.append(Source(s["route"], None, int(s["priority"]),
                                      float(s["capacity"]) if "capacity" in s else None))
            elif "share" in s:
                sources.append(Source(s["route"], s["share"] if isinstance(s["share"], list) else float(s["share"])))
            else:
                sources.append(Source(s["route"]))
        shares = []
        if c.get("shares"):  # a table of dated shares: date, then a column per source's route
            path = Path(c["shares"])
            path = path if path.is_absolute() or base is None else Path(base) / path
            with open(path, newline="") as f:
                rows = list(csv.DictReader(f))
            shares = sorted((dt.date.fromisoformat(r["date"]), {k: float(v) for k, v in r.items() if k != "date"})
                            for r in rows)
            if not shares:
                raise ValueError(f"utility {name}: {path} has no shares")
        leaks = c.get("leaks") or {}
        return cls(name=name, sources=sources, outfall=c["outfall"], losses=float(c.get("losses", 0.0)),
                   demand=float(c.get("demand", 0.0)), sewered=float(c.get("sewered", 1.0)),
                   infiltration=float(c.get("infiltration", 0.0)),
                   customers=[_customers(name, x, base) for x in c.get("customers", [])],
                   leaks=float(leaks.get("share", 0.0)), leaks_to=leaks.get("to"), shares=shares)

    # ---------- routes ----------
    def routes(self) -> list[str]:
        """The routes this system places water on (not another system's effluent)."""
        out = [s.route for s in self.sources if s.effluent is None] + [self.outfall]
        out += [c.septic_to for c in self.customers if c.septic > 0]
        if self.leaks > 0:
            out.append(self.leaks_to)
        return list(dict.fromkeys(out))

    # ---------- the day ----------
    def add_user(self, withdrawal: float, returned: float):
        """A user's day on this system: m3 withdrawn from the tap, m3 returned to the sewer (its share delivered,
        ``add_return``, after the day; this argument is kept for the first form's callers and adds to it)."""
        self._withdrawals += withdrawal
        self._returns += returned

    def add_return(self, returned: float):
        """A user's return to this sewer after the day: its withdrawal's delivered share of what it returned."""
        self._returns += returned

    def plan(self, ledger, day: dt.date, systems: dict | None = None) -> dict[str, float]:
        """Before the models' day: the demand, the production asked, and each source's request (m3). Another system's
        effluent is asked of it here (its ``release`` gives it)."""
        k = f"utility:{self.name}:"
        classes = {c.name: c.on(day) for c in self.customers}
        demand = sum(classes.values()) + self._withdrawals
        asked = demand / (1.0 - self.losses)
        for n, v in classes.items():
            ledger.enter(k + f"demand:{n}", v)
        ledger.enter(k + "users_withdrawal", self._withdrawals)
        ledger.enter(k + "demand", demand)
        ledger.enter(k + "asked", asked)
        requests: dict[str, float] = {}
        if self.shares:  # the last row dated on or before the day
            i = bisect.bisect_right([d for d, _ in self.shares], day) - 1
            if i < 0:
                raise ValueError(f"utility {self.name}: no source shares on or before {day} (the table starts "
                                 f"{self.shares[0][0]})")
            for route, share in self.shares[i][1].items():
                requests[route] = requests.get(route, 0.0) + asked * share
        elif self.sources[0].priority is None:
            for s in self.sources:
                requests[s.route] = requests.get(s.route, 0.0) + asked * _monthly(s.route, s.share)[day.month - 1]
        else:
            left = asked
            for s in sorted(self.sources, key=lambda s: s.priority):
                r = left if s.capacity is None else min(left, s.capacity * DAY)
                requests[s.route] = requests.get(s.route, 0.0) + r
                left -= r
        for route, v in requests.items():
            ledger.enter(k + f"source:{route}:requested", v)
            src = route[len(EFFLUENT):] if route.startswith(EFFLUENT) else None
            if src is not None:
                systems[src]._reuse_asked[self.name] = v
        self._day = {"classes": classes, "demand": demand, "asked": asked, "requests": requests}
        return {r: v for r, v in requests.items() if not r.startswith(EFFLUENT)}

    def release(self, ledger) -> dict[str, float]:
        """At the day's start, after every system's plan: yesterday's sewer is today's plant outflow, less what other
        systems take from it (in proportion to what they asked, up to all of it), to the outfall; yesterday's septic
        water and leaks to the ground. Returns the requests (m3)."""
        k = f"utility:{self.name}:"
        out = self.held["sewer"]
        asked = sum(self._reuse_asked.values())
        given = min(asked, out)
        self._reuse_given = {n: (v / asked * given if asked > 0 else 0.0) for n, v in self._reuse_asked.items()}
        for n, v in self._reuse_given.items():
            ledger.enter(k + f"reuse:{n}", v)
        ledger.enter(k + "plant_outflow", out)
        ledger.enter(k + "reuse", given)
        ledger.enter(k + "outfall", out - given)
        requests = {self.outfall: out - given}
        for route, v in self.held["ground"].items():
            requests[route] = requests.get(route, 0.0) + v
            ledger.enter(k + f"to_ground:{route}", v)
        self._reuse_asked = {}
        return requests

    def reuse_given(self, to: str) -> float:
        return self._reuse_given.get(to, 0.0)

    def share(self, ledger, systems: dict) -> float:
        """After the models' day: what each source delivered, and the share of production asked that came."""
        k = f"utility:{self.name}:"
        produced = 0.0
        for route, req in self._day["requests"].items():
            src = route[len(EFFLUENT):] if route.startswith(EFFLUENT) else None
            if src is not None:
                got = systems[src].reuse_given(self.name)
            else:  # a route shared with others: this system's part of what it delivered
                asked_r = ledger.today(f"route:{route}:requested")
                got = req * min(ledger.today(f"route:{route}:delivered") / asked_r, 1.0) if asked_r > 0 else 0.0
            ledger.enter(k + f"source:{route}:delivered", got)
            produced += got
        asked = self._day["asked"]
        share = min(produced / asked, 1.0) if asked > 0 else 1.0
        ledger.enter(k + "produced", produced)
        ledger.enter(k + "shortfall", asked - produced)
        ledger.enter(k + "share", share)
        self._day.update(produced=produced, share=share)
        return share

    def settle(self, ledger):
        """After ``share`` and the users' returns (``add_return``): the day's water by where it went, and tomorrow's
        held returns."""
        k = f"utility:{self.name}:"
        d, s = self._day, self._day["share"]
        delivered = s * d["demand"]
        losses = d["produced"] - delivered
        ledger.enter(k + "delivered", delivered)
        ledger.enter(k + "unmet", d["demand"] - delivered)
        ledger.enter(k + "users_delivered", s * self._withdrawals)
        ledger.enter(k + "losses", losses)
        leaks = self.leaks * losses
        ledger.enter(k + "leaks", leaks)
        ledger.enter(k + "unaccounted", losses - leaks)
        sewer, ground = 0.0, {}
        consumed = 0.0
        for c in self.customers:
            got = s * d["classes"][c.name]
            ledger.enter(k + f"delivered:{c.name}", got)
            ledger.enter(k + f"sewered:{c.name}", c.sewered * got)
            ledger.enter(k + f"septic:{c.name}", c.septic * got)
            ledger.enter(k + f"consumed:{c.name}", c.consumed * got)
            sewer += c.sewered * got
            consumed += c.consumed * got
            if c.septic > 0:
                ground[c.septic_to] = ground.get(c.septic_to, 0.0) + c.septic * got
        if leaks > 0:
            ground[self.leaks_to] = ground.get(self.leaks_to, 0.0) + leaks
        ledger.enter(k + "consumed", consumed)
        ledger.enter(k + "users_return", self._returns)
        ledger.enter(k + "infiltration", self.infiltration * DAY)
        sewer += self._returns + self.infiltration * DAY
        ledger.enter(k + "sewer", sewer)
        ledger.enter(k + "held_to_ground", sum(ground.values()))
        self.held = {"sewer": sewer, "ground": ground}
        self._withdrawals = self._returns = 0.0

    def close_day(self, ledger, day: dt.date | None = None) -> dict[str, float]:
        """One system alone, its sources delivering all they're asked (the first form's day, for tests): plan,
        release, share, settle. Returns the day's requests."""
        day = day or dt.date(2000, 1, 1)
        req = self.plan(ledger, day, {self.name: self})
        for r, v in self.release(ledger).items():
            req[r] = req.get(r, 0.0) + v
        for r, v in self._day["requests"].items():
            ledger.enter(f"route:{r}:requested", v)
            ledger.enter(f"route:{r}:delivered", v)
        self.share(ledger, {self.name: self})
        self.settle(ledger)
        return req

    def declare(self, ledger, rtol: float | None = None):
        k = f"utility:{self.name}:"
        ledger.pair(k + "asked", [k + "produced", k + "shortfall"], rtol=rtol)
        ledger.pair(k + "produced", [k + f"source:{r}:delivered" for r in dict.fromkeys(s.route for s in self.sources)],
                    rtol=rtol)
        ledger.pair(k + "produced", [k + "delivered", k + "losses"], rtol=rtol)
        ledger.pair(k + "demand", [k + f"demand:{c.name}" for c in self.customers] + [k + "users_withdrawal"])
        ledger.pair(k + "delivered", [k + f"delivered:{c.name}" for c in self.customers] + [k + "users_delivered"],
                    rtol=rtol)
        for c in self.customers:
            ledger.pair(k + f"delivered:{c.name}", [k + f"sewered:{c.name}", k + f"septic:{c.name}",
                                                    k + f"consumed:{c.name}"])
        ledger.pair(k + "losses", [k + "leaks", k + "unaccounted"])
        ledger.pair(k + "plant_outflow", [k + "outfall", k + "reuse"])
