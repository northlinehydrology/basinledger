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
    # with transport (docs/design.md, Mass through the human-system routes), kg m-3 by species:
    effluent = {nitrate = "wwtp_dmr.csv", chloride = "sewer"}  # the plant's outflow: a constant, 12 monthly values,
                               # a CSV (date, then a column per species, each row held to the next), or "sewer"
                               # (the sewer's mass passes through)
    septic_effluent = {nitrate = 0.04, chloride = 0.08}          # what septic water carries to the water table
                               # (the same forms, or "tap", the mixed water's); a class's own wins
    # a source's {route, ..., concentration = {chloride = 0.02}} states its water's [A]

The first form still reads: ``demand`` (m3 s-1) is a class of existing customers, ``sewered`` its sewered share
(the rest consumed), ``infiltration`` (m3 s-1) water from outside into the sewer.

Ledger accounts, all m3 a day, under ``utility:<name>:``; the pairs are in ``declare``. With transport, each
species' mass, kg, under ``utility:<name>:<species>:`` (``mass_release``, ``mass_settle``; pairs in
``declare_mass``): the sources' concentrations mixed by what each delivered, carried by everything the system
delivers, loses and leaks; each class's septic water at its septic effluent concentration; the sewer's mass the tap
water's in the sewered water and the users' returns; the plant's outflow at its effluent concentration the day it
leaves; outdoor use's and unaccounted water's mass dropped, booked as leaving.
"""

from __future__ import annotations

import bisect
import csv
import datetime as dt
from dataclasses import dataclass, field
from pathlib import Path

from .units import DAY

EFFLUENT = "effluent:"
TAP, SEWER = "tap", "sewer"  # a return at the system's mixed water's concentration; a plant passing its sewer through


class Quality:
    """A stated concentration through time, kg m-3: a constant, 12 monthly values, or dated values from a CSV
    (``date``, then a column per species, each row held until the next date; a day before the first is refused).
    ``keyword``: ``"tap"`` or ``"sewer"``, a concentration the system works out (``None`` for a stated one)."""

    _tables: dict = {}  # a CSV read once: path -> [(date, {column: value})]

    def __init__(self, where: str, species: str, spec, base: Path | None, keywords: tuple[str, ...]):
        self.where, self.keyword, self.monthly, self.dated = where, None, None, None
        if isinstance(spec, str) and spec in keywords:
            self.keyword = spec
        elif isinstance(spec, str):
            path = Path(spec) if Path(spec).is_absolute() or base is None else Path(base) / spec
            if path not in Quality._tables:
                with open(path, newline="") as f:
                    Quality._tables[path] = sorted((dt.date.fromisoformat(r["date"]), r) for r in csv.DictReader(f))
            rows = Quality._tables[path]
            if not rows or species not in rows[0][1]:
                raise ValueError(f"{where}: {path} has no column {species!r}")
            self.dated = ([d for d, _ in rows], [float(r[species]) for _, r in rows])
            if min(self.dated[1]) < 0:
                raise ValueError(f"{where}: {path} has a concentration below zero")
        else:
            self.monthly = _monthly(where, spec)
            if min(self.monthly) < 0:
                raise ValueError(f"{where}: a concentration below zero")

    def on(self, day: dt.date) -> float:
        if self.monthly is not None:
            return self.monthly[day.month - 1]
        dates, values = self.dated
        i = bisect.bisect_right(dates, day) - 1
        if i < 0:
            raise ValueError(f"{self.where}: no concentration on or before {day} (the table starts {dates[0]})")
        return values[i]


def _qualities(where: str, spec, base: Path | None, keywords: tuple[str, ...]) -> dict[str, Quality]:
    if not isinstance(spec, dict):
        raise ValueError(f"{where} is {{<species> = kg m-3, 12 monthly values, a CSV or {' or '.join(keywords)}}}; "
                         f"got {spec!r}")
    return {sp: Quality(f"{where} {sp}", sp, v, base, keywords) for sp, v in spec.items()}


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
    concentration: dict = field(default_factory=dict)  # species -> kg m-3 of its water, stated [A]

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
    septic_effluent: dict = field(default_factory=dict)  # species -> Quality, over the system's

    def on(self, day: dt.date) -> float:
        """The class's demand on the day, m3."""
        if isinstance(self.demand, dict):
            return self.demand.get(day, 0.0)
        return self.demand[day.month - 1] * DAY


def _customers(system: str, c: dict, base: Path | None) -> Customers:
    where = f"utility {system}: customers {c.get('name')!r}"
    if set(c) - {"name", "demand", "sewered", "septic", "septic_to", "consumed", "septic_effluent"} \
            or "name" not in c or "demand" not in c:
        raise ValueError(f"{where}: takes name, demand, sewered, septic, septic_to, consumed and septic_effluent; got "
                         f"{sorted(c)}")
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
    return Customers(c["name"], demand, sewered, septic, consumed, c.get("septic_to"),
                     _qualities(f"{where} septic_effluent", c.get("septic_effluent", {}), base, (TAP,)))


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
    effluent: dict = field(default_factory=dict)         # species -> Quality: the plant's outflow
    septic_effluent: dict = field(default_factory=dict)  # species -> Quality: septic water at the water table
    # the day, between plan and settle
    _withdrawals: float = field(default=0.0, repr=False)
    _returns: float = field(default=0.0, repr=False)          # users' returns to this sewer, scaled (settle)
    _day: dict = field(default_factory=dict, repr=False)
    # held overnight: tomorrow's plant outflow, septic water and leaks; reuse asked of this plant today
    held: dict = field(default_factory=lambda: {"sewer": 0.0, "ground": {}, "kg": {}}, repr=False)
    # each species' mixed concentration on the day last settled (its users' intake), and its plant outflow's today
    concentration: dict = field(default_factory=dict, repr=False)
    _effluent_c: dict = field(default_factory=dict, repr=False)
    _return_parts: list = field(default_factory=list, repr=False)  # users' returns: (m3, {species: kg reported})
    _mass_day: dict = field(default_factory=dict, repr=False)
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
            conc = s.get("concentration", {})
            if not isinstance(conc, dict) or any(isinstance(v, bool) or not isinstance(v, (int, float)) or v < 0
                                                 for v in conc.values()):
                raise ValueError(f"utility {name}: source {s['route']} concentration is {{<species> = kg m-3}}, none "
                                 f"below zero; got {conc!r}")
            conc = {k: float(v) for k, v in conc.items()}
            if "priority" in s:
                sources.append(Source(s["route"], None, int(s["priority"]),
                                      float(s["capacity"]) if "capacity" in s else None, conc))
            elif "share" in s:
                sources.append(Source(s["route"], s["share"] if isinstance(s["share"], list) else float(s["share"]),
                                      concentration=conc))
            else:
                sources.append(Source(s["route"], concentration=conc))
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
                   leaks=float(leaks.get("share", 0.0)), leaks_to=leaks.get("to"), shares=shares,
                   effluent=_qualities(f"utility {name}: effluent", c.get("effluent", {}), base, (SEWER,)),
                   septic_effluent=_qualities(f"utility {name}: septic_effluent", c.get("septic_effluent", {}), base,
                                              (TAP,)))

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

    def add_return(self, returned: float, kg: dict | None = None):
        """A user's return to this sewer after the day: its withdrawal's delivered share of what it returned, and the
        mass it reports by species (a species it doesn't report returns at the system's mixed concentration [A])."""
        self._returns += returned
        self._return_parts.append((returned, dict(kg or {})))

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
            self._day.setdefault("got", {})[route] = got
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

    # ---------- mass (with transport) ----------
    def needs(self, users_return: bool) -> dict[str, bool]:
        """What the system's mass needs stated: an effluent concentration if its plant can have an outflow (a
        sewered class, a user returning to its sewer, or inflow and infiltration), a septic one for each class with
        septic water (``(class, needed)``)."""
        plant = users_return or self.infiltration > 0 or any(c.sewered > 0 for c in self.customers)
        return {"effluent": plant, **{f"septic:{c.name}": c.septic > 0 for c in self.customers}}

    def missing(self, species: list[str], users_return: bool, sources: dict[str, set]) -> list[str]:
        """The concentrations the system's mass needs and doesn't have, by species: ``sources``: route -> the species
        Basin Ledger can work out its water's concentration for (the rest need one stated)."""
        out = []
        need = self.needs(users_return)
        for sp in species:
            if need["effluent"] and sp not in self.effluent:
                out.append(f"effluent {sp}")
            for c in self.customers:
                if need[f"septic:{c.name}"] and sp not in c.septic_effluent and sp not in self.septic_effluent:
                    out.append(f"septic_effluent {sp} (class {c.name})")
            for s in self.sources:
                if sp not in s.concentration and sp not in sources.get(s.route, set()):
                    out.append(f"source {s.route} concentration {sp}")
        given = set(self.effluent) | set(self.septic_effluent) | {sp for s in self.sources for sp in s.concentration}
        given |= {sp for c in self.customers for sp in c.septic_effluent}
        out += [f"{sorted(given - set(species))}, which aren't transported"] if given - set(species) else []
        return out

    def _septic(self, c: Customers, sp: str) -> Quality:
        return c.septic_effluent.get(sp) or self.septic_effluent[sp]

    def mass_release(self, ledger, day: dt.date, species: list[str], scale: dict | None = None) -> dict:
        """At the day's start, after ``release``: the plant's outflow at today's effluent concentration (or its sewer's
        mass passed through), to the outfall and the systems reusing it in proportion to their water, and yesterday's
        septic water's and leaks' mass to the ground. ``scale``: species -> {"effluent", "septic"} multipliers on the
        stated concentrations (salt's groups). Returns {"outfall": {species: kg}, "ground": {route: {species: kg}}}."""
        k = f"utility:{self.name}:"
        out = ledger.today(k + "plant_outflow")
        out_kg, ground = {}, {}
        self._effluent_c = {}
        self._mass_day = {"held_start": {}, "released": {}}
        for sp in species:
            ks = k + f"{sp}:"
            held = self.held.get("kg", {}).get(sp, {"sewer": 0.0, "ground": {}})
            self._mass_day["held_start"][sp] = held["sewer"] + sum(held["ground"].values())
            q = self.effluent.get(sp)
            if q is None or q.keyword == SEWER:
                kg = held["sewer"]
                c = kg / out if out > 0 else 0.0
            else:
                c = q.on(day) * (scale or {}).get(sp, {}).get("effluent", 1.0)
                kg = c * out
            self._effluent_c[sp] = c
            reuse = 0.0
            for n, v in self._reuse_given.items():
                ledger.enter(ks + f"reuse:{n}:kg", c * v)
                reuse += c * v
            ledger.enter(ks + "plant_outflow:kg", kg)
            ledger.enter(ks + "plant:added", kg - held["sewer"])
            ledger.enter(ks + "reuse:kg", reuse)
            ledger.enter(ks + "outfall:kg", kg - reuse)
            out_kg[sp] = kg - reuse
            to_ground = 0.0
            for route, m in held["ground"].items():
                ledger.enter(ks + f"to_ground:{route}:kg", m)
                ground.setdefault(route, {})[sp] = m
                to_ground += m
            self._mass_day["released"][sp] = kg + to_ground
        return {"outfall": out_kg, "ground": ground}

    def effluent_concentration(self, sp: str) -> float:
        """Today's plant outflow concentration (``mass_release``), for the systems reusing it."""
        return self._effluent_c.get(sp, 0.0)

    def mass_settle(self, ledger, day: dt.date, species: list[str], c_sources: dict[str, dict[str, float]],
                    scale: dict | None = None):
        """After ``settle``: each species' mass by where the day's water went, and what is held overnight.
        ``c_sources``: route -> species -> its water's concentration today (kg m-3), stated ones already in it."""
        k = f"utility:{self.name}:"
        d, s = self._day, self._day["share"]
        got = d.get("got", {})
        delivered = s * d["demand"]
        losses = d["produced"] - delivered
        leaks = self.leaks * losses
        users = ledger.today(k + "users_delivered")
        held = {}
        for sp in species:
            ks = k + f"{sp}:"
            produced_kg = 0.0
            for route, v in got.items():
                m = v * c_sources[route][sp]
                ledger.enter(ks + f"source:{route}:kg", m)
                produced_kg += m
            c = produced_kg / d["produced"] if d["produced"] > 0 else 0.0
            self.concentration[sp] = c
            ledger.enter(ks + "produced:kg", produced_kg)
            ledger.enter(ks + "concentration", c)
            ledger.enter(ks + "delivered:kg", c * delivered)
            ledger.enter(ks + "losses:kg", c * losses)
            ledger.enter(ks + "leaks:kg", c * leaks)
            ledger.enter(ks + "unaccounted:kg", c * (losses - leaks))
            sewer, consumed, added, ground = 0.0, 0.0, 0.0, {}
            for cl in self.customers:
                w = s * d["classes"][cl.name]
                ledger.enter(ks + f"delivered:{cl.name}:kg", c * w)
                ledger.enter(ks + f"sewered:{cl.name}:kg", c * cl.sewered * w)
                ledger.enter(ks + f"consumed:{cl.name}:kg", c * cl.consumed * w)
                sep = 0.0
                if cl.septic > 0:
                    q = self._septic(cl, sp)
                    cs = c if q.keyword == TAP else q.on(day) * (scale or {}).get(sp, {}).get("septic", 1.0)
                    sep = cs * cl.septic * w
                    ground[cl.septic_to] = ground.get(cl.septic_to, 0.0) + sep
                ledger.enter(ks + f"septic:{cl.name}:kg", sep)
                ledger.enter(ks + f"septic_tap:{cl.name}:kg", c * cl.septic * w)  # the tap water's in it
                ledger.enter(ks + f"septic_added:{cl.name}", sep - c * cl.septic * w)
                sewer += c * cl.sewered * w
                consumed += c * cl.consumed * w
                added += sep - c * cl.septic * w
            if leaks > 0:
                ground[self.leaks_to] = ground.get(self.leaks_to, 0.0) + c * leaks
            back = sum(kg[sp] if sp in kg else c * v for v, kg in self._return_parts)
            ledger.enter(ks + "users_delivered:kg", c * users)
            ledger.enter(ks + "users_return:kg", back)
            ledger.enter(ks + "consumed:kg", consumed)
            sewer += back  # inflow and infiltration, water from outside, carry none [A]
            ledger.enter(ks + "sewer:kg", sewer)
            ledger.enter(ks + "held_to_ground:kg", sum(ground.values()))
            held[sp] = {"sewer": sewer, "ground": ground}
            # the system's closure: what it holds overnight changes by what came in and what left today
            start = self._mass_day.get("held_start", {}).get(sp, 0.0)
            released = self._mass_day.get("released", {}).get(sp, 0.0)
            plant_added = ledger.today(ks + "plant:added")
            gain = produced_kg + back + added + plant_added
            loss = c * users + consumed + c * (losses - leaks) + released
            end = sewer + sum(ground.values())
            gross = abs(gain) + abs(loss) + start + end
            ledger.enter(ks + "closure_share", ((end - start) - (gain - loss)) / max(gross, 1e-12))
        self.held["kg"] = held
        self.held["concentration"] = dict(self.concentration)  # a restart's users' intake on its first day
        self._return_parts = []

    def declare_mass(self, ledger, species: list[str]):
        k = f"utility:{self.name}:"
        for sp in species:
            ks = k + f"{sp}:"
            routes = list(dict.fromkeys(s.route for s in self.sources))
            ledger.pair(ks + "produced:kg", [ks + f"source:{r}:kg" for r in routes])
            ledger.pair(ks + "produced:kg", [ks + "delivered:kg", ks + "losses:kg"])
            ledger.pair(ks + "delivered:kg", [ks + f"delivered:{c.name}:kg" for c in self.customers]
                        + [ks + "users_delivered:kg"])
            for c in self.customers:
                ledger.pair(ks + f"delivered:{c.name}:kg", [ks + f"sewered:{c.name}:kg", ks + f"septic_tap:{c.name}:kg",
                                                            ks + f"consumed:{c.name}:kg"])
                ledger.pair(ks + f"septic:{c.name}:kg", [ks + f"septic_tap:{c.name}:kg", ks + f"septic_added:{c.name}"])
            ledger.pair(ks + "losses:kg", [ks + "leaks:kg", ks + "unaccounted:kg"])
            ledger.pair(ks + "plant_outflow:kg", [ks + "outfall:kg", ks + "reuse:kg"])
            ledger.zero(ks + "closure_share", atol=1e-12)

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
