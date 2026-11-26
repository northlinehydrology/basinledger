"""Scenarios from a shared warm-up (design.md; testing G5).

One warm-up runs the basin with no water users, and its state is saved at the pause after its last day
(``Basin.save``). Each scenario then carries on from that state (``Basin(..., restart=...)``), with or without the
users and with its own changes to the config, so every case starts from the same SWAT+ and MODFLOW 6 states and
differs only from the end of the warm-up. A case with no change matches the unbroken run with no users, bit for bit
(``tests/test_scenarios.py``).

The config the warm-up uses must already hold every case's MODFLOW 6 wells and SWAT+ objects (at zero), since a
restart carries on the same models: a case changes users, rules and rates, not the models' structure.
"""

from __future__ import annotations

import copy
import datetime as dt
import hashlib
import json
import os
import shutil
import tomllib
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from . import config as C
from .attribution import BASELINE
from .balance import soil_water_trend
from .coupling.swatplus import StartWarning
from .driver import Basin
from .temperature import print_daily

SOIL_SETTLED_SHARE = 0.12  # [A] above what weather alone moves a settled year (10.4%, a study), below a fill (15%)


def warm_up(cfg: C.ExampleConfig, days: int, run: str = "warmup", prepare=None, wrap_surface=None,
            soil_settled_share: float = SOIL_SETTLED_SHARE, users: bool = False) -> Path:
    """Run ``days`` with no users (or, ``users``, with the example's) and save the state; returns the state's folder.

    With SWAT+, the warm-up also reports whether its soils have settled. ``warm_up.json`` in the state's folder holds
    the basin's soil water change over the warm-up's last year (all of it, if shorter) as a share of that time's
    precipitation (``balance.soil_water_trend``). A share above ``soil_settled_share`` is a warning
    (``StartWarning``), not a stop: a soil still filling or draining carries its trend into every case."""
    swat = cfg.surface == "swatplus"

    def prep(model, d):
        if model == "swatplus":
            print_daily(d, ["basin_wb"])
        if prepare:
            prepare(model, d)

    b = Basin(cfg, run=run, users=users, days=None, prepare=prep if swat else prepare, wrap_surface=wrap_surface)
    try:
        b.run(days)
        state = b.folder / "state"
        b.save(state)
    finally:
        b.finalize()
    if swat:
        trend = soil_water_trend(b.folder / "swatplus", min(days, 365))
        trend |= {"limit": soil_settled_share, "settled": bool(trend["share"] <= soil_settled_share)}
        (state / "warm_up.json").write_text(json.dumps(trend, indent=1))
        if not trend["settled"]:
            warnings.warn(f"the warm-up's soils haven't settled: soil water changed {trend['soil_water_change_mm']:+.0f}"
                          f" mm over its last {trend['days']} days, {trend['share']:.0%} of their precipitation "
                          f"(limit {soil_settled_share:.0%}); every case carries the trend, so warm up longer or "
                          "start the soils wetter", StartWarning, stacklevel=2)
    return state


def run_case(cfg: C.ExampleConfig, state: Path, name: str, users: bool = True,
             change: Callable[[C.ExampleConfig], None] | None = None, prepare=None, wrap_surface=None) -> Basin:
    """One scenario from the saved warm-up: ``change(cfg)`` edits a copy of the config first (users, rules, rates).
    Returns the finished Basin, with its outputs written; the caller finalizes it."""
    c = copy.deepcopy(cfg)
    if change:
        change(c)
    b = Basin(c, run=name, users=users, restart=state, prepare=prepare, wrap_surface=wrap_surface)
    b.run()
    b.write_outputs()
    return b


# ---------- scenarios of water users (docs/design.md, "Scenarios of water users") ----------

SET_KEYS = {"name", "description", "example", "folder", "warm_up", "run", "compare", "scenario", "combination"}
SCENARIO_KEYS = {"name", "description", "users", "wells", "ground", "places_of_use", "point_sources", "utilities",
                 "scale", "remove", "rules", "transfers"}
UTILITY_KEYS = {"name", "sources", "shares", "losses", "leaks", "customers", "outfall", "demand", "sewered",
                "infiltration"}


@dataclass
class Scenario:
    """One scenario's changes, as read from the set's file (paths made absolute)."""
    name: str
    description: str = ""
    users: list = field(default_factory=list)          # C.UserConfig
    wells: list = field(default_factory=list)          # C.WellConfig
    ground: list = field(default_factory=list)         # C.WellConfig
    places_of_use: list = field(default_factory=list)  # (name, channel or "res:<reservoir>")
    point_sources: list = field(default_factory=list)  # (name, channel or "res:<reservoir>")
    transfers: list = field(default_factory=list)      # [[transfers]] tables (operations.py)
    utilities: list = field(default_factory=list)      # raw tables
    scale: list = field(default_factory=list)          # {user|utility, schedule: [Phase], customers}
    remove: list = field(default_factory=list)         # "user:<name>", "utility:<name>"
    rules: list = field(default_factory=list)          # raw tables
    raw: dict = field(default_factory=dict)            # the table as written, for the fingerprint


@dataclass
class ScenarioSet:
    path: Path
    name: str
    description: str
    example: Path
    folder: Path
    warm_up_days: int
    warm_up_users: bool
    warm_up_state: Path | None
    days: int
    workers: int
    nice: int
    checkpoint_days: int
    compare: dict
    scenarios: dict[str, Scenario]
    combinations: dict[str, list[str]]

    @property
    def cases(self) -> list[str]:
        return [BASELINE, *self.scenarios, *self.combinations]

    def members(self, case: str) -> list[Scenario]:
        if case == BASELINE:
            return []
        if case in self.scenarios:
            return [self.scenarios[case]]
        if case in self.combinations:
            return [self.scenarios[n] for n in self.combinations[case]]
        raise KeyError(f"{self.path}: no case {case!r} (cases: {', '.join(self.cases)})")

    @property
    def state(self) -> Path:
        return self.warm_up_state or self.folder / "warm_up" / "state"


def _refuse(path, msg):
    raise ValueError(f"{path}: {msg}")


def _absolute_paths(base: Path, table: dict, keys: tuple[str, ...]) -> dict:
    t = copy.deepcopy(table)
    for k in keys:
        if isinstance(t.get(k), str):
            t[k] = str(C._path(base, t[k]))
    return t


def _utility(path: Path, u: dict) -> dict:
    if not isinstance(u, dict) or "name" not in u or set(u) - UTILITY_KEYS:
        _refuse(path, f"a scenario's utilities: each has name and any of {sorted(UTILITY_KEYS - {'name'})}; got "
                      f"{sorted(u) if isinstance(u, dict) else u!r}")
    t = _absolute_paths(path.parent, u, ("shares",))
    for c in t.get("customers", []):
        if isinstance(c, dict) and isinstance(c.get("demand"), str):
            c["demand"] = str(C._path(path.parent, c["demand"]))
    return t


def _scale(path: Path, where: str, x: dict) -> dict:
    from .phases import scale_schedule
    if not isinstance(x, dict) or ("user" in x) == ("utility" in x) or \
            set(x) - {"user", "utility", "multiplier", "from", "phases", "customers"}:
        _refuse(path, f"{where}: a scale names one user or utility and has multiplier (and from) or phases; got {x!r}")
    if "customers" in x and "utility" not in x:
        _refuse(path, f"{where}: customers are a system's classes; this scale is a user's")
    start = None
    if "from" in x:
        if "phases" in x:
            _refuse(path, f"{where}: from goes with multiplier; phases have their own dates")
        start = dt.date.fromisoformat(str(x["from"]))
    target = ("user", x["user"]) if "user" in x else ("utility", x["utility"])
    return {"target": target, "schedule": scale_schedule(f"{path}: {where}", x.get("multiplier"), start,
                                                         x.get("phases")),
            "customers": list(x["customers"]) if "customers" in x else None}


def _scenario(path: Path, s: dict) -> Scenario:
    if not isinstance(s, dict) or "name" not in s:
        _refuse(path, "each [[scenario]] has a name")
    where = f"scenario {s['name']!r}"
    if set(s) - SCENARIO_KEYS:
        _refuse(path, f"{where}: unknown keys {sorted(set(s) - SCENARIO_KEYS)} (it takes {sorted(SCENARIO_KEYS)})")
    if s["name"] in (BASELINE,):
        _refuse(path, f"{where}: {BASELINE!r} is the example itself")
    base = path.parent
    out = Scenario(s["name"], s.get("description", ""), raw=copy.deepcopy(s))
    out.users = [C._user(path, base, u) for u in s.get("users", [])]
    out.wells = [C._well(path, w) for w in s.get("wells", [])]
    out.ground = [C._ground(path, g) for g in s.get("ground", [])]
    out.places_of_use = [C._place_of_use(path, x) for x in s.get("places_of_use", [])]
    out.point_sources = [C._point_source(path, x) for x in s.get("point_sources", [])]
    out.transfers = [C._transfer(path, x) for x in s.get("transfers", [])]
    out.utilities = [_utility(path, u) for u in s.get("utilities", [])]
    out.scale = [_scale(path, f"{where} scale", x) for x in s.get("scale", [])]
    for r in s.get("remove", []):
        if not isinstance(r, str) or r.split(":", 1)[0] not in ("user", "utility") or ":" not in r:
            _refuse(path, f"{where}: remove lists \"user:<name>\" or \"utility:<name>\"; got {r!r}")
        out.remove.append(r)
    for r in s.get("rules", []):
        if r.get("kind") == "drawdown_limit":
            C._cell(path, f"{where} rule {r.get('name', r['kind'])!r}", r.get("cell"))
        out.rules.append(copy.deepcopy(r))
    return out


def load_set(path: str | Path) -> ScenarioSet:
    """A set of scenarios: its TOML file read and checked (the scenarios' changes aren't applied here)."""
    path = Path(path).resolve()
    with open(path, "rb") as f:
        c = tomllib.load(f)
    if set(c) - SET_KEYS:
        _refuse(path, f"unknown keys {sorted(set(c) - SET_KEYS)}")
    for k in ("name", "example"):
        if k not in c:
            _refuse(path, f"a scenario set needs {k}")
    example = C._path(path.parent, c["example"])
    wu = c.get("warm_up", {})
    if set(wu) - {"days", "users", "state"} or ("days" in wu) == ("state" in wu):
        _refuse(path, "[warm_up] has days (and users), or state: a saved warm-up (Basin.save) used as is")
    run = c.get("run", {})
    if set(run) - {"days", "workers", "nice", "checkpoint_days"}:
        _refuse(path, f"[run] unknown keys {sorted(set(run) - {'days', 'workers', 'nice', 'checkpoint_days'})}")
    cmp_ = dict(c.get("compare", {}))
    known = {"channels", "heads", "from", "to", "low_flow", "exceedance", "units"}
    if set(cmp_) - known:
        _refuse(path, f"[compare] unknown keys {sorted(set(cmp_) - known)}")
    for h in cmp_.get("heads", []):
        C._cell(path, f"[compare] heads {h.get('name')!r}", h.get("cell"))
    lf = cmp_.get("low_flow", {})
    if set(lf) - {"days", "return_period", "year_start"}:
        _refuse(path, "[compare] low_flow has days, return_period and year_start (MM-DD)")
    if cmp_.get("units", "SI") not in ("SI", "US"):
        _refuse(path, "[compare] units is \"SI\" or \"US\"")
    scenarios: dict[str, Scenario] = {}
    for s in c.get("scenario", []):
        sc = _scenario(path, s)
        if sc.name in scenarios:
            _refuse(path, f"two scenarios named {sc.name!r}")
        scenarios[sc.name] = sc
    combos: dict[str, list[str]] = {}
    for x in c.get("combination", []):
        if not isinstance(x, dict) or set(x) - {"name", "scenarios", "description"} or "name" not in x \
                or len(x.get("scenarios", [])) < 2:
            _refuse(path, "each [[combination]] has a name and scenarios, two or more")
        if x["name"] in scenarios or x["name"] in combos or x["name"] == BASELINE:
            _refuse(path, f"combination {x['name']!r}: the name is taken")
        missing = [n for n in x["scenarios"] if n not in scenarios]
        if missing:
            _refuse(path, f"combination {x['name']!r}: no scenarios {missing}")
        if len(set(x["scenarios"])) != len(x["scenarios"]):
            _refuse(path, f"combination {x['name']!r}: a scenario listed twice")
        combos[x["name"]] = list(x["scenarios"])
    ex_folder = C.load(example).run_folder
    return ScenarioSet(
        path=path, name=c["name"], description=c.get("description", ""), example=example,
        folder=C._path(path.parent, c["folder"]) if "folder" in c else ex_folder / c["name"],
        warm_up_days=int(wu.get("days", 0)), warm_up_users=C._flag(path, "[warm_up] users", wu.get("users", True)),
        warm_up_state=C._path(path.parent, wu["state"]) if "state" in wu else None,
        days=int(run.get("days", 0)), workers=int(run.get("workers", 1)), nice=int(run.get("nice", 10)),
        checkpoint_days=int(run.get("checkpoint_days", 365)), compare=cmp_, scenarios=scenarios,
        combinations=combos)


def base_config(s: ScenarioSet) -> C.ExampleConfig:
    """The example with every scenario's model parts (wells, returns to the ground, intakes, point sources) and the
    comparison's outputs: the models every case and the warm-up share."""
    cfg = C.load(s.example)
    cfg.run_folder = s.folder / "cases"
    if s.days:
        cfg.days = s.days
    cfg.outputs = copy.deepcopy(cfg.outputs)

    def union(kind, have, new, key):
        names = {key(x): x for x in have}
        for sc in s.scenarios.values():
            for x in getattr(sc, new):
                k = key(x)
                if k in names and names[k] != x:
                    raise ValueError(f"{s.path}: scenario {sc.name!r}: {kind} {k!r} is defined differently "
                                     "elsewhere (the example or another scenario)")
                if k not in names:
                    names[k] = x
                    have.append(x)

    union("well", cfg.wells, "wells", lambda w: w.name)
    union("return to the ground", cfg.ground, "ground", lambda w: w.name)
    union("place of use", cfg.places_of_use, "places_of_use", lambda x: x[0])
    union("point source", cfg.point_sources, "point_sources", lambda x: x[0])
    if s.compare.get("heads"):
        heads = cfg.outputs.setdefault("heads", [])
        have = {h["name"] for h in heads}
        for h in s.compare["heads"]:
            if h["name"] not in have:
                heads.append(dict(h))
            elif next(x for x in heads if x["name"] == h["name"]) != h:
                raise ValueError(f"{s.path}: [compare] heads {h['name']!r}: the example has a head output by that "
                                 "name elsewhere")
    if s.compare.get("channels") and "channels" in cfg.outputs:
        cfg.outputs["channels"] = list(dict.fromkeys([*cfg.outputs["channels"], *s.compare["channels"]]))
    return cfg


def case_config(s: ScenarioSet, case: str) -> tuple[C.ExampleConfig, dict]:
    """A case's config: the shared models (``base_config``) with its scenarios' changes, and each system's demand
    scale (system -> [(schedule, customers)]), which ``run_set``'s worker puts on the started systems."""
    cfg = base_config(s)
    members = s.members(case)
    where = f"{s.path}: case {case!r}"
    users = {u.name: u for u in cfg.users}
    utilities = {u["name"]: u for u in cfg.utilities}
    changed: dict[str, str] = {}  # what a scenario changes -> the scenario, so two don't change the same thing

    def claim(thing, sc):
        if thing in changed:
            raise ValueError(f"{where}: {changed[thing]!r} and {sc.name!r} both change {thing}")
        changed[thing] = sc.name

    for sc in members:
        for u in sc.users:
            if u.name in users:
                raise ValueError(f"{where}: scenario {sc.name!r} adds user {u.name!r}, which is already there")
            claim(f"user:{u.name}", sc)
            users[u.name] = u
            cfg.users.append(u)
        for u in sc.utilities:
            claim(f"utility:{u['name']}", sc)
            if u["name"] in utilities:
                old = utilities[u["name"]]
                if "sources" in u and "shares" not in u:
                    old.pop("shares", None)
                old.update(copy.deepcopy({k: v for k, v in u.items() if k != "name"}))
            else:
                if "sources" not in u or "outfall" not in u:
                    raise ValueError(f"{where}: scenario {sc.name!r} adds system {u['name']!r}: it needs sources "
                                     "and an outfall")
                new = copy.deepcopy(u)
                utilities[u["name"]] = new
                cfg.utilities.append(new)
        cfg.rules.extend(copy.deepcopy(sc.rules))
        for x in sc.transfers:
            if any(y["name"] == x["name"] for y in cfg.transfers):
                raise ValueError(f"{where}: scenario {sc.name!r} adds transfer {x['name']!r}, which is already there")
            claim(f"transfer:{x['name']}", sc)
            cfg.transfers.append(dict(x))
    from .phases import combine, scale_schedule
    system_scales: dict[str, list] = {}
    for sc in members:
        scales = list(sc.scale) + [{"target": tuple(r.split(":", 1)), "schedule": scale_schedule(r, 0.0),
                                    "customers": None} for r in sc.remove]
        for x in scales:
            k, name = x["target"]
            claim(f"{k}:{name}" + (":scale" if f"{k}:{name}" in changed and changed[f"{k}:{name}"] == sc.name
                                   else ""), sc)
            if k == "user":
                if name not in users:
                    raise ValueError(f"{where}: scenario {sc.name!r} scales user {name!r}, which isn't there")
                u = copy.copy(users[name])
                u.phases = combine(u.phases, x["schedule"])
                cfg.users[cfg.users.index(users[name])] = u
                users[name] = u
            else:
                if name not in utilities:
                    raise ValueError(f"{where}: scenario {sc.name!r} scales system {name!r}, which isn't there")
                system_scales.setdefault(name, []).append((x["schedule"], x["customers"]))
    for sc in members:  # a system removed takes its users with it, named
        for r in sc.remove:
            k, name = r.split(":", 1)
            if k != "utility":
                continue
            on = [u.name for u in cfg.users if u.withdrawal == f"utility:{name}"]
            kept = [n for n in on if f"user:{n}" not in {x for m in members for x in m.remove}]
            if kept:
                raise ValueError(f"{where}: scenario {sc.name!r} removes system {name!r}, but its users {kept} "
                                 "aren't removed with it")
    return cfg, system_scales


class _ScaledCustomers:
    """A system's class of customers with its demand scaled by a schedule of multipliers (``phases.Phase``)."""

    def __init__(self, base, schedules):
        self._base, self._schedules = base, schedules
        for k in ("name", "sewered", "septic", "consumed", "septic_to"):
            setattr(self, k, getattr(base, k))

    def on(self, day: dt.date) -> float:
        f = 1.0
        for sched in self._schedules:
            cur = sched[0]
            for p in sched:
                if p.start is None or p.start <= day:
                    cur = p
            f *= cur.multiplier
        return f * self._base.on(day)


def scale_systems(b: Basin, system_scales: dict):
    """Put a case's demand scales on its started systems' customers."""
    for name, scales in system_scales.items():
        chain = b.utilities[name]
        from .operations import Restricted
        for i, c in enumerate(chain.customers):
            mine = [sched for sched, only in scales if only is None or c.name in only]
            if mine and isinstance(c, Restricted):  # the demand scaled, then restricted
                c._base = _ScaledCustomers(c._base, mine)
            elif mine:
                chain.customers[i] = _ScaledCustomers(c, mine)
        for _, only in scales:
            missing = sorted(set(only or []) - {c.name for c in chain.customers})
            if missing:
                raise ValueError(f"system {name}: no customers {missing}")


# ---------- fingerprints: what a finished case was run from ----------

def _digest(x) -> str:
    return hashlib.sha256(json.dumps(x, sort_keys=True, default=str).encode()).hexdigest()


def _file_digest(p) -> str | None:
    p = Path(p) if p else None
    return hashlib.sha256(p.read_bytes()).hexdigest() if p and p.is_file() else None


def warm_up_fingerprint(s: ScenarioSet) -> str:
    cfg = C.load(s.example)
    return _digest({
        "example": _file_digest(s.example), "parameters": _file_digest(cfg.calibration_parameters),
        "users": [_file_digest(u.config) for u in cfg.users] if s.warm_up_users else None,
        "parts": [sc.raw.get(k) for sc in s.scenarios.values()
                  for k in ("wells", "ground", "places_of_use", "point_sources")],
        "outputs": {k: s.compare.get(k) for k in ("channels", "heads")},
        "warm_up": [s.warm_up_days, s.warm_up_users, str(s.warm_up_state)], "days": s.days})


def case_fingerprint(s: ScenarioSet, case: str) -> str:
    members = s.members(case)
    return _digest({"warm_up": warm_up_fingerprint(s), "case": case,
                    "changes": [m.raw for m in members],
                    "files": [_file_digest(u.config) for m in members for u in m.users]})


# ---------- the runner ----------

def _done(folder: Path, fingerprint: str) -> bool:
    f = folder / "done.json"
    return f.exists() and json.loads(f.read_text()).get("fingerprint") == fingerprint


def _worker_start(nice: int):
    if nice and hasattr(os, "nice"):  # POSIX; on Windows the workers run at the runner's priority
        os.nice(nice)
    warnings.simplefilter("ignore")  # each case's start-up notes are the example's, said once by the warm-up's log


def _warm_up_task(set_path: str) -> dict:
    s = load_set(set_path)
    cfg = base_config(s)
    cfg.run_folder = s.folder
    state = warm_up(cfg, s.warm_up_days, run="warm_up", users=s.warm_up_users)
    day = json.loads((state / "basin.json").read_text())["day"]
    info = {"fingerprint": warm_up_fingerprint(s), "day": day}
    (s.folder / "warm_up" / "done.json").write_text(json.dumps(info))
    return info


def _case_task(set_path: str, case: str) -> dict:
    s = load_set(set_path)
    cfg, system_scales = case_config(s, case)
    fp = case_fingerprint(s, case)
    ckpt = s.folder / "checkpoints" / case
    saved = json.loads((ckpt / "checkpoint.json").read_text()) if (ckpt / "checkpoint.json").exists() else {}
    resume = saved.get("fingerprint") == fp
    b = Basin(cfg, run=case, users=True, restart=ckpt if resume else s.state)
    try:
        scale_systems(b, system_scales)
        while b.day < b.days:
            n = b.days - b.day if not s.checkpoint_days else min(s.checkpoint_days, b.days - b.day)
            b.run(n)
            if b.day < b.days and s.checkpoint_days:
                new, old = ckpt.with_name(case + ".new"), ckpt.with_name(case + ".old")
                b.save(new)
                (new / "checkpoint.json").write_text(json.dumps({"fingerprint": fp, "day": b.day}))
                if ckpt.exists():
                    ckpt.rename(old)
                new.rename(ckpt)
                shutil.rmtree(old, ignore_errors=True)
        b.write_outputs()
    finally:
        b.finalize()
    info = {"fingerprint": fp, "days": b.day, "resumed_from": saved["day"] if resume else None}
    (b.folder / "done.json").write_text(json.dumps(info))
    shutil.rmtree(ckpt, ignore_errors=True)
    return info


def run_set(path: str | Path, workers: int | None = None, cases: list[str] | None = None, fresh: bool = False,
            log: Callable[[str], None] = print) -> Path:
    """Run a set of scenarios (``load_set``): the warm-up once, then the baseline, each scenario and each
    combination in parallel processes, each skipped if done from the same inputs, resumed from its checkpoint if
    stopped; then the tables (``attribution.write_tables``). Returns the tables' folder."""
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor, as_completed
    s = load_set(path)
    todo = cases or s.cases
    for c in todo:
        s.members(c)
    s.folder.mkdir(parents=True, exist_ok=True)
    if fresh:
        for sub in ("cases", "checkpoints", "tables") + (() if s.warm_up_state else ("warm_up",)):
            shutil.rmtree(s.folder / sub, ignore_errors=True)
    n = workers or s.workers
    ctx = mp.get_context("spawn")  # libmf6 and SWAT+ hold one model a process: a fresh process for each run
    failed = {}
    with ProcessPoolExecutor(max_workers=max(1, n), mp_context=ctx, initializer=_worker_start,
                             initargs=(s.nice,), max_tasks_per_child=1) as pool:
        if s.warm_up_state is None:
            fp = warm_up_fingerprint(s)
            if _done(s.folder / "warm_up", fp):
                log(f"warm-up: done ({s.warm_up_days} days)")
            else:
                if (s.folder / "warm_up" / "done.json").exists():
                    log("warm-up: the example, the scenarios' model parts or the points compared changed, so it "
                        "and every case run again")
                    shutil.rmtree(s.folder / "cases", ignore_errors=True)
                    shutil.rmtree(s.folder / "checkpoints", ignore_errors=True)
                log(f"warm-up: running {s.warm_up_days} days")
                pool.submit(_warm_up_task, str(s.path)).result()
        jobs = {}
        for c in todo:
            if _done(s.folder / "cases" / c, case_fingerprint(s, c)):
                log(f"{c}: done")
                continue
            log(f"{c}: running")
            jobs[pool.submit(_case_task, str(s.path), c)] = c
        for f in as_completed(jobs):
            c = jobs[f]
            try:
                info = f.result()
                log(f"{c}: finished" + (f" (resumed from day {info['resumed_from']})" if info["resumed_from"] else ""))
            except Exception as e:  # one case's failure doesn't stop the others
                failed[c] = f"{type(e).__name__}: {e}"
                log(f"{c}: FAILED: {failed[c]}")
    (s.folder / "status.json").write_text(json.dumps({"failed": failed}, indent=1))
    from .attribution import write_tables
    out = write_tables(s)
    if failed:
        raise RuntimeError(f"{len(failed)} case(s) failed: " + "; ".join(f"{c}: {e}" for c, e in failed.items()))
    return out
