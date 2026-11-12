"""Calibration: parameters, observations, the forward model, statistics, and the calibration record.

Built:
- the record and its status (no record, or a status other than "calibrated", means a screening run);
- parameters and observations as tables (CSV);
- simulated values at observations, from a run's series (channel flow, heads at named cells);
- statistics: NSE, KGE, percent bias, a low-flow NSE (on log flows), mean error and RMSE;
- the forward model: run an example and write simulated values at every observation, for PEST++ to read;
- applying a parameter set to a run's copy of the models' files before they start (``apply_parameters``): SWAT+
  through its own ``calibration.cal`` (swatplus-bmi 0.9.0, which refuses changes SWAT+ 62.0.1 wouldn't make as
  asked), MODFLOW 6 as multipliers on package values, everywhere, by layer or by zone (flopy);
- the example's own parameter set, applied on every run (``[calibration] parameters``), with the forward model's
  values in place of any it names;
- PEST++ IES files (``write_pest_files``): the control file, a template for the parameter table and an instruction
  file for the simulated values; the forward model as the model command (``python -m basinledger forward``).

Not built yet (docs/design.md, Calibration):
- water users' parameters (set from the applicant's data in each user's config, not calibrated);
- MODFLOW 6 values by pilot point (pyEMU);
- fetching USGS groundwater levels, matched to MODFLOW 6 cells.

Tables:

    parameters.csv    name, model, target, value, lower, upper, transform
                      model: swatplus | modflow6 | user; transform: none | log; target: what the value changes
                      swatplus: <calibration.cal name>:<change>[:objects=1,2,...][:layers=1-2][:<condition>=<value>]
                        change absval | abschg | pctchg | relchg; conditions hsg, texture, landuse (HRUs only);
                        e.g. cn2:pctchg, awc:absval:layers=1-2, cn2:pctchg:hsg=D
                      modflow6: <package>.<value>:<where>, the value multiplied by the parameter's value
                        a griddata array (npf.k, npf.k22, npf.k33, sto.sy, sto.ss) or a column of a list package (riv.cond,
                        drn.cond, ghb.cond, sfr.rhk, sfr.man); where all | layer=<n> (0-based) | zone=<id>
                        (``[calibration] zones``: a CSV of cells, layer, row, col or layer, icell2d, and zone)
    observations.csv  name, kind, object, date, value, weight, source
                      kind: flow (m3 s-1, object = channel or hydrofabric id) | head (m, object = a named head
                      output); date ISO; source: e.g. "USGS 04119055 dv"
"""

from __future__ import annotations

import csv
import tomllib
from dataclasses import dataclass
from pathlib import Path

import numpy as np

STATUSES = ("calibrated", "screening")
SCREENING = "screening"


def model_statuses(record: Path | None, models: list[str]) -> dict[str, str]:
    """Each model's status from an example's record.

    A record gives either one status for the example (top-level ``status``, ``[calibration]``, ``[validation]``)
    or one per model (``[models.<name>]`` with the same keys), e.g. streamflow from NextGen calibrated and MODFLOW 6
    a screening model. A model the record doesn't name is screening. A calibrated model needs its calibration and
    validation statistics.
    """
    if record is None:
        return {m: SCREENING for m in models}
    with open(record, "rb") as f:
        rec = tomllib.load(f)
    per = rec.get("models")
    out = {}
    for m in models:
        entry = per.get(m) if per is not None else rec
        status = (entry or {}).get("status", SCREENING)
        if status not in STATUSES:
            raise ValueError(f"{record}: status must be one of {STATUSES}, not {status!r}")
        if status == "calibrated":
            missing = [k for k in ("calibration", "validation") if k not in entry]
            if missing:
                raise ValueError(f"{record}: a calibrated record needs its {' and '.join(missing)} statistics")
        out[m] = status
    return out


def record_targets(record: Path | None) -> dict[str, str]:
    """What each model in a record was calibrated to (its ``target``), e.g. streamflow at a gage."""
    if record is None:
        return {}
    with open(record, "rb") as f:
        rec = tomllib.load(f)
    return {m: e["target"] for m, e in (rec.get("models") or {}).items() if "target" in e}


def record_status(record: Path | None, models: list[str] | None = None) -> str:
    """'calibrated' only if every model the example runs is calibrated; otherwise a screening run."""
    st = model_statuses(record, models or ["all"])
    return "calibrated" if all(v == "calibrated" for v in st.values()) else SCREENING


@dataclass
class Parameter:
    name: str
    model: str
    target: str
    value: float
    lower: float
    upper: float
    transform: str = "none"

    def __post_init__(self):
        if self.model not in ("swatplus", "modflow6", "user"):
            raise ValueError(f"parameter {self.name}: model must be swatplus, modflow6 or user")
        if not self.lower <= self.value <= self.upper:
            raise ValueError(f"parameter {self.name}: {self.value} outside [{self.lower}, {self.upper}]")


@dataclass
class Observation:
    name: str
    kind: str       # flow | head
    object: str     # channel (or hydrofabric id) for flow; named head output for head
    date: str
    value: float
    weight: float = 1.0
    source: str = ""


def read_parameters(path: Path) -> list[Parameter]:
    with open(path, newline="") as f:
        return [Parameter(r["name"], r["model"], r["target"], float(r["value"]), float(r["lower"]),
                          float(r["upper"]), r.get("transform") or "none") for r in csv.DictReader(f)]


def read_observations(path: Path) -> list[Observation]:
    with open(path, newline="") as f:
        obs = [Observation(r["name"], r["kind"], r["object"], r["date"], float(r["value"]),
                           float(r.get("weight") or 1.0), r.get("source", "")) for r in csv.DictReader(f)]
    for o in obs:
        if o.kind not in ("flow", "head"):
            raise ValueError(f"observation {o.name}: kind must be flow or head")
    return obs


def merged(base: list[Parameter], over: list[Parameter]) -> list[Parameter]:
    """``base`` with ``over``'s values in place of any it names, and ``over``'s others after."""
    names = {p.name for p in over}
    return [p for p in base if p.name not in names] + list(over)


def apply_parameters(model: str, folder: Path, params: list[Parameter], zones: Path | None = None,
                     gwf: str = "gwf"):
    """Write a parameter set into a run's copy of one model's inputs, before it starts (the driver's ``prepare``
    hook): ``swatplus`` into the run's ``calibration.cal``, ``modflow6`` into the packages of model ``gwf``.
    Water users' parameters aren't built."""
    mine = [p for p in params if p.model == model]
    if not mine:
        return
    if model == "swatplus":
        from swatplus_bmi.calibration import write_changes
        write_changes(folder, [swatplus_change(p) for p in mine])
    elif model == "modflow6":
        apply_modflow6(Path(folder), mine, zones, gwf)
    else:
        raise NotImplementedError(f"{model} parameters: set them in the user's config (from the applicant's data); "
                                  "calibrating them isn't built")


def swatplus_change(p: Parameter):
    """A SWAT+ parameter's ``target`` (module docstring) as a calibration.cal row."""
    from swatplus_bmi.calibration import CONDITIONS, Change
    name, *rest = p.target.split(":")
    if not rest or "=" in rest[0]:
        raise ValueError(f"parameter {p.name}: a SWAT+ target is <name>:<change>[:...], not {p.target!r}")
    change, objects, layers, conditions = rest[0], (), (0, 0), []
    for part in rest[1:]:
        key, _, val = part.partition("=")
        if key == "objects":
            objects = tuple(int(i) for i in val.split(","))
        elif key == "layers":
            first, _, last = val.partition("-")
            layers = (int(first), int(last or first))
        elif key in CONDITIONS:
            conditions.append((key, val))
        else:
            raise ValueError(f"parameter {p.name}: {key!r} in {p.target!r}; a SWAT+ target takes objects, layers "
                             f"and the conditions {CONDITIONS}")
    return Change(name, p.value, change, objects, layers, tuple(conditions))


MF6_ARRAYS = {"npf": ("k", "k22", "k33"), "sto": ("sy", "ss")}
MF6_COLUMNS = {"riv": ("cond",), "drn": ("cond",), "ghb": ("cond",), "sfr": ("rhk", "man")}


def read_zones(path: Path, shape: tuple[int, ...]) -> np.ndarray:
    """A zone per cell (-1 where none is given), in the grid's user shape: (nlay, nrow, ncol) or (nlay, ncpl)."""
    z = np.full(shape, -1, dtype=int)
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            cols = ("layer", "row", "col") if len(shape) == 3 else ("layer", "icell2d")
            if any(not (r.get(c) or "").strip() for c in cols):
                raise ValueError(f"{path}: a zone row gives {', '.join(cols)} and zone on this grid")
            z[tuple(int(r[c]) for c in cols)] = int(r["zone"])
    return z


def apply_modflow6(folder: Path, params: list[Parameter], zones: Path | None, gwf: str):
    """Each parameter's value multiplies a package value (module docstring) where its target says. Only the
    packages changed are written again, at full precision."""
    import flopy

    sim = flopy.mf6.MFSimulation.load(sim_ws=str(folder), verbosity_level=0)
    sim.simulation_data.float_precision, sim.simulation_data.float_characters = 17, 26
    m = sim.get_model(gwf)
    shape = tuple(int(n) for n in m.modelgrid.shape)
    zone = read_zones(zones, shape) if zones else None
    changed = {}
    for p in params:
        what, _, where = p.target.partition(":")
        pkg_name, _, value = what.partition(".")
        pkg = m.get_package(pkg_name)
        if pkg is None:
            raise ValueError(f"parameter {p.name}: model {gwf} has no package {pkg_name!r}")
        ptype = pkg.package_type
        key, _, arg = (where or "all").partition("=")
        if key == "all":
            mask = np.ones(shape, bool)
        elif key == "layer":
            if not 0 <= int(arg) < shape[0]:
                raise ValueError(f"parameter {p.name}: layer {arg} isn't in the grid's {shape[0]} (from 0)")
            mask = np.zeros(shape, bool)
            mask[int(arg)] = True
        elif key == "zone":
            if zone is None:
                raise ValueError(f"parameter {p.name}: zone={arg} needs [calibration] zones")
            mask = zone == int(arg)
        else:
            raise ValueError(f"parameter {p.name}: where is all, layer=<n> or zone=<id>, not {where!r}")
        if not mask.any():
            raise ValueError(f"parameter {p.name}: {where} has no cells")
        if value in MF6_ARRAYS.get(ptype, ()):
            arr = getattr(pkg, value)
            data = np.array(arr.array, dtype=float).reshape(shape)
            data[mask] *= p.value
            setattr(pkg, value, data)
        elif value in MF6_COLUMNS.get(ptype, ()):
            hit = 0
            blocks = ([pkg.packagedata] if ptype == "sfr" else [pkg.stress_period_data])
            for block in blocks:
                data = block.get_data()
                periods = data if isinstance(data, dict) else {None: data}
                for per, rows in periods.items():
                    if rows is None:
                        continue
                    rows = rows.copy()
                    sel = np.array([isinstance(c, tuple) and bool(mask[c]) for c in rows["cellid"]])
                    rows[value][sel] = rows[value][sel] * p.value
                    hit += int(sel.sum())
                    if per is None:
                        block.set_data(rows)
                    else:
                        block.set_data(rows, key=per)
            if not hit:
                raise ValueError(f"parameter {p.name}: {pkg_name} has no entries in {where or 'all'}")
        else:
            known = MF6_ARRAYS.get(ptype, ()) + MF6_COLUMNS.get(ptype, ())
            raise ValueError(f"parameter {p.name}: {ptype}.{value} isn't a value Basin Ledger multiplies"
                             + (f" ({ptype}: {', '.join(known)})" if known else ""))
        changed[pkg_name] = pkg
    for pkg in changed.values():
        pkg.write()


def simulated(series: dict[str, list[float]], dates: list[str], observations: list[Observation]) -> np.ndarray:
    """The run's value at each observation: channel flow (``flow:<object>``) or head (``head:<object>``) on its
    date."""
    index = {d: i for i, d in enumerate(dates)}
    out = np.empty(len(observations))
    for j, o in enumerate(observations):
        key = f"{o.kind}:{o.object}"
        if key not in series:
            raise KeyError(f"observation {o.name}: the run has no {key!r} output (add it under [outputs])")
        if o.date not in index:
            raise KeyError(f"observation {o.name}: {o.date} is outside the run")
        out[j] = series[key][index[o.date]]
    return out


def nse(sim, obs) -> float:
    sim, obs = np.asarray(sim, float), np.asarray(obs, float)
    return float(1 - np.sum((sim - obs) ** 2) / np.sum((obs - obs.mean()) ** 2))


def log_nse(sim, obs, floor: float = 1e-6) -> float:
    """NSE on log flows: weights low flows, where a large user shows."""
    return nse(np.log(np.maximum(sim, floor)), np.log(np.maximum(obs, floor)))


def kge(sim, obs) -> float:
    """Kling-Gupta efficiency (Gupta et al. 2009)."""
    sim, obs = np.asarray(sim, float), np.asarray(obs, float)
    r = np.corrcoef(sim, obs)[0, 1]
    return float(1 - np.sqrt((r - 1) ** 2 + (sim.std() / obs.std() - 1) ** 2 + (sim.mean() / obs.mean() - 1) ** 2))


def pbias(sim, obs) -> float:
    """Percent bias, positive when the model overestimates."""
    sim, obs = np.asarray(sim, float), np.asarray(obs, float)
    return float(100 * np.sum(sim - obs) / np.sum(obs))


def mean_error(sim, obs) -> float:
    return float(np.mean(np.asarray(sim, float) - np.asarray(obs, float)))


def rmse(sim, obs) -> float:
    return float(np.sqrt(np.mean((np.asarray(sim, float) - np.asarray(obs, float)) ** 2)))


def statistics(observations: list[Observation], sim: np.ndarray) -> dict[str, dict[str, float]]:
    """Statistics per site (kind and object): flow gets NSE, log NSE, KGE and percent bias; heads mean error, RMSE
    and RMSE over the observed range."""
    out = {}
    for key in dict.fromkeys((o.kind, o.object) for o in observations):
        idx = [i for i, o in enumerate(observations) if (o.kind, o.object) == key]
        s = sim[idx]
        o = np.array([observations[i].value for i in idx])
        if key[0] == "flow":
            out[f"flow:{key[1]}"] = {"n": len(idx), "nse": nse(s, o), "log_nse": log_nse(s, o), "kge": kge(s, o),
                                     "pbias": pbias(s, o)}
        else:
            rng = float(o.max() - o.min())
            out[f"head:{key[1]}"] = {"n": len(idx), "mean_error": mean_error(s, o), "rmse": rmse(s, o),
                                     "rmse_over_range": rmse(s, o) / rng if rng > 0 else float("nan")}
    return out


def run_forward_model(config: str | Path, parameters: Path | None, observations: Path, out: Path,
                      run: str = "calibration", run_folder: str | Path | None = None) -> Path:
    """The forward model PEST++ calls: apply a parameter set (the example's own, with these values in place of any
    it names), run the example, and write the simulated value at every observation (name, simulated) to ``out``.
    ``run_folder``: where the run's files go, instead of the example's, so parallel workers don't share one."""
    from . import config as C
    from .driver import Basin

    cfg = C.load(config)
    if run_folder is not None:
        cfg.run_folder = Path(run_folder).resolve()
    obs = read_observations(observations)
    params = read_parameters(parameters) if parameters else []
    b = Basin(cfg, run=run, parameters=params)
    try:
        b.run()
        sim = simulated(b.series, b.dates(), obs)
    finally:
        b.finalize()
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "simulated"])
        for o, s in zip(obs, sim):
            w.writerow([o.name, repr(float(s))])
    return out


PEST_NAME = 200  # PEST++'s longest parameter or observation name


def _pest_name(kind: str, name: str) -> str:
    if not name or len(name) > PEST_NAME or any(c.isspace() for c in name):
        raise ValueError(f"{kind} {name!r}: PEST++ names are 1 to {PEST_NAME} characters with no spaces")
    return name.lower()


def write_pest_files(config: str | Path, parameters: Path, observations: Path, folder: Path, reals: int = 50,
                     iterations: int = 3, python: str | None = None) -> Path:
    """PEST++ IES files for an example, in ``folder``: ``parameters.csv.tpl`` (the parameter table with each value
    a PEST++ marker), ``simulated.csv.ins`` (where each observation's value is in the forward model's output), a copy
    of the observations, and ``basin.pst``, whose model command runs the forward model in the worker's own folder
    (``python -m basinledger forward``). Parameters start at their values in ``parameters``; observations carry their
    weights, grouped by kind and object. ``reals``: the ensemble's size; ``iterations``: IES's iterations
    (NOPTMAX). Returns the control file."""
    import shutil
    import sys

    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    params, obs = read_parameters(parameters), read_observations(observations)
    pnames = [_pest_name("parameter", p.name) for p in params]
    onames = [_pest_name("observation", o.name) for o in obs]
    for kind, names in (("parameter", pnames), ("observation", onames)):
        dup = sorted({n for n in names if names.count(n) > 1})
        if dup:
            raise ValueError(f"{kind} names {dup} repeat (PEST++ doesn't tell case apart)")
    for p in params:
        if p.transform not in ("none", "log"):
            raise ValueError(f"parameter {p.name}: transform is none or log, not {p.transform!r}")
        if p.transform == "log" and p.lower <= 0:
            raise ValueError(f"parameter {p.name}: a log parameter's lower bound must be above 0")

    with open(folder / "parameters.csv.tpl", "w", newline="") as f:
        f.write("ptf ~\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["name", "model", "target", "value", "lower", "upper", "transform"])
        for p, n in zip(params, pnames):
            w.writerow([p.name, p.model, p.target, f"~{n:^24}~", p.lower, p.upper, p.transform])
    with open(folder / "simulated.csv.ins", "w") as f:
        f.write("pif ~\nl1\n" + "".join(f"l1 ~,~ !{n}!\n" for n in onames))
    shutil.copyfile(observations, folder / "observations.csv")

    groups = list(dict.fromkeys(p.model for p in params))
    ogroups = [f"{o.kind}_{o.object}".lower() for o in obs]
    py = python or sys.executable
    lines = [
        "pcf",
        "* control data",
        "restart estimation",
        f"{len(params)} {len(obs)} {len(groups)} 0 {len(set(ogroups))}",
        "1 1 single point 1 0 0",
        "10.0 -3.0 0.3 0.03 10",
        "10.0 10.0 0.001",
        "0.1",
        f"{iterations} 0.005 4 4 0.005 4",
        "1 1 1",
        "* parameter groups",
        *[f"{g} relative 0.01 0.0 switch 2.0 parabolic" for g in groups],
        "* parameter data",
        *[f"{n} {p.transform} {'factor' if p.transform == 'log' else 'relative'} {p.value!r} {p.lower!r} "
          f"{p.upper!r} {p.model} 1.0 0.0 1" for p, n in zip(params, pnames)],
        "* observation groups",
        *list(dict.fromkeys(ogroups)),
        "* observation data",
        *[f"{n} {o.value!r} {o.weight!r} {g}" for o, n, g in zip(obs, onames, ogroups)],
        "* model command line",
        f"{py} -m basinledger forward {Path(config).resolve()} parameters.csv observations.csv simulated.csv "
        "--run-folder runs",
        "* model input/output",
        "parameters.csv.tpl parameters.csv",
        "simulated.csv.ins simulated.csv",
        f"++ ies_num_reals({reals})",
    ]
    pst = folder / "basin.pst"
    pst.write_text("\n".join(lines) + "\n")
    return pst


USGS_DAILY = "https://api.waterdata.usgs.gov/ogcapi/v0/collections/daily/items"
CFS_TO_CMS = 0.0283168466


def usgs_daily(site: str, start: str, end: str, cache: Path | None = None) -> list[tuple[str, float, str]]:
    """USGS daily mean discharge (m3 s-1) with its approval status, from the USGS Water Data API's daily values
    (parameter 00060, statistic 00003), cached as CSV: [(date, flow, approval), ...]."""
    import json
    import urllib.request

    if cache and Path(cache).exists():
        with open(cache, newline="") as f:
            return [(r["date"], float(r["flow"]), r["approval"]) for r in csv.DictReader(f)]
    url = (f"{USGS_DAILY}?f=json&monitoring_location_id=USGS-{site}&parameter_code=00060&statistic_id=00003"
           f"&time={start}/{end}&limit=10000")
    rows = []
    while url:
        with urllib.request.urlopen(url, timeout=60) as r:
            page = json.load(r)
        rows += [(f["properties"]["time"][:10], float(f["properties"]["value"]) * CFS_TO_CMS,
                  f["properties"].get("approval_status") or "") for f in page["features"]
                 if f["properties"]["value"] is not None]
        url = next((link["href"] for link in page.get("links", []) if link["rel"] == "next"), None)
    rows.sort()
    if cache:
        Path(cache).parent.mkdir(parents=True, exist_ok=True)
        with open(cache, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["date", "flow", "approval"])
            w.writerows(rows)
    return rows


def recession_ratios(q: np.ndarray, run_days: int = 5) -> np.ndarray:
    """Ratios of successive days' flow on recession days: days in runs of at least ``run_days`` falling days, with
    flow below the median (as nextgen-michigan's scripts/recession.py)."""
    q = np.asarray(q, float)
    falling = np.r_[False, q[1:] < q[:-1]]
    run = np.zeros(len(q), int)
    for i in range(1, len(q)):
        run[i] = run[i - 1] + 1 if falling[i] else 0
    sel = (run >= run_days) & (q < np.median(q))
    return q[1:][sel[1:]] / q[:-1][sel[1:]]


def eckhardt(q: np.ndarray, alpha: float, bfi_max: float = 0.8) -> np.ndarray:
    """Base flow by the Eckhardt (2005) recursive filter. BFImax 0.8 is Eckhardt's value for perennial streams on
    porous aquifers [A]; alpha is the recession constant (here the 90th percentile of the recession ratios)."""
    q = np.asarray(q, float)
    b = np.zeros_like(q)
    b[0] = q[0]
    for i in range(1, len(q)):
        b[i] = min(((1 - bfi_max) * alpha * b[i - 1] + (1 - alpha) * bfi_max * q[i]) / (1 - alpha * bfi_max), q[i])
    return b


def fetch_usgs(*args, **kwargs):
    raise NotImplementedError("USGS groundwater level downloads are planned; see docs/design.md")
