"""Lakes: SWAT+'s reservoir holds a lake's water, MODFLOW 6's lake package (LAK) is its face to the aquifer, as SFR is
a stream's (docs/design.md, "Extending an existing MODFLOW 6 model", lakes; swatplus-bmi 0.11.0, request 6).

Each day, after SWAT+'s day, the lake in MODFLOW 6 is set from the reservoir and solved with the aquifer, and its
exchange with the aquifer is set on the reservoir for the next day, as a channel's is. Two forms, chosen per lake:

- **A, active within the day** (``form = "active"``): the lake starts the day at the stage of the reservoir's volume
  at the day's start (``XOLDPAK``), takes the day's inflow and precipitation (``RUNOFF``) and the day's outflow,
  evaporation and SWAT+'s own seepage (``WITHDRAWAL``), all in m3 as SWAT+ reported them, and MODFLOW 6 solves its
  stage and leakage within the day on the lake's own bathymetry and lakebed. What its stage did within the day isn't
  kept: tomorrow starts from SWAT+'s volume again, so the lake's water is held once, in SWAT+.
- **B, at a fixed stage** (``form = "fixed_stage"``): status constant, the stage that of the reservoir's volume at the
  end of the day; LAK computes the leakage at that stage.

SWAT+ keeps a reservoir's volume, not its level; the stage is LAK's own at that volume, from the lake's table if it
has one, else from its connections as MODFLOW 6 computes a lake's volume (``lak_calculate_vol`` in 6.8.2: each
connection's surface area, smoothed at its top, times the depth over its bottom), inverted by bisection. The lake's
own rainfall, evaporation, runoff, inflow and withdrawal from its files are replaced each day (SWAT+ owns that water),
and its leakage is SWAT+'s reservoir's exchange: positive into the reservoir, a loss limited by SWAT+ to what the
reservoir holds (the rest booked as a shortfall). Refused: a lake with outlets (water between water bodies, which
SWAT+ routes; cutting them isn't built), and a reservoir whose own seepage (``k`` in hydrology.res) isn't 0.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

VOLUME = "reservoir_water__volume"                                  # m3, the end of the day
INFLOW = "reservoir_water__volume_inflow_rate"                      # m3 s-1, the day just done
OUTFLOW = "reservoir_water__volume_outflow_rate"
PRECIP = "reservoir_water__precipitation_volume_flow_rate"
EVAP = "reservoir_water__evaporation_volume_flow_rate"
SEEP = "reservoir_water__seepage_volume_flow_rate"                  # SWAT+'s own, from k
GW = "reservoir_water__groundwater_volume_exchange_rate"            # m3 s-1, set for the next day, + into it
FORMS = ("active", "fixed_stage")


@dataclass
class LakeConfig:
    reservoir: str          # the SWAT+ reservoir that holds the lake's water
    package: str            # MODFLOW 6's lake package
    lake: int = 1           # the lake's number in the package, 1-based as in its file
    form: str = "active"    # "active" (A) or "fixed_stage" (B)


def eps_quadratic(top: np.ndarray, bot: np.ndarray, x: float, eps: float = 1e-6) -> np.ndarray:
    """MODFLOW 6's sQuadraticSaturation (SmoothingFunctions.f90, 6.8.2), elementwise."""
    b = top - bot
    out = np.where(x < bot, 0.0, 1.0)
    pos = b > 0
    br = np.clip((x - bot[pos]) / b[pos], 0.0, 1.0)
    av = 1.0 / (1.0 - eps)
    bri = 1.0 - br
    y = np.where(br < eps, av * 0.5 * br * br / eps,
                 np.where(br < 1.0 - eps, av * br + 0.5 * (1.0 - av),
                          np.where(br < 1.0, 1.0 - av * 0.5 * bri * bri / eps, 1.0)))
    out[pos] = y
    return out


class LakeCurve:
    """A lake's volume at a stage, as MODFLOW 6 computes it, and the stage at a volume."""

    def __init__(self, top=None, bot=None, sarea=None, table=None):
        self.top, self.bot, self.sarea = top, bot, sarea
        self.table = table  # (stage, volume, sarea) arrays, from the lake's table file
        if table is not None:
            self.lo, self.hi = float(table[0][0]), float(table[0][-1])
        else:
            self.lo, self.hi = float(np.min(bot)), float(np.max(top))

    def volume(self, stage: float) -> float:
        if self.table is not None:
            s, v, a = self.table
            if stage <= s[0]:
                return float(v[0])
            if stage >= s[-1]:
                return float(v[-1] + (stage - s[-1]) * a[-1])
            return float(np.interp(stage, s, v))
        sa = eps_quadratic(self.top, self.bot, stage) * self.sarea
        depth = np.where(stage < self.bot, 0.0, stage - self.bot)  # below or between top and bottom, or above both
        return float(np.sum(sa * depth))

    def stage(self, volume: float) -> float:
        """The stage holding ``volume`` (bisection; the volume rises with the stage). At no volume, the lake's bottom."""
        lo, hi = self.lo, self.hi
        if volume <= 0:
            return lo
        while self.volume(hi) < volume:  # above the highest top: the volume rises linearly with the whole area
            hi = lo + 2 * (hi - lo) + 1.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if self.volume(mid) < volume:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-10 * max(1.0, abs(hi)):
                break
        return 0.5 * (lo + hi)


def _table_file(mf6, package: str, lake: int) -> Path | None:
    """The lake's table file, from the LAK file's TABLES block (``<lake> TAB6 FILEIN <file>``), or None. The LAK file
    is the one the model's name file gives the package."""
    from .existing import _model_nam
    from .mover import _blocks
    nam, _, _ = _model_nam(mf6.ws, mf6.model)
    lak = next((mf6.ws / w[1].strip("'\"") for name, _, lines in _blocks(nam, mf6.ws) if name == "PACKAGES"
                for w in lines if len(w) > 2 and w[2].lower() == package.lower()), None)
    if lak is None:
        return None
    for name, _, lines in _blocks(lak, mf6.ws):
        if name == "TABLES":
            for w in lines:
                if len(w) >= 4 and w[0].isdigit() and int(w[0]) == lake and w[1].upper() == "TAB6":
                    return mf6.ws / w[3].strip("'\"")
    return None


def read_table(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A LAK table file: its TABLE block's stage, volume and surface area columns."""
    from .mover import _blocks
    rows = [w for name, _, lines in _blocks(path, path.parent) if name == "TABLE" for w in lines]
    a = np.array([[float(x) for x in w[:3]] for w in rows])
    return a[:, 0], a[:, 1], a[:, 2]


def curve(mf6, package: str, lake: int) -> LakeCurve:
    """The lake's volume-stage curve, from MODFLOW 6's memory (its connections), or its table if it has one."""
    p = package.upper()
    ntab = mf6.get_model_units("NTABROW", mf6.model, p)
    if int(ntab[lake - 1]) > 0:
        f = _table_file(mf6, package, lake)
        if f is None or not Path(f).exists():
            raise ValueError(f"{package}: lake {lake} has a table, but its file wasn't found in the LAK file's TABLES")
        s, v, a = read_table(Path(f))
        f_len = mf6.factor((1, 0))
        return LakeCurve(table=(s * f_len, v * f_len ** 3, a * f_len ** 2))
    idx = mf6.get_model_units("IDXLAKECONN", mf6.model, p).astype(int)
    first, last = idx[lake - 1] - 1, idx[lake] - 1
    top = mf6.get("TELEV", mf6.model, p)[first:last]
    bot = mf6.get("BELEV", mf6.model, p)[first:last]
    sarea = mf6.get("SAREA", mf6.model, p)[first:last]
    return LakeCurve(top, bot, sarea)


def connections(mf6, package: str, lake: int) -> slice:
    """The lake's connections among the package's entries (SIMVALS is one value per connection)."""
    idx = mf6.get_model_units("IDXLAKECONN", mf6.model, package.upper()).astype(int)
    return slice(idx[lake - 1] - 1, idx[lake] - 1)


def own_seepage(project: str | Path) -> dict[str, float]:
    """Each reservoir's bottom conductivity ``k`` (mm/hr, SWAT+'s own seepage), by its name in reservoir.con."""
    project = Path(project)

    def rows(name):
        f = project / name
        return [ln.split() for ln in f.read_text(errors="replace").splitlines()[2:] if ln.split()] if f.exists() else []

    res = {r[0]: r[3] for r in rows("reservoir.res") if len(r) > 3}            # id -> its hydrology's name
    hyd = {r[0]: float(r[7]) for r in rows("hydrology.res") if len(r) > 7}     # name -> k
    return {r[1]: hyd.get(res.get(r[7], ""), 0.0) for r in rows("reservoir.con") if len(r) > 7}


def problems(mf6, lakes: list[LakeConfig], reservoirs: list[str], project) -> list[str]:
    """What the start-up refuses in the lakes as configured."""
    bad = []
    k = own_seepage(project) if project else {}
    seen = set()
    for lk in lakes:
        where = f"[[lakes]] {lk.reservoir}"
        if lk.form not in FORMS:
            bad.append(f"{where}: form is one of {', '.join(FORMS)}")
        if lk.reservoir not in reservoirs:
            bad.append(f"{where}: no such reservoir in the SWAT+ project")
        elif k.get(lk.reservoir, 0.0) > 0:
            bad.append(f"{where}: its bottom conductivity k in hydrology.res is {k[lk.reservoir]:g} mm/hr, so SWAT+ "
                       "seeps the lake's water itself and LAK's exchange would take it again; set k to 0")
        if (lk.package.lower(), lk.lake) in seen:
            bad.append(f"{where}: lake {lk.lake} of {lk.package} is coupled twice")
        seen.add((lk.package.lower(), lk.lake))
        try:
            t = str(mf6.get("PACKAGE_TYPE", mf6.model, lk.package.upper())[0]).strip().lower()
        except Exception:  # noqa: BLE001  (no such package in the model)
            bad.append(f"{where}: the model has no package {lk.package!r}")
            continue
        if t != "lak":
            bad.append(f"{where}: {lk.package} is a {t.upper()} package, not a lake (LAK)")
            continue
        n = int(mf6.get_model_units("NLAKES", mf6.model, lk.package.upper())[0])
        if not 1 <= lk.lake <= n:
            bad.append(f"{where}: {lk.package} has lakes 1 to {n}; lake {lk.lake} isn't one")
            continue
        if int(mf6.get_model_units("NOUTLETS", mf6.model, lk.package.upper())[0]) > 0:
            bad.append(f"{where}: {lk.package} has outlets, water between water bodies, which SWAT+ routes; cutting "
                       "a coupled lake's outlets isn't built")
    return bad


class LakeLink:
    """The lakes coupled, day by day: the reservoirs' water read after SWAT+'s day, each lake set and solved with the
    aquifer, its leakage set on the reservoir for the next day."""

    def __init__(self, mf6, sw, lakes: list[LakeConfig]):
        self.mf6, self.sw, self.lakes = mf6, sw, lakes
        names = list(sw.names(VOLUME))
        self.index = [names.index(lk.reservoir) for lk in lakes]
        self.curves = [curve(mf6, lk.package, lk.lake) for lk in lakes]
        self.conn = [connections(mf6, lk.package, lk.lake) for lk in lakes]
        self.due = {lk.reservoir: 0.0 for lk in lakes}         # m3 for the next day, positive into the reservoir
        self.start_volume = sw.get(VOLUME)[self.index].astype(float)  # the reservoirs at the first day's start
        self.day: dict[str, dict[str, float]] = {}

    def before_swat(self):
        """Yesterday's leakage, on the reservoirs for today (SWAT+ takes a loss up to what a reservoir holds)."""
        names = list(self.sw.names(GW))
        q = np.zeros(len(names))
        for lk, i in zip(self.lakes, self.index):
            q[i] = self.due[lk.reservoir] / 86400.0
        self.sw.set(GW, q)

    def after_swat(self):
        """The reservoirs' day, as SWAT+ reports it, m3."""
        got = {v: self.sw.get(v) for v in (VOLUME, INFLOW, OUTFLOW, PRECIP, EVAP, SEEP, GW)}
        self.day = {}
        for lk, i, v0 in zip(self.lakes, self.index, self.start_volume):
            self.day[lk.reservoir] = {"start": float(v0), "volume": float(got[VOLUME][i]),
                                      **{k: float(got[v][i]) * 86400.0 for k, v in
                                         (("inflow", INFLOW), ("outflow", OUTFLOW), ("precipitation", PRECIP),
                                          ("evaporation", EVAP), ("seepage", SEEP), ("applied", GW))}}
        self.start_volume = got[VOLUME][self.index].astype(float)  # tomorrow's start

    def before_solve(self):
        """Each lake set from its reservoir's day, after MODFLOW 6 prepares the day."""
        mf = self.mf6
        for lk, c in zip(self.lakes, self.curves):
            d, p, j = self.day[lk.reservoir], lk.package, lk.lake - 1
            for var in ("RAINFALL", "EVAPORATION", "RUNOFF", "INFLOW", "WITHDRAWAL"):  # SWAT+ owns the lake's water
                x = mf.get(var, mf.model, p.upper())
                x[j] = 0.0
                mf.set(var, p, x)
            status = mf.get_model_units("IBOUND", mf.model, p.upper())
            if lk.form == "fixed_stage":
                d["stage"] = c.stage(d["volume"])
                status[j] = -1
                mf.api.set_value(mf.address("IBOUND", mf.model, p.upper()), status)
                x = mf.get("STAGE", mf.model, p.upper())
                x[j] = d["stage"]
                mf.set("STAGE", p, x)
            else:
                d["stage"] = c.stage(d["start"])
                status[j] = 1
                mf.api.set_value(mf.address("IBOUND", mf.model, p.upper()), status)
                for var in ("XOLDPAK", "XNEWPAK"):
                    x = mf.get(var, mf.model, p.upper())
                    x[j] = d["stage"]
                    mf.set(var, p, x)
                for var, v in (("RUNOFF", d["inflow"] + d["precipitation"]),
                               ("WITHDRAWAL", d["outflow"] + d["evaporation"] + d["seepage"])):
                    x = mf.get(var, mf.model, p.upper())
                    x[j] = v
                    mf.set(var, p, x)

    def after_solve(self) -> dict[str, dict[str, float]]:
        """Each lake's leakage today (positive into the aquifer), due on its reservoir tomorrow; and its day."""
        mf = self.mf6
        for lk, s in zip(self.lakes, self.conn):
            d = self.day[lk.reservoir]
            d["to_aquifer"] = float(mf.simvals(lk.package)[s].sum())
            d["mf6_end_stage"] = float(mf.get("XNEWPAK", mf.model, lk.package.upper())[lk.lake - 1])
            d["due"] = self.due[lk.reservoir]
            self.due[lk.reservoir] = -d["to_aquifer"]
        return self.day
