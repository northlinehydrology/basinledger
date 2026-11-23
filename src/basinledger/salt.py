"""Salt: conservative species (chloride; weathering's ions; a mass, kg) through the land, the aquifer and the
streams of a coupled SWAT+ / MODFLOW 6 run, each day (docs/design.md, Salt). Built Nov 22, 2026 for the Virginia
study's chloride and conductance. Each species (``[salt.<species>]``) has its own GWT model, balance, ledger and
closure; specific conductance is a linear combination of them.

SWAT+ 62.0.1's salt module can't run with aquifers or point sources, so Basin Ledger carries the mass itself:

1. **Land** (``land_day``): a well-mixed store per HRU, its water the profile's (wilting point included) and SWAT+'s
   surface and lateral lags and wetland at the end of the day [A]. In: deposition (a concentration in precipitation
   [A]), road salt (a rate per m2 of road on snow and ice days, a series, or an annual amount per calendar year
   spread over the year's snow and ice days [A]), land sources (septic, fertilizer: kg m-2 a year, by HRU, land use
   or calendar year, timed evenly, by month or on the first fertilizer day [A]), and yesterday's groundwater into the
   soil at its cells' concentration. Out, at the store's concentration: runoff, lateral and tile flow to the streams,
   and percolation to the aquifer. Optionally a surface store (snowpack and road) in front of it (``surface_store``),
   whose melt and rain carry its mass to the stream fast [A]. Each input group's mass can be scaled by a multiplier
   (``Settings.multipliers``; calibration's ``salt`` parameters), applied where it enters.
2. **Groundwater**: the GWT species (coupling/transport.py), its recharge concentration each day each cell's
   share of the HRUs' percolation mass over their water, placed by the recharge map exactly as the water is;
   weathering as GWT's own mass source (``add_weathering``) [A].
3. **Streams** (``water_bodies_day``): each channel and reservoir a well-mixed store, upstream first. In: upstream
   outflow, the land's mass through the routing units, rain on the water, point sources (a user's returned mass, or
   a stated concentration [A]) and the aquifer's gains at their cells' concentration, a day after GWT gave them. Out,
   at the mixed concentration: the outflow, withdrawals, a reservoir's seepage, and the water the aquifer took the
   day before (which SWAT+ takes from the channel today). That water entered the aquifer with no chloride; its mass
   enters the same cells today, by a GWT source (``add_stream_sources``), so the link conserves exactly with the
   water's own one-day lag.
4. **Conductance** (``conductance``): a constant plus each species' slope times its mixed concentration (mg/L), per
   study [A], an output only.

Every source, store, link and outflow is booked in kg, and the land, the streams and the whole system (with GWT's
storage) close every day (``SaltBalance``; the ledger's ``salt:<species>:`` accounts).

SWAT+'s water is read through swatplus-bmi after each day: its printed daily columns (``hru_wb:...``,
``hru_pw:tmn``, ``channel_sd:precip``, read without printing, but only after ``print.prt``'s ``nyskip`` years, so
the run's copy starts at year 0) and the curated variables the other links use (percolation, channel and reservoir
volumes and flows, point sources, places of use).
"""

from __future__ import annotations

import csv
import datetime as dt
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

DAY = 86400.0
MG_L = 1000.0  # kg m-3 to mg L-1

# SWAT+'s day just done, per HRU, in mm unless named (swatplus-bmi's printed columns)
HRU_READS = {"precip": "hru_wb:precip", "snofall": "hru_wb:snofall", "snomlt": "hru_wb:snomlt",
             "surq": "hru_wb:surq_gen",        # surface runoff reaching the stream today (after SWAT+'s surface lag)
             "surq_cont": "hru_wb:surq_cont",  # surface runoff made today (before the lag)
             "latq": "hru_wb:latq", "qtile": "hru_wb:qtile", "sw": "hru_wb:sw_final", "lagsurf": "hru_wb:lagsurf",
             "laglatq": "hru_wb:laglatq", "wet": "hru_wb:wet_stor", "sno": "hru_wb:sno_final", "tmn": "hru_pw:tmn"}
PERCOLATION = "soil_water__percolation_rate"          # m s-1, the recharge link's own
SOIL_WATER = "soil_water__depth"                      # m above the wilting point, for the starting store
WILTING = "soil_layer_water__wilting_point_depth"     # m per layer
CH_VOLUME = "channel_water__volume"                   # m3 at the end of the day
CH_OUT = "channel_water__volume_outflow_rate"         # m3 s-1
CH_PRECIP = "channel_sd:precip"                       # m3, rain on the channel
RES_VOLUME = "reservoir_water__volume"
RES_OUT = "reservoir_water__volume_outflow_rate"
RES_SEEP = "reservoir_water__seepage_volume_flow_rate"
RES_PRECIP = "reservoir_water__precipitation_volume_flow_rate"
POINT = "point_source_water__volume_flow_rate"
WITHDRAWN = "place_of_use_water__volume_withdrawal_rate"
FERTILIZER = ("hru_nb:fertn", "hru_nb:fertp")         # kg/ha of N and P fertilizer applied today
COMPONENTS = ("sur", "lat", "til")  # the land's water to the streams, as SWAT+'s routing names its hydrographs


# ---------- settings ----------

@dataclass
class Annual:
    """An amount per calendar year (kg per m2 of the HRU), read at start-up in the run's HRU order
    (``read_annual``): a table ``hru, year, kg_m2``, or ``hru, kg_m2`` times a ``factor`` table ``year, factor``.
    A year after the table's last takes each HRU's mean of its last ``repeat_last`` years, and one before its first
    the mean of its first ``repeat_first`` [A]; without them such a year is refused."""
    path: Path
    factor: Path | None = None
    repeat_last: int | None = None
    repeat_first: int | None = None


@dataclass
class RoadSalt:
    """Road salt [A], one of three ways:

    - **a rate**: ``rate`` kg per m2 of road on each application day, on each HRU's road area (``share`` of its area,
      a number or {hru: share}); an application day has SWAT+'s snowfall of at least ``snowfall_mm``, or
      precipitation of at least ``precipitation_mm`` with the HRU's minimum temperature at or below ``tmin_c``;
    - **a series**: ``series`` {date: kg m-2 of road}, on those days only, on the road area;
    - **an annual amount** (``annual``, an ``Annual``): kg per m2 of the HRU each calendar year (as the USGS grid of
      Bock and others 2018 gives it), spread over that year's application days in proportion to their
      precipitation (``spread = "precipitation"``) or equally (``"equal"``). A winter's two halves come from two
      calendar years' amounts. The days come from the weather SWAT+ reads (the HRU's station's ``.pcp`` and ``.tmp``
      records, a day's snowfall being its precipitation at a mean temperature at or below the HRU's ``fall_tmp`` in
      ``snow.sno``), so a year's days are known before it starts and its amount is applied whole
      (``road_calendar``).

    ``group``: the multiplier group its mass is scaled by (``Settings.multipliers``), "road_salt"."""
    rate: float = 0.0
    share: float | dict = 0.0
    snowfall_mm: float = 1.0
    precipitation_mm: float = 1.0
    tmin_c: float | None = 0.0
    series: dict | None = None
    annual: Annual | None = None
    spread: str = "precipitation"
    group: str = "road_salt"

    @property
    def mode(self) -> str:
        return "series" if self.series is not None else "annual" if self.annual is not None else "rate"


@dataclass
class LandSource:
    """A chloride source on the land [A], into each HRU's soil store, kg per m2 of the HRU a year: ``rate`` one
    value, {hru: rate}, {land use: rate} (SWAT+'s ``lu_mgt`` in ``hru-data.hru``; ``by = "landuse"``) or an
    ``Annual`` (``hru, year, kg_m2`` by calendar year). ``timing``: "uniform" (the year's amount spread evenly over
    its days: septic effluent), 12 monthly weights, or "fertilizer" (the year's amount on the HRU's first fertilizer
    day of the calendar year, as SWAT+ reports it in ``hru_nb``'s ``fertn`` and ``fertp``; on 31 December if it had
    none: potash with the year's first fertilizer). ``group``: its multiplier group (its name unless given)."""
    name: str
    rate: float | dict | Annual
    by: str = "hru"
    timing: str | list = "uniform"
    group: str = ""

    def __post_init__(self):
        self.group = self.group or self.name


TIMINGS = ("uniform", "fertilizer")


@dataclass
class Settings:
    """``[salt.<species>]``. Concentrations in kg m-3 (= g L-1). ``multipliers``: {group: factor} on the inputs'
    mass, set per study or by calibration (a parameter of model ``salt``, target ``<species>:<group>``), applied
    where the mass enters, so the ledger books the scaled mass. Groups: ``precipitation``, ``road_salt``,
    ``point_sources`` (the stated concentrations; a user's reported mass is its own) and each land source's
    ``group``."""
    precipitation: float                   # in precipitation (bulk deposition) [A]
    road_salt: RoadSalt | None = None
    sources: list = field(default_factory=list)  # LandSource
    point_sources: dict = field(default_factory=dict)  # point source -> concentration of the water no user reports [A]
    soil_start: float = 0.0                # the land's starting concentration [A]
    stream_start: float = 0.0              # the channels' and reservoirs' [A]
    surface_store: bool = False            # a snowpack and road store in front of the soil's [A]
    multipliers: dict = field(default_factory=dict)

    @property
    def groups(self) -> set[str]:
        return {"precipitation", "road_salt", "point_sources"} | {x.group for x in self.sources}

    @classmethod
    def from_config(cls, sp: str, s: dict, base: Path) -> "Settings":
        known = {"precipitation", "road_salt", "sources", "point_sources", "soil_start", "stream_start",
                 "surface_store", "multipliers"}
        if set(s) - known:
            raise ValueError(f"[salt.{sp}] unknown settings {sorted(set(s) - known)}"
                             + (" (conductance is [salt] conductance, over the species)" if "conductance" in s else ""))
        if "precipitation" not in s:
            raise ValueError(f"[salt.{sp}] needs precipitation: the concentration in precipitation, kg m-3 [A]")
        road = _road_salt(sp, s["road_salt"], base) if s.get("road_salt") is not None else None
        sources = [_land_source(sp, x, base) for x in s.get("sources", [])]
        names = [x.name for x in sources]
        if len(set(names)) != len(names):
            raise ValueError(f"[salt.{sp}] sources: each needs a name of its own")
        surface = s.get("surface_store", False)
        if not isinstance(surface, bool):
            raise ValueError(f"[salt.{sp}] surface_store is true or false")
        out = cls(precipitation=float(s["precipitation"]), road_salt=road, sources=sources,
                  point_sources={k: float(v) for k, v in s.get("point_sources", {}).items()},
                  soil_start=float(s.get("soil_start", 0.0)), stream_start=float(s.get("stream_start", 0.0)),
                  surface_store=surface, multipliers={k: float(v) for k, v in s.get("multipliers", {}).items()})
        out.check_groups(out.multipliers, f"[salt.{sp}] multipliers")
        return out

    def check_groups(self, groups, what: str):
        unknown = sorted(set(groups) - self.groups)
        if unknown:
            raise ValueError(f"{what}: no input group {unknown} (groups: {sorted(self.groups)})")

    def multiplier(self, group: str) -> float:
        return self.multipliers.get(group, 1.0)


def _annual(what: str, x: dict, base: Path, path_key: str) -> Annual:
    rl, rf = x.get("repeat_last"), x.get("repeat_first")
    for k, v in (("repeat_last", rl), ("repeat_first", rf)):
        if v is not None and (not isinstance(v, int) or v < 1):
            raise ValueError(f"{what}: {k} is a number of years, 1 or more")
    return Annual(base / x[path_key], base / x["factor"] if x.get("factor") else None, rl, rf)


def _road_salt(sp: str, r: dict, base: Path) -> RoadSalt:
    r = dict(r)
    what = f"[salt.{sp}] road_salt"
    rk = {"rate", "share", "snowfall_mm", "precipitation_mm", "tmin_c", "series", "annual", "factor", "spread",
          "repeat_last", "repeat_first", "group"}
    if set(r) - rk:
        raise ValueError(f"{what}: unknown settings {sorted(set(r) - rk)}")
    ways = [k for k in ("rate", "series", "annual") if k in r]
    if len(ways) != 1:
        raise ValueError(f"{what} takes one of rate (kg m-2 of road a day, with share), series (a CSV of date, kg_m2, "
                         "with share) or annual (a CSV of hru, year, kg_m2 by calendar year, or of hru, kg_m2 with "
                         "factor)")
    if ways[0] != "annual" and "share" not in r:
        raise ValueError(f"{what} needs share (a number or a CSV of hru, share) with {ways[0]}")
    if ways[0] == "annual" and "share" in r:
        raise ValueError(f"{what}: an annual amount is per m2 of the HRU; share is for a rate or a series")
    if ways[0] != "annual" and set(r) & {"factor", "repeat_last", "repeat_first"}:
        raise ValueError(f"{what}: factor, repeat_last and repeat_first are for an annual amount")
    share = r.get("share", 0.0)
    if isinstance(share, str):
        share = _read_table(base / share, "hru", "share", f"{what} share")
    series = None
    if r.get("series"):
        series = {dt.date.fromisoformat(k): v for k, v in _read_table(base / r["series"], "date", "kg_m2",
                                                                      f"{what} series").items()}
    spread = r.get("spread", "precipitation")
    if spread not in ("precipitation", "equal"):
        raise ValueError(f"{what}: spread is \"precipitation\" or \"equal\"")
    return RoadSalt(rate=float(r.get("rate", 0.0)), share=share, snowfall_mm=float(r.get("snowfall_mm", 1.0)),
                    precipitation_mm=float(r.get("precipitation_mm", 1.0)),
                    tmin_c=None if r.get("tmin_c", 0.0) is None else float(r.get("tmin_c", 0.0)),
                    series=series, annual=_annual(what, r, base, "annual") if "annual" in r else None,
                    spread=spread, group=r.get("group", "road_salt"))


def _land_source(sp: str, x: dict, base: Path) -> LandSource:
    what = f"[salt.{sp}] sources"
    keys = {"name", "rate", "by", "timing", "group", "factor", "repeat_last", "repeat_first"}
    if not isinstance(x, dict) or set(x) - keys or not {"name", "rate"} <= set(x):
        raise ValueError(f"{what}: each is {{name, rate, by, timing, group, repeat_last, repeat_first}}, name and "
                         f"rate required; got {x!r}")
    what = f"{what} {x['name']}"
    by = x.get("by", "hru")
    if by not in ("hru", "landuse"):
        raise ValueError(f"{what}: by is \"hru\" or \"landuse\"")
    rate = x["rate"]
    if isinstance(rate, str):
        with open(base / rate, newline="") as f:
            head = [h.strip() for h in next(csv.reader(f))]
        if "year" in head or "factor" in x:  # hru, year, kg_m2 by calendar year (or hru, kg_m2 with a factor)
            if by != "hru":
                raise ValueError(f"{what}: a table by year is by HRU")
            rate = _annual(what, x, base, "rate")
        else:
            rate = _read_table(base / rate, by, "kg_m2_yr", what)
    elif by == "landuse":
        raise ValueError(f"{what}: by = \"landuse\" takes a CSV of landuse, kg_m2_yr")
    else:
        rate = float(rate)
    if not isinstance(rate, Annual) and set(x) & {"factor", "repeat_last", "repeat_first"}:
        raise ValueError(f"{what}: factor, repeat_last and repeat_first are for a table by year")
    timing = x.get("timing", "uniform")
    if isinstance(timing, list):
        if len(timing) != 12 or min(timing) < 0 or sum(timing) <= 0:
            raise ValueError(f"{what}: timing by month is 12 weights, none negative")
        timing = [float(v) for v in timing]
    elif timing not in TIMINGS:
        raise ValueError(f"{what}: timing is \"uniform\", \"fertilizer\" or 12 monthly weights")
    return LandSource(x["name"], rate, by, timing, x.get("group", ""))


def read_annual(a: Annual, hrus: list[str], years: range, what: str) -> dict[int, np.ndarray]:
    """kg per m2 of each HRU (in ``hrus``' order) for each calendar year in ``years``. HRUs the table doesn't name
    take zero; a name the model doesn't have is refused."""
    index = {h: i for i, h in enumerate(hrus)}
    with open(a.path, newline="") as f:
        rd = csv.reader(f)
        head = [h.strip() for h in next(rd)]
        rows = list(rd)
    if a.factor is None:
        if head[:3] != ["hru", "year", "kg_m2"] and set(head) != {"hru", "year", "kg_m2"}:
            raise ValueError(f"{what}: {a.path} needs columns hru, year, kg_m2 (or hru, kg_m2 with factor)")
        ih, iy, iv = (head.index(c) for c in ("hru", "year", "kg_m2"))
        ys = sorted({int(r[iy]) for r in rows if r})
        col = {y: j for j, y in enumerate(ys)}
        table = np.zeros((len(hrus), len(ys)))
        unknown = set()
        for r in rows:
            if not r:
                continue
            i = index.get(r[ih].strip())
            if i is None:
                unknown.add(r[ih].strip())
                continue
            table[i, col[int(r[iy])]] = float(r[iv])
    else:
        if "hru" not in head or "kg_m2" not in head:
            raise ValueError(f"{what}: {a.path} needs columns hru, kg_m2 with a factor")
        ih, iv = head.index("hru"), head.index("kg_m2")
        base = np.zeros(len(hrus))
        unknown = set()
        for r in rows:
            if r:
                i = index.get(r[ih].strip())
                if i is None:
                    unknown.add(r[ih].strip())
                else:
                    base[i] = float(r[iv])
        fac = {int(k): v for k, v in _read_table(a.factor, "year", "factor", f"{what} factor").items()}
        ys = sorted(fac)
        table = base[:, None] * np.array([fac[y] for y in ys])[None, :]
    if unknown:
        raise ValueError(f"{what}: {a.path} names HRUs SWAT+ doesn't have: {sorted(unknown)[:5]}")
    out = {}
    for y in years:
        if y in ys:
            out[y] = table[:, ys.index(y)]
        elif y > ys[-1] and a.repeat_last:
            out[y] = table[:, -a.repeat_last:].mean(axis=1)
        elif y < ys[0] and a.repeat_first:
            out[y] = table[:, :a.repeat_first].mean(axis=1)
        else:
            raise ValueError(f"{what}: {a.path} has no year {y} ({ys[0]} to {ys[-1]}); repeat_last or repeat_first "
                             "take the mean of the table's last or first years for years outside it")
    return out


def _read_table(path: Path, key: str, value: str, what: str) -> dict:
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows or key not in rows[0] or value not in rows[0]:
        raise ValueError(f"{what}: {path} needs columns {key}, {value}")
    return {r[key].strip(): float(r[value]) for r in rows}


# ---------- the land ----------

def land_day(m_soil, m_surface, area, w: dict, perc_mm, store_mm, deposition, road, from_groundwater,
             surface_store: bool = False, to_soil=0.0) -> dict:
    """One day of every HRU's store (arrays, one value per HRU). ``m_soil``, ``m_surface``: kg at the start of the
    day; ``area`` m2; ``w``: SWAT+'s day in mm (``HRU_READS``' keys); ``perc_mm`` percolation; ``store_mm`` the
    store's water at the end of the day; ``deposition``, ``road``, ``from_groundwater``: kg in today; ``to_soil``:
    kg in today straight into the soil store either way (the land sources: septic effluent, fertilizer) [A].

    Single store: the inputs join the soil store, and runoff, lateral and tile flow and percolation leave at
    (mass) / (end-of-day water + the day's outflows), so the mass closes by construction. With ``surface_store``,
    deposition and road salt land on the surface (the snowpack and the road) instead; the share of its water that
    leaves as liquid (rain + melt over that + the snowpack left) takes the same share of its mass, of which the share
    that ran off today (``surq_cont`` over rain + melt) goes to the stream as surface runoff's mass and the rest into
    the soil store; the soil then gives only lateral, tile and percolation.

    Returns the new masses and kg out by path (``sur``, ``lat``, ``til``, ``perc``), and the soil store's
    concentration (kg m-3)."""
    a = np.asarray(area, float) / 1000.0  # mm -> m3
    m_soil = np.asarray(m_soil, float).copy()
    m_surface = np.asarray(m_surface, float).copy()
    sur_kg = np.zeros_like(m_soil)
    if surface_store:
        m_surface += deposition + road
        liquid = np.clip(w["precip"] - w["snofall"], 0.0, None) + np.clip(w["snomlt"], 0.0, None)
        held = liquid + np.clip(w["sno"], 0.0, None)
        release = m_surface * np.divide(liquid, held, out=np.zeros_like(liquid), where=held > 0)
        m_surface -= release
        ran = np.minimum(1.0, np.divide(np.clip(w["surq_cont"], 0.0, None), liquid, out=np.zeros_like(liquid),
                                        where=liquid > 0))
        sur_kg = release * ran
        m_soil += release - sur_kg + from_groundwater + to_soil
        q = {"lat": w["latq"], "til": w["qtile"], "perc": perc_mm}
    else:
        m_soil += deposition + road + from_groundwater + to_soil
        q = {"sur": w["surq"], "lat": w["latq"], "til": w["qtile"], "perc": perc_mm}
    q = {k: np.clip(np.asarray(v, float), 0.0, None) * a for k, v in q.items()}
    water = np.clip(np.asarray(store_mm, float), 0.0, None) * a + sum(q.values())
    c = np.divide(m_soil, water, out=np.zeros_like(m_soil), where=water > 0)
    out = {k: c * v for k, v in q.items()}
    m_soil -= sum(out.values())
    if surface_store:
        out["sur"] = sur_kg
    return {"m_soil": m_soil, "m_surface": m_surface, "out": out, "concentration": c}


def cell_concentration(amap, perc_kg, perc_m, area) -> np.ndarray:
    """Each cell's recharge concentration (kg m-3, per node): the HRUs' percolation mass (kg) and water (m) placed by
    the recharge map (``amap``, each HRU's rows scaled to its ``area``, m2) exactly as the water is, so a cell under
    two HRUs takes their flow-weighted mean and the cells' mass adds to the HRUs'."""
    area = np.asarray(area, float)
    per_m2 = np.divide(perc_kg, area, out=np.zeros_like(area), where=area > 0)
    kg = np.bincount(amap.node, weights=per_m2[amap.obj] * amap.area)
    m3 = np.bincount(amap.node, weights=np.asarray(perc_m, float)[amap.obj] * amap.area, minlength=len(kg))
    return np.divide(kg, m3, out=np.zeros_like(kg), where=m3 > 0)


def application_days(w: dict, road: RoadSalt, day: dt.date) -> np.ndarray:
    """kg of chloride per m2 of road on each HRU today [A]: the series' rate on its days, or the rate on a snow or ice
    day (snowfall of ``snowfall_mm``, or precipitation of ``precipitation_mm`` at or below ``tmin_c``)."""
    if road.series is not None:
        return np.full(len(w["precip"]), float(road.series.get(day, 0.0)))
    snow = w["snofall"] >= road.snowfall_mm
    ice = (w["precip"] >= road.precipitation_mm) & (w["tmn"] <= road.tmin_c) if road.tmin_c is not None else False
    return np.where(snow | ice, road.rate, 0.0)


def source_share(timing, day: dt.date) -> float:
    """The share of a year's amount a land source puts on ``day``: uniform, or by 12 monthly weights (each month's
    share spread evenly over its days). The fertilizer timing is per HRU (``SaltBalance``)."""
    import calendar
    if timing == "uniform":
        return 1.0 / (366 if calendar.isleap(day.year) else 365)
    w = np.asarray(timing, float)
    return float(w[day.month - 1] / w.sum()) / calendar.monthrange(day.year, day.month)[1]


def read_station_weather(swat: Path) -> dict[str, dict[dt.date, tuple[float, float, float]]]:
    """Each station's daily (precipitation mm, tmax, tmin) from the files ``weather-sta.cli`` names."""
    swat = Path(swat)
    lines = (swat / "weather-sta.cli").read_text().splitlines()
    head = lines[1].split()
    out = {}
    cache = {}

    def series(fname, ncol):
        if fname not in cache:
            d = {}
            for ln in (swat / fname).read_text().splitlines()[3:]:
                f = ln.split()
                if len(f) >= 2 + ncol:
                    day = dt.date(int(f[0]), 1, 1) + dt.timedelta(days=int(f[1]) - 1)
                    d[day] = tuple(float(v) for v in f[2:2 + ncol])
            cache[fname] = d
        return cache[fname]

    for ln in lines[2:]:
        f = ln.split()
        if not f:
            continue
        row = dict(zip(head, f))
        if row.get("pcp", "sim") in ("sim", "null") or row.get("tmp", "sim") in ("sim", "null"):
            raise ValueError(f"salt: station {row['name']}: an annual road salt amount needs its precipitation and "
                             "temperature records (SWAT+ generates them)")
        pcp, tmp = series(row["pcp"], 1), series(row["tmp"], 2)
        out[row["name"]] = {d: (pcp[d][0], *tmp[d]) for d in pcp if d in tmp}
    return out


def road_calendar(road: RoadSalt, hru_station: list[str], hru_fall: np.ndarray, weather: dict, first: dt.date,
                  last: dt.date):
    """The annual road salt's application days, worked out before the run from the weather SWAT+ reads: for each
    group of HRUs sharing a station and a snowfall temperature, each day's weight (0 off an application day) and
    each calendar year's total weight. An application day: precipitation of at least ``snowfall_mm`` at a mean
    temperature at or below the HRU's ``fall_tmp`` (snow), or of at least ``precipitation_mm`` with the minimum at or
    below ``tmin_c`` (ice). Its weight: its precipitation, or 1 (``spread``). A year with no application day in the
    records takes its amount equally over its days with the minimum at or below ``tmin_c`` (0 C if unset) [A]; with
    none of those either it isn't applied, and is listed (``unapplied``). The years are ``first``'s to ``last``'s,
    whole, as far as the records go.

    Returns (group of each HRU, {date: weight per group}, {(group, year): total weight}, frost years used,
    unapplied [(group, year)])."""
    keys = sorted(set(zip(hru_station, np.round(hru_fall, 6).tolist())))
    gi = {k: i for i, k in enumerate(keys)}
    group = np.array([gi[(s, round(float(f), 6))] for s, f in zip(hru_station, hru_fall)], int)
    y0, y1 = first.year, last.year
    tmin_c = 0.0 if road.tmin_c is None else road.tmin_c
    weights: dict[dt.date, np.ndarray] = {}
    total: dict[tuple[int, int], float] = {}
    frost: dict[tuple[int, int], list[dt.date]] = {}
    for g, (station, fall) in enumerate(keys):
        for d, (p, tmax, tmin) in weather[station].items():
            if not y0 <= d.year <= y1:
                continue
            snow = p >= road.snowfall_mm and 0.5 * (tmax + tmin) <= fall
            ice = p >= road.precipitation_mm and tmin <= tmin_c
            key = (g, d.year)
            total.setdefault(key, 0.0)
            if tmin <= tmin_c:
                frost.setdefault(key, []).append(d)
            if snow or ice:
                wt = p if road.spread == "precipitation" else 1.0
                weights.setdefault(d, np.zeros(len(keys)))[g] += wt
                total[key] += wt
    used, unapplied = [], []
    for key, t in total.items():
        if t > 0:
            continue
        days = frost.get(key, [])
        if not days:
            unapplied.append(key)
            continue
        used.append(key)
        for d in days:
            weights.setdefault(d, np.zeros(len(keys)))[key[0]] += 1.0
        total[key] = float(len(days))
    return group, weights, total, used, unapplied


def conductance(base: float, slopes: dict[str, float], mg_l: dict[str, float]) -> float:
    """Specific conductance (uS/cm): ``base`` plus each species' slope (uS/cm per mg/L) times its concentration
    (mg/L) [A], fitted per study."""
    return float(base + sum(slopes[sp] * mg_l[sp] for sp in slopes))


def _hru_data(swat: Path, hrus: list[str]) -> dict:
    """Each HRU's weather station (``hru.con``'s wst), land use (``hru-data.hru``'s lu_mgt) and snowfall temperature
    (``snow.sno``'s fall_tmp for its snow parameters)."""
    swat = Path(swat)

    def table(path):
        lines = path.read_text().splitlines()
        head = lines[1].split()
        return [dict(zip(head, ln.split())) for ln in lines[2:] if ln.split()]
    con = {r["name"]: r for r in table(swat / "hru.con")}
    data = {r["id"]: r for r in table(swat / "hru-data.hru")}
    snow = {r["name"]: float(r["fall_tmp"]) for r in table(swat / "snow.sno")} if (swat / "snow.sno").exists() else {}
    out = {"station": [], "landuse": [], "fall_tmp": []}
    for h in hrus:
        c = con[h]
        d = data[c["hru"]]
        out["station"].append(c["wst"])
        out["landuse"].append(d["lu_mgt"])
        out["fall_tmp"].append(snow.get(d.get("snow", ""), 1.0))  # SWAT+'s default snowfall temperature, 1 C
    out["fall_tmp"] = np.array(out["fall_tmp"])
    return out


# ---------- the streams ----------

def water_bodies_day(order, down, m_store, inflow, v_end, v_out, v_withdrawn, v_seep, v_loss) -> dict:
    """One day of every channel's and reservoir's store (arrays by body), upstream first (``order``).

    ``down[b]``: [(body, share), ...] downstream of b (none: the basin's outlet). ``inflow``: kg in today besides the
    upstream bodies' (the land, rain, point sources, the aquifer's gains). Volumes, m3: at the end of the day, out
    downstream, withdrawn, seeping (a reservoir's own), and lost to the aquifer (the water MODFLOW 6 took the day
    before, which SWAT+ took from the channel today). Each leaves at (mass) / (end volume + outflows); a body with no
    water keeps its mass."""
    n = len(m_store)
    m = np.asarray(m_store, float) + np.asarray(inflow, float)
    out, wd, seep, loss, outlet, c = (np.zeros(n) for _ in range(6))
    for b in order:
        vol = v_end[b] + v_out[b] + v_withdrawn[b] + v_seep[b] + v_loss[b]
        if vol > 0.0:
            c[b] = m[b] / vol
            out[b], wd[b], seep[b], loss[b] = (c[b] * v for v in (v_out[b], v_withdrawn[b], v_seep[b], v_loss[b]))
            m[b] -= out[b] + wd[b] + seep[b] + loss[b]
        if down[b]:
            for d, share in down[b]:
                m[d] += out[b] * share
        else:
            outlet[b] = out[b]
    return {"m": m, "out": out, "withdrawn": wd, "seepage": seep, "loss": loss, "outlet": outlet, "concentration": c}


# ---------- SWAT+'s routing ----------

def _con(path: Path) -> list[dict]:
    """A SWAT+ ``.con`` file's rows: id, name, area (ha) and outputs [(type, id, hydrograph, share), ...]."""
    lines = Path(path).read_text().splitlines()
    head = [h.lower() for h in lines[1].split()]  # Basin Ledger's own recall.con heads its columns in capitals
    i_out = head.index("out_tot")
    i_area = next((head.index(a) for a in ("area", "area_ha") if a in head), None)
    rows = []
    for ln in lines[2:]:
        f = ln.split()
        if not f:
            continue
        n = int(f[i_out])
        outs = [(f[i_out + 1 + 4 * k], int(f[i_out + 2 + 4 * k]), f[i_out + 3 + 4 * k], float(f[i_out + 4 + 4 * k]))
                for k in range(n)]
        rows.append({"id": int(f[0]), "name": f[1], "area": float(f[i_area]) if i_area is not None else 0.0,
                     "out": outs})
    return rows


def _elements(tokens: list[str]) -> list[int]:
    """A routing unit's element numbers; ``a -b`` is the range a to b (SWAT+'s define_unit_elements)."""
    out, i = [], 0
    while i < len(tokens):
        a = int(tokens[i])
        if i + 1 < len(tokens) and int(tokens[i + 1]) < 0:
            out.extend(range(a, -int(tokens[i + 1]) + 1))
            i += 2
        else:
            out.append(a)
            i += 1
    return out


def _cio_files(swat: Path) -> dict[str, list[str]]:
    out = {}
    for ln in (swat / "file.cio").read_text().splitlines()[1:]:
        f = ln.split()
        if f:
            out[f[0]] = f[1:]
    return out


@dataclass
class Routing:
    """Where SWAT+ sends each HRU's water and each water body's outflow. ``bodies``: the channels, then the
    reservoirs; ``order``: upstream first; ``down[b]``: [(body, share)]; ``land[comp]``: (hru index, body index,
    share) rows for each of the land's paths to the streams (surface, lateral, tile)."""
    bodies: list[str]
    n_channels: int
    order: list[int]
    down: list[list[tuple[int, float]]]
    land: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]
    points: dict[str, int]  # recall object -> body


def routing(swat: Path, hrus: list[str], channels: list[str], reservoirs: list[str]) -> Routing:
    """Read the routing from a SWAT+ project (the run's copy). An HRU's water reaches the streams through its routing
    units (``rout_unit.def``, ``.ele``, ``.con``) or its own outputs in ``hru.con``; each HRU's mass is sent in full,
    split among its units by their share of its water (``ru_control.f90``'s expansion factor) and then as each unit
    splits its hydrographs (``tot`` all, ``sur``, ``lat``, ``til`` each its own; ``rhg``, recharge, is percolation and
    goes through the recharge map). Refused: water from the land or a body to anything but channels, reservoirs and
    routing units (runon to HRUs, an aquifer's surface or lateral water), and shares that don't add to 1."""
    swat = Path(swat)
    cio = _cio_files(swat)
    connect = cio.get("connect", [])
    bodies = list(channels) + list(reservoirs)
    index = {b: i for i, b in enumerate(bodies)}
    cha = {r["id"]: r["name"] for r in _con(swat / "chandeg.con")}
    res = {r["id"]: r["name"] for r in _con(swat / "reservoir.con")} if (swat / "reservoir.con").exists() else {}

    def body(typ, i, what):
        name = {"sdc": cha, "res": res}[typ].get(i)
        if name is None or name not in index:
            raise ValueError(f"salt: {what} sends water to {typ} {i}, which isn't among the model's channels and "
                             "reservoirs")
        return index[name]

    # routing units: elements, their HRUs' shares, and outputs
    rus = {}
    if (swat / "rout_unit.con").exists() and "rout_unit.con" in connect:
        con = {r["id"]: r for r in _con(swat / "rout_unit.con")}
        ele = {}
        for ln in (swat / "rout_unit.ele").read_text().splitlines()[2:]:
            f = ln.split()
            if f:
                ele[int(f[0])] = (f[2], int(f[3]), float(f[4]))
        for ln in (swat / "rout_unit.def").read_text().splitlines()[2:]:
            f = ln.split()
            if f:
                rus[int(f[0])] = {"elements": [ele[e] for e in _elements(f[3:])], **con[int(f[0])]}
    hru_con = _con(swat / "hru.con")
    hru_by_id = {r["id"]: r for r in hru_con}
    hidx = {h: i for i, h in enumerate(hrus)}

    def unit_shares(outs, comp, what, depth=0) -> dict[int, float]:
        """body -> share of a unit's ``comp`` water."""
        if depth > 50:
            raise ValueError(f"salt: {what}: routing units route in a loop")
        got: dict[int, float] = {}
        for typ, i, hyd, frac in outs:
            if hyd == "rhg":
                continue  # recharge: percolation, through the recharge map
            if hyd not in ("tot", comp):
                continue
            if typ in ("sdc", "res"):
                b = body(typ, i, what)
                got[b] = got.get(b, 0.0) + frac
            elif typ == "ru":
                for b, s in unit_shares(rus[i]["out"], comp, f"routing unit {rus[i]['name']}", depth + 1).items():
                    got[b] = got.get(b, 0.0) + frac * s
            else:
                raise ValueError(f"salt: {what} sends its {hyd} water to {typ} {i}; only channels, reservoirs and "
                                 "routing units are built")
        return got

    # each HRU's units, by their share of its water
    units: dict[int, list[tuple[list, float, str]]] = {i: [] for i in range(len(hrus))}
    for ru in rus.values():
        for typ, oid, frac in ru["elements"]:
            if typ != "hru":
                raise ValueError(f"salt: routing unit {ru['name']} has a {typ} element; HRUs only are built")
            name = hru_by_id[oid]["name"]
            if name not in hidx:
                continue
            h = hidx[name]
            ef = 1.0 if frac >= 0.99999 else frac * ru["area"] / max(hru_by_id[oid]["area"], 1e-12)
            units[h].append((ru["out"], ef, f"routing unit {ru['name']}"))
    for r in hru_con:
        if r["out"] and r["name"] in hidx:
            if units[hidx[r["name"]]]:
                raise ValueError(f"salt: {r['name']} routes by its own outputs in hru.con and through a routing unit")
            units[hidx[r["name"]]].append((r["out"], 1.0, r["name"]))
    land = {}
    for comp in COMPONENTS:
        hh, bb, ss = [], [], []
        for h, us in units.items():
            total = sum(ef for _, ef, _ in us)
            got: dict[int, float] = {}
            for outs, ef, what in us:
                for b, s in unit_shares(outs, comp, what).items():
                    got[b] = got.get(b, 0.0) + s * ef / total
            total_share = sum(got.values())
            if abs(total_share - 1.0) > 1e-4:
                raise ValueError(f"salt: {hrus[h]}'s {comp} water reaches the streams in shares adding to "
                                 f"{total_share:.6f}, not 1")
            for b, s in got.items():  # held to 1 exactly, so the land's mass reaches the streams whole
                hh.append(h); bb.append(b); ss.append(s / total_share)
        land[comp] = (np.array(hh, int), np.array(bb, int), np.array(ss))

    # the water bodies' outflow
    down: list[list[tuple[int, float]]] = [[] for _ in bodies]
    for f, kind in (("chandeg.con", "channel"), ("reservoir.con", "reservoir")):
        if not (swat / f).exists() or f not in connect:
            continue
        for r in _con(swat / f):
            if r["name"] not in index:
                continue
            b = index[r["name"]]
            for typ, i, hyd, frac in r["out"]:
                if typ not in ("sdc", "res"):
                    raise ValueError(f"salt: {kind} {r['name']} sends water to {typ} {i}; only channels and "
                                     "reservoirs are built downstream")
                down[b].append((body(typ, i, f"{kind} {r['name']}"), frac))
            if down[b] and abs(sum(s for _, s in down[b]) - 1.0) > 1e-4:
                raise ValueError(f"salt: {kind} {r['name']}'s outflow shares add to {sum(s for _, s in down[b])}")
    # upstream first (Kahn)
    n_up = np.zeros(len(bodies), int)
    for b in range(len(bodies)):
        for d, _ in down[b]:
            n_up[d] += 1
    ready = [b for b in range(len(bodies)) if n_up[b] == 0]
    order = []
    while ready:
        b = ready.pop()
        order.append(b)
        for d, _ in down[b]:
            n_up[d] -= 1
            if n_up[d] == 0:
                ready.append(d)
    if len(order) != len(bodies):
        raise ValueError("salt: the channels and reservoirs route in a loop")
    points = {}
    if (swat / "recall.con").exists() and "recall.con" in connect:
        for r in _con(swat / "recall.con"):
            outs = [(t, i) for t, i, _, _ in r["out"] if t in ("sdc", "res")]
            if len(outs) != 1 or len(r["out"]) != 1:
                raise ValueError(f"salt: point source {r['name']} must send its water to one channel or reservoir")
            points[r["name"]] = body(*outs[0], f"point source {r['name']}")
    return Routing(bodies, len(channels), order, down, land, points)


# ---------- weathering in the aquifer ----------

WEATHERING = "weathering"  # the GWT source weathering's mass enters by (add_weathering)


def add_weathering(gwt, rate, zones=None, pname: str = WEATHERING):
    """Weathering as a mass source in a GWT model (an SRC package, ``pname``): ``rate`` kg per m3 of aquifer (the
    cell's bulk volume) a year [A], one value, or {zone: rate} with ``zones`` a zone number per cell in the grid's
    shape (layers first). Zero-rate cells get no entry. Call from an example's builder, after ``add_transport``.

    The rates written are the cells' at full saturation. In a run with ``[salt]`` Basin Ledger scales each entry every
    day by its cell's saturated share at the start of the day (``start_saturation``, then MODFLOW 6's own), so a dry
    cell weathers nothing and a partly saturated one in proportion to its water: a convertible layer above the water
    table took the full rate into almost no water and the transport solve diverged (the Virginia study, Nov 22)."""
    import flopy
    dis = gwt.get_package("dis") or gwt.get_package("disv")
    top = np.asarray(dis.top.get_data(), float)
    botm = np.asarray(dis.botm.get_data(), float)
    thick = np.concatenate([top[None], botm[:-1]]) - botm
    if dis.package_type == "dis":
        area = np.outer(np.asarray(dis.delc.get_data(), float), np.asarray(dis.delr.get_data(), float))
    else:
        verts = {int(v[0]): (float(v[1]), float(v[2])) for v in dis.vertices.get_data()}
        area = np.zeros(int(dis.ncpl.get_data()))
        for c in dis.cell2d.get_data():
            ids = [int(i) for i in list(c)[4:4 + int(c[3])]]
            x, y = np.array([verts[i] for i in ids]).T
            area[int(c[0])] = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    volume = thick * area[None]
    if zones is None:
        r = np.full(volume.shape, float(rate))
    else:
        z = np.asarray(zones).reshape(volume.shape)
        r = np.vectorize(lambda k: float(rate.get(int(k), 0.0)))(z)
    idomain = dis.idomain.get_data()
    active = np.ones(volume.shape, bool) if idomain is None else np.asarray(idomain).reshape(volume.shape) > 0
    kg_d = r * volume / 365.25
    cells = np.argwhere((kg_d > 0) & active)
    spd = [(*map(int, c), float(kg_d[tuple(c)])) for c in cells]
    return flopy.mf6.ModflowGwtsrc(gwt, pname=pname, filename=f"{gwt.name}.{pname}.src",
                                   stress_period_data={0: [((*c[:-1],), c[-1]) for c in spd]})


def start_saturation(heads, top, bot, icelltype) -> np.ndarray:
    """Each cell's saturated share from its heads, where MODFLOW 6 hasn't solved a day yet (a run's first day, or a
    restart without it saved): 1 in a confined cell (ICELLTYPE 0), else (min(h, top) - bot) / (top - bot) clipped to
    0..1 (MODFLOW 6's own without Newton's smoothing [A]), 0 in a dry cell (MODFLOW 6's HDRY head, -1e30)."""
    h, top, bot = (np.asarray(v, float) for v in (heads, top, bot))
    thick = top - bot
    s = np.divide(np.minimum(h, top) - bot, thick, out=np.zeros_like(h), where=thick > 0)
    s = np.clip(s, 0.0, 1.0)
    return np.where(np.asarray(icelltype) == 0, 1.0, s)


STREAM_SOURCES = "stream_losses"  # the GWT source the streams' losses enter by


def add_stream_sources(gwt, gwf, packages: list[str], pname: str = STREAM_SOURCES):
    """The GWT mass source the streams' losses carry their chloride into the aquifer by (an SRC package, ``pname``):
    one entry on each entry's cell of each coupled river (and lakes') package of the flow model, in order, at zero.
    Basin Ledger sets each entry's rate every day. Call from an example's builder, after ``add_transport``."""
    import flopy
    cells = []
    for pname_ in packages:
        pkg = gwf.get_package(pname_)
        if pkg is None:
            raise KeyError(f"add_stream_sources: the flow model has no package {pname_!r}")
        rec = pkg.stress_period_data.get_data(0)
        cells += [tuple(r[0]) if not isinstance(r[0], (int, np.integer)) else (int(r[0]),) for r in rec]
    return flopy.mf6.ModflowGwtsrc(gwt, pname=pname, filename=f"{gwt.name}.{pname}.src", maxbound=len(cells),
                                   stress_period_data={0: [(c, 0.0) for c in cells]})


# ---------- the balance in a run ----------

class SaltBalance:
    """One species' salt in a running Basin (driver.py), each day:

    - ``after_surface(L, day)``, once SWAT+ has run the day: the land, then the streams (with yesterday's gains in
      and yesterday's losses out); the recharge's concentration per cell, and the losses' mass for the aquifer;
    - ``recharge_concentration(entry_nodes)``, ``stream_sources(n)``: for the transport model, after MODFLOW 6
      prepares the day;
    - ``after_groundwater(L, budget, ...)``, after the solve: the gains, the lost water and the soil's water due
      tomorrow, and the system's closure."""

    def __init__(self, species: str, settings: Settings, sw, swat_dir: Path, recharge_map, hru_area, places,
                 direct_returns: dict[str, list[tuple[str, set]]], gw_soil=None, saved: dict | None = None,
                 start: dt.date | None = None, days: int = 0):
        from .coupling import swatplus as S
        self.sp, self.s, self.sw, self.gw_soil = species, settings, sw, gw_soil
        self.k = f"salt:{species}:"
        self.hrus = sw.names(PERCOLATION)
        self.area = np.asarray(hru_area, float)
        missing = [v for v in [*HRU_READS.values(), CH_PRECIP] if not sw.has(v)]
        if missing:
            raise ValueError(f"salt needs swatplus-bmi's printed daily columns ({', '.join(missing[:3])}); "
                             "swatplus-bmi 0.9.0 or later")
        self.map = recharge_map  # each HRU's rows scaled to its area: the recharge water's map
        if list(self.map.objects) != list(self.hrus):
            raise ValueError("salt: the recharge map's objects aren't SWAT+'s HRUs ([recharge] source = \"hru\")")
        self.channels = sw.names(S.QOUT)
        self.reservoirs = sw.names(RES_VOLUME) if sw.has(RES_VOLUME) else []
        self.route = routing(swat_dir, self.hrus, self.channels, self.reservoirs)
        nb = len(self.route.bodies)
        self.ch_index = {c: i for i, c in enumerate(self.channels)}
        # wilting point water per HRU, m (constant through a run)
        at = {h: i for i, h in enumerate(self.hrus)}
        lay = np.array([at[n.split(":")[0]] for n in sw.names(WILTING)])
        self.wp_mm = np.bincount(lay, weights=sw.get(WILTING), minlength=len(self.hrus)) * 1000.0
        hru_data = _hru_data(swat_dir, self.hrus)  # each HRU's station, land use and snowfall temperature
        # road area per HRU, m2 (a rate or a series), or the years' amounts and their days (an annual amount)
        self.road_m2 = np.zeros(len(self.hrus))
        self.road_annual = None
        rs = settings.road_salt
        if rs is not None and rs.mode != "annual":
            sh = rs.share
            if isinstance(sh, dict):
                self._known(sh, "road_salt share")
                self.road_m2 = np.array([sh.get(h, 0.0) for h in self.hrus]) * self.area
            else:
                self.road_m2 = float(sh) * self.area
        elif rs is not None:
            self._road_annual(rs, swat_dir, hru_data, start, days)
        # land sources: kg a year per HRU, scaled by their groups' multipliers (a table by year: {year: kg})
        self.sources = []
        last = start + dt.timedelta(days=max(days, 1) - 1) if start is not None else None
        for src in settings.sources:
            m = settings.multiplier(src.group)
            if isinstance(src.rate, Annual):
                if start is None:
                    raise ValueError(f"salt: source {src.name}'s table by year needs the run's start and length")
                by_year = read_annual(src.rate, self.hrus, range(start.year, last.year + 1), f"source {src.name}")
                self.sources.append((src, {y: v * self.area * m for y, v in by_year.items()}))
                continue
            if isinstance(src.rate, dict):
                keys = self.hrus if src.by == "hru" else hru_data["landuse"]
                if src.by == "hru":
                    self._known(src.rate, f"source {src.name}")
                else:
                    unknown = sorted(set(src.rate) - set(hru_data["landuse"]))
                    if unknown:
                        raise ValueError(f"salt: source {src.name} names land uses no HRU has: {unknown[:5]}")
                rate = np.array([src.rate.get(x, 0.0) for x in keys])
            else:
                rate = np.full(len(self.hrus), float(src.rate))
            self.sources.append((src, rate * self.area * m))
        if any(src.timing == "fertilizer" for src in settings.sources):
            for v in FERTILIZER:
                if not sw.has(v):
                    raise ValueError(f"salt: a source timed by fertilizer reads SWAT+'s {v} (swatplus-bmi 0.9.0)")
        # point sources and withdrawals
        self.points = sw.names(POINT) if sw.has(POINT) else []
        unknown = sorted(set(settings.point_sources) - set(self.points))
        if unknown:
            raise ValueError(f"salt: point_sources names {unknown}, which SWAT+ doesn't have")
        missing_route = [p for p in self.points if p not in self.route.points]
        if missing_route:
            raise ValueError(f"salt: point sources {missing_route} aren't in recall.con")
        self.direct_returns = direct_returns  # point source -> (user, its constituents) returning to it directly
        self.places = sw.names(WITHDRAWN) if sw.has(WITHDRAWN) else []
        to_ch = dict(places)
        for p in self.places:
            if p not in to_ch:
                raise ValueError(f"salt: place of use {p} isn't Basin Ledger's; SWAT+'s own water allocation isn't "
                                 "built")
        self.place_body = np.array([self.ch_index[to_ch[p]] for p in self.places], int)
        # state
        if saved is not None:
            self.m_soil = np.asarray(saved["m_soil"], float)
            self.m_surface = np.asarray(saved["m_surface"], float)
            self.m_body = np.asarray(saved["m_body"], float)
            self.gain_due = np.asarray(saved["gain_due"], float)
            self.loss_q = np.asarray(saved["loss_q"], float)
            self.loss_body = np.asarray(saved["loss_body"], int)
            self.soil_due = np.asarray(saved["soil_due"], float)
            # each cell's saturated share at the end of the day saved, for weathering (saved since Nov 22, 2026)
            self.saturation = np.asarray(saved["saturation"], float) if "saturation" in saved else None
            self.total = float(saved["total"][0])
            self.fertilized = {src.name: (int(saved[f"fertilized_year_{src.name}"][0]),
                                          np.asarray(saved[f"fertilized_{src.name}"], bool))
                               for src, _ in self.sources if src.timing == "fertilizer"}
        else:
            store0 = (sw.get(SOIL_WATER) * 1000.0 + self.wp_mm) * self.area / 1000.0
            self.m_soil = settings.soil_start * store0
            self.m_surface = np.zeros(len(self.hrus))
            v0 = np.concatenate([sw.get(CH_VOLUME), sw.get(RES_VOLUME) if self.reservoirs else np.zeros(0)])
            self.m_body = settings.stream_start * np.clip(v0, 0.0, None)
            self.gain_due = np.zeros(nb)
            self.loss_q = np.zeros(0)   # m3 each coupled river and lake entry lost to the aquifer yesterday
            self.loss_body = np.zeros(0, int)  # and its channel (-1: not coupled)
            self.soil_due = np.zeros(len(self.hrus))
            self.total = None  # the system's mass outside the aquifer, once the first day sets it
            # a source timed by fertilizer: the year, and the HRUs that have had it that year
            self.fertilized = {src.name: (0, np.zeros(len(self.hrus), bool))
                               for src, _ in self.sources if src.timing == "fertilizer"}
            self.saturation = None  # each cell's saturated share at the end of the day just solved (weathering)
        self.weathering_entries = None  # (rate at full saturation, kg d-1; node) per weathering entry, read on day 1
        self.c_body = np.zeros(nb)
        self.c_cell = None
        self.perc_kg = 0.0
        self.src_rates = np.zeros(len(self.loss_q))

    def state(self) -> dict:
        out = {"m_soil": self.m_soil, "m_surface": self.m_surface, "m_body": self.m_body, "gain_due": self.gain_due,
               "loss_q": self.loss_q, "loss_body": self.loss_body, "soil_due": self.soil_due,
               "total": np.array([self._outside()])}
        if self.saturation is not None:
            out["saturation"] = self.saturation
        for name, (year, done) in self.fertilized.items():
            out[f"fertilized_year_{name}"] = np.array([year])
            out[f"fertilized_{name}"] = done
        return out

    def _known(self, table: dict, what: str):
        unknown = sorted(set(table) - set(self.hrus))
        if unknown:
            raise ValueError(f"salt: {what} names HRUs SWAT+ doesn't have: {unknown[:5]}")

    def _road_annual(self, rs: RoadSalt, swat_dir: Path, hru_data: dict, start: dt.date | None, days: int):
        """The annual amounts by calendar year (kg per HRU, scaled) and their days (``road_calendar``)."""
        if start is None or days <= 0:
            raise ValueError("salt: an annual road salt amount needs the run's start and length")
        last = start + dt.timedelta(days=days - 1)
        years = range(start.year, last.year + 1)
        m = self.s.multiplier(rs.group)
        amount = {y: v * self.area * m for y, v in read_annual(rs.annual, self.hrus, years, "road_salt").items()}
        weather = read_station_weather(swat_dir)
        group, weights, total, frost, unapplied = road_calendar(rs, hru_data["station"], hru_data["fall_tmp"],
                                                                weather, start, last)
        unapplied = [(g, y) for g, y in unapplied if amount[y][group == g].sum() > 0]
        import warnings
        if unapplied:
            warnings.warn(f"salt: road salt for {len(unapplied)} station-years has no application day or frost day "
                          f"in the weather records and isn't applied (years {sorted({y for _, y in unapplied})})",
                          stacklevel=3)
        if frost:
            warnings.warn(f"salt: {len(frost)} station-years have no snow or ice day; their road salt is spread over "
                          "their frost days", stacklevel=3)
        self.road_annual = {"amount": amount, "group": group, "weights": weights, "total": total}

    def road_today(self, w: dict, day: dt.date) -> np.ndarray:
        """kg of road salt on each HRU today."""
        rs = self.s.road_salt
        if rs is None:
            return np.zeros(len(self.hrus))
        if rs.mode != "annual":
            return application_days(w, rs, day) * self.road_m2 * self.s.multiplier(rs.group)
        a = self.road_annual
        wt = a["weights"].get(day)
        y = day.year
        if wt is None or y not in a["amount"]:
            return np.zeros(len(self.hrus))
        g = a["group"]
        tot = np.array([a["total"].get((i, y), 0.0) for i in range(len(wt))])
        share = np.divide(wt, tot, out=np.zeros_like(wt), where=tot > 0)
        return a["amount"][y] * share[g]

    def sources_today(self, day: dt.date, reads) -> dict[str, np.ndarray]:
        """kg from each land source on each HRU today."""
        out = {}
        for src, kg in self.sources:
            kg_yr = kg[day.year] if isinstance(kg, dict) else kg
            if src.timing != "fertilizer":
                out[src.name] = kg_yr * source_share(src.timing, day)
                continue
            year, done = self.fertilized[src.name]
            if year != day.year:
                year, done = day.year, np.zeros(len(self.hrus), bool)
            fert = (reads(FERTILIZER[0]) + reads(FERTILIZER[1])) > 0
            now = ~done & (fert | (day.month == 12 and day.day == 31))
            out[src.name] = np.where(now, kg_yr, 0.0)
            self.fertilized[src.name] = (year, done | now)
        return out

    def _outside(self) -> float:
        """The mass outside the aquifer, kg: land, water bodies, and the links in transit (the gains to the streams
        and the soil's water, due tomorrow)."""
        return float(self.m_soil.sum() + self.m_surface.sum() + self.m_body.sum() + self.gain_due.sum()
                     + self.soil_due.sum())

    # ---------- after SWAT+'s day ----------
    def after_surface(self, L, day: dt.date, cache: dict | None = None):
        """``cache``: the day's reads, shared by the species (a dict the driver empties each day)."""
        from .coupling import swatplus as S
        k = self.k
        cache = {} if cache is None else cache

        def reads(v):
            if v not in cache:
                cache[v] = self.sw.get(v)
            return cache[v]
        if self.total is None:
            self.total = self._outside()
        start_total = self.total
        w = {n: reads(v) for n, v in HRU_READS.items()}
        if not np.isfinite(w["precip"]).all():
            raise ValueError("salt: SWAT+ has no daily record for the day (print.prt's nyskip years)")
        perc_m = np.clip(reads(PERCOLATION), 0.0, None) * DAY
        store = w["sw"] + self.wp_mm + w["lagsurf"] + w["laglatq"] + w["wet"]
        c_rain = self.s.precipitation * self.s.multiplier("precipitation")
        dep = c_rain * np.clip(w["precip"], 0.0, None) / 1000.0 * self.area
        road = self.road_today(w, day)
        sources = self.sources_today(day, reads)
        to_soil = sum(sources.values()) if sources else 0.0
        land0 = self.m_soil.sum() + self.m_surface.sum()
        r = land_day(self.m_soil, self.m_surface, self.area, w, perc_m * 1000.0, store, dep, road, self.soil_due,
                     self.s.surface_store, to_soil)
        self.m_soil, self.m_surface = r["m_soil"], r["m_surface"]
        out = r["out"]
        received = float(self.soil_due.sum())
        self.soil_due = np.zeros(len(self.hrus))
        to_streams = sum(float(out[c].sum()) for c in ("sur", "lat", "til") if c in out)
        land_sources = float(sum(float(v.sum()) for v in sources.values()))
        land_in = float(dep.sum() + road.sum()) + land_sources + received
        land_out = to_streams + float(out["perc"].sum())
        land1 = float(self.m_soil.sum() + self.m_surface.sum())
        L.enter(k + "land:precipitation", float(dep.sum()))
        L.enter(k + "land:road_salt", float(road.sum()))
        for name, v in sources.items():
            L.enter(k + f"land:source:{name}", float(v.sum()))
        L.enter(k + "land:from_groundwater", received)
        L.enter(k + "land:to_streams", to_streams)
        L.enter(k + "land:recharge", float(out["perc"].sum()))
        L.enter(k + "land:storage", land1)
        L.enter(k + "land:closure_share", (land1 - land0 - land_in + land_out) / max(land_in + land_out + land0, 1e-6))
        # the recharge's concentration per cell: the HRUs' percolation mass and water placed by the recharge map
        self.perc_kg = float(out["perc"].sum())
        self.c_cell = cell_concentration(self.map, out["perc"], perc_m, self.area)

        # the streams
        nb, nc = len(self.route.bodies), self.route.n_channels
        inflow = np.zeros(nb)
        from_land = 0.0
        for comp in ("sur", "lat", "til"):
            if comp not in out:
                continue
            h, b, s = self.route.land[comp]
            add = np.bincount(b, weights=out[comp][h] * s, minlength=nb)
            inflow += add
            from_land += float(add.sum())
        rain = np.zeros(nb)
        rain[:nc] = c_rain * np.clip(reads(CH_PRECIP), 0.0, None)
        v_end, v_out, v_seep = np.zeros(nb), np.zeros(nb), np.zeros(nb)
        v_end[:nc] = reads(CH_VOLUME)
        v_out[:nc] = reads(S.QOUT) * DAY
        if self.reservoirs:
            rain[nc:] = c_rain * np.clip(reads(RES_PRECIP), 0.0, None) * DAY
            v_end[nc:] = reads(RES_VOLUME)
            v_out[nc:] = reads(RES_OUT) * DAY
            v_seep[nc:] = reads(RES_SEEP) * DAY
        inflow += rain
        points = np.zeros(nb)
        if self.points:
            vol = reads(POINT) * DAY
            for p, v in zip(self.points, vol):
                if v <= 0:
                    continue
                users = [u for u, reported in self.direct_returns.get(p, []) if self.sp in reported]
                kg_u = sum(L.today(f"user:{u}:{self.sp}:kg") for u in users)
                v_u = sum(L.today(f"user:{u}:return") for u in users)
                rest = v - v_u
                rest = rest if rest > 1e-6 * max(v, 1.0) else 0.0  # SWAT+'s 4-byte rounding of the users' water
                if rest > 0.0 and p not in self.s.point_sources:
                    raise ValueError(f"salt: point source {p} carries {rest:.1f} m3 of water with no {self.sp} "
                                     f"stated for it ([salt.{self.sp}] point_sources) and no user reporting it")
                points[self.route.points[p]] += kg_u + (self.s.point_sources.get(p, 0.0) * rest
                                                        * self.s.multiplier("point_sources"))
        inflow += points
        inflow += self.gain_due
        v_wd = np.zeros(nb)
        if self.places:
            np.add.at(v_wd, self.place_body, np.clip(reads(WITHDRAWN), 0.0, None) * DAY)
        lost = self.loss_body >= 0
        v_loss = np.bincount(self.loss_body[lost], weights=self.loss_q[lost], minlength=nb)
        streams0 = float(self.m_body.sum())
        gains = float(self.gain_due.sum())
        q = water_bodies_day(self.route.order, self.route.down, self.m_body, inflow, v_end, np.clip(v_out, 0.0, None),
                             v_wd, np.clip(v_seep, 0.0, None), v_loss)
        self.m_body = q["m"]
        self.c_body = q["concentration"]
        # the losses' mass, each entry's water at its channel's mixed concentration: the aquifer's source today
        self.src_rates = np.where(lost, self.c_body[np.maximum(self.loss_body, 0)] * self.loss_q, 0.0)
        losses = float(q["loss"].sum())
        s_in = from_land + float(rain.sum() + points.sum()) + gains
        s_out = losses + float(q["outlet"].sum() + q["withdrawn"].sum() + q["seepage"].sum())
        streams1 = float(self.m_body.sum())
        L.enter(k + "streams:from_land", from_land)
        L.enter(k + "streams:precipitation", float(rain.sum()))
        L.enter(k + "streams:point_sources", float(points.sum()))
        L.enter(k + "streams:gains", gains)
        L.enter(k + "streams:losses", losses)
        L.enter(k + "streams:withdrawn", float(q["withdrawn"].sum()))
        L.enter(k + "streams:outlet", float(q["outlet"].sum()))
        L.enter(k + "streams:reservoir_seepage", float(q["seepage"].sum()))
        L.enter(k + "streams:storage", streams1)
        L.enter(k + "streams:closure_share", (streams1 - streams0 - s_in + s_out) / max(s_in + s_out + streams0, 1e-6))
        self.gain_due = np.zeros(nb)
        # the day's external flows outside the aquifer, for the system's closure after the solve
        ext_in = float(dep.sum() + road.sum() + rain.sum() + points.sum()) + land_sources
        self._day = {"start": start_total, "in": ext_in,
                     "out": float(q["outlet"].sum() + q["withdrawn"].sum() + q["seepage"].sum()),
                     "gross": ext_in + land_out + s_out}

    def concentration(self, channel: str) -> float:
        """A channel's outflow concentration today, kg m-3."""
        return float(self.c_body[self.ch_index[channel]])

    # ---------- MODFLOW 6's sources ----------
    def recharge_concentration(self, entry_nodes: np.ndarray) -> np.ndarray:
        """Per recharge entry (its 0-based node, -1 for none), today's concentration, kg m-3."""
        n = np.asarray(entry_nodes, int)
        c = np.zeros(len(n))
        ok = (n >= 0) & (n < len(self.c_cell))
        c[ok] = self.c_cell[n[ok]]
        return c

    def stream_sources(self, n: int) -> np.ndarray:
        """Today's rate (kg d-1) for each of the aquifer's ``n`` stream-loss source entries: yesterday's water lost
        there at its channel's mixed concentration today; zero on the first day."""
        if len(self.src_rates) == 0:
            return np.zeros(n)
        if len(self.src_rates) != n:
            raise ValueError(f"salt: {len(self.src_rates)} stream-loss rates for {n} source entries")
        return self.src_rates

    def weathering_rates(self, base: np.ndarray, nodes: np.ndarray, first) -> np.ndarray:
        """Today's rate (kg d-1) for each weathering entry (``add_weathering``): its rate at full saturation, the
        builder's, times its cell's saturated share at the start of the day [A]: yesterday's at the end of its solve
        (MODFLOW 6's NPF saturation, kept by ``after_groundwater``) or, on a run's first day, ``first()`` (from the
        starting heads, ``start_saturation``). A dry cell weathers nothing. The day's own saturation isn't known
        until the flow is solved, which is after the transport's sources are set: a one-day lag, as the streams'
        losses have."""
        if self.saturation is None:
            self.saturation = np.asarray(first(), float)
        return np.asarray(base, float) * self.saturation[np.asarray(nodes, int)]

    # ---------- after MODFLOW 6's day ----------
    def after_groundwater(self, L, b, mf6=None, soil_package: str | None = None, x: np.ndarray | None = None,
                          losses: list[tuple[np.ndarray, np.ndarray, list[str]]] = ()):
        """``b``: the transport model's MassBudget for the day. ``losses``: for each coupled river and lakes package,
        in the stream sources' order, each entry's flow (m3, positive into the aquifer), its channel's index in
        ``channels`` (-1: not coupled) and ``channels``. The gains due to the streams tomorrow, the water lost to the
        aquifer (its mass the aquifer's source tomorrow), the soil's water due to the land, and the system's
        closure."""
        k = self.k
        if mf6 is not None:  # each cell's saturated share at the end of the day, for tomorrow's weathering
            self.saturation = mf6.get("SAT", mf6.model, "NPF").astype(float)
        for c, m in b.to_stream.items():
            self.gain_due[self.ch_index[c]] += m
        qq, bb = [], []
        for q, entry_channel, channels in losses:
            body = np.array([self.ch_index[c] for c in channels] + [-1])
            e = np.asarray(entry_channel, int)[: len(q)]
            on = (e >= 0) & (q > 0)
            qq.append(np.where(on, q, 0.0))
            bb.append(np.where(on, body[np.where(e >= 0, e, len(channels))], -1))
        self.loss_q = np.concatenate(qq) if qq else np.zeros(0)
        self.loss_body = np.concatenate(bb).astype(int) if bb else np.zeros(0, int)
        src = b.boundary.get(STREAM_SOURCES, 0.0)
        L.enter(k + "aquifer:from_streams", src)
        soil_sent = 0.0
        if self.gw_soil is not None:
            gs = self.gw_soil
            out = np.clip(-mf6.simvals(gs.package), 0.0, None)
            cell = np.zeros(mf6.nodes)
            cell[gs.nodes] = out
            if gs.rule == "drain":
                row = cell[gs.map.node] * gs.shares
            else:
                planned = np.zeros(mf6.nodes)
                planned[gs.nodes] = gs.planned
                frac = np.divide(cell, planned, out=np.zeros_like(cell), where=planned > 0)
                row = gs.row_planned * frac[gs.map.node]
            self.soil_due = np.bincount(gs.map.obj, weights=row * x[gs.map.node], minlength=len(self.hrus))
            soil_sent = -b.boundary.get(soil_package.lower(), 0.0)
            L.enter(k + "soil:sent", soil_sent)
            L.enter(k + "soil:mapped", float(self.soil_due.sum()))
        L.enter(k + "transit:gains", float(self.gain_due.sum()))
        L.enter(k + "transit:soil", float(self.soil_due.sum()))
        L.enter(k + "transit:losses_m3", float(self.loss_q.sum()))
        # the system: land + water bodies + in transit + the aquifer change by what enters and leaves from outside
        internal = {STREAM_SOURCES} | ({soil_package.lower()} if soil_package else set())
        boundary = sum(m for p, m in b.boundary.items() if p not in internal)
        wells = sum(b.wells.values())
        total = self._outside()
        d = self._day
        change = (total - d["start"]) - b.storage_release
        ext = d["in"] - d["out"] + boundary - wells
        gross = d["gross"] + b.gross + abs(total - d["start"])
        L.enter(k + "closure_share", (change - ext) / max(gross, 1e-6))
        self.total = total
