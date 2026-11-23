"""Salt: chloride (a conservative mass, kg) through the land, the aquifer and the streams of a coupled SWAT+ / MODFLOW 6
run, each day (docs/design.md, Salt). Built Nov 22, 2026 for the Virginia study's chloride and conductance.

SWAT+ 62.0.1's salt module can't run with aquifers or point sources, so Basin Ledger carries the mass itself:

1. **Land** (``land_day``): a well-mixed store per HRU, its water the profile's (wilting point included) and SWAT+'s
   surface and lateral lags and wetland at the end of the day [A]. In: deposition (a concentration in precipitation
   [A]), road salt (a rate per m2 of road on snow and ice days, or a stated series [A]), and yesterday's groundwater
   into the soil at its cells' concentration. Out, at the store's concentration: runoff, lateral and tile flow to the
   streams, and percolation to the aquifer. Optionally a surface store (snowpack and road) in front of it
   (``surface_store``), whose melt and rain carry its mass to the stream fast [A].
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
4. **Conductance**: background + slope x chloride (mg/L), per study [A], an output only.

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
COMPONENTS = ("sur", "lat", "til")  # the land's water to the streams, as SWAT+'s routing names its hydrographs


# ---------- settings ----------

@dataclass
class RoadSalt:
    """Road salt [A]: ``rate`` kg of chloride per m2 of road on an application day; ``share`` each HRU's road area as a
    share of its area (a number, or {hru: share}); an application day has snowfall of at least ``snowfall_mm`` or
    precipitation of at least ``precipitation_mm`` with the HRU's minimum temperature at or below ``tmin_c`` (a snow
    or ice event). ``series``, if given, replaces the rule: {date: kg m-2 of road}, applied on those days only."""
    rate: float = 0.0
    share: float | dict = 0.0
    snowfall_mm: float = 1.0
    precipitation_mm: float = 1.0
    tmin_c: float | None = 0.0
    series: dict | None = None


@dataclass
class Settings:
    """``[salt.<species>]``. Concentrations in kg m-3 (= g L-1)."""
    precipitation: float                   # in precipitation (bulk deposition) [A]
    road_salt: RoadSalt | None = None
    point_sources: dict = field(default_factory=dict)  # point source -> concentration of the water no user reports [A]
    soil_start: float = 0.0                # the land's starting concentration [A]
    stream_start: float = 0.0              # the channels' and reservoirs' [A]
    surface_store: bool = False            # a snowpack and road store in front of the soil's [A]
    conductance: dict | None = None        # {background (uS/cm), slope (uS/cm per mg/L)} [A]

    @classmethod
    def from_config(cls, sp: str, s: dict, base: Path) -> "Settings":
        known = {"precipitation", "road_salt", "point_sources", "soil_start", "stream_start", "surface_store",
                 "conductance"}
        if set(s) - known:
            raise ValueError(f"[salt.{sp}] unknown settings {sorted(set(s) - known)}")
        if "precipitation" not in s:
            raise ValueError(f"[salt.{sp}] needs precipitation: the concentration in precipitation, kg m-3 [A]")
        road = None
        if s.get("road_salt") is not None:
            r = dict(s["road_salt"])
            rk = {"rate", "share", "snowfall_mm", "precipitation_mm", "tmin_c", "series"}
            if set(r) - rk:
                raise ValueError(f"[salt.{sp}] road_salt: unknown settings {sorted(set(r) - rk)}")
            if "share" not in r or ("rate" not in r and "series" not in r):
                raise ValueError(f"[salt.{sp}] road_salt needs share (a number or a CSV of hru, share) and rate "
                                 "(kg m-2 of road a day) or series (a CSV of date, kg m-2)")
            share = r["share"]
            if isinstance(share, str):
                share = _read_table(base / share, "hru", "share", f"[salt.{sp}] road_salt share")
            series = None
            if r.get("series"):
                series = {dt.date.fromisoformat(k): v for k, v in
                          _read_table(base / r["series"], "date", "kg_m2", f"[salt.{sp}] road_salt series").items()}
            road = RoadSalt(rate=float(r.get("rate", 0.0)), share=share,
                            snowfall_mm=float(r.get("snowfall_mm", 1.0)),
                            precipitation_mm=float(r.get("precipitation_mm", 1.0)),
                            tmin_c=None if r.get("tmin_c", 0.0) is None else float(r.get("tmin_c", 0.0)),
                            series=series)
        cond = s.get("conductance")
        if cond is not None and set(cond) != {"background", "slope"}:
            raise ValueError(f"[salt.{sp}] conductance is {{background = uS/cm, slope = uS/cm per mg/L}}")
        surface = s.get("surface_store", False)
        if not isinstance(surface, bool):
            raise ValueError(f"[salt.{sp}] surface_store is true or false")
        return cls(precipitation=float(s["precipitation"]), road_salt=road,
                   point_sources={k: float(v) for k, v in s.get("point_sources", {}).items()},
                   soil_start=float(s.get("soil_start", 0.0)), stream_start=float(s.get("stream_start", 0.0)),
                   surface_store=surface, conductance=dict(cond) if cond else None)


def _read_table(path: Path, key: str, value: str, what: str) -> dict:
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows or key not in rows[0] or value not in rows[0]:
        raise ValueError(f"{what}: {path} needs columns {key}, {value}")
    return {r[key].strip(): float(r[value]) for r in rows}


# ---------- the land ----------

def land_day(m_soil, m_surface, area, w: dict, perc_mm, store_mm, deposition, road, from_groundwater,
             surface_store: bool = False) -> dict:
    """One day of every HRU's store (arrays, one value per HRU). ``m_soil``, ``m_surface``: kg at the start of the
    day; ``area`` m2; ``w``: SWAT+'s day in mm (``HRU_READS``' keys); ``perc_mm`` percolation; ``store_mm`` the
    store's water at the end of the day; ``deposition``, ``road``, ``from_groundwater``: kg in today.

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
        m_soil += release - sur_kg + from_groundwater
        q = {"lat": w["latq"], "til": w["qtile"], "perc": perc_mm}
    else:
        m_soil += deposition + road + from_groundwater
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

def add_weathering(gwt, rate, zones=None, pname: str = "weathering"):
    """Weathering as a mass source in a GWT model (an SRC package, ``pname``): ``rate`` kg per m3 of aquifer (the
    cell's bulk volume) a year [A], one value, or {zone: rate} with ``zones`` a zone number per cell in the grid's
    shape (layers first). Zero-rate cells get no entry. Call from an example's builder, after ``add_transport``."""
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
                 direct_returns: dict[str, list[str]], gw_soil=None, saved: dict | None = None):
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
        # road area per HRU, m2
        self.road_m2 = np.zeros(len(self.hrus))
        if settings.road_salt is not None:
            sh = settings.road_salt.share
            if isinstance(sh, dict):
                unknown = sorted(set(sh) - set(self.hrus))
                if unknown:
                    raise ValueError(f"salt: road_salt share names HRUs SWAT+ doesn't have: {unknown[:5]}")
                self.road_m2 = np.array([sh.get(h, 0.0) for h in self.hrus]) * self.area
            else:
                self.road_m2 = float(sh) * self.area
        # point sources and withdrawals
        self.points = sw.names(POINT) if sw.has(POINT) else []
        unknown = sorted(set(settings.point_sources) - set(self.points))
        if unknown:
            raise ValueError(f"salt: point_sources names {unknown}, which SWAT+ doesn't have")
        missing_route = [p for p in self.points if p not in self.route.points]
        if missing_route:
            raise ValueError(f"salt: point sources {missing_route} aren't in recall.con")
        self.direct_returns = direct_returns  # point source -> users returning to it directly
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
            self.total = float(saved["total"][0])
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
        self.c_body = np.zeros(nb)
        self.c_cell = None
        self.perc_kg = 0.0
        self.src_rates = np.zeros(len(self.loss_q))

    def state(self) -> dict:
        return {"m_soil": self.m_soil, "m_surface": self.m_surface, "m_body": self.m_body, "gain_due": self.gain_due,
                "loss_q": self.loss_q, "loss_body": self.loss_body, "soil_due": self.soil_due,
                "total": np.array([self._outside()])}

    def _outside(self) -> float:
        """The mass outside the aquifer, kg: land, water bodies, and the links in transit (the gains to the streams
        and the soil's water, due tomorrow)."""
        return float(self.m_soil.sum() + self.m_surface.sum() + self.m_body.sum() + self.gain_due.sum()
                     + self.soil_due.sum())

    # ---------- after SWAT+'s day ----------
    def after_surface(self, L, day: dt.date):
        from .coupling import swatplus as S
        k = self.k
        if self.total is None:
            self.total = self._outside()
        start_total = self.total
        w = {n: self.sw.get(v) for n, v in HRU_READS.items()}
        if not np.isfinite(w["precip"]).all():
            raise ValueError("salt: SWAT+ has no daily record for the day (print.prt's nyskip years)")
        perc_m = np.clip(self.sw.get(PERCOLATION), 0.0, None) * DAY
        store = w["sw"] + self.wp_mm + w["lagsurf"] + w["laglatq"] + w["wet"]
        dep = self.s.precipitation * np.clip(w["precip"], 0.0, None) / 1000.0 * self.area
        road = (application_days(w, self.s.road_salt, day) * self.road_m2 if self.s.road_salt is not None
                else np.zeros(len(self.hrus)))
        land0 = self.m_soil.sum() + self.m_surface.sum()
        r = land_day(self.m_soil, self.m_surface, self.area, w, perc_m * 1000.0, store, dep, road, self.soil_due,
                     self.s.surface_store)
        self.m_soil, self.m_surface = r["m_soil"], r["m_surface"]
        out = r["out"]
        received = float(self.soil_due.sum())
        self.soil_due = np.zeros(len(self.hrus))
        to_streams = sum(float(out[c].sum()) for c in ("sur", "lat", "til") if c in out)
        land_in = float(dep.sum() + road.sum()) + received
        land_out = to_streams + float(out["perc"].sum())
        land1 = float(self.m_soil.sum() + self.m_surface.sum())
        L.enter(k + "land:precipitation", float(dep.sum()))
        L.enter(k + "land:road_salt", float(road.sum()))
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
        rain[:nc] = self.s.precipitation * np.clip(self.sw.get(CH_PRECIP), 0.0, None)
        v_end, v_out, v_seep = np.zeros(nb), np.zeros(nb), np.zeros(nb)
        v_end[:nc] = self.sw.get(CH_VOLUME)
        v_out[:nc] = self.sw.get(S.QOUT) * DAY
        if self.reservoirs:
            rain[nc:] = self.s.precipitation * np.clip(self.sw.get(RES_PRECIP), 0.0, None) * DAY
            v_end[nc:] = self.sw.get(RES_VOLUME)
            v_out[nc:] = self.sw.get(RES_OUT) * DAY
            v_seep[nc:] = self.sw.get(RES_SEEP) * DAY
        inflow += rain
        points = np.zeros(nb)
        if self.points:
            vol = self.sw.get(POINT) * DAY
            for p, v in zip(self.points, vol):
                if v <= 0:
                    continue
                kg_u = sum(L.today(f"user:{u}:{self.sp}:kg") for u in self.direct_returns.get(p, []))
                v_u = sum(L.today(f"user:{u}:return") for u in self.direct_returns.get(p, []))
                rest = v - v_u
                rest = rest if rest > 1e-6 * max(v, 1.0) else 0.0  # SWAT+'s 4-byte rounding of the users' water
                if rest > 0.0 and p not in self.s.point_sources:
                    raise ValueError(f"salt: point source {p} carries {rest:.1f} m3 of water with no {self.sp} "
                                     f"stated for it ([salt.{self.sp}] point_sources) and no user reporting it")
                points[self.route.points[p]] += kg_u + self.s.point_sources.get(p, 0.0) * rest
        inflow += points
        inflow += self.gain_due
        v_wd = np.zeros(nb)
        if self.places:
            np.add.at(v_wd, self.place_body, np.clip(self.sw.get(WITHDRAWN), 0.0, None) * DAY)
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
        self._day = {"start": start_total, "in": float(dep.sum() + road.sum() + rain.sum() + points.sum()),
                     "out": float(q["outlet"].sum() + q["withdrawn"].sum() + q["seepage"].sum()),
                     "gross": float(dep.sum() + road.sum() + rain.sum() + points.sum()) + land_out + s_out}

    def concentration(self, channel: str) -> float:
        """A channel's outflow concentration today, kg m-3."""
        return float(self.c_body[self.ch_index[channel]])

    def conductance(self, channel: str) -> float | None:
        """Specific conductance, uS/cm, from the channel's chloride [A]; None without the study's relation."""
        if self.s.conductance is None:
            return None
        return self.s.conductance["background"] + self.s.conductance["slope"] * self.concentration(channel) * MG_L

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

    # ---------- after MODFLOW 6's day ----------
    def after_groundwater(self, L, b, mf6=None, soil_package: str | None = None, x: np.ndarray | None = None,
                          losses: list[tuple[np.ndarray, np.ndarray, list[str]]] = ()):
        """``b``: the transport model's MassBudget for the day. ``losses``: for each coupled river and lakes package,
        in the stream sources' order, each entry's flow (m3, positive into the aquifer), its channel's index in
        ``channels`` (-1: not coupled) and ``channels``. The gains due to the streams tomorrow, the water lost to the
        aquifer (its mass the aquifer's source tomorrow), the soil's water due to the land, and the system's
        closure."""
        k = self.k
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
