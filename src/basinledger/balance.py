"""The whole basin over a coupled run (testing G3): SWAT+'s own budgets, the links and MODFLOW 6's, in m3.

Over the area the surface model covers:

- **land** (SWAT+'s ``basin_wb_day``, summed over the run): precipitation + irrigation + groundwater into the soil
  (``gwsoilq``) = ET + runoff to the channels + percolation + the gain in soil water and snow + a residual, which is SWAT+'s own (wetland storage, for one, isn't
  in it). Irrigation is SWAT+'s ``irr``, what reached the soil, from any source. What SWAT+'s own irrigation took
  from its channels and reservoirs (swatplus-bmi 0.14.0, booked by the driver) is reported beside it:
  ``irrigation_from_streams``, ``irrigation_from_lakes`` and the application's ``irrigation_loss``, which leaves
  the basin; the rest of ``irr`` is water from outside, pumped from MODFLOW 6's wells (``[irrigation]``), or put on
  the fields by ``land:`` routes (``irrigation_from_routes``), whose runoff share SWAT+ prints in its runoff, not in
  ``irr``, and the land takes as an input (``irrigation_runoff_from_routes``);
- **the link** (the ledger): what SWAT+'s aquifers recharged, what MODFLOW 6 received, and what the map added;
- **groundwater** (the ledger): recharge + boundary inflow + users' returns to the ground (``route:ground:``) +
  release from storage = discharge to the coupled streams + pumping + water up into the soil (``gw_soil:sent``),
  which closes by MODFLOW 6's budget;
- **channels and lakes**: runoff + groundwater SWAT+'s channels received (the exchange as applied, gains less
  losses) + precipitation on lakes and channels + water users' returns to the streams + inflow from elsewhere =
  evaporation from lakes and channels + withdrawals from the streams + irrigation from the streams and lakes +
  outflow at the outlet + the gain in lake and channel storage + a residual. Precipitation, evaporation and storage come from the lakes' and channels' daily
  files (else the average annual ones times the run's years, with no storage term). Inflow from elsewhere is what
  SWAT+'s channels report taking in beyond their upstream channels, the lakes' releases, the land's runoff, the
  groundwater they gained and Basin Ledger's returns: the project's own point sources (recall objects), for one,
  which are water from outside. A channel's loss to the aquifer comes out of its outflow after routing (SWAT+'s
  seepage), not its inflow.

An example asks for it with ``[outputs] balance = "<outlet channel>"``: SWAT+ then prints ``basin_wb``,
``channel_sd`` and ``reservoir`` daily.

**Recirculation** (with ``[groundwater_to_soil]``). Recharge is gross: some of the water MODFLOW 6 sends up into the
soil percolates again, the same day or later, and is counted again. The balance reports ``recharge_gross`` (what
MODFLOW 6 received from the land, ``recharge_received``), ``groundwater_to_soil`` (what it sent up,
``groundwater_to_soil_sent``), ``recharge_net`` (the difference) and ``recirculation_share`` (sent up over gross).
With ``rule = "water_table"``, the water the soils didn't take and returned the same day (``gw_soil:returned``) is
taken out of the first two and booked apart as ``recharge_returned``: it never entered the soil. Bookkeeping only:
nothing it reads changes a flow.

``recirculation`` splits it by HRU group (land use, by default) from the run's ``outputs/recirculation_hru.csv``,
which the driver writes with ``[outputs] recirculation = true``: each HRU's percolation, groundwater into the soil,
soil ET, lateral flow and saturation excess, summed from SWAT+'s day (``hru_wb``) over the run. Where the water that
came up went is an estimate [A]: SWAT+ mixes the soil's water, so it is taken to leave the soil in the proportions
all the soil's water left it over the run (well mixed). An exact split would need a tracer carried through SWAT+'s
soil.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .analysis import _sum, read_ledger


def _table(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = path.read_text().splitlines()
    head = lines[1].split()
    rows = [ln.split() for ln in lines[2:] if ln.split() and ln.split()[0].lstrip("-").isdigit()]
    return head, rows


def _aa(path: Path) -> dict[str, float]:
    if not path.exists():
        return {}
    head, rows = _table(path)
    return {k: float(v) for k, v in zip(head, rows[-1]) if _num(v)} if rows else {}


def _num(v: str) -> bool:
    try:
        float(v)
        return True
    except ValueError:
        return False


def basin_balance(run: str | Path, outlet: str) -> dict[str, float]:
    """The run's balance in m3 (``run`` holds ``swatplus/`` and ``outputs/``), with the residual of each part."""
    run = Path(run)
    sw = run / "swatplus"
    land = sum(float(r[2]) for r in _table_skip(sw / "ls_unit.def")) * 1e4
    head, rows = _table(sw / "basin_wb_day.txt")
    col = {k: head.index(k) for k in head}
    mm = lambda k: sum(float(r[col[k]]) for r in rows)  # noqa: E731
    d_sw = sum(float(r[col["sw_final"]]) - float(r[col["sw_init"]]) for r in rows)
    d_sno = sum(float(r[col["sno_final"]]) - float(r[col["sno_init"]]) for r in rows)
    aa = _aa(sw / "basin_wb_aa.txt")
    years = mm("precip") / aa["precip"] if aa.get("precip") else len(rows) / 365.25
    res, cha = _aa(sw / "basin_res_aa.txt"), _aa(sw / "basin_sd_cha_aa.txt")
    if not (sw / "reservoir.con").exists():
        res = {}
    t = read_ledger(run / "outputs" / "ledger.csv")
    import csv
    with open(run / "outputs" / "series.csv", newline="") as f:
        out_gage = sum(float(r[f"flow:{outlet}"]) for r in csv.DictReader(f)) * 86400.0
    to_m3 = lambda x: x / 1000.0 * land  # noqa: E731  mm over the land
    b = {
        "years": years, "land_m2": land,
        "precipitation_land": to_m3(mm("precip")),
        "et_land": to_m3(mm("et")),
        "runoff_to_channels": to_m3(mm("surq_gen") + mm("latq")),
        "percolation": to_m3(mm("perc")),
        "soil_water_gain": to_m3(d_sw), "snow_gain": to_m3(d_sno),
        "precipitation_lakes": _daily_sum(sw / "reservoir_day.txt", "precip", res.get("precip", 0.0) * years),
        "evaporation_lakes": _daily_sum(sw / "reservoir_day.txt", "evap", res.get("evap", 0.0) * years),
        "precipitation_channels": _daily_sum(sw / "channel_sd_day.txt", "precip", cha.get("precip", 0.0) * years),
        "evaporation_channels": _daily_sum(sw / "channel_sd_day.txt", "evap", cha.get("evap", 0.0) * years),
        "recharge_made": float(np.nansum(_sum(t, "recharge:made"))),
        "recharge_received": float(np.nansum(_sum(t, "recharge:received"))),
        "recharge_by_the_map": float(np.nansum(_sum(t, "recharge:by_the_map"))),
        "groundwater_boundary_in": float(np.nansum(_sum(t, "groundwater:boundary:"))),
        "groundwater_storage_release": float(np.nansum(_sum(t, "groundwater:storage_release"))),
        "groundwater_to_streams": float(np.nansum(_sum(t, "exchange:", ":to_stream"))),
        "groundwater_into_channels": float(np.nansum(_sum(t, "exchange:", ":applied"))),
        "pumping": float(np.nansum(_sum(t, "route:well:", ":delivered"))),
        "returns_to_ground": float(np.nansum(_sum(t, "route:ground:", ":delivered"))),
        "returns_to_streams": float(np.nansum(_sum(t, "route:swat:point_source:", ":delivered"))),
        "withdrawals_from_streams": float(np.nansum(_sum(t, "route:swat:place_of_use:", ":delivered"))),
        "irrigation": to_m3(mm("irr")) if "irr" in col else 0.0,
        "irrigation_from_streams": float(np.nansum(_sum(t, "irrigation:native:from_channels"))),
        "irrigation_from_lakes": float(np.nansum(_sum(t, "irrigation:native:from_reservoirs"))),
        "irrigation_loss": float(np.nansum(_sum(t, "irrigation:native:loss"))),
        # water put on the fields by land: routes: into the soil (in irr) and running off (SWAT+ prints it in
        # surq_cont, not irr, so the land takes it here)
        "irrigation_from_routes": float(np.nansum(_sum(t, "route:land:", ":delivered"))),
        "irrigation_runoff_from_routes": float(np.nansum(_sum(t, "route:land:", ":runoff"))),
        "groundwater_into_soil": to_m3(mm("gwsoilq")) if "gwsoilq" in col else 0.0,
        "groundwater_to_soil_sent": float(np.nansum(_sum(t, "gw_soil:sent"))),
        # rule = "water_table": what no soil layer took, back in the recharge received (a report), and what full
        # profiles spilled at the surface: a part of runoff_to_channels already, or on a paddy or wetland, of its
        # ponded water (reports, not added)
        "groundwater_to_soil_not_taken": float(np.nansum(_sum(t, "gw_soil:not_taken"))),
        "runoff_saturation_excess": float(np.nansum(_sum(t, "soil:saturation_excess", "excess"))),
        "saturation_excess_to_ponds": float(np.nansum(_sum(t, "soil:saturation_excess_to_ponds"))),
        "inflow_from_elsewhere": _other_inflow(sw, t),
        "outflow": out_gage,
        "max_daily_closure_share": float(np.nanmax(np.abs(t.get("groundwater:closure_share", np.zeros(1))))),
    }
    # with the exchange's point sources (swatplus-bmi before 0.2.0), the exchange is in the routes; don't count twice
    b["returns_to_streams"] -= float(np.nansum(_sum(t, "route:swat:point_source:gw_", ":delivered")))
    b["withdrawals_from_streams"] -= float(np.nansum(_sum(t, "route:swat:place_of_use:seep_", ":delivered")))
    # recirculation (bookkeeping only): gross recharge, what came back up into the soil, and the net; what
    # water_table's soils didn't take and returned the same day is neither, booked apart
    b["recharge_returned"] = float(np.nansum(_sum(t, "gw_soil:returned")))
    b["recharge_gross"] = b["recharge_received"] - b["recharge_returned"]
    b["groundwater_to_soil"] = b["groundwater_to_soil_sent"] - b["recharge_returned"]
    b["recharge_net"] = b["recharge_gross"] - b["groundwater_to_soil"]
    b["recirculation_share"] = (b["groundwater_to_soil"] / b["recharge_gross"] if b["recharge_gross"]
                                else float("nan"))
    b["land_residual"] = (b["precipitation_land"] + b["irrigation"] + b["irrigation_runoff_from_routes"]
                          + b["groundwater_into_soil"] - b["et_land"] - b["runoff_to_channels"] - b["percolation"]
                          - b["soil_water_gain"] - b["snow_gain"])
    b["groundwater_residual"] = (b["recharge_received"] + b["groundwater_boundary_in"] + b["returns_to_ground"]
                                 + b["groundwater_storage_release"] - b["groundwater_to_streams"] - b["pumping"]
                                 - b["groundwater_to_soil_sent"])
    b["lake_storage_gain"] = _storage_gain(sw / "reservoir_day.txt")
    b["channel_storage_gain"] = _storage_gain(sw / "channel_sd_day.txt")
    b["channels_residual"] = (- b["lake_storage_gain"] - b["channel_storage_gain"]
                              + b["runoff_to_channels"] + b["groundwater_into_channels"] + b["precipitation_lakes"]
                              + b["precipitation_channels"] + b["returns_to_streams"]
                              + b["inflow_from_elsewhere"] - b["evaporation_lakes"]
                              - b["evaporation_channels"] - b["withdrawals_from_streams"]
                              - b["irrigation_from_streams"] - b["irrigation_from_lakes"] - b["outflow"])
    return b


def _storage_gain(path: Path) -> float:
    """The change in storage (``flo_stor``, m3) over the run, summed over the objects of a SWAT+ daily file: each
    object's last day less its first (the first day's start isn't printed) [A]."""
    if not path.exists():
        return 0.0
    head, rows = _table(path)
    i, n = head.index("flo_stor"), head.index("name")
    first, last = {}, {}
    for r in rows:
        first.setdefault(r[n], float(r[i]))
        last[r[n]] = float(r[i])
    return sum(last[k] - first[k] for k in last)


def _daily_sum(path: Path, column: str, otherwise: float) -> float:
    """A column of a SWAT+ daily file summed over the run (m3), or ``otherwise`` (from the average annual file)."""
    if not path.exists():
        return otherwise
    head, rows = _table(path)
    return sum(float(r[head.index(column)]) for r in rows) if rows else otherwise


def _other_inflow(sw: Path, t: dict) -> float:
    """m3 SWAT+'s channels took in over the run from anything the balance doesn't otherwise count (channel_sd_day and
    reservoir_day, printed daily for it): each channel's inflow less its upstream channels' and lakes' outflow, less
    the land's runoff into it, the groundwater it received and Basin Ledger's point sources into it."""
    path = sw / "channel_sd_day.txt"
    if not path.exists():
        return float("nan")
    head, rows = _table(path)
    i_name, i_in, i_out = head.index("name"), head.index("flo_in"), head.index("flo_out")
    flo_in = sum(float(r[i_in]) for r in rows) * 86400.0
    out = {}
    for r in rows:
        out[r[i_name]] = out.get(r[i_name], 0.0) + float(r[i_out]) * 86400.0
    down = {}
    for f in ("chandeg.con", "reservoir.con"):
        p = sw / f
        if p.exists():
            for r in [ln.split() for ln in p.read_text().splitlines()[2:] if ln.split()]:
                down[r[1]] = {"sdc": f"cha{int(r[-3]):03d}", "res": f"res{int(r[-3]):03d}"}.get(r[-4])
    from_channels = sum(v for c, v in out.items() if down.get(c, "") and down[c].startswith("cha"))
    from_lakes = 0.0
    rp = sw / "reservoir_day.txt"
    if rp.exists():
        h2, r2 = _table(rp)
        j_name, j_out = h2.index("name"), h2.index("flo_out")
        from_lakes = sum(float(r[j_out]) for r in r2 if down.get(r[j_name], "").startswith("cha"))
    # a channel's inflow holds the groundwater it gained; a loss comes out of its outflow after routing (SWAT+'s
    # seepage), so only the gains come off here
    gw = float(sum(np.nansum(np.clip(v, 0.0, None)) for k, v in t.items()
                   if k.startswith("exchange:") and k.endswith(":applied")))
    returns = float(np.nansum(_sum(t, "route:swat:point_source:", ":delivered")))
    wb_head, wb = _table(sw / "basin_wb_day.txt")
    lsu = sum(float(r[2]) for r in _table_skip(sw / "ls_unit.def")) * 1e4
    runoff = sum(float(r[wb_head.index("surq_cha")]) + float(r[wb_head.index("latq_cha")]) for r in wb) / 1000 * lsu
    return flo_in - from_channels - from_lakes - runoff - gw - returns


def _table_skip(path: Path) -> list[list[str]]:
    """ls_unit.def: title, count, header, then rows (id, name, area ha, ...)."""
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    return [ln.split() for ln in lines[3:] if ln.split() and ln.split()[0].isdigit()]


def soil_water_trend(swat_dir: str | Path, days: int = 365) -> dict[str, float]:
    """SWAT+'s basin soil water over the last ``days`` days its ``basin_wb_day.txt`` prints (all of them if fewer):
    the change (mm, final less initial), the precipitation (mm) and the change as a share of the precipitation. A
    warm-up whose soil water still changes by a sizable share of its last year's precipitation hasn't settled
    (``scenarios.warm_up``)."""
    head, rows = _table(Path(swat_dir) / "basin_wb_day.txt")
    col = {k: head.index(k) for k in head}
    rows = rows[-days:]
    change = float(rows[-1][col["sw_final"]]) - float(rows[0][col["sw_init"]])
    precip = sum(float(r[col["precip"]]) for r in rows)
    return {"days": len(rows), "soil_water_change_mm": change, "precipitation_mm": precip,
            "share": abs(change) / precip if precip > 0 else float("nan"),
            "soil_water_final_mm": float(rows[-1][col["sw_final"]])}


# ---------- recirculation by HRU group ----------

# SWAT+'s day per HRU (mm, swatplus-bmi's printed hru_wb columns), summed over a run by HruWater
HRU_TERMS = {"precipitation": "hru_wb:precip", "irrigation": "hru_wb:irr", "percolation": "hru_wb:perc",
             "groundwater_to_soil": "hru_wb:gwsoilq", "et": "hru_wb:et", "et_plant": "hru_wb:eplant",
             "et_soil": "hru_wb:esoil", "lateral": "hru_wb:latq", "tile": "hru_wb:qtile",
             "saturation_excess": "hru_wb:satex", "surface_runoff": "hru_wb:surq_cont"}
# where the soil's water leaves it, for the well-mixed split [A]
LEAVES = {"percolation": ("percolation",), "et": ("et_plant", "et_soil"), "lateral": ("lateral", "tile"),
          "surface": ("saturation_excess",)}
# beside them, booked apart: what water_table's held layers didn't take (m s-1 per HRU), returned to MODFLOW 6
APART = {"groundwater_not_taken": "soil_water__groundwater_inflow_not_taken_rate"}
HRU_FILE = "recirculation_hru.csv"


class HruWater:
    """Each HRU's water over a run (m3), from SWAT+'s day as swatplus-bmi reads it (``HRU_TERMS``) and the HRUs'
    areas (``hru.con``, ``areas`` in m2): the driver's ``[outputs] recirculation``. Reads only. SWAT+ makes the records it reads from
    the first year (``records_every_year``); a NaN read is refused."""

    def __init__(self, sw, areas: dict[str, float], saved: dict | None = None):
        self.sw = sw
        self.hrus = sw.names(HRU_TERMS["percolation"])
        missing = [v for v in HRU_TERMS.values() if not sw.has(v)]
        if missing:
            raise ValueError(f"[outputs] recirculation reads SWAT+'s {missing}, which swatplus-bmi doesn't have")
        self.area = np.array([areas[h] for h in self.hrus])  # m2 (swatplus.hru_areas)
        self.total = {k: np.array(saved["total"][k], float) for k in HRU_TERMS} if saved else \
            {k: np.zeros(len(self.hrus)) for k in HRU_TERMS}
        self.apart = {}  # each APART term swatplus-bmi has: its HRUs in this order
        for k, v in APART.items():
            if sw.has(v):
                where = {h: i for i, h in enumerate(sw.names(v))}
                self.apart[k] = np.array([where[h] for h in self.hrus])
                kept = (saved or {}).get("total", {}).get(k)
                self.total[k] = np.zeros(len(self.hrus)) if kept is None else np.array(kept, float)
        self.days = int(saved["days"]) if saved else 0

    def day(self):
        """After SWAT+'s day."""
        for k, v in HRU_TERMS.items():
            mm = self.sw.get(v)
            if not np.isfinite(mm).all():
                raise ValueError(f"[outputs] recirculation: SWAT+ made no record of {v} for the day (print.prt's "
                                 "nyskip years)")
            self.total[k] += mm / 1000.0 * self.area
        for k, idx in self.apart.items():
            self.total[k] += np.clip(self.sw.get(APART[k]).astype(float)[idx], 0.0, None) * 86400.0 * self.area
        self.days += 1

    def state(self) -> dict:
        return {"days": self.days, "total": {k: v.tolist() for k, v in self.total.items()}}

    def write(self, path: Path):
        import csv
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            terms = [*HRU_TERMS, *self.apart]
            w.writerow(["hru", "area_m2", "days", *terms])
            for i, h in enumerate(self.hrus):
                w.writerow([h, repr(float(self.area[i])), self.days, *(repr(float(self.total[k][i]))
                                                                        for k in terms)])


def recirculation(run: str | Path, groups: str | dict = "landuse") -> dict[str, dict[str, float]]:
    """Recharge and the water that came back up into the soil, by HRU group, from a run's
    ``outputs/recirculation_hru.csv`` (``[outputs] recirculation = true``). ``groups``: ``"landuse"`` (each HRU's
    ``lu_mgt``), ``"hru"`` (each HRU alone) or a dict HRU -> group naming every HRU.

    Per group, and for the whole basin (``"basin"``): ``area_m2``, ``years``; ``percolation``, ``groundwater_to_soil``
    and ``net_recharge`` (the difference), each in m3 and in mm a year (``_mm_yr``); ``recirculation_share``
    (groundwater to soil over percolation); and where the water that came up went, an estimate [A]: SWAT+ mixes the
    soil's water, so it is taken to leave the soil in the proportions all the soil's water left it over the run
    (percolation, soil and plant ET, lateral and tile flow, saturation excess to the surface).
    ``from_groundwater_<path>`` (m3) and ``from_groundwater_<path>_share`` for each path, which add up to the
    groundwater to soil. ``groundwater_not_taken`` (m3 and mm a year) is booked apart, in none of them: what
    ``rule = "water_table"``'s held layers had no room for, returned to MODFLOW 6 as recharge without entering the
    soil (0 without it, or from a file written before it). Saturation excess is the part SWAT+ records as ``satex``; where it doesn't (its own gwflow
    off), what the full profile sheds is in its surface runoff with the rain's and isn't split out. An exact split
    would need a tracer through SWAT+'s soil."""
    import csv
    run = Path(run)
    path = run / "outputs" / HRU_FILE
    if not path.exists():
        raise FileNotFoundError(f"{path}: the run didn't write it; recirculation needs [outputs] recirculation = true")
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    hrus = [r["hru"] for r in rows]
    if groups == "landuse":
        from .salt import _hru_data
        names = _hru_data(run / "swatplus", hrus)["landuse"]
    elif groups == "hru":
        names = hrus
    elif isinstance(groups, dict):
        missing = sorted(set(hrus) - set(groups))
        if missing:
            raise ValueError(f"recirculation: groups names no group for HRUs {missing[:5]}")
        names = [groups[h] for h in hrus]
    else:
        raise ValueError(f"recirculation: groups is 'landuse', 'hru' or a dict HRU -> group; got {groups!r}")
    years = int(rows[0]["days"]) / 365.25 if rows else float("nan")
    sums: dict[str, dict[str, float]] = {}
    for g, r in zip(names, rows):
        for key in ("basin", g):
            s = sums.setdefault(key, dict.fromkeys(["area_m2", *HRU_TERMS, *APART], 0.0))
            for k in s:
                s[k] += float(r.get(k) or 0.0) if k in APART else float(r[k])  # an older file has no APART
    out = {}
    for g in sorted(sums, key=lambda k: (k == "basin", k)):
        s = sums[g]
        o = {"area_m2": s["area_m2"], "years": years, "percolation": s["percolation"],
             "groundwater_to_soil": s["groundwater_to_soil"],
             "net_recharge": s["percolation"] - s["groundwater_to_soil"],
             "groundwater_not_taken": s["groundwater_not_taken"]}
        for k in ("percolation", "groundwater_to_soil", "net_recharge", "groundwater_not_taken"):
            o[k + "_mm_yr"] = o[k] / s["area_m2"] / years * 1000.0 if s["area_m2"] and years else float("nan")
        o["recirculation_share"] = o["groundwater_to_soil"] / o["percolation"] if o["percolation"] else float("nan")
        leaving = {p: sum(max(s[k], 0.0) for k in ks) for p, ks in LEAVES.items()}
        total = sum(leaving.values())
        for p, v in leaving.items():
            share = v / total if total else float("nan")
            o[f"from_groundwater_{p}_share"] = share
            o[f"from_groundwater_{p}"] = s["groundwater_to_soil"] * share if total else 0.0
        out[g] = o
    return out
