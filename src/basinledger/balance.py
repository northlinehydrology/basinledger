"""The whole basin over a coupled run (testing G3): SWAT+'s own budgets, the links and MODFLOW 6's, in m3.

Over the area the surface model covers:

- **land** (SWAT+'s ``basin_wb_day``, summed over the run): precipitation + irrigation = ET + runoff to the channels
  + percolation + the gain in soil water and snow + a residual, which is SWAT+'s own (wetland storage, for one, isn't
  in it). Irrigation is counted as water from outside [A]: where a project irrigates from its own channels or
  aquifers, that is a transfer inside the basin, and the example should say so;
- **the link** (the ledger): what SWAT+'s aquifers recharged, what MODFLOW 6 received, and what the map added;
- **groundwater** (the ledger): recharge + boundary inflow + release from storage = discharge to the coupled streams
  + pumping, which closes by MODFLOW 6's budget;
- **channels and lakes**: runoff + groundwater SWAT+'s channels received (the exchange as applied, gains less
  losses) + precipitation on lakes and channels + water users' returns to the streams + inflow from elsewhere =
  evaporation from lakes and channels + withdrawals from the streams + outflow at the outlet + the gain in lake and
  channel storage + a residual. Precipitation, evaporation and storage come from the lakes' and channels' daily
  files (else the average annual ones times the run's years, with no storage term). Inflow from elsewhere is what
  SWAT+'s channels report taking in beyond their upstream channels, the lakes' releases, the land's runoff, the
  groundwater they gained and Basin Ledger's returns: the project's own point sources (recall objects), for one,
  which are water from outside. A channel's loss to the aquifer comes out of its outflow after routing (SWAT+'s
  seepage), not its inflow.

An example asks for it with ``[outputs] balance = "<outlet channel>"``: SWAT+ then prints ``basin_wb``,
``channel_sd`` and ``reservoir`` daily.
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
        "returns_to_streams": float(np.nansum(_sum(t, "route:swat:point_source:", ":delivered"))),
        "withdrawals_from_streams": float(np.nansum(_sum(t, "route:swat:place_of_use:", ":delivered"))),
        "irrigation": to_m3(mm("irr")) if "irr" in col else 0.0,
        "inflow_from_elsewhere": _other_inflow(sw, t),
        "outflow": out_gage,
        "max_daily_closure_share": float(np.nanmax(np.abs(t.get("groundwater:closure_share", np.zeros(1))))),
    }
    # with the exchange's point sources (swatplus-bmi before 0.2.0), the exchange is in the routes; don't count twice
    b["returns_to_streams"] -= float(np.nansum(_sum(t, "route:swat:point_source:gw_", ":delivered")))
    b["withdrawals_from_streams"] -= float(np.nansum(_sum(t, "route:swat:place_of_use:seep_", ":delivered")))
    b["land_residual"] = (b["precipitation_land"] + b["irrigation"] - b["et_land"] - b["runoff_to_channels"] - b["percolation"]
                          - b["soil_water_gain"] - b["snow_gain"])
    b["groundwater_residual"] = (b["recharge_received"] + b["groundwater_boundary_in"]
                                 + b["groundwater_storage_release"] - b["groundwater_to_streams"] - b["pumping"])
    b["lake_storage_gain"] = _storage_gain(sw / "reservoir_day.txt")
    b["channel_storage_gain"] = _storage_gain(sw / "channel_sd_day.txt")
    b["channels_residual"] = (- b["lake_storage_gain"] - b["channel_storage_gain"]
                              + b["runoff_to_channels"] + b["groundwater_into_channels"] + b["precipitation_lakes"]
                              + b["precipitation_channels"] + b["returns_to_streams"]
                              + b["inflow_from_elsewhere"] - b["evaporation_lakes"]
                              - b["evaporation_channels"] - b["withdrawals_from_streams"] - b["outflow"])
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
