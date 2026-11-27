"""Recirculation (balance.py, docs/design.md, Recirculation): recharge is gross where groundwater comes back up into the
soil. Bookkeeping only:
- the whole-basin balance and the summary give gross recharge, groundwater to soil and net recharge;
- ``[outputs] recirculation`` sums each HRU's water over the run, and ``recirculation`` splits it by HRU group: the
  groups add up to the basin, the HRUs' sums agree with the ledger, and the well-mixed split [A] adds up to what came
  up;
- it changes no flow: the ledger is the same to the bit with it on, and a restart carries its sums."""

from __future__ import annotations

import csv
import json

import numpy as np
import pytest
from test_gwsoil import gw_config, gw_run, need_gw_soil

from basinledger.balance import HRU_FILE, HRU_TERMS, basin_balance, recirculation


def _write(folder, rows, days=365):
    (folder / "outputs").mkdir(parents=True)
    with open(folder / "outputs" / HRU_FILE, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["hru", "area_m2", "days", *HRU_TERMS])
        for h, area, terms in rows:
            w.writerow([h, area, days, *(terms.get(k, 0.0) for k in HRU_TERMS)])


def test_the_groups_add_up_and_the_split_follows_the_soils_water(tmp_path):
    a = {"percolation": 300.0, "groundwater_to_soil": 100.0, "et_plant": 200.0, "et_soil": 100.0, "lateral": 50.0,
         "tile": 50.0}
    b = {"percolation": 100.0, "groundwater_to_soil": 0.0, "et_plant": 100.0}
    c = {"percolation": 50.0, "groundwater_to_soil": 60.0, "et_soil": 150.0}
    _write(tmp_path, [("hru1", 1e4, a), ("hru2", 2e4, b), ("hru3", 1e4, c)])
    r = recirculation(tmp_path, {"hru1": "wet", "hru2": "dry", "hru3": "wet"})
    assert list(r) == ["dry", "wet", "basin"]
    wet, basin = r["wet"], r["basin"]
    assert wet["area_m2"] == 2e4 and basin["area_m2"] == 4e4
    assert wet["percolation"] == 350.0 and wet["groundwater_to_soil"] == 160.0 and wet["net_recharge"] == 190.0
    assert wet["percolation_mm_yr"] == pytest.approx(350.0 / 2e4 * 1000.0 / (365 / 365.25))   # 365 days
    assert wet["recirculation_share"] == pytest.approx(160.0 / 350.0)
    # the wet group's soil water left by percolation 350, ET 450, lateral 100, the surface 0 (900 in all)
    for path, share in (("percolation", 350 / 900), ("et", 450 / 900), ("lateral", 100 / 900), ("surface", 0.0)):
        assert wet[f"from_groundwater_{path}_share"] == pytest.approx(share)
        assert wet[f"from_groundwater_{path}"] == pytest.approx(160.0 * share)
    for g in r.values():
        parts = sum(g[f"from_groundwater_{p}"] for p in ("percolation", "et", "lateral", "surface"))
        assert parts == pytest.approx(g["groundwater_to_soil"], abs=1e-12)
    assert basin["percolation"] == r["wet"]["percolation"] + r["dry"]["percolation"]
    assert recirculation(tmp_path, "hru")["hru2"]["recirculation_share"] == 0.0


def test_a_run_without_the_table_or_bad_groups_is_refused(tmp_path):
    with pytest.raises(FileNotFoundError, match="recirculation = true"):
        recirculation(tmp_path)
    _write(tmp_path, [("hru1", 1e4, {"percolation": 1.0})])
    with pytest.raises(ValueError, match="no group for HRUs"):
        recirculation(tmp_path, {"hru9": "x"})
    with pytest.raises(ValueError, match="groups is"):
        recirculation(tmp_path, "soil")


# ---------- on SWAT+ and MODFLOW 6: Osu_1hru with its water table in the soil (test_gwsoil.py) ----------

def _cfg(link=True, on=True):
    cfg = gw_config("drain", "none", link=link)
    if on:
        cfg.outputs = dict(cfg.outputs, recirculation=True)
    return cfg


@pytest.fixture(scope="module")
def run(swatplus):
    need_gw_soil()
    return gw_run(_cfg(), "recirculation")


@pytest.mark.slow
def test_gross_groundwater_to_soil_and_net(run):
    r = basin_balance(run["folder"], "cha01")
    t = run["table"]
    assert r["recharge_gross"] == r["recharge_received"] > 0
    assert r["groundwater_to_soil"] == r["groundwater_to_soil_sent"] > 0
    assert r["recharge_net"] == r["recharge_gross"] - r["groundwater_to_soil"]
    assert r["recirculation_share"] == pytest.approx(r["groundwater_to_soil"] / r["recharge_gross"])
    s = json.loads((run["folder"] / "outputs" / "summary.json").read_text())["recirculation"]
    assert s["recharge_gross"] == pytest.approx(np.nansum(t["recharge:received"]), rel=1e-12)
    assert s["groundwater_to_soil"] == pytest.approx(np.nansum(t["gw_soil:sent"]), rel=1e-12)
    assert s["recharge_net"] == pytest.approx(s["recharge_gross"] - s["groundwater_to_soil"], rel=1e-12)


@pytest.mark.slow
def test_each_hrus_sums_agree_with_the_ledger(run):
    t = run["table"]
    r = recirculation(run["folder"], "hru")
    hru = r["hru0001"]
    # SWAT+'s day read as 4-byte reals, summed: what came up is what SWAT+'s layers received, and the percolation
    # is what the HRU sent down (recharge:made)
    assert hru["groundwater_to_soil"] == pytest.approx(np.nansum(t["gw_soil:received"]), rel=1e-5)
    assert hru["percolation"] == pytest.approx(np.nansum(t["recharge:made"]), rel=1e-5)
    parts = sum(hru[f"from_groundwater_{p}"] for p in ("percolation", "et", "lateral", "surface"))
    assert parts == pytest.approx(hru["groundwater_to_soil"], rel=1e-12)
    assert hru["from_groundwater_percolation"] > 0 and hru["from_groundwater_et"] > 0
    landuse = recirculation(run["folder"])
    assert len(landuse) == 2 and all(v["percolation"] == hru["percolation"] for v in landuse.values())


@pytest.mark.slow
def test_it_changes_no_flow(run):
    off = gw_run(_cfg(on=False), "recirculation-off")
    a, b = run["table"], off["table"]
    assert a.keys() == b.keys() and all(a[k].tobytes() == b[k].tobytes() for k in a)
    assert "recirculation" in json.loads((off["folder"] / "outputs" / "summary.json").read_text())  # the link is on
    assert not (off["folder"] / "outputs" / HRU_FILE).exists()


@pytest.mark.slow
def test_without_the_link_the_recharge_is_net(swatplus):
    need_gw_soil()
    alone = gw_run(_cfg(link=False), "recirculation-alone", days=60)
    r = basin_balance(alone["folder"], "cha01")
    assert r["groundwater_to_soil"] == 0.0 and r["recharge_net"] == r["recharge_gross"]
    assert r["recirculation_share"] == 0.0 or r["recharge_gross"] == 0.0
    assert "recirculation" not in json.loads((alone["folder"] / "outputs" / "summary.json").read_text())
    assert recirculation(alone["folder"], "hru")["hru0001"]["groundwater_to_soil"] == 0.0


@pytest.mark.slow
def test_a_restart_carries_the_sums(swatplus):
    need_gw_soil()
    whole = gw_run(_cfg(), "recirculation-whole", days=120)
    first = gw_run(_cfg(), "recirculation-first", days=120, save_at=60)
    rest = gw_run(_cfg(), "recirculation-rest", days=120, restart=first["folder"] / "state")
    assert (rest["folder"] / "outputs" / HRU_FILE).read_text() == (whole["folder"] / "outputs" / HRU_FILE).read_text()
