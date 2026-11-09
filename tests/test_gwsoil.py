"""Groundwater up into the soil, the parts built before the rule is chosen (coupling/gwsoil.py): a cell's water split
among the HRUs over it by their shares of its ground, and among each HRU's layers below the water table by their
thickness below it; the profile's bottom under each cell. Volume is kept exactly."""

from __future__ import annotations

import numpy as np

from basinledger.coupling import gwsoil as G
from basinledger.mapping import AreaMap


def two_hrus():
    # cell 0: HRU a 60 m2, HRU b 40 m2; cell 1: HRU b 100 m2; cell 2: nothing
    return AreaMap(["a", "b"], np.array([0, 1, 1]), np.array([0, 0, 1]), np.array([60.0, 40.0, 100.0]))


def test_a_cells_water_goes_to_its_hrus_by_their_share_of_its_ground():
    s = G.cell_shares(two_hrus(), 3)
    np.testing.assert_allclose(s, [0.6, 0.4, 1.0])


def test_the_layers_below_the_water_table_take_it_by_thickness():
    tops, bottoms = np.array([0.0, 0.3, 1.0]), np.array([0.3, 1.0, 2.0])
    np.testing.assert_allclose(G.layer_split(tops, bottoms, 0.5), np.array([0.0, 0.5, 1.0]) / 1.5)
    np.testing.assert_allclose(G.layer_split(tops, bottoms, 0.0), [0.15, 0.35, 0.5])   # at the surface: all of it
    np.testing.assert_array_equal(G.layer_split(tops, bottoms, 2.5), [0.0, 0.0, 0.0])  # below the profile


def test_volume_is_kept_from_cells_to_layers():
    m = two_hrus()
    layer_hru = np.array([0, 0, 1, 1, 1])
    layer_bottom = np.array([0.4, 1.5, 0.2, 0.8, 1.2])
    cell_volume = np.array([10.0, 5.0, 0.0])
    row_volume = cell_volume[m.node] * G.cell_shares(m, 3)
    head_depth = np.array([0.5, 1.0, 9.0])           # depth to water per cell
    vol, at_bottom = G.to_layers(row_volume, head_depth[m.node], m, layer_hru, layer_bottom)
    assert at_bottom == 0.0
    assert vol.sum() == np.float64(row_volume.sum())
    # HRU a, 6 m3 at 0.5 m: its first layer is above the water table
    np.testing.assert_allclose(vol[:2], [0.0, 6.0])
    # HRU b: 4 m3 at 0.5 m (layers 0.3 and 0.4 m below it), 5 m3 at 1.0 m (only the 0.2 m of its last layer)
    np.testing.assert_allclose(vol[2:], [0.0, 4.0 * 0.3 / 0.7, 4.0 * 0.4 / 0.7 + 5.0])
    rates = G.layer_rates(vol, layer_hru, np.array([60.0, 140.0]))
    np.testing.assert_allclose((rates * 86400.0 * np.array([60.0, 140.0])[layer_hru]).sum(), 15.0)


def test_water_with_its_table_below_the_profile_goes_to_the_bottom_layer():
    m = two_hrus()
    vol, at_bottom = G.to_layers(np.array([1.0, 2.0, 3.0]), np.array([5.0, 5.0, 0.1]), m, np.array([0, 1, 1]),
                                 np.array([1.0, 0.5, 1.0]))
    # rows 0 and 1 (a and b in cell 0, water 5 m down) go to each HRU's bottom layer; b's 3 m3 in cell 1, water
    # 0.1 m down, splits 0.4 to 0.5 between its layers
    assert at_bottom == 3.0
    np.testing.assert_allclose(vol, [1.0, 3.0 * 0.4 / 0.9, 2.0 + 3.0 * 0.5 / 0.9])


def test_the_drain_sits_at_the_profiles_bottom_under_each_cell():
    m = two_hrus()
    depth = G.profile_depth(np.array([0, 0, 1, 1, 1]), np.array([0.4, 1.5, 0.2, 0.8, 1.2]), 2)
    np.testing.assert_allclose(depth, [1.5, 1.2])
    e = G.drain_elevation(m, np.array([100.0, 101.0, 102.0]), depth, 3)
    np.testing.assert_allclose(e[:2], [100.0 - (0.6 * 1.5 + 0.4 * 1.2), 101.0 - 1.2])
    assert np.isnan(e[2])


# --- on the models: Osu_1hru over the small aquifer, the water table inside the soil profile ----------------------

import warnings  # noqa: E402

import flopy  # noqa: E402
import pytest  # noqa: E402

from basinledger import config as C  # noqa: E402
from conftest import ROOT, RUNS  # noqa: E402

TOP = 88.5        # m: the small aquifer's heads start at the river's 88 m, half a metre down; Osu_1hru's profile is 0.75 m
DAYS = 365
CASES = [("drain", "none"), ("drain", "room"), ("gwflow", "none"), ("gwflow", "room")]


def need_gw_soil():
    import shutil
    import tempfile
    from swatplus_bmi import SwatPlusBmi
    from basinledger.config import swatplus_refdata
    with tempfile.TemporaryDirectory() as d:
        shutil.copytree(swatplus_refdata() / "Osu_1hru", f"{d}/p")
        m = SwatPlusBmi()
        m.initialize(f"{d}/p")
        has = G.INFLOW in m.get_output_var_names() and G.BOTTOM in m.get_output_var_names()
        m.finalize()
    if not has:
        pytest.skip(f"swatplus-bmi has no {G.INFLOW} (0.7.0)")


def gw_config(rule, cap, top=TOP, link=True):
    cfg = C.load(ROOT / "examples" / "synthetic" / "two_way.toml")
    cfg.run_folder = RUNS
    cfg.recharge_source, cfg.recharge_map = "hru", "all:hru0001"
    cfg.mf6_builder_args = {"top": top, "soil": "drn" if rule == "drain" else "wel"}
    cfg.groundwater_to_soil = {"rule": rule, "cap": cap, "package": "soil"} if link else None
    cfg.outputs = dict(cfg.outputs, balance="cha01")
    return cfg


def gw_run(cfg, name, days=DAYS, restart=None, save_at=None):
    from basinledger.driver import Basin
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run=name, users=False, days=days, restart=restart)
    if save_at:
        b.run(save_at - b.day)
        b.save(b.folder / "state")
    else:
        b.run()
    b.write_outputs()
    out = {"table": b.ledger.table(), "ws": b.mf6.ws, "folder": b.folder, "series": dict(b.series),
           "heads": b.mf6.head()}
    b.finalize()
    return out


def soil_budget(ws, text):
    cbc = flopy.utils.CellBudgetFile(ws / "gwf.cbc", precision="double")
    return np.array([float(d["q"].sum()) for d in cbc.get_data(text=text, paknam2="SOIL")])


@pytest.fixture(scope="module")
def runs(swatplus):
    need_gw_soil()
    return {(r, c): gw_run(gw_config(r, c), f"gw-soil-{r}-{c}") for r, c in CASES}


@pytest.mark.slow
@pytest.mark.parametrize("rule,cap", CASES)
def test_what_leaves_the_cells_is_what_the_layers_receive(runs, rule, cap):
    t = runs[(rule, cap)]["table"]
    sent, due, got = t["gw_soil:sent"], t["gw_soil:due"], t["gw_soil:received"]
    assert sent.sum() > 0
    assert due[0] == 0.0
    np.testing.assert_allclose(due[1:], sent[:-1], rtol=1e-12, atol=1e-9)      # the next day, whole
    np.testing.assert_allclose(got, due, rtol=1e-6, atol=1e-6)                  # SWAT+ keeps 4-byte reals
    # MODFLOW 6's budget: the soil package took out what was sent, every day
    out = -soil_budget(runs[(rule, cap)]["ws"], "DRN" if rule == "drain" else "WEL")
    np.testing.assert_allclose(out, sent, rtol=1e-12, atol=1e-9)


@pytest.mark.slow
@pytest.mark.parametrize("rule,cap", CASES)
def test_the_whole_basin_closes_with_the_link(runs, rule, cap):
    from basinledger.balance import basin_balance
    r = basin_balance(runs[(rule, cap)]["folder"], "cha01")
    t = runs[(rule, cap)]["table"]
    assert abs(r["groundwater_residual"]) <= 1e-6 * (r["recharge_received"] + r["groundwater_to_soil_sent"])
    # SWAT+'s own land budget counts what came up (gwsoilq), to the digits it prints
    assert r["groundwater_into_soil"] == pytest.approx(t["gw_soil:received"].sum(), rel=1e-3)
    assert r["groundwater_to_soil_sent"] == pytest.approx(t["gw_soil:sent"].sum(), rel=1e-12)


@pytest.mark.slow
def test_the_cap_holds_water_back(runs):
    for rule in ("drain", "gwflow"):
        free, capped = runs[(rule, "none")]["table"], runs[(rule, "room")]["table"]
        assert capped["gw_soil:sent"].sum() < free["gw_soil:sent"].sum()
    assert runs[("drain", "room")]["table"]["gw_soil:drains_raised"].max() > 0


@pytest.mark.slow
def test_a_restart_carries_what_is_due(swatplus):
    need_gw_soil()
    cfg = gw_config("drain", "none")
    whole = gw_run(cfg, "gw-soil-whole", days=120)
    first = gw_run(cfg, "gw-soil-first", days=120, save_at=60)
    rest = gw_run(cfg, "gw-soil-rest", days=120, restart=first["folder"] / "state")
    assert rest["table"]["gw_soil:due"][60] > 0
    np.testing.assert_array_equal(rest["table"]["gw_soil:sent"], whole["table"]["gw_soil:sent"])
    np.testing.assert_array_equal(rest["series"]["flow:cha01"], whole["series"]["flow:cha01"])
    np.testing.assert_array_equal(rest["heads"], whole["heads"])


@pytest.mark.slow
def test_a_water_table_below_every_profile_moves_nothing(swatplus):
    need_gw_soil()
    linked = gw_run(gw_config("drain", "none", top=100.0), "gw-soil-deep")
    alone = gw_run(gw_config("drain", "none", top=100.0, link=False), "gw-soil-deep-alone")
    assert np.all(linked["table"]["gw_soil:sent"] == 0.0)
    np.testing.assert_array_equal(linked["series"]["flow:cha01"], alone["series"]["flow:cha01"])
    np.testing.assert_array_equal(linked["heads"], alone["heads"])


def test_the_config_names_a_rule_and_a_cap(tmp_path):
    base = ('name = "x"\n[swatplus]\nproject = "p"\n[modflow6]\nbuilder = "m:b"\n'
            '[recharge]\nsource = "hru"\nmap = "all:hru0001"\n[exchange]\nmap = "all:cha01"\n')
    (tmp_path / "a.toml").write_text(base + '[groundwater_to_soil]\nrule = "drain"\npackage = "soil"\n')
    with pytest.raises(ValueError, match="no default"):
        C.load(tmp_path / "a.toml")
    (tmp_path / "b.toml").write_text(base.replace('source = "hru"', 'source = "aquifer"')
                                     + '[groundwater_to_soil]\nrule = "drain"\ncap = "none"\npackage = "soil"\n')
    with pytest.raises(ValueError, match="HRU map"):
        C.load(tmp_path / "b.toml")
    (tmp_path / "c.toml").write_text(base + '[groundwater_to_soil]\nrule = "uzf"\ncap = "none"\npackage = "soil"\n')
    with pytest.raises(ValueError, match="rule is one of"):
        C.load(tmp_path / "c.toml")
    (tmp_path / "d.toml").write_text(base + '[groundwater_to_soil]\nrule = "gwflow"\ncap = "room"\npackage = "soil"\n')
    assert C.load(tmp_path / "d.toml").groundwater_to_soil == {"rule": "gwflow", "cap": "room", "package": "soil"}
