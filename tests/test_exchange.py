"""Two-way stream-aquifer exchange on the synthetic example (examples/synthetic/two_way.toml).

Osu_1hru over a 10 ha MODFLOW 6 aquifer of the same area. SWAT+'s aquifer passes its recharge through; MODFLOW 6's
river exchange goes back to SWAT+'s channel the next day. Three years.
"""

from __future__ import annotations

import flopy
import numpy as np
import pytest

from basinledger import config as C
from basinledger.coupling import swatplus as S
from basinledger.driver import Basin
from conftest import ROOT, RUNS

pytestmark = pytest.mark.slow

DAYS = 1096


@pytest.fixture(scope="module")
def two_way(swatplus):
    cfg = C.load(ROOT / "examples" / "synthetic" / "two_way.toml")
    cfg.run_folder = RUNS
    b = Basin(cfg, run="two-way", days=DAYS)
    storage = []
    for _ in range(DAYS):
        b.step()
        storage.append(float(b.sw.get("aquifer_water__volume").sum()))
    out = dict(table=b.ledger.table(), due=dict(b.due), ws=b.mf6.ws, storage=np.array(storage))
    b.finalize()
    return out


def test_swat_aquifer_returns_nothing(two_way):
    t = two_way["table"]
    assert (t["swatplus:aquifer_outflow"] == 0).all()


def test_swat_aquifer_keeps_exactly_its_recharge(two_way):
    """No double counting: SWAT+'s aquifer only stores what it sent to MODFLOW 6."""
    t, v = two_way["table"], two_way["storage"]
    eps = np.finfo(np.float32).eps
    # SWAT+ holds storage as a 4-byte real, so a day's change is known to the rounding of the storage itself
    np.testing.assert_allclose(np.diff(v), t["recharge:sent"][1:], rtol=0, atol=4 * eps * v.max())
    assert v[-1] - v[0] == pytest.approx(t["recharge:sent"][1:].sum(), rel=1e-5)


def test_exchange_arrives_the_next_day(two_way):
    t = two_way["table"]
    k = "exchange:cha01:"
    np.testing.assert_allclose(t[k + "due"][1:], t[k + "to_stream"][:-1], rtol=0)
    assert t[k + "due"][0] == 0
    np.testing.assert_allclose(t[k + "due"], t[k + "applied"] + t[k + "shortfall"], rtol=1e-6, atol=1e-6)
    total = t[k + "to_stream"].sum()
    accounted = (t[k + "applied"] + t[k + "shortfall"]).sum() + two_way["due"]["cha01"]
    assert accounted == pytest.approx(total, rel=1e-6)


def test_the_river_gains_and_loses(two_way):
    t = two_way["table"]
    x = t["exchange:cha01:to_stream"]
    assert (x > 0).sum() > 100 and (x < 0).sum() > 10
    applied = t["exchange:cha01:applied"]
    assert (applied > 0).sum() > 100 and (applied < 0).sum() > 0


def test_shortfalls_are_reported_not_hidden(two_way):
    """A losing reach SWAT+ can't supply is entered as a shortfall; SWAT+ never takes more than was due."""
    t = two_way["table"]
    due, applied, short = (t[f"exchange:cha01:{k}"] for k in ("due", "applied", "shortfall"))
    assert (short <= 1e-9).all()
    loss = due < 0
    assert (applied[loss] >= due[loss] * (1 + 1e-6) - 1e-6).all()
    np.testing.assert_allclose(due, applied + short, rtol=1e-6, atol=1e-6)


def test_modflow_budget_closes_and_matches_the_ledger(two_way):
    ws = two_way["ws"]
    cbc = flopy.utils.CellBudgetFile(ws / "gwf.cbc", precision="double")
    terms = {}
    for text in ("STO-SS", "STO-SY", "RCH", "WEL", "RIV"):
        terms[text] = np.array([float(d["q"].sum()) if d.dtype.names else float(np.sum(d))
                                for d in cbc.get_data(text=text)])
    net = sum(terms.values())
    gross = sum(np.abs(v) for v in terms.values())
    assert np.max(np.abs(net) / gross) < 1e-7
    np.testing.assert_allclose(-terms["RIV"], two_way["table"]["exchange:cha01:to_stream"], rtol=1e-12)
    np.testing.assert_allclose(terms["RCH"], two_way["table"]["recharge:sent"], rtol=1e-12, atol=1e-9)


def test_pass_through_settings(tmp_path, swatplus):
    import shutil
    d = tmp_path / "p"
    shutil.copytree(swatplus / "Osu_1hru", d)
    S.aquifer_pass_through(d)
    lines = (d / "aquifer.aqu").read_text().splitlines()
    cols = lines[1].split()
    for l in lines[2:]:
        f = l.split()
        assert float(f[cols.index("flo_min")]) == -1e6 and float(f[cols.index("rchg_dp")]) == 0


def test_stage_from_swat_depth(swatplus):
    """With swatplus-bmi 0.2.0, river stage = RBOT + bed thickness + SWAT+'s channel depth (the pause's, yesterday's)."""
    if not S.has_groundwater_exchange():
        pytest.skip("needs swatplus-bmi 0.2.0 (channel depth)")
    cfg = C.load(ROOT / "examples" / "synthetic" / "two_way.toml")
    cfg.run_folder = RUNS
    cfg.rivers = {"package": "riv", "map": "all:cha01", "stage": "surface_depth", "bed_thickness": 1.0}
    cfg.exchange = dict(cfg.exchange, limit=False)  # the stage as set from SWAT+'s depth, before any lowering
    b = Basin(cfg, run="two-way-stage", days=200)
    stages, depths = [], []
    for _ in range(200):
        b.step()
        stages.append(b.mf6.get("STAGE", "GWF", "RIV").copy())
        depths.append(float(b.sw.get(S.DEPTH)[0]))  # day d's depth, which set day d's stage
    rbot = b.mf6.get("RBOT", "GWF", "RIV")
    t = b.ledger.table()
    b.finalize()
    # the stage set on day d is RBOT + 1 m + the depth SWAT+ reported for day d (4-byte reals)
    st = np.array(stages)[:, 0]
    np.testing.assert_allclose(st, rbot[0] + 1.0 + np.array(depths), rtol=0, atol=1e-5)
    assert np.ptp(st) > 0
    assert (t["exchange:cha01:due"] != 0).sum() > 100


def test_limit_to_supply_fits_the_loss():
    """The stage is lowered only where a channel's net exchange would take more than it has, never below the bottom,
    and then takes no more than it has; a channel's gaining entries count against its losses."""
    from basinledger.coupling.exchange import limit_to_supply
    stage = np.array([10.0, 10.0, 10.0, 5.0])
    rbot = np.array([9.0, 9.0, 9.0, 4.0])
    cond = np.array([100.0, 100.0, 100.0, 50.0])
    head = np.array([8.0, 9.5, 11.0, 3.0])  # below the bottom, between, above (gaining); channel 2 below
    ch = np.array([0, 0, 0, 1])
    net = lambda st: (cond * (np.clip(st - np.maximum(head, rbot), 0, None) - np.clip(head - st, 0, None)))  # noqa
    assert net(stage)[:3].sum() == pytest.approx(50.0)  # losses of 150 less a gain of 100
    new, d = limit_to_supply(stage, rbot, cond, head, ch, np.array([20.0, 1000.0]))
    assert net(new)[:3].sum() == pytest.approx(20.0, rel=1e-9)
    assert d[0] > 0 and d[1] == 0 and new[3] == stage[3]
    assert (new >= rbot).all()
    new, d = limit_to_supply(stage, rbot, cond, head, ch, np.array([0.0, 0.0]))
    assert net(new)[:3].sum() <= 1e-9 and net(new)[3] == pytest.approx(0.0, abs=1e-9)
    assert (new >= rbot).all()


def test_the_limit_removes_the_shortfall(swatplus):
    """Lowering the stage where the channel can't supply the loss leaves a smaller shortfall than booking it."""
    out = {}
    for limit in (False, True):
        cfg = C.load(ROOT / "examples" / "synthetic" / "two_way.toml")
        cfg.run_folder = RUNS
        cfg.exchange = dict(cfg.exchange, limit=limit)
        b = Basin(cfg, run=f"two-way-limit-{limit}", days=DAYS)
        b.run()
        t = b.ledger.table()
        b.finalize()
        out[limit] = (-t["exchange:cha01:shortfall"].sum(), t["exchange:cha01:stage_lowered_m"] if limit else None)
    assert out[False][0] > 100.0  # m3: the synthetic example's losing days
    # what remains is the lag: today's loss is set against today's flow and taken from tomorrow's
    assert out[True][0] < 0.5 * out[False][0]
    assert (out[True][1] > 0).any()


def test_recharge_made_is_recharge_sent(two_way):
    """The two-way grid is exactly the aquifer's area: the map adds nothing and places everything."""
    t = two_way["table"]
    np.testing.assert_allclose(t["recharge:made"], t["recharge:sent"], rtol=1e-9, atol=1e-9)
    assert np.abs(t["recharge:by_the_map"]).max() < 1e-6


def test_the_area_check_is_on_by_default(swatplus):
    """A map whose area isn't the aquifer's stops the run unless the example turns the check off."""
    cfg = C.load(ROOT / "examples" / "synthetic" / "basin.toml")
    cfg.run_folder = RUNS
    assert cfg.area_rtol is None  # basin.toml spreads the aquifer over the whole grid on purpose
    cfg.area_rtol = 1e-6
    with pytest.raises(ValueError, match="mapped area"):
        Basin(cfg, run="area-check", days=1)


def test_no_water_from_outside_in_the_two_way_example(two_way):
    """Every MODFLOW 6 flow in the two-way example is a link (recharge, the well, the exchange), so the boundary
    accounts are zero, and the aquifer's storage closes with the links."""
    t = two_way["table"]
    for k in [k for k in t if k.startswith("groundwater:boundary:")]:
        assert np.abs(t[k]).max() < 1e-9, k
    assert np.abs(t["groundwater:closure_share"]).max() < 1e-6
    links = t["recharge:received"] - t["exchange:cha01:to_stream"] + t["groundwater:storage_release"]
    wells = [k for k in t if k.startswith("route:well:") and k.endswith(":delivered")]
    assert wells
    np.testing.assert_allclose(links - sum(t[k] for k in wells), 0.0, atol=1e-6 * np.abs(t["recharge:received"]).max())


def test_swat_settings_that_would_fight_are_refused(tmp_path, swatplus):
    import shutil
    d = tmp_path / "p"
    shutil.copytree(swatplus / "Osu_1hru", d)
    S.aquifer_pass_through(d)
    assert any("transmission loss" in p for p in S.coupling_problems(d, two_way=True))  # Osu_1hru's channels lose
    S.no_transmission_loss(d)
    assert S.coupling_problems(d, two_way=True) == []
    (d / "file.cio").write_text((d / "file.cio").read_text().replace("cal_parms.cal     null",
                                                                       "cal_parms.cal     calibration.cal"))
    (d / "calibration.cal").write_text("calibration.cal\n 1\n NAME CHG_TYP VAL\nflo_min absval 5.0 0 0 0 0 0 0 0 0\n")
    assert any("flo_min" in p for p in S.coupling_problems(d, two_way=True))
    assert S.coupling_problems(d, two_way=False) == []


def test_drains_reach_the_stream(swatplus):
    """Seepage at land surface (a drain) is coupled to the channel like the river's gains: it arrives the next day,
    the drain is no boundary flow, and the aquifer closes with the links."""
    cfg = C.load(ROOT / "examples" / "synthetic" / "two_way.toml")
    cfg.run_folder = RUNS
    cfg.mf6_builder = "model:build_small"
    cfg.mf6_builder_args = {"drain_elevation": 87.5}  # half a meter below the river's stage (88 m)
    cfg.exchange = dict(cfg.exchange, drains={"package": "drn", "map": "all:cha01"})
    b = Basin(cfg, run="two-way-drains", days=400)
    b.run()
    t = b.ledger.table()
    b.finalize()
    seep = t["exchange:cha01:seepage"]
    assert (seep > 0).sum() > 50 and (seep >= -1e-9).all()
    assert np.abs(t["groundwater:boundary:drn"]).max() < 1e-9
    assert np.abs(t["groundwater:closure_share"]).max() < 1e-6
    np.testing.assert_allclose(t["exchange:cha01:due"][1:], t["exchange:cha01:to_stream"][:-1], rtol=0)


def test_spin_up_settles_storage(swatplus):
    """Cycling the first year from the last cycle's heads brings the year's release from storage toward zero, and a
    run from those heads starts where the spin-up ended."""
    cfg = C.load(ROOT / "examples" / "synthetic" / "two_way.toml")
    cfg.run_folder = RUNS
    heads, report = Basin.spin_up(cfg, days=365, cycles=4, tol=1e-4, run="two-way-spinup")
    shares = [abs(r["share"]) for r in report]
    assert len(report) >= 2 and shares[-1] < shares[0]
    b = Basin(cfg, run="two-way-spun", days=2, heads=heads)
    np.testing.assert_array_equal(b.mf6.head(), heads)
    b.finalize()
