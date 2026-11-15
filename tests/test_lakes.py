"""Lakes: SWAT+'s reservoir holds the lake's water, MODFLOW 6's LAK is its face to the aquifer (coupling/lakes.py;
docs/design.md, "Extending an existing MODFLOW 6 model", lakes).

- the lake's volume at a stage as MODFLOW 6 computes it, from its connections: over a day the change in that volume
  is MODFLOW 6's own lake storage (QSTO), and the stage at a volume inverts it (tier 1, MODFLOW 6 alone);
- a lake's table read; a reservoir's own seepage read from the project; what the start-up refuses (tier 1);
- coupled on the synthetic foreign model with a pond on its land surface and Osu_1hru with a reservoir below its
  channel, both forms, 75 days: the lake set from the reservoir each day, its leakage the reservoir's exchange the
  next day, as SWAT+ took it (a loss it couldn't supply booked as a shortfall), the LAK package Basin Ledger's in
  the start-up check, and the ledger closing every day (tier 2).
"""

from __future__ import annotations

import shutil
import warnings

import numpy as np
import pytest
from conftest import RUNS
from test_existing import E, config
from test_existing_api import _build

from basinledger import config as C
from basinledger.coupling import lakes as LK
from basinledger.coupling.modflow6 import Modflow6


@pytest.fixture
def api_model(libmf6):
    ws = RUNS / "lakes-api"
    if ws.exists():
        shutil.rmtree(ws)
    ws.mkdir(parents=True)
    _build(ws, days=3)
    mf = Modflow6(ws, "gwf")
    yield mf
    mf.finalize()


def test_the_lakes_volume_is_modflow_6s(api_model):
    mf = api_model
    for day in range(2):
        mf.prepare_day()
        c = LK.curve(mf, "lake", 1)
        for var in ("XOLDPAK", "XNEWPAK"):
            x = mf.get(var, "GWF", "LAKE")
            x[0] = 15.2
            mf.set(var, "lake", x)
        x = mf.get("INFLOW", "GWF", "LAKE")
        x[0] = 1000.0
        mf.set("INFLOW", "lake", x)
        mf.solve_day()
        end = float(mf.get("XNEWPAK", "GWF", "LAKE")[0])
        qsto = float(mf.get_model_units("QSTO", "GWF", "LAKE")[0])  # MODFLOW 6's: (old - new volume) a day
        assert c.volume(end) - c.volume(15.2) == pytest.approx(-qsto, rel=1e-9)
    for stage in (0.5, 15.2, 20.0, 31.0):
        assert c.stage(c.volume(stage)) == pytest.approx(stage, abs=1e-8)
    assert c.stage(0.0) == c.lo


def test_the_smoothing_is_modflow_6s():
    top, bot = np.array([1.0, 1.0, 2.0]), np.array([0.0, 1.0, 0.0])
    y = LK.eps_quadratic(top, bot, 0.5)
    assert y[0] == pytest.approx(0.5) and y[1] == 0.0 and y[2] == pytest.approx(0.25, rel=1e-5)  # eps 1e-6
    assert (LK.eps_quadratic(top, bot, 3.0) == 1.0).all() and (LK.eps_quadratic(top, bot, -1.0) == 0.0).all()


def test_a_lake_table_read(tmp_path):
    (tmp_path / "pond.tab").write_text("BEGIN dimensions\n  NROW 3\n  NCOL 3\nEND dimensions\n\nBEGIN table\n"
                                       "  10.0  0.0  0.0\n  11.0  500.0  1000.0\n  12.0  2000.0  2000.0\nEND table\n")
    c = LK.LakeCurve(table=LK.read_table(tmp_path / "pond.tab"))
    assert c.volume(9.0) == 0.0 and c.volume(11.5) == pytest.approx(1250.0) and c.volume(13.0) == 4000.0
    assert c.stage(1250.0) == pytest.approx(11.5)


def test_a_reservoirs_own_seepage_read(swatplus, tmp_path):
    p = E.reservoir_project(tmp_path / "osu", swatplus)
    assert LK.own_seepage(p) == {"res01": 0.0}
    hyd = p / "hydrology.res"
    hyd.write_text(hyd.read_text().replace(" 2.0 0.0 0.6", " 2.0 0.05 0.6"))
    assert LK.own_seepage(p) == {"res01": 0.05}


def test_what_the_start_up_refuses(api_model, swatplus, tmp_path):
    p = E.reservoir_project(tmp_path / "osu", swatplus)
    hyd = p / "hydrology.res"
    hyd.write_text(hyd.read_text().replace(" 2.0 0.0 0.6", " 2.0 0.05 0.6"))
    bad = LK.problems(api_model, [LK.LakeConfig("res01", "lake", 1, "drained"), LK.LakeConfig("res02", "lake"),
                                  LK.LakeConfig("res01", "river"), LK.LakeConfig("res01", "lake", 2)],
                      ["res01"], p)
    text = "\n".join(bad)
    assert "[[lakes]] res01: form is one of active, fixed_stage" in text
    assert "[[lakes]] res02: no such reservoir in the SWAT+ project" in text
    assert "its bottom conductivity k in hydrology.res is 0.05 mm/hr" in text
    assert "[[lakes]] res01: river is a SFR package, not a lake (LAK)" in text
    assert "[[lakes]] res01: lake has lakes 1 to 1; lake 2 isn't one" in text


def test_the_lakes_settings_read(tmp_path):
    from test_existing import CONFIG, SYNTHETIC
    base = CONFIG.read_text().replace('"maps/', f'"{SYNTHETIC}/maps/')
    (tmp_path / "x.toml").write_text(base + '\n[[lakes]]\nreservoir = "res01"\npackage = "lak_pond"\n'
                                     'form = "fixed_stage"\n')
    assert C.load(tmp_path / "x.toml").lakes == [LK.LakeConfig("res01", "lak_pond", 1, "fixed_stage")]
    (tmp_path / "x.toml").write_text(base + '\n[[lakes]]\nreservoir = "res01"\npackage = "lak_pond"\nform = "B"\n')
    with pytest.raises(ValueError, match="form is one of active, fixed_stage"):
        C.load(tmp_path / "x.toml")


def with_the_pond(form: str, swatplus):
    cfg = config(lake=True)
    cfg.swat_project = E.reservoir_project(RUNS / "existing-lake-osu", swatplus)
    cfg.lakes = [LK.LakeConfig("res01", "lak_pond", 1, form)]
    return cfg


@pytest.mark.parametrize("form", ["active", "fixed_stage"])
def test_a_lake_coupled(swatplus, form):
    from basinledger.driver import Basin
    probe = Basin.check(with_the_pond(form, swatplus), run=f"existing-check-lake-{form}", takeover=(rows := []))
    if any("swatplus-bmi 0.11.0" in p for p in probe[0]):
        pytest.skip("swatplus-bmi before 0.11.0: no reservoirs")
    assert probe[0] == []
    lak = {r.package: r for r in rows}["lak_pond"]
    assert lak.role == "Basin Ledger's" and "SWAT+'s reservoir res01" in lak.what
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(with_the_pond(form, swatplus), run=f"existing-run-lake-{form}", days=75)
        try:
            b.run()
            t = b.ledger.table()
        finally:
            b.finalize()
    k = "lake:res01:"
    leak, due, applied, short = t[k + "to_aquifer"], t[k + "due"], t[k + "applied"], t[k + "shortfall"]
    assert t[k + "volume"].max() > 100.0                       # the reservoir fills from the channel
    assert leak.min() < 0.0                                     # below the aquifer's heads, the pond gains
    np.testing.assert_allclose(due[1:], -leak[:-1], rtol=1e-12)  # today's leakage is due on the reservoir tomorrow
    np.testing.assert_allclose(applied + short, due, rtol=1e-6, atol=1e-6)
    assert (short <= 0).all()
    assert (t[k + "stage_m"] >= E.MIDDLE - 1e-9).all()          # the pond's bed is the bottom layer's top
    if form == "active":
        assert (t[k + "mf6_end_stage_m"] >= E.MIDDLE - 1e-9).all()
    assert abs(t["groundwater:closure_share"]).max() < 1e-6
    assert np.nansum(t["groundwater:boundary:lak_pond"]) == 0.0  # none of it the model's


def test_a_restart_carries_the_lake_and_the_irrigation(swatplus):
    """The pond (A) and the paddy's irrigation from wells together, saved on the paddy's third irrigated day and
    carried on: the lake's leakage due and the reservoir's water come back with the state, and the run is the unbroken
    run's, to the last digit."""
    from test_irrigation import FIRST, irrigated

    from basinledger.driver import Basin
    days, cut = FIRST + 10, FIRST + 3

    def cfg():
        return irrigated(with_the_pond("active", swatplus))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = Basin(cfg(), run="existing-restart-lake-unbroken", days=days)
        if a.lakes is None:
            a.finalize()
            pytest.skip("swatplus-bmi before 0.11.0: no reservoirs")
        a.run()
        ta, ha = a.ledger.table(), a.mf6.head()
        a.finalize()
        w = Basin(cfg(), run="existing-restart-lake-warm", days=days)
        w.run(cut)
        w.save(RUNS / "existing-restart-lake-state")
        w.finalize()
        b = Basin(cfg(), run="existing-restart-lake-restarted", days=days, restart=RUNS / "existing-restart-lake-state")
        b.run()
        tb, hb = b.ledger.table(), b.mf6.head()
        b.finalize()
    assert ta["irrigation:pumped"][cut:].sum() > 0 and np.nanmax(np.abs(ta["lake:res01:due"][cut:])) > 0
    assert np.array_equal(ha, hb)
    assert set(ta) == set(tb)
    for key in ta:
        assert np.array_equal(ta[key], tb[key], equal_nan=True), key


def outlet_model(ws):
    """Two lakes on a small model: lake 1 with a weir to outside and a specified outlet into lake 2."""
    import flopy
    if ws.exists():
        shutil.rmtree(ws)
    ws.mkdir(parents=True)
    sim = flopy.mf6.MFSimulation(sim_name="p", sim_ws=str(ws), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, time_units="days", nper=1, perioddata=[(3.0, 3, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="moderate")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf")
    idomain = np.ones((2, 10, 10), int)
    idomain[0, 4:6, 4:6] = 0
    idomain[0, 1, 1] = 0
    flopy.mf6.ModflowGwfdis(gwf, length_units="meters", nlay=2, nrow=10, ncol=10, delr=100., delc=100., top=20.,
                            botm=[0., -20.], idomain=idomain)
    flopy.mf6.ModflowGwfnpf(gwf, k=10., icelltype=1)
    flopy.mf6.ModflowGwfsto(gwf, sy=0.2, ss=1e-5, iconvert=1, transient={0: True})
    flopy.mf6.ModflowGwfic(gwf, strt=15.)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[((0, r, 0), 15.) for r in range(2, 10)])
    conn = [(0, k, (1, r, c), "vertical", 0.1, 0., 0., 0., 0.) for k, (r, c) in enumerate([(4, 4), (4, 5), (5, 4),
                                                                                          (5, 5)])]
    conn.append((1, 0, (1, 1, 1), "vertical", 0.1, 0., 0., 0., 0.))
    flopy.mf6.ModflowGwflak(gwf, pname="lake", nlakes=2, packagedata=[(0, 16., 4), (1, 15., 1)], connectiondata=conn,
                            noutlets=2, outlets=[(0, 0, -1, "weir", 15.5, 1.0, 1.0, 0.0),
                                                 (1, 0, 1, "specified", 0.0, 0.0, 0.0, 0.0)],
                            perioddata={0: [(1, "rate", -50.0)]}, length_conversion=1.0, time_conversion=86400.0)
    sim.write_simulation(silent=True)
    return ws


def test_a_coupled_lakes_outlets_are_cut(libmf6):
    mf = Modflow6(outlet_model(RUNS / "lakes-outlets"), "gwf")
    try:
        assert LK.outlets(mf, "lake") == [(1, 0), (1, 2)]
        mf.prepare_day()
        mf.solve_day()
        assert LK.outlet_flow(mf, "lake", 1) > 1000.0  # the weir and the 50 m3 a day into lake 2, as the model has them
        mf.prepare_day()
        LK.cut_outlets(mf, "lake", 1)
        mf.solve_day()
        assert LK.outlet_flow(mf, "lake", 1) == 0.0
        assert LK.problems(mf, [LK.LakeConfig("res01", "lake", 1)], ["res01"], None) == []  # its own outlets: cut
        bad = LK.problems(mf, [LK.LakeConfig("res02", "lake", 2)], ["res02"], None)
        assert bad == ["[[lakes]] res02: lake's outlet 2 brings lake 1's water into lake 2, water SWAT+'s reservoir "
                       "doesn't see; couple lake 1 too"]
        both = [LK.LakeConfig("res01", "lake", 1), LK.LakeConfig("res02", "lake", 2)]
        assert LK.problems(mf, both, ["res01", "res02"], None) == []
    finally:
        mf.finalize()
