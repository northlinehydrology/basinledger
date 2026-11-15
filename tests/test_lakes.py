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
