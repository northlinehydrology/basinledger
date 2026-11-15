"""Several groundwater models coupled together (``[modflow6] models``; coupling/modflow6.py ``Models``; docs/design.md,
"Several groundwater models coupled together"): the synthetic foreign model and its neighbor east of it, joined by a
GWF-GWF exchange, both coupled.

- the models' cells in one sequence, each model's after the one before; a cell named by its model; a package in
  both models taken as one; the settings read and refused (tier 1);
- Osu_1hru's aquifer recharging both models through one map with a model column, a user's well added in the
  neighbor, 75 days: each model's recharge received is what the map sent it, the well delivers, the exchange
  between the two is booked inside the groundwater, each model's closure and the whole's hold, the start-up check
  lists both models' packages; a restart continues the unbroken run to the last digit (tier 2).
"""

from __future__ import annotations

import copy
import csv
import shutil
import warnings

import numpy as np
import pytest
from conftest import RUNS
from test_existing import CONFIG, START, SYNTHETIC, E, config

from basinledger import config as C
from basinledger.coupling.modflow6 import Models, Modflow6

MODELS = ["kzoo", "kzoo_east"]


@pytest.fixture
def models(libmf6):
    ws = RUNS / "models-unit"
    if ws.exists():
        shutil.rmtree(ws)
    E.build(ws, 31, START, neighbor=True, neighbor_recharge=True)
    m = Models(Modflow6(ws, "kzoo"), MODELS)
    yield m
    m.finalize()


def test_the_cells_in_one_sequence(models):
    kzoo, east = models.views
    assert models.offset.tolist() == [0, kzoo.nodes, kzoo.nodes + east.nodes]
    assert models.node(1, 5, 6) == kzoo.node(1, 5, 6)                       # the first model's, as before
    assert models.node(0, 3, 1, model="kzoo_east") == kzoo.nodes + east.node(0, 3, 1)
    assert models.node(model="kzoo_east", node=5) == kzoo.nodes + 5
    assert models.model_of(kzoo.nodes + 7) == (east, 7)
    nodes = models.package_nodes("rch_swb")                                 # in both models, one after the other
    assert len(nodes) == kzoo.nbound("rch_swb") + east.nbound("rch_swb") == models.nbound("rch_swb")
    assert (nodes[-east.nbound("rch_swb"):] >= kzoo.nodes).all()
    with pytest.raises(ValueError, match="'kzoo_west' isn't one of the coupled models"):
        models.node(0, 1, 1, model="kzoo_west")


def test_a_model_not_in_the_simulation_refused(libmf6):
    ws = RUNS / "models-unit-bad"
    if ws.exists():
        shutil.rmtree(ws)
    E.build(ws, 31, START, neighbor=True)
    mf = Modflow6(ws, "kzoo")
    try:
        with pytest.raises(ValueError, match="'kzoo_west' isn't a groundwater model of the simulation"):
            Models(mf, ["kzoo", "kzoo_west"])
        with pytest.raises(ValueError, match="the first is kzoo_east, the model started is kzoo"):
            Models(mf, ["kzoo_east", "kzoo"])
    finally:
        mf.finalize()


def test_the_models_settings_read(tmp_path):
    base = CONFIG.read_text().replace('"maps/', f'"{SYNTHETIC}/maps/')
    (tmp_path / "x.toml").write_text(base.replace('model = "kzoo"', 'models = ["kzoo", "kzoo_east"]'))
    cfg = C.load(tmp_path / "x.toml")
    assert cfg.mf6_models == MODELS and cfg.mf6_model == "kzoo"
    (tmp_path / "x.toml").write_text(base.replace('model = "kzoo"', 'model = "kzoo"\nmodels = ["kzoo"]'))
    with pytest.raises(ValueError, match="takes model .one. or models .several coupled together., not both"):
        C.load(tmp_path / "x.toml")
    (tmp_path / "x.toml").write_text(base.replace('model = "kzoo"', 'models = "kzoo"'))
    with pytest.raises(ValueError, match="models is a list"):
        C.load(tmp_path / "x.toml")


def both(cfg=None):
    """Both models coupled: Osu_1hru's 10 ha aquifer over 64 of kzoo's top cells and 36 of its neighbor's (all but
    its column at the regional head), and a second plant on a well added in the neighbor."""
    cfg = cfg or config(neighbor=True, neighbor_recharge=True)
    cfg.mf6_models = list(MODELS)
    path = RUNS / "models-recharge.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["object", "model", "layer", "row", "col", "area_m2"])
        w.writerows(["aqu001", "kzoo", 0, r, c, 1000.0] for r in range(1, 9) for c in range(1, 9))
        w.writerows(["aqu001", "kzoo_east", 0, r, c, 1000.0] for r in range(12) for c in range(3))
    cfg.recharge_map = str(path)
    cfg.wells.append(C.WellConfig("east1", "wel_east", (1, 5, 1), add=True, model="kzoo_east"))
    plant = copy.deepcopy(cfg.users[0])
    plant.name, plant.withdrawal = "plant2", "well:east1"
    cfg.users.append(plant)
    return cfg


def test_the_start_up_check_lists_both_models(swatplus):
    from basinledger.driver import Basin
    rows = []
    found, notes = Basin.check(both(), run="models-check", takeover=rows)
    assert found == []
    roles = {r.package: r.role for r in rows}
    assert roles["kzoo/rch_swb"] == roles["kzoo_east/rch_swb"] == "Basin Ledger's"
    assert roles["kzoo_east/wel_east"] == "Basin Ledger's" and roles["kzoo_east/chd_regional"] == "the model's"
    assert roles["kzoo/riv_creek"] == "partly Basin Ledger's"
    assert any(n.startswith("gwf-gwf_1: groundwater between kzoo and kzoo_east, both coupled") for n in notes)


def run(cfg, name, days=75, restart=None, until=None):
    from basinledger.driver import Basin
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run=name, days=days, restart=restart)
        try:
            if until is not None:
                b.run(until)
                b.save(RUNS / f"{name}-state")
            else:
                b.run()
            return b.ledger.table(), b.mf6.head(), b
        finally:
            b.finalize()


def test_both_models_coupled(swatplus):
    t, _, b = run(both(), "models-run")
    kzoo, east = b.mf6.views
    tot = {k: np.nansum(v) for k, v in t.items()}
    assert tot["recharge:received"] == pytest.approx(tot["recharge:made"], rel=1e-9)  # the whole 10 ha, both models
    rk, re = t["recharge:received:kzoo"], t["recharge:received:kzoo_east"]
    np.testing.assert_allclose(re, rk * 36 / 64, rtol=1e-9)  # by the map: 36 of the neighbor's cells, 64 of kzoo's
    assert re.sum() > 0
    assert tot["route:well:east1:requested"] > 0
    assert tot["route:well:east1:delivered"] == pytest.approx(tot["route:well:east1:requested"])
    q = t["groundwater:exchange:gwf-gwf_1"]
    assert np.isfinite(q).all() and np.abs(q).max() > 0  # across the shared boundary, inside the groundwater
    assert "groundwater:boundary:kzoo/gwf-gwf_1" not in t and "groundwater:boundary:kzoo_east/gwf-gwf_1" not in t
    for m in ("kzoo", "kzoo_east"):
        assert np.abs(t[f"groundwater:{m}:closure_share"]).max() < 1e-6
    assert np.abs(t["groundwater:closure_share"]).max() < 1e-6
    assert tot["groundwater:boundary:kzoo_east/chd_regional"] != 0.0      # the regional head, from outside
    assert tot["groundwater:boundary:kzoo_east/rch_swb"] == 0.0           # all of it SWAT+'s


def test_a_restart_of_both_models_continues_the_unbroken_run(swatplus):
    ta, ha, _ = run(both(), "models-unbroken")
    run(both(), "models-warm", until=40)
    tb, hb, b = run(both(), "models-restarted", restart=RUNS / "models-warm-state")
    assert np.array_equal(ha, hb)
    assert set(ta) == set(tb)
    for key in ta:
        assert np.array_equal(ta[key], tb[key], equal_nan=True), key


def test_streams_in_both_models(libmf6):
    """A channel's reaches in each model (kzoo's creek, the neighbor's losing river), one exchange through both
    (``StreamExchanges``), MODFLOW 6 alone: the supply limit holds each model's channel to what it carries, and the
    day's exchange comes back by channel from both."""
    from basinledger.coupling.exchange import (
        ChannelMap,
        StreamExchange,
        StreamExchanges,
    )
    ws = RUNS / "models-streams"
    if ws.exists():
        shutil.rmtree(ws)
    E.build(ws, 31, START, neighbor=True, neighbor_river=True)
    m = Models(Modflow6(ws, "kzoo"), MODELS)
    try:
        kz, ea = m.views
        m.prepare_day()  # the packages have their entries once the first day is prepared
        creek = ChannelMap.from_spec("maps/existing_river.csv", SYNTHETIC, kz, "riv_creek", ["cha01", "cha02"],
                                     channels=["cha01", "cha02"], first_model="KZOO")
        assert creek.channels == ["cha01"]
        river = ChannelMap(["cha02"], np.zeros(ea.nbound("riv_east"), int))
        ex = StreamExchanges([StreamExchange(kz, "riv_creek", creek), StreamExchange(ea, "riv_east", river)])
        assert ex.channels == ["cha01", "cha02"] and ex.kind == "riv"
        free = {}
        for day, have in enumerate(({"cha01": 1e6, "cha02": 1e6}, {"cha01": 1e6, "cha02": 0.0})):
            if day:
                m.prepare_day()
            m.solve_day(ex.limit(have))
            free[day] = ex.after_solve()
            lowered = ex.lowered
            ex.end_day()
        assert free[0]["cha02"] < -1.0                       # the neighbor's river loses, with water to lose
        # with none, the limit holds it to nothing, within its stage tolerance: 1 mm on 8 entries of 10 m2 d-1
        from basinledger.coupling.exchange import STAGE_TOLERANCE_M
        assert free[1]["cha02"] >= -8 * 10.0 * STAGE_TOLERANCE_M and free[1]["cha02"] > free[0]["cha02"] / 100
        assert lowered["cha02"] > 0 and lowered["cha01"] == 0.0
        assert np.isfinite(free[1]["cha01"])
    finally:
        m.finalize()
