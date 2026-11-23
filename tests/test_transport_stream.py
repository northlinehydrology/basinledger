"""Groundwater transport, stage 2: the stream side (coupling/surface_mass.py), decay, sorbed storage in the ledger
and the steady-state start (coupling/transport.py).

- decay (first order, dissolved and sorbed) is booked, and the mass still closes every day (tier 1);
- sorbed mass is storage in the ledger: the change in the dissolved and sorbed mass the cells hold is the storage
  release, every day (tier 1);
- a stream's mass source is booked to its channel (tier 1);
- HRU by HRU: recharge mass is placed on the cells the water goes to, in the same shares (tier 1);
- the steady state: uniform recharge gives its concentration everywhere, and with decay and no dispersion it is the
  transport model's own long-run answer (tier 3);
- the coupled synthetic run with nitrate and soluble P (examples/synthetic/nitrate.toml): the mass closes every day,
  a gain's mass arrives in the channel, a loss takes the channel's concentration, the soil link and recharge carry
  SWAT+'s mass, a restart carries it all, and the steady start (tier 2).
"""

from __future__ import annotations

import datetime as dt
import importlib.util
import warnings

import numpy as np
import pytest

from basinledger.coupling.modflow6 import Modflow6
from basinledger.coupling.surface_mass import SurfaceMass, spread
from basinledger.coupling.transport import Transport, add_transport
from basinledger.mapping import AreaMap
from conftest import ROOT, RUNS


def synthetic_model():
    spec = importlib.util.spec_from_file_location("synthetic_model", ROOT / "examples" / "synthetic" / "model.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ---------- MODFLOW 6 alone ----------

def small_run(ws, days, species="no3", **settings):
    """build_small with drains and a transport model, driven as the driver does, the river losing on alternate weeks
    (the well pumps hard), the mass sources on the river's cells put to work. Returns each day's budget and the
    model's dissolved and sorbed mass before and after."""
    m = synthetic_model()
    m.build_small(ws, days, dt.date(2010, 1, 1), drain_elevation=89.0,
                  transport={species: dict(porosity=0.25, dispersivity=5.0, strt=0.01, stream_side=True, **settings)})
    mf = Modflow6(ws)
    t = Transport(mf, species)
    river = (np.zeros(10, int), ["cha01"])
    vol = mf.cell_area * 30.0  # 100 - 70 m
    out = []
    for d in range(days):
        mf.prepare_day()
        mf.set_recharge("rch", np.full(mf.nodes, 0.002 + 0.002 * np.sin(d / 5.0)))
        mf.set_well("wel", 0, -900.0 if d % 14 < 7 else 0.0)
        t.set_inflow("rch", 0.03)
        t.set_inflow("riv", 0.0)
        src = np.zeros(10)
        src[d % 10] = 0.05 * (d % 3)  # kg d-1 on one river cell
        t.set_stream_mass("riv", src)
        t.set_stream_mass("drn", np.zeros(10))
        x0 = t.concentration().copy()
        sat0 = mf.get("SAT", "GWF", "NPF").copy()
        mf.solve_day(lambda: False)
        b = t.budget("rch", {"pw1": ("wel", 0)}, [("riv", *river), ("drn", *river)])
        out.append((b, x0, t.concentration().copy(), sat0, mf.get("SAT", "GWF", "NPF").copy(), src.sum()))
    mf.finalize()
    return out, vol


@pytest.mark.parametrize("sorbs", [False, True])
def test_decay_is_booked_and_the_mass_closes(libmf6, sorbs):
    sorption = {"distcoef": 2e-3, "bulk_density": 1700.0} if sorbs else None
    days, vol = small_run(RUNS / f"stream-decay-{sorbs}", 60, decay=0.002, sorption=sorption,
                          decay_sorbed=0.001 if sorbs else None)
    for b, x0, x1, s0, s1, _ in days:
        assert abs(b.closure_share) < 1e-9
        assert b.decay > 0
        # MODFLOW 6's implicit step: decay at the end of the day's concentrations, dissolved and sorbed
        want = 0.002 * np.sum(0.25 * vol * s1 * x1) + (0.001 * 1700.0 * 2e-3 * np.sum(vol * x1) if sorbs else 0.0)
        assert b.decay == pytest.approx(want, rel=1e-9)


def test_sorbed_mass_is_storage_in_the_ledger(libmf6):
    """With linear sorption, what the cells hold dissolved (porosity x volume x C) and sorbed (bulk density x Kd x
    volume x C) changes by exactly the storage release, mobile and sorbed, every day, and the sorbed part is booked on
    its own. Steady flow, so the water each cell holds doesn't change (with transient flow MODFLOW 6 moves the mobile
    store's volume with the flow model's storage, and the closure is what's checked: test_decay...)."""
    import flopy

    from test_transport import flow_model
    kd, rho, theta, days = 2e-3, 1700.0, 0.3, 200
    ws = RUNS / "stream-sorbed"
    sim, gwf = flow_model(ws, days, steady=True, nrow=10, ncol=10, cell=20.0, top=13.0, botm=10.0, k=20.0)
    flopy.mf6.ModflowGwfic(gwf, strt=12.0)
    flopy.mf6.ModflowGwfrcha(gwf, recharge=0.005, pname="rch")
    flopy.mf6.ModflowGwfriv(gwf, stress_period_data={0: [((0, r, 0), 12.0, 50.0, 11.0) for r in range(10)]},
                            pname="riv")
    add_transport(sim, gwf, "solp", theta, 2.0, strt=0.0, sources=[("rch", "solp")],
                  sorption={"distcoef": kd, "bulk_density": rho})
    sim.write_simulation(silent=True)
    mf = Modflow6(ws)
    t = Transport(mf, "solp")
    vol = 400.0 * 3.0
    budgets = []
    try:
        for d in range(days):
            mf.prepare_day()
            t.set_inflow("rch", 0.025 if d < 100 else 0.0)  # loading, then clean water: the sorbed store releases
            x0 = t.concentration().copy()
            mf.solve_day()
            b = t.budget("rch", {}, [("riv", np.zeros(10, int), ["r"])])
            x1 = t.concentration()
            assert abs(b.closure_share) < 1e-9
            dissolved, sorbed = theta * vol * np.sum(x1 - x0), rho * kd * vol * np.sum(x1 - x0)
            assert -b.storage_release == pytest.approx(dissolved + sorbed, rel=1e-9, abs=1e-14)
            assert -b.sorbed_release == pytest.approx(sorbed, rel=1e-9, abs=1e-14)
            budgets.append(b)
    finally:
        mf.finalize()
    assert sum(b.sorbed_release < 0 for b in budgets[:100]) == 100   # taking up P while loaded
    assert sum(b.sorbed_release > 0 for b in budgets[100:]) == 100   # and giving it back after


def test_a_streams_mass_source_is_booked_to_its_channel(libmf6):
    days, _ = small_run(RUNS / "stream-source", 30)
    for b, *_, src in days:
        assert abs(b.closure_share) < 1e-9
        assert b.returned["cha01"] == pytest.approx(src, rel=1e-12, abs=1e-15)
        assert not b.boundary  # the mass sources on coupled cells are the stream's, not a boundary


def test_the_mass_source_must_mirror_its_package(tmp_path, libmf6):
    import flopy
    m = synthetic_model()
    ws = tmp_path / "w"
    m.build_small(ws, 3, dt.date(2010, 1, 1), transport={"no3": dict(porosity=0.25, dispersivity=5.0)})
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(ws), verbosity_level=0)
    gwt = sim.get_model("gwt_no3")
    flopy.mf6.ModflowGwtsrc(gwt, pname="src_riv", filename="gwt_no3.src_riv", maxbound=10,
                            stress_period_data={0: [((0, r, 1), 0.0) for r in range(10)]})  # the wrong column
    sim.write_simulation(silent=True)
    mf = Modflow6(ws)
    try:
        t = Transport(mf, "no3")
        mf.prepare_day()
        with pytest.raises(ValueError, match="aren't riv's"):
            t.set_stream_mass("riv", np.zeros(10))
    finally:
        mf.finalize()


def test_kd_by_zone(tmp_path):
    import flopy
    m = synthetic_model()
    zones = np.ones((1, 10, 10), int)
    zones[0, :, 5:] = 2
    m.build_small(tmp_path / "z", 3, dt.date(2010, 1, 1),
                  transport={"solp": dict(porosity=0.25, dispersivity=5.0,
                                          sorption={"distcoef": {1: 1e-3, 2: 5e-3}, "zones": zones,
                                                    "bulk_density": 1700.0})})
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(tmp_path / "z"), verbosity_level=0)
    kd = sim.get_model("gwt_solp").get_package("mst").distcoef.array
    assert kd[0, 0, 0] == 1e-3 and kd[0, 0, 9] == 5e-3
    with pytest.raises(ValueError, match="zones \\[3\\] have no value"):
        m.build_small(tmp_path / "z2", 3, dt.date(2010, 1, 1),
                      transport={"solp": dict(porosity=0.25, dispersivity=5.0,
                                              sorption={"distcoef": {1: 1e-3, 2: 5e-3}, "zones": zones + 1,
                                                        "bulk_density": 1700.0})})


def test_spread():
    np.testing.assert_allclose(spread(6.0, np.array([1.0, 2.0, 0.0]), np.array([1.0, 1.0, 1.0])), [2.0, 4.0, 0.0])
    np.testing.assert_allclose(spread(6.0, np.zeros(3), np.array([0.0, 0.0, 3.0])), [0.0, 0.0, 6.0])
    np.testing.assert_allclose(spread(6.0, np.zeros(2), np.zeros(2)), [3.0, 3.0])


# ---------- recharge mass, HRU by HRU ----------

class _Sw:
    def __init__(self, values):
        self.values = values

    def has(self, var):
        return var in self.values

    def names(self, var):
        return ["h1", "h2"]

    def get(self, var):
        return np.asarray(self.values[var], float)


class _T:
    def entry_nodes(self, package):
        return np.array([0, 1, 2, -1])  # the last an inactive column


class _Link:
    def __init__(self, amap, rate):
        self.map, self.rate, self.package = amap, rate, "rch"


class _Mf:
    cell_area = np.array([400.0, 400.0, 400.0])


def test_recharge_mass_is_placed_where_the_water_goes():
    """Two HRUs over three cells (h1 on cells 0 and 1, h2 on 1 and 2): each HRU's mass goes to the cells its water
    goes to, in the same shares, so a cell's concentration is the mix of the HRUs' water and mass it gets."""
    amap = AreaMap(["h1", "h2"], np.array([0, 0, 1, 1]), np.array([0, 1, 1, 2]),
                   np.array([300.0, 100.0, 200.0, 400.0])).conserving(np.array([500.0, 600.0]))
    perc = np.array([0.002, 0.004])          # m d-1
    leached = np.array([0.5, 0.0])           # kg ha-1: h2 leaches none
    sw = _Sw({"hru_pw:percn": leached, "soil_water__percolation_rate": perc})
    sm = SurfaceMass(sw, _T(), "no3", recharge=True, stream=False, soil=False, source_area=np.array([500.0, 600.0]))
    c, made, sent, dry = sm.recharge_mass(_Link(amap, perc), _Mf())
    assert made == pytest.approx(0.5 / 1e4 * 500.0, rel=1e-14) and sent == pytest.approx(made) and dry == 0.0
    water = amap.to_cells(perc, _Mf.cell_area) * _Mf.cell_area
    h1_water = amap.to_cells(np.array([0.002, 0.0]), _Mf.cell_area) * _Mf.cell_area
    mass = c[:3] * water
    # the mass goes where h1's water goes, in its shares
    np.testing.assert_allclose(mass / mass.sum(), h1_water / h1_water.sum(), rtol=1e-12)
    # cell 0: h1 only, so h1's own concentration; cell 1: h1's mass in both HRUs' water; cell 2: none
    assert c[0] == pytest.approx(0.5 / 1e4 / 0.002, rel=1e-12)
    assert 0 < c[1] < c[0] and c[2] == 0.0 and c[3] == 0.0


def test_soluble_p_has_no_aquifer_recharge_path():
    sw = _Sw({"aquifer_water_solp__recharge_mass_rate": [0.0]})
    with pytest.raises(ValueError, match="source = \"hru\""):
        SurfaceMass(sw, _T(), "solp", recharge=True, stream=False, soil=False, source="aquifer")


def test_aquifer_recharge_mass_says_so_loudly():
    from basinledger.coupling.surface_mass import AquiferRechargeMass
    sw = _Sw({"aquifer_water_no3__recharge_mass_rate": [0.0]})
    with pytest.warns(AquiferRechargeMass, match="NOT HRU by HRU"):
        SurfaceMass(sw, _T(), "no3", recharge=True, stream=False, soil=False, source="aquifer",
                    source_area=np.array([1.0]))


# ---------- the steady state ----------

def steady_model(ws, days, decay=None, dispersivity=2.0, scheme="tvd"):
    import flopy

    from test_transport import flow_model
    sim, gwf = flow_model(ws, days, steady=True, nrow=10, ncol=10, cell=20.0, top=13.0, botm=10.0, k=20.0)
    flopy.mf6.ModflowGwfic(gwf, strt=12.0)
    flopy.mf6.ModflowGwfrcha(gwf, recharge=0.005, pname="rch")
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data={0: [((0, 5, 7), -2.0)]}, pname="wel")
    flopy.mf6.ModflowGwfriv(gwf, stress_period_data={0: [((0, r, 0), 12.0, 50.0, 11.0) for r in range(10)]},
                            pname="riv")
    add_transport(sim, gwf, "no3", 0.3, dispersivity, strt=0.0, sources=[("rch", "no3")], decay=decay,
                  scheme=scheme)
    sim.write_simulation(silent=True)
    return Modflow6(ws)


@pytest.mark.parametrize("decay", [None, 0.002])
def test_steady_state_is_the_long_run(libmf6, decay):
    """Steady flow, recharge at 0.025 kg m-3 leaving through a gaining river and a well (test_transport's steady
    case). Without decay every cell's steady state is the recharge's concentration; with decay, and with no
    dispersion and upstream advection (the steady state's own assumptions), it is where the transport model itself
    settles after 22 turnovers."""
    c0, days = 0.025, 4000
    mf = steady_model(RUNS / f"steady-{decay}", days, decay=decay, dispersivity=0.0, scheme="upstream")
    t = Transport(mf, "no3")
    try:
        for d in range(days):
            mf.prepare_day()
            t.set_inflow("rch", c0)
            mf.solve_day()
            b = t.budget("rch", {"w": ("wel", 0)}, [("riv", np.zeros(10, int), ["r"])])
            assert abs(b.closure_share) < 1e-8  # the solve's residual criterion, 1e-8 kg d-1, on 10 kg d-1 gross
            if d < 30:
                t.track()  # the steady state needs only the mean flows, which are steady
        x = t.concentration()
        s = t.steady_state()
    finally:
        mf.finalize()
    if decay is None:
        np.testing.assert_allclose(s, c0, rtol=1e-12)
    else:
        assert s.min() < 0.9 * c0  # decay lowers it downstream
    np.testing.assert_allclose(s, x, rtol=1e-6)


# ---------- the coupled synthetic run (SWAT+ two-way, nitrate and soluble P) ----------

NITRATE = ROOT / "examples" / "synthetic" / "nitrate.toml"
SPECIES = ("nitrate", "phosphorus")


def nitrate_config():
    from basinledger import config as C
    cfg = C.load(NITRATE)
    cfg.run_folder = RUNS
    return cfg


def run(cfg, name, days, **kw):
    from basinledger.driver import Basin
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run=name, days=days, **kw)
    return b


@pytest.fixture(scope="module")
def coupled(swatplus):
    b = run(nitrate_config(), "nitrate", 730)
    b.run()
    t = b.ledger.table()
    t["series:flow:cha01"] = np.asarray(b.series["flow:cha01"])
    b.finalize()
    return t


@pytest.mark.slow
def test_coupled_mass_closes_every_day(coupled):
    t = coupled
    for sp in SPECIES:
        k = f"transport:{sp}:"
        assert (np.abs(t[k + "closure_share"]) < 1e-9).all()
        # the aquifer's mass over the run, from the links' accounts: in, less out, plus storage, less decay, is zero
        total = (t[k + "recharge:received"] + t[k + "stream:cha01:returned_received"] - t[k + "stream:cha01:to_stream"]
                 - t[k + "well:pw1:withdrawn"] - t[k + "gw_soil:sent"] + t[k + "storage_release"]
                 - np.nan_to_num(t.get(k + "decay", np.zeros(1)))).sum()
        assert abs(total) < 1e-9 * t[k + "recharge:received"].sum()
        assert not any(a.startswith(k + "boundary:") for a in t)
        # every link: what the aquifer sent a channel is what SWAT+ kept plus what came back, the next day
        np.testing.assert_allclose(t[k + "stream:cha01:due"][1:], t[k + "stream:cha01:to_stream"][:-1], rtol=0)
        np.testing.assert_allclose(t[k + "stream:cha01:due"],
                                   t[k + "stream:cha01:in_swat"] + t[k + "stream:cha01:returned"], rtol=1e-9, atol=1e-12)
        np.testing.assert_allclose(t[k + "stream:cha01:returned"], t[k + "stream:cha01:returned_received"],
                                   rtol=1e-9, atol=1e-12)
        assert t[k + "recharge:received"].sum() > 0 and (t[k + "gw_soil:sent"] > 0).sum() > 100
    assert t["transport:nitrate:decay"].sum() > 0
    assert (t["transport:phosphorus:sorbed_release"] != 0).all()
    for x in ("no3", "solp"):  # SWAT+'s aquifers pass their recharge mass through, as their water
        assert not np.abs(t[f"swatplus:aquifer_{x}_outflow"]).any()


@pytest.mark.slow
def test_recharge_carries_swats_leached_mass(coupled):
    t = coupled
    for sp in SPECIES:
        k = f"transport:{sp}:"
        np.testing.assert_allclose(t[k + "recharge:received"], t[k + "recharge:made"], rtol=1e-9, atol=1e-12)
        assert not t[k + "recharge:no_water"].any()
    # nitrate: the leached mass HRU by HRU is SWAT+'s own aquifer recharge mass, day by day
    np.testing.assert_allclose(t["transport:nitrate:recharge:made"], t["transport:nitrate:recharge:swat_aquifers"],
                               rtol=1e-5, atol=1e-9)


@pytest.mark.slow
def test_a_gains_mass_arrives_in_the_channel(coupled):
    t = coupled
    gain = t["exchange:cha01:applied"] > 0
    assert gain.sum() > 100
    for sp in SPECIES:
        k = f"transport:{sp}:"
        mass_set = np.nan_to_num(t[k + "stream:cha01:mass_set"])
        np.testing.assert_allclose(t[k + "stream:cha01:in_swat"][gain], mass_set[gain], rtol=1e-6, atol=1e-12)
        assert (mass_set[gain] > 0).sum() > 100


@pytest.mark.slow
def test_a_loss_takes_the_channels_concentration(coupled):
    """On a day the channel loses (its net exchange), SWAT+'s loss takes the outflow's concentration for that day,
    which the aquifer gets on its loss cells the same day, through the mass sources. A loss that takes the whole of
    the day's outflow leaves no outflow to have a concentration, and takes all its mass."""
    t = coupled
    loss = t["exchange:cha01:applied"] < 0
    flowing = loss & (t["series:flow:cha01"] > 0)
    assert flowing.sum() > 10 and (loss & ~flowing).sum() > 10
    for sp in SPECIES:
        k = f"transport:{sp}:"
        np.testing.assert_allclose(t[k + "stream:cha01:in_swat"][flowing],
                                   t[k + "stream:cha01:concentration"][flowing] * t["exchange:cha01:applied"][flowing],
                                   rtol=1e-5, atol=1e-12)
        assert (t[k + "stream:cha01:in_swat"][loss] <= 0).all()
        assert (t[k + "stream:cha01:returned"][loss] > 0).all()


@pytest.mark.slow
def test_the_soil_link_carries_the_cells_mass(coupled):
    t = coupled
    for sp in SPECIES:
        k = f"transport:{sp}:"
        np.testing.assert_allclose(t[k + "gw_soil:due"][1:], t[k + "gw_soil:split"][:-1], rtol=1e-6, atol=1e-12)
        np.testing.assert_allclose(t[k + "gw_soil:received"], t[k + "gw_soil:due"], rtol=1e-6, atol=1e-12)


@pytest.mark.slow
def test_restart_carries_the_mass_links(swatplus):
    n, k = 160, 70
    a = run(nitrate_config(), "nitrate-unbroken", n)
    a.run()
    ta = a.ledger.table()
    a.finalize()
    w = run(nitrate_config(), "nitrate-warm", n)
    w.run(k)
    w.save(RUNS / "nitrate-state")
    w.finalize()
    b = run(nitrate_config(), "nitrate-restarted", n, restart=RUNS / "nitrate-state")
    b.run()
    tb = b.ledger.table()
    b.finalize()
    for sp in SPECIES:
        for acct in ("well:pw1:concentration", "stream:cha01:in_swat", "stream:cha01:returned", "gw_soil:received",
                     "recharge:received", "storage_release"):
            np.testing.assert_array_equal(tb[f"transport:{sp}:{acct}"], ta[f"transport:{sp}:{acct}"])


@pytest.mark.slow
def test_steady_start(swatplus):
    from basinledger.driver import Basin
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x, mean = Basin.transport_steady_state(nitrate_config(), days=120, run="nitrate-steady")
    assert set(x) == set(SPECIES) and all(np.isfinite(v).all() and (v >= 0).all() for v in x.values())
    assert mean["nitrate"] > 0
    b = run(nitrate_config(), "nitrate-from-steady", 5, concentrations=x)
    np.testing.assert_array_equal(b.transport["nitrate"].concentration(), x["nitrate"])
    b.run()
    assert (np.abs(b.ledger.table()["transport:nitrate:closure_share"]) < 1e-9).all()
    b.finalize()


def test_config_checks(tmp_path):
    from basinledger import config as C
    text = NITRATE.read_text()
    bad = tmp_path / "bad.toml"
    bad.write_text(text.replace('surface = "solp"', 'surface = "no3"'))
    with pytest.raises(ValueError, match="two species carry"):
        C.load(bad)
    bad.write_text(text.replace('surface = "no3"', 'surface = "nh4"'))
    with pytest.raises(ValueError, match="surface is SWAT\\+'s species"):
        C.load(bad)
    bad.write_text(text.replace('[transport.nitrate]\nsurface = "no3"\nrecharge = "swatplus"',
                                '[transport.nitrate]\nsurface = "no3"\nrecharge = "swat"'))
    with pytest.raises(ValueError, match="a concentration \\(kg m-3\\) or \"swatplus\""):
        C.load(bad)


@pytest.mark.slow
def test_check_finds_a_missing_mass_source(swatplus):
    from basinledger.driver import Basin
    problems, _ = Basin.check(nitrate_config(), run="nitrate-check")
    assert problems == []
    cfg = nitrate_config()
    tr = {k: dict(v, stream_side=False) for k, v in cfg.mf6_builder_args["transport"].items()}
    cfg.mf6_builder_args = dict(cfg.mf6_builder_args, transport=tr)
    problems, _ = Basin.check(cfg, run="nitrate-check-bad")
    assert len(problems) == 1 and "stream_mass" in problems[0] and "['riv', 'drn']" in problems[0]


@pytest.mark.slow
def test_a_lake_at_a_fixed_stage_is_its_outlet_channels(swatplus):
    """[exchange] lakes with the stream side: the lake's cells exchange with its outlet channel, so its gains carry the
    cells' concentrations into the channel and its losses the channel's mass back, through the lake's own mass sources;
    the mass still closes every day and no lake cell is a boundary."""
    cfg = nitrate_config()
    cfg.mf6_builder_args = dict(cfg.mf6_builder_args, lake_stage=88.5)
    cfg.exchange = dict(cfg.exchange, lakes={"package": "lakes", "map": "maps/lakes.csv"})
    b = run(cfg, "nitrate-lakes", 365)
    b.run()
    t = b.ledger.table()
    src = b.transport["nitrate"].last["SRC_LAKES"][1]
    b.finalize()
    assert (t["exchange:cha01:lakes"] < 0).sum() > 100  # the lake loses to the aquifer
    for sp in SPECIES:
        k = f"transport:{sp}:"
        assert (np.abs(t[k + "closure_share"]) < 1e-9).all()
        assert not any(a.startswith(k + "boundary:") for a in t)
        np.testing.assert_allclose(t[k + "stream:cha01:due"],
                                   t[k + "stream:cha01:in_swat"] + t[k + "stream:cha01:returned"], rtol=1e-9, atol=1e-12)
    assert src.sum() > 0  # the lake's losses took their mass back onto its cells
