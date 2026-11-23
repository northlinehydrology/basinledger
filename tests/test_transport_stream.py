"""Groundwater transport, stage 2: decay, sorbed storage in the ledger, the stream side's mass sources and the
steady-state start (coupling/transport.py).

- decay (first order, dissolved and sorbed) is booked, and the mass still closes every day (tier 1);
- sorbed mass is storage in the ledger: the change in the dissolved and sorbed mass the cells hold is the storage
  release, every day (tier 1);
- a stream's mass source is booked to its channel (tier 1);
- the steady state: uniform recharge gives its concentration everywhere, and with decay and no dispersion it is the
  transport model's own long-run answer (tier 3).
"""

from __future__ import annotations

import datetime as dt
import importlib.util
import warnings

import numpy as np
import pytest

from basinledger.coupling.modflow6 import Modflow6
from basinledger.coupling.transport import Transport, add_transport
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
