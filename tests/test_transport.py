"""Groundwater transport (coupling/transport.py): MODFLOW 6 GWT through Basin Ledger's coupling, its mass ledger, and
the intake it gives a well's user.

- the builder: the transport model, its own solution, the exchange, and the auxiliary variables it adds;
- mass closes every day, with recharge, a well, a river that gains and loses, drains and sorption (tier 1 and 2);
- a steady flow under recharge of one concentration ends with that concentration everywhere (tier 3);
- a one-dimensional column against Ogata and Banks (1961) (tier 3);
- the coupled synthetic run (examples/synthetic/transport.toml): the ledger, and the user's intake following the
  aquifer's concentration at its well a day later (tier 2).
"""

from __future__ import annotations

import datetime as dt
import importlib.util

import numpy as np
import pytest

from basinledger.analytic import ogata_banks
from basinledger.coupling.modflow6 import Modflow6
from basinledger.coupling.transport import Transport, add_transport
from conftest import ROOT, RUNS


def synthetic_model():
    spec = importlib.util.spec_from_file_location("synthetic_model", ROOT / "examples" / "synthetic" / "model.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def flow_model(ws, days, steady=False, nrow=5, ncol=5, cell=10.0, top=100.0, botm=70.0, k=2.0):
    import flopy
    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=str(ws), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, time_units="days", nper=1, perioddata=[(float(days), days, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="simple", outer_dvclose=1e-9, inner_dvclose=1e-11,
                         rcloserecord=[1e-9, "strict"])
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf")
    flopy.mf6.ModflowGwfdis(gwf, length_units="meters", nlay=1, nrow=nrow, ncol=ncol, delr=cell, delc=cell, top=top,
                            botm=botm)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k)
    if not steady:
        flopy.mf6.ModflowGwfsto(gwf, iconvert=0, sy=0.15, ss=1e-5, transient={0: True})
    return sim, gwf


# ---------- the builder (no MODFLOW 6 run) ----------

def test_builder_adds_the_model_its_solution_and_the_auxiliaries(tmp_path):
    import flopy
    sim, gwf = flow_model(tmp_path, 5)
    flopy.mf6.ModflowGwfic(gwf, strt=88.0)
    flopy.mf6.ModflowGwfrcha(gwf, recharge=0.001, pname="rch")
    flopy.mf6.ModflowGwfriv(gwf, stress_period_data={0: [((0, r, 0), 88.0, 20.0, 86.0) for r in range(5)]},
                            pname="riv")
    flopy.mf6.ModflowGwfghb(gwf, stress_period_data={0: [((0, 2, 4), 88.0, 5.0, 0.4)]}, pname="ghb",
                            auxiliary=["salt"])
    add_transport(sim, gwf, "chloride", 0.3, 2.0, strt=0.01, sources=[("rch", "chloride"), ("riv", "chloride"),
                                                                     ("ghb", "chloride")])
    add_transport(sim, gwf, "phosphorus", 0.3, 2.0, sorption={"distcoef": 1e-3, "bulk_density": 1800.0},
                  sources=[("rch", "phosphorus")])
    sim.write_simulation(silent=True)
    nam = (tmp_path / "mfsim.nam").read_text()
    for line in ("gwt6  gwt_chloride.nam  gwt_chloride", "GWF6-GWT6  gwf_gwt_chloride.gwfgwt  gwf  gwt_chloride",
                 "ims6  gwt_chloride.ims  gwt_chloride", "ims6  gwt_phosphorus.ims  gwt_phosphorus"):
        assert line in nam
    assert "ims6  t.ims  gwf\n" in nam  # the flow's solution keeps only the flow model
    rch = (tmp_path / "gwf.rcha").read_text()
    assert "auxiliary  chloride  phosphorus" in rch
    riv = (tmp_path / "gwf.riv").read_text()
    assert "auxiliary  chloride" in riv and riv.count("0.00000000E+00") == 5  # a zero column for each entry
    ghb = (tmp_path / "gwf.ghb").read_text()
    assert "auxiliary  salt  chloride" in ghb and "4.00000000E-01 0.00000000E+00" in ghb  # the old column kept
    mst = (tmp_path / "gwt_phosphorus.mst").read_text()
    assert "sorption  linear" in mst.lower() and "distcoef" in mst.lower()
    dsp = (tmp_path / "gwt_chloride.dsp").read_text().lower()
    assert "ath1" in dsp and "0.20000000" in dsp  # a tenth of the longitudinal [A]


def test_builder_refuses_an_unknown_package(tmp_path):
    sim, gwf = flow_model(tmp_path, 5)
    with pytest.raises(KeyError, match="no package 'rch'"):
        add_transport(sim, gwf, "chloride", 0.3, 1.0, sources=[("rch", "chloride")])


# ---------- MODFLOW 6 with transport, through the coupling ----------

def run_small(ws, days, sorption=None, pump=-900.0):
    """The synthetic two-way grid with chloride (examples/synthetic/model.py, build_small, with drains), driven as
    the driver does: varying recharge, a well that pumps hard enough on alternate weeks that the river loses."""
    m = synthetic_model()
    m.build_small(ws, days, dt.date(2010, 1, 1), drain_elevation=89.0,
                  transport={"chloride": dict(porosity=0.25, dispersivity=5.0, strt=0.01, sorption=sorption)})
    mf = Modflow6(ws)
    t = Transport(mf, "chloride")
    well = mf.node(*m.SMALL_WELL)
    river = (np.zeros(10, int), ["cha01"])
    out = []
    for d in range(days):
        mf.prepare_day()
        mf.set_recharge("rch", np.full(mf.nodes, 0.002 + 0.002 * np.sin(d / 5.0)))
        mf.set_well("wel", 0, pump if d % 14 < 7 else 0.0)
        t.set_inflow("rch", 0.03)
        t.set_inflow("riv", 0.05)
        mf.solve_day(lambda: False)  # the callback path, as the supply limit's
        b = t.budget("rch", {"pw1": ("wel", 0)}, [("riv", *river), ("drn", *river)])
        out.append((b, float(mf.simvals("rch").sum()), float(t.concentration()[well])))
    mf.finalize()
    return out


@pytest.mark.parametrize("sorption", [None, {"distcoef": 1e-4, "bulk_density": 1800.0}])
def test_mass_closes_every_day(libmf6, sorption):
    tag = "sorbed" if sorption else "conservative"
    days = run_small(RUNS / f"transport-closure-{tag}", 90, sorption)
    for b, rch_m3, _ in days:
        assert abs(b.closure_share) < 1e-9
        assert b.recharge == pytest.approx(0.03 * rch_m3, rel=1e-12)              # what was set is what it took in
        assert b.from_stream["cha01"] == pytest.approx(0.05 * b.leakage_m3["cha01"], rel=1e-12, abs=1e-15)
    budgets = [b for b, _, _ in days]
    assert sum(b.leakage_m3["cha01"] > 0 for b in budgets) > 10                  # the river lost on some days
    assert sum(b.to_stream["cha01"] > 0 for b in budgets) > 10                   # and gained on others
    assert sum(b.wells["pw1"] for b in budgets) > 0
    assert not any(b.boundary for b in budgets)                                  # everything here is coupled


def test_the_callback_path_solves_the_transport(libmf6):
    """Modflow6.solve_day with a callback (the supply limit's path) solves the transport after the flow, as
    do_time_step does: both give the same concentrations."""
    m = synthetic_model()
    x = {}
    for cb in (False, True):
        ws = RUNS / f"transport-path-{cb}"
        m.build_small(ws, 10, dt.date(2010, 1, 1), transport={"chloride": dict(porosity=0.25, dispersivity=5.0)})
        mf = Modflow6(ws)
        t = Transport(mf, "chloride")
        for _ in range(10):
            mf.prepare_day()
            mf.set_recharge("rch", np.full(mf.nodes, 0.003))
            t.set_inflow("rch", 0.03)
            t.set_inflow("riv", 0.0)
            mf.solve_day((lambda: False) if cb else None)
        x[cb] = t.concentration()
        mf.finalize()
    assert x[False].max() > 0
    np.testing.assert_allclose(x[True], x[False], rtol=1e-9, atol=1e-15)


def test_uniform_recharge_concentration_reaches_everywhere(libmf6):
    """Steady flow, all of it recharge leaving through a gaining river and a well: the aquifer, starting free of the
    species, ends at the recharge's concentration in every cell. A thin aquifer (3 m, porosity 0.3) under 5 mm/d
    turns its water over in 180 days; 4,000 days is 22 turnovers."""
    import flopy
    ws = RUNS / "transport-steady"
    days, c0 = 4000, 0.025
    sim, gwf = flow_model(ws, days, steady=True, nrow=10, ncol=10, cell=20.0, top=13.0, botm=10.0, k=20.0)
    flopy.mf6.ModflowGwfic(gwf, strt=12.0)
    flopy.mf6.ModflowGwfrcha(gwf, recharge=0.005, pname="rch")
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data={0: [((0, 5, 7), -2.0)]}, pname="wel")
    flopy.mf6.ModflowGwfriv(gwf, stress_period_data={0: [((0, r, 0), 12.0, 50.0, 11.0) for r in range(10)]},
                            pname="riv")
    add_transport(sim, gwf, "chloride", 0.3, 2.0, strt=0.0, sources=[("rch", "chloride")])
    sim.write_simulation(silent=True)
    mf = Modflow6(ws)
    t = Transport(mf, "chloride")
    for _ in range(days):
        mf.prepare_day()
        t.set_inflow("rch", c0)
        mf.solve_day()
        assert (mf.simvals("riv") <= 0).all()  # the river only gains
        b = t.budget("rch", {"w": ("wel", 0)}, [("riv", np.zeros(10, int), ["r"])])
        assert abs(b.closure_share) < 1e-9
    x = t.concentration()
    mf.finalize()
    np.testing.assert_allclose(x, c0, rtol=1e-8)
    # and at steady state the mass in equals the mass out, at that concentration
    assert b.recharge == pytest.approx(b.wells["w"] + b.to_stream["r"], rel=1e-7)
    assert b.wells["w"] == pytest.approx(2.0 * c0, rel=1e-8)


# Ogata and Banks: a column of 1 m cells, a uniform pore velocity of 0.5 m/d, longitudinal dispersivity 10 m (D = 5
# m2/d), the first cell held at C0 = 1 (GWT's constant concentration) from the first day. Daily steps, as the driver
# takes. MODFLOW 6's implicit step adds a numerical dispersion of about v^2 dt / 2 = 0.125 m2/d (2.5% of D), and TVD
# advection keeps the spatial part small.
NCOL, V, ALPHA, THETA = 501, 0.5, 10.0, 0.3
OB_TOL = 0.01  # of C0, at every cell


def test_column_matches_ogata_banks(libmf6):
    import flopy
    ws = RUNS / "transport-ogata-banks"
    days = 300
    q = V * THETA                       # Darcy flux, m/d
    k = 10.0
    grad = q / k
    sim, gwf = flow_model(ws, days, steady=True, nrow=1, ncol=NCOL, cell=1.0, top=1.0, botm=0.0, k=k)
    h0 = 10.0
    flopy.mf6.ModflowGwfic(gwf, strt=h0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [((0, 0, 0), h0), ((0, 0, NCOL - 1),
                                                                            h0 - grad * (NCOL - 1))]}, pname="chd")
    gwt = add_transport(sim, gwf, "tracer", THETA, (ALPHA, 0.0, 0.0), strt=0.0, sources=[])
    flopy.mf6.ModflowGwtcnc(gwt, stress_period_data={0: [((0, 0, 0), 1.0)]}, pname="cnc")
    sim.write_simulation(silent=True)
    mf = Modflow6(ws)
    t = Transport(mf, "tracer")
    assert t.boundaries == ["CNC"]
    x = np.arange(NCOL, dtype=float)    # distance from the held cell's center, m
    errors = {}
    for d in range(1, days + 1):
        mf.prepare_day()
        mf.solve_day()
        b = t.budget()
        assert abs(b.closure_share) < 1e-9
        if d in (100, 200, 300):
            c = t.concentration()
            errors[d] = np.abs(c - ogata_banks(x, d, V, ALPHA * V)).max()
    mf.finalize()
    print({d: round(e, 4) for d, e in errors.items()})
    assert max(errors.values()) < OB_TOL


# ---------- the coupled synthetic run (SWAT+ two-way) ----------

TRANSPORT = ROOT / "examples" / "synthetic" / "transport.toml"


def transport_config():
    from basinledger import config as C
    cfg = C.load(TRANSPORT)
    cfg.run_folder = RUNS
    return cfg


@pytest.fixture(scope="module")
def coupled(swatplus):
    from basinledger.driver import Basin
    b = Basin(transport_config(), run="transport", days=730)
    b.run()
    t = b.ledger.table()
    b.finalize()
    return t


@pytest.mark.slow
def test_coupled_mass_ledger(coupled):
    t = coupled
    k = "transport:chloride:"
    assert (np.abs(t[k + "closure_share"]) < 1e-9).all()
    np.testing.assert_allclose(t[k + "recharge:received"], t[k + "recharge:sent"], rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(t[k + "recharge:sent"], 0.03 * t["recharge:sent"], rtol=1e-12)
    lost = t[k + "stream:cha01:leakage_received"]
    assert (lost > 0).sum() > 10 and (t[k + "stream:cha01:to_stream"] > 0).sum() > 100  # the river loses and gains
    # the run's mass, from the accounts: in, less out, plus the release from storage, is zero
    total = (t[k + "recharge:received"] + lost - t[k + "stream:cha01:to_stream"] - t[k + "well:pw1:withdrawn"]
             + t[k + "storage_release"]).sum()
    assert abs(total) < 1e-9 * t[k + "recharge:received"].sum()
    assert not any(a.startswith(k + "boundary:") for a in t)


@pytest.mark.slow
def test_intake_follows_the_aquifer_a_day_later(coupled):
    t = coupled
    c_well = t["transport:chloride:well:pw1:concentration"]
    intake = t["user:plant:intake:chloride"]
    assert intake[0] == pytest.approx(0.01, rel=1e-12)       # the starting concentration on the first day
    np.testing.assert_array_equal(intake[1:], c_well[:-1])    # then the aquifer's at the end of the day before
    assert c_well[-1] > c_well[0] * 1.05                       # recharge at 0.03 raises the aquifer's
    # the user takes in the well's water at that concentration, and returns all of it (the synthetic user)
    w = t["user:plant:withdrawal"]
    np.testing.assert_allclose(t["user:plant:intake:chloride:kg"], w * intake, rtol=1e-12)
    np.testing.assert_allclose(t["user:plant:chloride:kg"], t["user:plant:intake:chloride:kg"], rtol=1e-9)
    # what the well removed and what the user took in differ only by the day's lag
    removed, taken = t["transport:chloride:well:pw1:withdrawn"].sum(), t["user:plant:intake:chloride:kg"].sum()
    assert removed == pytest.approx(taken, rel=1e-3)
    np.testing.assert_allclose(t["transport:chloride:well:pw1:withdrawn"], t["route:well:pw1:delivered"] * c_well,
                               rtol=1e-9)


@pytest.mark.slow
def test_restart_carries_the_concentrations(swatplus):
    from basinledger.driver import Basin
    n, k = 200, 80
    a = Basin(transport_config(), run="transport-unbroken", days=n)
    a.run()
    ta = a.ledger.table()
    a.finalize()
    w = Basin(transport_config(), run="transport-warm", days=n)
    w.run(k)
    w.save(RUNS / "transport-state")
    w.finalize()
    b = Basin(transport_config(), run="transport-restarted", days=n, restart=RUNS / "transport-state")
    b.run()
    tb = b.ledger.table()
    b.finalize()
    for acct in ("transport:chloride:well:pw1:concentration", "user:plant:intake:chloride",
                 "transport:chloride:storage_release"):
        np.testing.assert_array_equal(tb[acct], ta[acct])


def test_settings_checked_at_start_up(tmp_path):
    from basinledger import config as C
    text = TRANSPORT.read_text()
    bad = tmp_path / "bad.toml"
    bad.write_text(text.replace("stream = 0.02             # kg m-3 [A]\n", ""))
    with pytest.raises(ValueError, match="needs stream"):
        C.load(bad)
    bad.write_text(text.replace("recharge = 0.03", "recharg = 0.03"))
    with pytest.raises(ValueError, match="needs recharge"):
        C.load(bad)


@pytest.mark.slow
def test_driver_refuses_what_it_cannot_book(swatplus):
    """A transport model with no [transport] section, and a constant intake the aquifer would replace."""
    from basinledger.driver import Basin
    cfg = transport_config()
    cfg.transport = {}
    with pytest.raises(ValueError, match="no \\[transport\\] section"):
        Basin(cfg, run="transport-unbooked", days=5)
    cfg = transport_config()
    cfg.users[0].inputs = {"water_user_intake_water_chloride__mass_concentration": 0.012}
    with pytest.raises(ValueError, match="set from the aquifer"):
        Basin(cfg, run="transport-constant-intake", days=5)
