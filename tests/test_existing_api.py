"""What MODFLOW 6's API lets a coupler set each step in LAK, MAW and SFR (docs/design.md, "Extending an existing
MODFLOW 6 model", the open questions answered). Each answer the design rests on is a test here, so a MODFLOW 6 that
changes one fails here first. A small two-layer model: a lake, a multi-aquifer well over both layers and three SFR
reaches in a chain; a few one-day steps (tier 1).
"""

from __future__ import annotations

import flopy
import numpy as np
import pytest
from modflowapi import ModflowApi

from basinledger.coupling.modflow6 import find_library
from conftest import RUNS


def _build(ws, days=3):
    nlay, nrow, ncol = 2, 10, 10
    sim = flopy.mf6.MFSimulation(sim_name="p", sim_ws=str(ws), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, time_units="days", nper=1, perioddata=[(float(days), days, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="moderate", outer_maximum=200, inner_maximum=200)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    idomain = np.ones((nlay, nrow, ncol), int)
    idomain[0, 4:6, 4:6] = 0  # the lake's cells
    flopy.mf6.ModflowGwfdis(gwf, length_units="meters", nlay=nlay, nrow=nrow, ncol=ncol, delr=100., delc=100.,
                            top=20., botm=[0., -20.], idomain=idomain)
    flopy.mf6.ModflowGwfnpf(gwf, k=10., icelltype=1)
    flopy.mf6.ModflowGwfsto(gwf, sy=0.2, ss=1e-5, iconvert=1, transient={0: True})
    flopy.mf6.ModflowGwfic(gwf, strt=15.)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[((0, r, 0), 15.) for r in range(nrow)], pname="edge")
    conn = [(0, k, (1, r, c), "vertical", 0.1, 0., 0., 0., 0.)
            for k, (r, c) in enumerate([(4, 4), (4, 5), (5, 4), (5, 5)])]
    flopy.mf6.ModflowGwflak(gwf, pname="lake", nlakes=1, packagedata=[(0, 15., len(conn))], connectiondata=conn,
                            perioddata={0: [(0, "rainfall", 0.0), (0, "evaporation", 0.0), (0, "inflow", 0.0)]},
                            noutlets=0, length_conversion=1.0, time_conversion=86400.0)
    flopy.mf6.ModflowGwfmaw(gwf, pname="mw", nmawwells=1, packagedata=[(0, 0.1, -20., 15., "thiem", 2)],
                            connectiondata=[(0, 0, (0, 8, 8), 20., 0., 0., 0.), (0, 1, (1, 8, 8), 0., -20., 0., 0.)],
                            perioddata={0: [(0, "rate", 0.0)]})
    pd = [(i, (0, 2, 3 + i), 100., 3., 1e-3, 14., 1., 0.5, 0.035, 1 if i in (0, 2) else 2, 1.0, 0) for i in range(3)]
    flopy.mf6.ModflowGwfsfr(gwf, pname="river", nreaches=3, packagedata=pd,
                            connectiondata=[(0, -1), (1, 0, -2), (2, 1)],
                            perioddata={0: [(0, "inflow", 0.0)]}, unit_conversion=86400.)
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="gwf.cbc", printrecord=[("budget", "last")])
    sim.write_simulation(silent=True)


class Api:
    """libmf6 on the test model, one day at a time: set inputs after prepare, read after the solve."""

    def __init__(self, ws):
        self.ws = ws
        self.closed = False
        self.api = ModflowApi(str(find_library()), working_directory=str(ws))
        self.api.initialize()

    def a(self, var, pkg):
        return self.api.get_var_address(var, "GWF", pkg)

    def get(self, var, pkg):
        return self.api.get_value(self.a(var, pkg)).copy()

    def set(self, var, pkg, idx, val):
        x = self.get(var, pkg)
        x[idx] = val
        self.api.set_value(self.a(var, pkg), x)

    def day(self, before):
        self.api.prepare_time_step(1.0)
        before(self)
        self.api.prepare_solve(1)
        for _ in range(500):
            if self.api.solve(1):
                break
        self.api.finalize_solve(1)
        out = {k: self.get(*k) for k in [("XNEWPAK", "LAKE"), ("SIMVALS", "LAKE"), ("RATESIM", "MW"),
                                         ("SIMVALS", "MW"), ("STAGE", "RIVER"), ("SIMVALS", "RIVER"),
                                         ("QOUTFLOW", "RIVER"), ("QEXTOUTFLOW", "RIVER"), ("INFLOW", "RIVER")]}
        self.api.finalize_time_step()
        return out

    def close(self):
        if not self.closed:
            self.api.finalize()
            self.closed = True


@pytest.fixture
def model(libmf6, request):
    ws = RUNS / f"existing-api-{request.node.name}"
    _build(ws)
    m = Api(ws)
    yield m
    m.close()


def test_lake_rain_evaporation_runoff_inflow_and_withdrawal_set_each_day(model):
    def lake(m):
        m.set("RAINFALL", "LAKE", 0, 0.01)       # m d-1 over the lake's surface
        m.set("EVAPORATION", "LAKE", 0, 0.002)
        m.set("RUNOFF", "LAKE", 0, 300.)         # m3 d-1
        m.set("INFLOW", "LAKE", 0, 2000.)
        m.set("WITHDRAWAL", "LAKE", 0, 100.)
    stages = [float(model.day(lake)[("XNEWPAK", "LAKE")][0]) for _ in range(3)]
    # the lake's area is four 100 m cells: rain 400 m3, evaporation 80 m3 a day, at a stage above the cells' top
    assert model.get("PRECIP", "LAKE")[0] == pytest.approx(400.0, rel=1e-9)
    assert model.get("EVAP", "LAKE")[0] == pytest.approx(-80.0, rel=1e-9)  # out of the lake
    assert np.all(np.diff([15.0] + stages) > 0.05)  # 2,620 m3 a day in, a little out to the aquifer: it rises


def test_lake_left_alone_barely_moves(model):
    s = [float(model.day(lambda m: None)[("XNEWPAK", "LAKE")][0]) for _ in range(3)]
    assert abs(s[-1] - 15.0) < 1e-2


def test_maw_rate_is_split_among_its_connections(model):
    for d in range(1, 4):
        out = model.day(lambda m, d=d: m.set("RATE", "MW", 0, -200. * d))
        assert out[("RATESIM", "MW")][0] == pytest.approx(-200. * d)
        q = out[("SIMVALS", "MW")]
        assert len(q) == 2 and np.all(q < 0)
        # the connections add to the rate less what the wellbore's water level gave up: small, not zero
        assert q.sum() == pytest.approx(-200. * d, rel=1e-3)


def test_sfr_inflow_routes_down_the_reaches(model):
    out = model.day(lambda m: m.set("INFLOW", "RIVER", 0, 1000.))
    leak = out[("SIMVALS", "RIVER")]  # + into the aquifer
    q = out[("QOUTFLOW", "RIVER")]    # negative: leaving the reach
    assert -q[0] == pytest.approx(1000. - leak[0])
    assert -q[1] == pytest.approx(-q[0] - leak[1])
    # the last reach's outflow leaves the network
    assert -out[("QEXTOUTFLOW", "RIVER")][2] == pytest.approx(-q[1] - leak[2])


def test_sfr_connection_cut_at_run_time_books_the_outflow_as_leaving(model):
    def cut(m):
        m.set("INFLOW", "RIVER", 0, 1000.)
        m.set("INFLOW", "RIVER", 1, 500.)
        m.set("USTRF", "RIVER", 1, 0.0)  # reach 1 takes nothing from reach 0
    for _ in range(3):  # the period's three days: MODFLOW 6 prints SFR's budget at the period's end
        out = model.day(cut)
    leak, q = out[("SIMVALS", "RIVER")], out[("QOUTFLOW", "RIVER")]
    assert -q[1] == pytest.approx(500. - leak[1])  # only its own inflow
    model.close()
    lst = (model.ws / "gwf.lst").read_text()
    block = lst[lst.rindex("RIVER BUDGET FOR ENTIRE MODEL"):]
    rates = lambda name: [float(x.split("=")[-1]) for x in block.splitlines() if x.strip().startswith(name)]  # noqa: E731
    ext_in, ext_out = rates("EXT-INFLOW")[0], rates("EXT-OUTFLOW")[1]  # the day's rates: in, then out
    # reach 0's outflow leaves the network there, with the last reach's: the 1,500 m3 set, less the net leakage
    assert ext_in == pytest.approx(1500.0, rel=1e-6)
    assert ext_out == pytest.approx(1500.0 - leak.sum(), rel=1e-6)
    assert ext_out > -q[1]  # more than the last reach could carry: reach 0's water is in it
    disc = [float(x.split("=")[1].split()[0]) for x in block.splitlines() if "PERCENT DISCREPANCY" in x]
    assert abs(disc[0]) < 0.01


def test_sfr_simple_status_takes_the_stage_given_and_holds_its_loss_to_its_inflow(model):
    def simple(m, stage, inflow):
        ib = m.get("IBOUNDPAK", "RIVER")
        ib[:] = -1
        m.api.set_value(m.a("IBOUNDPAK", "RIVER"), ib)
        m.set("SSTAGE", "RIVER", slice(None), stage)
        m.set("INFLOW", "RIVER", 0, inflow)
    out = model.day(lambda m: simple(m, 14.3, 1000.))
    assert np.allclose(out[("STAGE", "RIVER")], 14.3)  # below the 15 m head: the reaches gain
    assert np.all(out[("SIMVALS", "RIVER")] < 0)
    out = model.day(lambda m: simple(m, 18.0, 100.))   # above the head: it would lose far more than 100 m3
    assert out[("SIMVALS", "RIVER")][0] == pytest.approx(100.0)
    assert np.allclose(out[("SIMVALS", "RIVER")][1:], 0.0)
    assert out[("STAGE", "RIVER")][0] == pytest.approx(14.0)  # put back to the streambed's top


def test_lake_at_a_constant_stage_set_each_day(model):
    """The lake as a fixed-stage face (status constant, -1): its stage is the one set, and its leakage is computed
    at that stage; the lake holds no water of its own from day to day."""
    for stage in (15.5, 16.0):
        def constant(m, stage=stage):
            m.set("IBOUND", "LAKE", 0, -1)
            m.set("STAGE", "LAKE", 0, stage)
        out = model.day(constant)
        assert out[("XNEWPAK", "LAKE")][0] == pytest.approx(stage)
        assert out[("SIMVALS", "LAKE")].sum() > 0  # above the 15 m heads: it loses to the aquifer
    higher = out[("SIMVALS", "LAKE")].sum()
    out = model.day(lambda m: (m.set("IBOUND", "LAKE", 0, -1), m.set("STAGE", "LAKE", 0, 15.5)))
    assert out[("SIMVALS", "LAKE")].sum() < higher


def test_lake_started_each_day_from_a_stage_given(model):
    """The lake active, its day started from a stage given (XOLDPAK, the stage at the step's start): it takes its
    inflows and leaks within the day from there."""
    for _ in range(2):
        def start(m):
            m.set("XOLDPAK", "LAKE", 0, 15.2)
            m.set("XNEWPAK", "LAKE", 0, 15.2)
            m.set("INFLOW", "LAKE", 0, 1000.)
        out = model.day(start)
        # 1,000 m3 into 4 ha of lake above the cells' top: about 2.5 cm, less what it leaks
        assert 15.2 < out[("XNEWPAK", "LAKE")][0] < 15.2 + 1000. / 40000. + 1e-9
