"""MODFLOW 6 through libmf6: values set between prepare and solve are the ones it uses."""

from __future__ import annotations

import datetime as dt
import importlib.util

import numpy as np
import pytest

from basinledger.coupling.modflow6 import Modflow6
from conftest import ROOT, RUNS


def build(ws, days):
    spec = importlib.util.spec_from_file_location("synthetic_model", ROOT / "examples" / "synthetic" / "model.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.build(ws, days, dt.date(2010, 1, 1))
    return m


def test_recharge_and_well_set_each_day(libmf6):
    ws = RUNS / "mf6-alone"
    model = build(ws, 5)
    mf = Modflow6(ws)
    well = mf.node(*model.WELL)
    for d in range(5):
        mf.prepare_day()
        rate = np.full(mf.nodes, 0.001 * (d + 1))  # m d-1
        mf.set_recharge("rch", rate)
        mf.set_well("wel", mf.entry("wel", well), -300.0 * (d + 1))
        mf.solve_day()
        assert mf.simvals("rch").sum() == pytest.approx(0.001 * (d + 1) * mf.cell_area.sum(), rel=1e-12)
        assert mf.simvals("wel")[0] == pytest.approx(-300.0 * (d + 1), rel=1e-12)
    h = mf.head()
    assert np.isfinite(h).all() and h.argmin() == well
    assert mf.time == pytest.approx(5.0)
    mf.finalize()


def test_a_solve_that_breaks_down_stops_at_once(libmf6):
    """MODFLOW 6 can return NaN heads from a linear solver that broke down and still report convergence (found in a
    Rappahannock calibration, Nov 22, 2026): the day stops there, with and without a callback, naming the iteration."""
    ws = RUNS / "mf6-nan"
    build(ws, 3)
    mf = Modflow6(ws)
    real = mf.head
    mf.prepare_day()
    mf.solve_day()  # a good day first
    mf.head = lambda: np.full_like(real(), np.nan)
    mf.prepare_day()
    with pytest.raises(FloatingPointError, match="day 1 .the day's solve"):
        mf.solve_day()
    mf.finalize()
    ws2 = RUNS / "mf6-nan-callback"
    build(ws2, 3)
    mf = Modflow6(ws2)
    mf.head = lambda: np.full(mf.nodes, np.nan)
    mf.prepare_day()
    with pytest.raises(FloatingPointError, match="outer iteration 1"):
        mf.solve_day(lambda: False)
    mf.finalize()


def test_runs_follow_one_another_in_a_process(libmf6):
    ws = RUNS / "mf6-again"
    build(ws, 2)
    for _ in range(2):
        mf = Modflow6(ws)
        mf.prepare_day()
        mf.solve_day()
        mf.finalize()


def test_model_must_step_daily(libmf6):
    ws = RUNS / "mf6-weekly"
    build(ws, 14)
    import flopy
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws, verbosity_level=0)
    sim.tdis.perioddata = [(14.0, 2, 1.0)]
    sim.write_simulation(silent=True)
    mf = Modflow6(ws)
    with pytest.raises(ValueError, match="one day at a time"):
        mf.prepare_day()
    mf.finalize()
