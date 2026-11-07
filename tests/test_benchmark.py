"""Tier 3: MODFLOW 6, driven through Basin Ledger's coupling, against Glover and Balmer (1954).

A uniform aquifer (transmissivity 500 m2/d, storativity 0.1) with one straight river along its west edge, held at a
fixed stage with a streambed conductance (1e5 m2/d a cell, 200 times the transmissivity) high enough to be no
resistance, and a well 500 m from it pumping 1,000
m3/d from the first day. The other edges are 5 km away, beyond the response over the year. The river's loss of
discharge is the depletion; Glover and Balmer's solution assumes the same, but an infinite aquifer and a line river.
"""

from __future__ import annotations

import datetime as dt

import numpy as np
import pytest

from basinledger.analytic import glover_depletion

CELL, NROW, NCOL = 50.0, 101, 101
T_M2D, S, Q, DIST = 500.0, 0.1, 1000.0, 500.0
DAYS = 365


def build(ws, days):
    import flopy
    sim = flopy.mf6.MFSimulation(sim_name="glover", sim_ws=str(ws), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, time_units="days", start_date_time="2020-01-01",
                          nper=1, perioddata=[(float(days), days, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="moderate", outer_dvclose=1e-7, inner_dvclose=1e-9,
                         rcloserecord=[1e-3, "strict"], outer_maximum=100, inner_maximum=300)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf")
    b = 50.0
    flopy.mf6.ModflowGwfdis(gwf, length_units="meters", nlay=1, nrow=NROW, ncol=NCOL, delr=CELL, delc=CELL,
                            top=100.0, botm=100.0 - b)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=T_M2D / b)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=S / b, sy=S, transient={0: True})
    flopy.mf6.ModflowGwfic(gwf, strt=90.0)
    flopy.mf6.ModflowGwfriv(gwf, pname="riv", stress_period_data={0: [((0, r, 0), 90.0, 1e5, 80.0)
                                                                      for r in range(NROW)]})
    col = int(round(DIST / CELL))
    flopy.mf6.ModflowGwfwel(gwf, pname="wel", stress_period_data={0: [((0, NROW // 2, col), 0.0)]})
    flopy.mf6.ModflowGwfoc(gwf)
    sim.write_simulation(silent=True)
    return (0, NROW // 2, col)


def test_depletion_matches_glover(libmf6, tmp_path):
    from basinledger.coupling.modflow6 import Modflow6
    cell = build(tmp_path, DAYS)
    m = Modflow6(tmp_path)
    try:
        dep = []
        for d in range(DAYS):
            m.prepare_day()
            if d == 0:  # list packages are filled at the first prepare
                well = m.entry("wel", m.node(*cell))
            m.set_well("wel", well, -Q)
            m.solve_day()
            q = -m.simvals("riv").sum()  # discharge to the river, m3/d (it was zero before pumping)
            dep.append(-q)
    finally:
        m.finalize()
    dep = np.array(dep)
    glover = glover_depletion(np.full(DAYS, Q), DIST, T_M2D, S)
    # a 50 m grid and a river a column of cells wide: within 2% of Q from day 10 (0.5% when written, Nov 6, 2026)
    after = slice(10, None)
    assert np.max(np.abs(dep[after] - glover[after])) <= 0.02 * Q
    assert dep[-1] == pytest.approx(glover[-1], abs=0.005 * Q)
