"""The synthetic example's MODFLOW 6 model: a small aquifer drained by a river, with one pumping well.

One layer, 41 x 41 cells of 50 m (about 4.2 km2), 100 m thick. Transmissivity is held constant (ICELLTYPE 0) and
storage uses a specific yield of 0.15 below the top, a common simplification for a water-table aquifer that stays far
from drying. A river runs along the west edge (stage 88 m). Other edges are no-flow. Recharge and the well's rate are
zero in the files; the driver sets them every day.
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import flopy

NROW = NCOL = 41
CELL = 50.0
RIVER_STAGE = 88.0
WELL = (0, 20, 28)


def build(workspace: str | Path, days: int, start: dt.date) -> Path:
    ws = Path(workspace)
    sim = flopy.mf6.MFSimulation(sim_name="synthetic", sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, time_units="days", start_date_time=start.isoformat(),
                          nper=1, perioddata=[(float(days), days, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="simple", outer_dvclose=1e-8, inner_dvclose=1e-10,
                         rcloserecord=[1e-8, "strict"])
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    flopy.mf6.ModflowGwfdis(gwf, length_units="meters", nlay=1, nrow=NROW, ncol=NCOL, delr=CELL, delc=CELL,
                            top=100.0, botm=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=5.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, sy=0.15, ss=1e-5, transient={0: True})
    flopy.mf6.ModflowGwfic(gwf, strt=RIVER_STAGE + 0.5)
    flopy.mf6.ModflowGwfrcha(gwf, recharge=0.0, pname="rch")
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data={0: [(WELL, 0.0)]}, pname="wel")
    # conductance: K (1 m/d) x width (5 m) x length (50 m) / bed thickness (1 m)
    flopy.mf6.ModflowGwfriv(gwf, stress_period_data={0: [((0, r, 0), RIVER_STAGE, 250.0, RIVER_STAGE - 3.0)
                                                         for r in range(NROW)]}, pname="riv")
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="gwf.cbc", head_filerecord="gwf.hds",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")], printrecord=[("BUDGET", "ALL")])
    sim.write_simulation(silent=True)
    return ws


SMALL_CELL = 1000.0 ** 0.5  # m: 10 x 10 cells of 1,000 m2 make 10 ha, the area of Osu_1hru's aquifer
SMALL_WELL = (0, 5, 6)


def build_small(workspace: str | Path, days: int, start: dt.date, drain_elevation: float | None = None,
                transport: dict | None = None) -> Path:
    """A 10 ha aquifer under Osu_1hru for the two-way exchange: 10 x 10 cells, the river along the west edge.
    Thinner and less permeable than the large grid so the river can lose water when the well pumps hard.
    ``drain_elevation``: seepage cells (a drain) along the second column at that elevation, for the coupled drains'
    tests. ``transport``: species -> ``add_transport`` settings (porosity, dispersivity, strt, sorption), a transport
    model per species with recharge and the river as its sources."""
    ws = Path(workspace)
    sim = flopy.mf6.MFSimulation(sim_name="synthetic_small", sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, time_units="days", start_date_time=start.isoformat(),
                          nper=1, perioddata=[(float(days), days, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="simple", outer_dvclose=1e-8, inner_dvclose=1e-10,
                         rcloserecord=[1e-8, "strict"])
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    flopy.mf6.ModflowGwfdis(gwf, length_units="meters", nlay=1, nrow=10, ncol=10, delr=SMALL_CELL, delc=SMALL_CELL,
                            top=100.0, botm=70.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=2.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, sy=0.15, ss=1e-5, transient={0: True})
    flopy.mf6.ModflowGwfic(gwf, strt=RIVER_STAGE)
    flopy.mf6.ModflowGwfrcha(gwf, recharge=0.0, pname="rch")
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data={0: [(SMALL_WELL, 0.0)]}, pname="wel")
    flopy.mf6.ModflowGwfriv(gwf, stress_period_data={0: [((0, r, 0), RIVER_STAGE, 20.0, RIVER_STAGE - 2.0)
                                                         for r in range(10)]}, pname="riv")
    if drain_elevation is not None:
        flopy.mf6.ModflowGwfdrn(gwf, stress_period_data={0: [((0, r, 1), drain_elevation, 5.0) for r in range(10)]},
                                pname="drn")
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="gwf.cbc", head_filerecord="gwf.hds",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")], printrecord=[("BUDGET", "ALL")])
    if transport:
        from basinledger.coupling.transport import add_transport
        for species, settings in transport.items():
            add_transport(sim, gwf, species, sources=[("rch", species), ("riv", species)], **settings)
    sim.write_simulation(silent=True)
    return ws
