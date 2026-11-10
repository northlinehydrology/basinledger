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
                transport: dict | None = None, top: float = 100.0, soil: str | None = None,
                soil_conductance: float = 10.0) -> Path:
    """A 10 ha aquifer under Osu_1hru for the two-way exchange: 10 x 10 cells, the river along the west edge.
    Thinner and less permeable than the large grid so the river can lose water when the well pumps hard.
    ``drain_elevation``: seepage cells (a drain) along the second column at that elevation, for the coupled drains'
    tests. ``transport``: species -> ``add_transport`` settings (porosity, dispersivity, strt, sorption), a transport
    model per species with recharge and the river as its sources. ``top`` and ``soil``: land surface (m) and a
    package named "soil" on every cell for groundwater up into the soil, a drain ("drn", conductance
    ``soil_conductance`` m2 d-1 [A], elevation set by Basin Ledger each day) or a well ("wel", gwflow's rule)."""
    ws = Path(workspace)
    sim = flopy.mf6.MFSimulation(sim_name="synthetic_small", sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, time_units="days", start_date_time=start.isoformat(),
                          nper=1, perioddata=[(float(days), days, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="simple", outer_dvclose=1e-8, inner_dvclose=1e-10,
                         rcloserecord=[1e-8, "strict"])
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    flopy.mf6.ModflowGwfdis(gwf, length_units="meters", nlay=1, nrow=10, ncol=10, delr=SMALL_CELL, delc=SMALL_CELL,
                            top=top, botm=70.0)
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
    cells = [(0, r, c) for r in range(10) for c in range(10)]
    if soil == "drn":
        flopy.mf6.ModflowGwfdrn(gwf, stress_period_data={0: [(x, top, soil_conductance) for x in cells]}, pname="soil",
                                filename="gwf.soil.drn")
    elif soil == "wel":
        flopy.mf6.ModflowGwfwel(gwf, stress_period_data={0: [(x, 0.0) for x in cells]}, pname="soil",
                                filename="gwf.soil.wel")
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="gwf.cbc", head_filerecord="gwf.hds",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")], printrecord=[("BUDGET", "ALL")])
    if transport:
        from basinledger.coupling.transport import add_transport
        for species, settings in transport.items():
            add_transport(sim, gwf, species, sources=[("rch", species), ("riv", species)], **settings)
    sim.write_simulation(silent=True)
    return ws


SMALL_SIDE = 10 * SMALL_CELL  # m: the small aquifer's side


def small_disv_grid(refine: bool = True) -> list[dict]:
    """The small aquifer as a vertex grid, cells in icell2d order: each a dict with ``icell2d``, ``x``, ``y`` (center),
    ``side`` (m), ``west`` (on the west edge, where the river is) and ``vertices`` (clockwise, (x, y)).

    ``refine``: 4 x 4 cells of 79 m with the 2 x 2 block east of center (rows 1-2, columns 2-3, row 0 at the top)
    split into 40 m cells, 28 cells in all, a two-level quadtree; a large cell beside two small ones lists the
    vertex mid-edge, so MODFLOW 6 connects it to both. Without, 10 x 10 cells of 31.6 m, the cells of ``build_small``
    numbered row by row (icell2d = row * 10 + col)."""
    n, coarse = (8, 2) if refine else (10, 1)   # fine lattice steps per side, a large cell's side in steps
    step = SMALL_SIDE / n
    split = {(1, 2), (1, 3), (2, 2), (2, 3)} if refine else set()
    boxes = []  # (i, j, size) in lattice steps, i east from the west edge, j south from the north edge
    for r in range(n // coarse):
        for c in range(n // coarse):
            if (r, c) in split:
                boxes += [(c * 2 + a, r * 2 + b, 1) for b in (0, 1) for a in (0, 1)]
            else:
                boxes.append((c * coarse, r * coarse, coarse))
    corners = {(i + a * s, j + b * s) for i, j, s in boxes for a in (0, 1) for b in (0, 1)}
    cells = []
    for k, (i, j, s) in enumerate(boxes):
        ring = ([(i + t, j) for t in range(s)] + [(i + s, j + t) for t in range(s)] +
                [(i + s - t, j + s) for t in range(s)] + [(i, j + s - t) for t in range(s)])
        ring = [p for p in ring if p in corners]  # clockwise from the north-west corner
        cells.append(dict(icell2d=k, x=(i + s / 2) * step, y=SMALL_SIDE - (j + s / 2) * step, side=s * step,
                          west=i == 0, vertices=[(a * step, SMALL_SIDE - b * step) for a, b in ring]))
    return cells


def small_disv_well(refine: bool = True) -> tuple[int, int]:
    """(layer, icell2d) of the pumping well: a small cell of the refined block, or ``SMALL_WELL``'s cell."""
    cells = small_disv_grid(refine)
    if not refine:
        return (0, SMALL_WELL[1] * 10 + SMALL_WELL[2])
    small = [c for c in cells if c["side"] < SMALL_SIDE / 4]
    near = min(small, key=lambda c: (c["x"] - 0.6 * SMALL_SIDE) ** 2 + (c["y"] - 0.5 * SMALL_SIDE) ** 2)
    return (0, near["icell2d"])


def build_small_disv(workspace: str | Path, days: int, start: dt.date, refine: bool = True,
                     transport: dict | None = None) -> Path:
    """``build_small`` on a vertex grid (DISV, ``small_disv_grid``), for the coupling's tests on refined grids: the
    same aquifer, area, river stage and bottom, the river's conductance per metre of river (20 m2 d-1 per 31.6 m
    cell) and the well in ``small_disv_well``. With ``refine = False`` it is ``build_small``'s grid cell for cell.
    ``transport`` as in ``build_small``."""
    ws = Path(workspace)
    cells = small_disv_grid(refine)
    points = sorted({v for c in cells for v in c["vertices"]})
    index = {v: k for k, v in enumerate(points)}
    sim = flopy.mf6.MFSimulation(sim_name="synthetic_disv", sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, time_units="days", start_date_time=start.isoformat(),
                          nper=1, perioddata=[(float(days), days, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="simple", outer_dvclose=1e-8, inner_dvclose=1e-10,
                         rcloserecord=[1e-8, "strict"])
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    flopy.mf6.ModflowGwfdisv(gwf, length_units="meters", nlay=1, ncpl=len(cells), nvert=len(points),
                             vertices=[[k, x, y] for k, (x, y) in enumerate(points)],
                             cell2d=[[c["icell2d"], c["x"], c["y"], len(c["vertices"])] +
                                     [index[v] for v in c["vertices"]] for c in cells],
                             top=100.0, botm=70.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=2.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, sy=0.15, ss=1e-5, transient={0: True})
    flopy.mf6.ModflowGwfic(gwf, strt=RIVER_STAGE)
    flopy.mf6.ModflowGwfrcha(gwf, recharge=0.0, pname="rch")
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data={0: [(small_disv_well(refine), 0.0)]}, pname="wel")
    flopy.mf6.ModflowGwfriv(gwf, stress_period_data={0: [((0, c["icell2d"]), RIVER_STAGE, 20.0 * c["side"] / SMALL_CELL,
                                                          RIVER_STAGE - 2.0) for c in cells if c["west"]]}, pname="riv")
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="gwf.cbc", head_filerecord="gwf.hds",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")], printrecord=[("BUDGET", "ALL")])
    if transport:
        from basinledger.coupling.transport import add_transport
        for species, settings in transport.items():
            add_transport(sim, gwf, species, sources=[("rch", species), ("riv", species)], **settings)
    sim.write_simulation(silent=True)
    return ws


def write_small_disv_maps(folder: str | Path, refine: bool = True):
    """The maps the DISV example reads (``maps/disv_recharge.csv``, ``maps/disv_river.csv``), in the vertex grid's
    form (object, layer, icell2d, ...): every cell's whole area to Osu_1hru's aquifer, and the river cells to its
    channel."""
    import csv
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    cells = small_disv_grid(refine)
    with open(folder / "disv_recharge.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["object", "layer", "icell2d", "area_m2"])
        for c in cells:
            w.writerow(["aqu001", 0, c["icell2d"], repr(c["side"] ** 2)])
    with open(folder / "disv_river.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["object", "layer", "icell2d"])
        for c in cells:
            if c["west"]:
                w.writerow(["cha01", 0, c["icell2d"]])
