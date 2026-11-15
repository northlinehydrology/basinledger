"""A MODFLOW 6 model written in someone else's style, for the second way in (docs/design.md, "Extending an existing
MODFLOW 6 model"): Basin Ledger takes it as given and couples SWAT+ to it. The tests use it as a stand-in for a
consultant's or USGS's model until a published one is chosen. Its places and numbers mean nothing.

What makes it foreign, against the models ``model.py`` writes for the studies' own builds:

- **Georeferenced and rotated:** the grid's origin is at UTM-like coordinates (611,000 E, 4,703,000 N) and it is
  rotated 15 degrees, so maps from geometry must go through the model's own grid.
- **Two layers, with inactive cells:** 12 x 12 cells of 31.6 m in each layer. The top layer's outer ring is inactive,
  leaving 10 x 10 active cells, Osu_1hru's 10 ha; the bottom layer is active everywhere.
- **Monthly stress periods,** each with daily steps (one step per month with ``daily=False``, which the driver
  refuses until it can rewrite the time discretization).
- **Its own package names and its own boundaries:** the recharge it was calibrated with (``rch_swb``, a rate for each
  month), a creek (``riv_creek``) down the second column with a tributary along row 5 that no SWAT+ channel reaches,
  tile drains (``drn_tile``), a regional boundary on the east (``ghb_regional``), a fixed head on the south
  (``chd_lake``), two city wells (``wel_city``), one of which a coupled user takes over, and riparian
  evapotranspiration (``evt_riparian``). With ``septic=True``, a second recharge package (``rch_septic``), which
  the start-up check refuses: recharge would be counted twice.

``grid="disv"`` writes the same cells as a vertex grid (icell2d = row x 12 + col).

``length_units`` and ``time_units`` write the same model in other units (feet, centimeters; seconds, minutes,
hours): every length, area, rate, conductance and solver tolerance converted, so it is the meters-and-days model's
twin and a coupled run of it gives the same ledger (Basin Ledger converts at each link). ``state_units=False`` leaves
LENGTH_UNITS and TIME_UNITS undefined, as mf5to6 does, for the config to state them. The twins take IMS's simple
settings (``complexity="simple"``): the moderate settings' backtracking limits are absolute, so they act differently
in other units and the solves part. In time units other than days the residual criterion is 1e-7 m3 d-1 rather than
1e-8: in seconds 1e-8 m3 d-1 is 4e-12 ft3 s-1, tighter than the linear solver gets on the first step.
"""

from __future__ import annotations

import calendar
import datetime as dt
from pathlib import Path

import flopy
import numpy as np

N = 12                       # rows and columns per layer
CELL = 1000.0 ** 0.5         # m: 10 x 10 active top cells of 1,000 m2 make 10 ha, Osu_1hru's aquifer
ORIGIN = (611000.0, 4703000.0)
ANGLE = 15.0                 # degrees, counterclockwise
TOP, MIDDLE, BOTTOM = 100.0, 75.0, 40.0
CREEK = [(0, r, 1) for r in range(1, 11)]       # (layer, row, col) down the second column
TRIBUTARY = [(0, 5, 2), (0, 5, 3)]              # along row 5: in the model, but not a SWAT+ channel
TILES = [(0, 9, c) for c in range(3, 9)]
EAST = [(1, r, N - 1) for r in range(1, N - 1)]
SOUTH = [(1, N - 1, c) for c in range(2, 6)]
CITY_WELLS = [(1, 5, 6), (1, 3, 8)]             # the first is the coupled user's
RIPARIAN = [(0, r, 2) for r in (2, 3, 7, 8)]
STAGE = 88.0


def months(start: dt.date, days: int) -> list[int]:
    """Lengths of the calendar months the run covers, the first and last cut to the run."""
    out, d, end = [], start, start + dt.timedelta(days=days)
    while d < end:
        nxt = dt.date(d.year + d.month // 12, d.month % 12 + 1, 1)
        out.append((min(nxt, end) - d).days)
        d = nxt
    return out


def units(length_units: str = "meters", time_units: str = "days"):
    """A value in meters and days, of a (length, time) dimension, in the model's units."""
    from basinledger.coupling.modflow6 import DAYS, METERS
    lf, tf = 1.0 / METERS[length_units], 1.0 / DAYS[time_units]
    return lambda v, length=1, time=0: np.asarray(v, float) * lf ** length * tf ** time if np.ndim(v) else \
        float(v) * lf ** length * tf ** time


def cellid(cell, grid: str):
    """(layer, row, col) in the grid's own form."""
    return cell if grid == "dis" else (cell[0], cell[1] * N + cell[2])


def build(workspace: str | Path, days: int, start: dt.date, grid: str = "dis", daily: bool = True,
          septic: bool = False, length_units: str = "meters", time_units: str = "days",
          state_units: bool = True, complexity: str = "moderate") -> Path:
    if grid not in ("dis", "disv"):
        raise ValueError("grid is dis or disv")
    u = units(length_units, time_units)
    Q, LT, COND = (3, -1), (1, -1), (2, -1)  # (length, time) exponents; lengths are u's default
    ws = Path(workspace)
    lengths = months(start, days)
    nper = len(lengths)
    sim = flopy.mf6.MFSimulation(sim_name="kzoo_regional", sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, time_units=time_units if state_units else None, start_date_time=start.isoformat(),
                          nper=nper, perioddata=[(u(n, 0, 1), n if daily else 1, 1.0) for n in lengths])
    rclose = 1e-8 if time_units == "days" else 1e-7  # m3 d-1: in seconds, 1e-8 is 4e-12 ft3 s-1, past the solver
    flopy.mf6.ModflowIms(sim, complexity=complexity, outer_dvclose=u(1e-8), inner_dvclose=u(1e-10),
                         rcloserecord=[u(rclose, *Q), "strict"], pname="solver")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="kzoo", save_flows=True)
    idomain = np.ones((2, N, N), int)
    idomain[0, [0, -1], :] = 0
    idomain[0, :, [0, -1]] = 0
    c = lambda x: cellid(x, grid)  # noqa: E731
    lu = length_units if state_units else None
    cell, top, botm, (x0, y0) = u(CELL), u(TOP), [u(MIDDLE), u(BOTTOM)], (u(ORIGIN[0]), u(ORIGIN[1]))
    if grid == "dis":
        flopy.mf6.ModflowGwfdis(gwf, length_units=lu, nlay=2, nrow=N, ncol=N, delr=cell, delc=cell, top=top,
                                botm=botm, idomain=idomain, xorigin=x0, yorigin=y0, angrot=ANGLE, pname="dis_kzoo")
    else:
        verts = [(i * (N + 1) + j, j * cell, (N - i) * cell) for i in range(N + 1) for j in range(N + 1)]
        cell2d = []
        for r in range(N):
            for col in range(N):
                v = [r * (N + 1) + col, r * (N + 1) + col + 1, (r + 1) * (N + 1) + col + 1, (r + 1) * (N + 1) + col]
                cell2d.append((r * N + col, (col + 0.5) * cell, (N - r - 0.5) * cell, 4, *v))
        flopy.mf6.ModflowGwfdisv(gwf, length_units=lu, nlay=2, ncpl=N * N, nvert=len(verts), top=top, botm=botm,
                                 vertices=verts, cell2d=cell2d, idomain=idomain.reshape(2, N * N), xorigin=x0,
                                 yorigin=y0, angrot=ANGLE, pname="disv_kzoo")
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=[1, 0], k=[u(8.0, *LT), u(2.0, *LT)], k33=[u(0.8, *LT), u(0.2, *LT)],
                            pname="k_glacial")
    flopy.mf6.ModflowGwfsto(gwf, iconvert=[1, 0], sy=0.2, ss=[u(1e-5, -1), u(5e-6, -1)], transient={0: True},
                            pname="sto_glacial")
    flopy.mf6.ModflowGwfic(gwf, strt=u(STAGE + 0.5))
    seasonal = [0.0008, 0.0006, 0.0012, 0.0015, 0.0007, 0.0003, 0.0001, 0.0001, 0.0002, 0.0005, 0.0009, 0.0010]
    month = lambda i: (start.month - 1 + i) % 12  # noqa: E731
    flopy.mf6.ModflowGwfrcha(gwf, recharge={i: u(seasonal[month(i)], *LT) for i in range(nper)}, pname="rch_swb",
                             filename="kzoo.rch_swb.rcha")
    if septic:
        flopy.mf6.ModflowGwfrch(gwf, stress_period_data={0: [(c((0, 6, 6)), u(0.002, *LT))]}, pname="rch_septic",
                                filename="kzoo.rch_septic.rch")
    flopy.mf6.ModflowGwfriv(gwf, stress_period_data={
        i: [(c(x), u(STAGE + 0.2 * np.sin(2 * np.pi * month(i) / 12)), u(20.0, *COND), u(STAGE - 2.0)) for x in CREEK]
        + [(c(x), u(STAGE + 0.6), u(5.0, *COND), u(STAGE - 1.0)) for x in TRIBUTARY] for i in range(nper)},
        pname="riv_creek", filename="kzoo.riv_creek.riv")
    flopy.mf6.ModflowGwfdrn(gwf, stress_period_data={0: [(c(x), u(STAGE + 0.8), u(2.0, *COND)) for x in TILES]},
                            pname="drn_tile", filename="kzoo.drn_tile.drn")
    flopy.mf6.ModflowGwfghb(gwf, stress_period_data={0: [(c(x), u(STAGE + 1.0), u(10.0, *COND)) for x in EAST]},
                            pname="ghb_regional", filename="kzoo.ghb_regional.ghb")
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [(c(x), u(STAGE - 0.5)) for x in SOUTH]}, pname="chd_lake",
                            filename="kzoo.chd_lake.chd")
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data={i: [(c(CITY_WELLS[0]), u(-40.0, *Q)),
                                                         (c(CITY_WELLS[1]), u(-25.0, *Q))] for i in range(nper)},
                            pname="wel_city", filename="kzoo.wel_city.wel")
    flopy.mf6.ModflowGwfevt(gwf, stress_period_data={0: [(c(x), u(TOP - 11.0), u(0.002, *LT), u(2.0))
                                                         for x in RIPARIAN]},
                            pname="evt_riparian", filename="kzoo.evt_riparian.evt")
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="kzoo.cbc", head_filerecord="kzoo.hds",
                           saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")], printrecord=[("BUDGET", "LAST")])
    sim.write_simulation(silent=True)
    return ws


def world(row: float, col: float, length_units: str = "meters") -> tuple[float, float]:
    """A point given in cell units (row down, col across, from the grid's top-left corner) in the grid's world
    coordinates, rotated about its origin as MODFLOW 6's XORIGIN, YORIGIN and ANGROT place it, in the grid's
    ``length_units``."""
    u = units(length_units)
    x, y = col * CELL, (N - row) * CELL
    a = np.radians(ANGLE)
    return (u(ORIGIN[0] + x * np.cos(a) - y * np.sin(a)), u(ORIGIN[1] + x * np.sin(a) + y * np.cos(a)))


def creek_line(length_units: str = "meters"):
    """The SWAT+ channel the creek stands for, drawn downstream (north to south) down the second column's middle,
    from the first active row's top to the last's bottom: 10 cells, 316 m."""
    import shapely
    return shapely.LineString([world(1.0, 1.5, length_units), world(11.0, 1.5, length_units)])


def hru_polygons(length_units: str = "meters"):
    """Two HRUs over the active top layer, split at row 6: hru01 to the north (6 rows of 10 cells, 6 ha) and hru02 to
    the south (4 ha)."""
    import shapely
    w = lambda r, c: world(r, c, length_units)  # noqa: E731
    poly = lambda r0, r1: shapely.Polygon([w(r0, 1), w(r0, 11), w(r1, 11), w(r1, 1)])  # noqa: E731
    return {"hru01": poly(1, 7), "hru02": poly(7, 11)}


def aquifer_polygon(length_units: str = "meters"):
    """Osu_1hru's aquifer, aqu001, over the active top layer: 10 ha."""
    import shapely
    w = lambda r, c: world(r, c, length_units)  # noqa: E731
    return shapely.Polygon([w(1, 1), w(1, 11), w(11, 11), w(11, 1)])


def write_maps(folder: str | Path, workspace: str | Path, length_units: str | None = None) -> tuple:
    """The example's maps on the model as built in ``workspace`` (``existing.toml`` reads them from ``maps/``):
    aqu001 to cells by area, and cha01 to the creek's RIV entries by length, through gridmaps, in square meters and
    meters whatever the model's units. ``length_units``: the model's, where it leaves them undefined. Returns the HRU
    map and the stream map, for their reports."""
    from basinledger import gridmaps as G
    _, gwf = G.load(workspace)
    grid = gwf.modelgrid
    lu = length_units or grid.units
    tag = "" if grid.grid_type == "structured" else "_disv"
    hm = G.hru_map(grid, {"aqu001": aquifer_polygon(lu)}, length_units=length_units)
    hm.write(Path(folder) / f"existing_recharge{tag}.csv")
    sm = G.channel_map(grid, G.stream_entries(gwf, "riv_creek"), {"cha01": creek_line(lu)},
                       length_units=length_units)
    sm.write(Path(folder) / f"existing_river{tag}.csv")
    return hm, sm
