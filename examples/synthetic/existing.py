"""A MODFLOW 6 model written in someone else's style, for the second way in (docs/design.md, "Extending an existing
MODFLOW 6 model"): Basin Ledger takes it as given and couples SWAT+ to it. The tests use it as a stand-in for a
consultant's or USGS's model until a published one is chosen. Its places and numbers mean nothing.

What makes it foreign, against the models ``model.py`` writes for the studies' own builds:

- **Georeferenced and rotated:** the grid's origin is at UTM-like coordinates (611,000 E, 4,703,000 N) and it is
  rotated 15 degrees, so maps from geometry must go through the model's own grid.
- **Two layers, with inactive cells:** 12 x 12 cells of 31.6 m in each layer. The top layer's outer ring is inactive,
  leaving 10 x 10 active cells, Osu_1hru's 10 ha; the bottom layer is active everywhere.
- **Monthly stress periods,** each with daily steps; or, with ``daily=False``, ``steps`` steps a month (one by
  default) growing by ``tsmult``, as a model built for a study usually has. A run rewrites those to daily steps in
  its copy (coupling/timesteps.py), so both give the same run.
- **Its own package names and its own boundaries:** the recharge it was calibrated with (``rch_swb``, a rate for each
  month), a creek (``riv_creek``) down the second column with a tributary along row 5 that no SWAT+ channel reaches,
  tile drains (``drn_tile``), a regional boundary on the east (``ghb_regional``), a fixed head on the south
  (``chd_lake``), two city wells (``wel_city``), one of which a coupled user takes over, and riparian
  evapotranspiration (``evt_riparian``). With ``septic=True``, a second recharge package (``rch_septic``), which
  the start-up check refuses: recharge would be counted twice.

``grid="disv"`` writes the same cells as a vertex grid (icell2d = row x 12 + col).

``neighbor=True`` (DIS only) puts a second groundwater model in the simulation, ``kzoo_east``: 12 x 4 cells of the
same aquifer east of the grid, held at the regional head on its far side, joined to ``kzoo``'s east column of the
bottom layer by a GWF-GWF exchange and solved with it, as an inset or a neighbor in someone's simulation would be.
Basin Ledger couples ``kzoo`` and books what crosses to ``kzoo_east`` as water from or to outside.

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
POND = [(1, 0, 5), (1, 0, 6)]  # a lake where the top layer is inactive, its bed on the bottom layer (lake=True)
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
          steps: int = 1, tsmult: float = 1.0, septic: bool | tuple = False, length_units: str = "meters",
          time_units: str = "days", state_units: bool = True, complexity: str = "moderate",
          neighbor: bool = False, lake: bool = False, neighbor_recharge: bool = False,
          neighbor_river: bool | str = False, maw: bool = False, steady: str | None = None,
          neighbor_drains: bool = False, neighbor_lake: bool = False, soil: bool = False,
          riparian: bool = True, land: float | None = None) -> Path:
    if grid not in ("dis", "disv"):
        raise ValueError("grid is dis or disv")
    if neighbor and grid != "dis":
        raise ValueError("neighbor=True is built on DIS")
    u = units(length_units, time_units)
    Q, LT, COND = (3, -1), (1, -1), (2, -1)  # (length, time) exponents; lengths are u's default
    ws = Path(workspace)
    lengths = months(start, days)
    nper = len(lengths)
    sim = flopy.mf6.MFSimulation(sim_name="kzoo_regional", sim_ws=ws, exe_name="mf6")
    periods = [(u(n, 0, 1), n, 1.0) if daily else (u(n, 0, 1), min(steps, n), tsmult) for n in lengths]
    flopy.mf6.ModflowTdis(sim, time_units=time_units if state_units else None, start_date_time=start.isoformat(),
                          nper=nper, perioddata=periods)
    rclose = 1e-8 if time_units == "days" else 1e-7  # m3 d-1: in seconds, 1e-8 is 4e-12 ft3 s-1, past the solver
    ims = flopy.mf6.ModflowIms(sim, complexity=complexity, outer_dvclose=u(1e-8), inner_dvclose=u(1e-10),
                         rcloserecord=[u(rclose, *Q), "strict"], pname="solver")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="kzoo", save_flows=True)
    idomain = np.ones((2, N, N), int)
    idomain[0, [0, -1], :] = 0
    idomain[0, :, [0, -1]] = 0
    c = lambda x: cellid(x, grid)  # noqa: E731
    lu = length_units if state_units else None
    land = TOP if land is None else land  # land surface, m: the top layer's top, in both models
    cell, top, botm, (x0, y0) = u(CELL), u(land), [u(MIDDLE), u(BOTTOM)], (u(ORIGIN[0]), u(ORIGIN[1]))
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
    if steady is None:
        flopy.mf6.ModflowGwfsto(gwf, iconvert=[1, 0], sy=0.2, ss=[u(1e-5, -1), u(5e-6, -1)], transient={0: True},
                                pname="sto_glacial")
    elif steady == "sto":  # handed over steady, its storage package's periods steady
        flopy.mf6.ModflowGwfsto(gwf, iconvert=[1, 0], sy=0.2, ss=[u(1e-5, -1), u(5e-6, -1)], steady_state={0: True},
                                pname="sto_glacial")
    elif steady != "no_sto":  # "no_sto": handed over steady with no storage package, as Kalamazoo's
        raise ValueError("steady is None, 'sto' or 'no_sto'")
    flopy.mf6.ModflowGwfic(gwf, strt=u(STAGE + 0.5))
    seasonal = [0.0008, 0.0006, 0.0012, 0.0015, 0.0007, 0.0003, 0.0001, 0.0001, 0.0002, 0.0005, 0.0009, 0.0010]
    month = lambda i: (start.month - 1 + i) % 12  # noqa: E731
    flopy.mf6.ModflowGwfrcha(gwf, recharge={i: u(seasonal[month(i)], *LT) for i in range(nper)}, pname="rch_swb",
                             filename="kzoo.rch_swb.rcha")
    if septic:  # True: on a cell of SWAT+'s ground; or the cell, such as one outside it
        cell_septic = (0, 6, 6) if septic is True else tuple(septic)
        flopy.mf6.ModflowGwfrch(gwf, stress_period_data={0: [(c(cell_septic), u(0.002, *LT))]}, pname="rch_septic",
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
    if riparian:  # left out for groundwater up into the soil, which SWAT+'s plants take ET from
        flopy.mf6.ModflowGwfevt(gwf, stress_period_data={0: [(c(x), u(TOP - 11.0), u(0.002, *LT), u(2.0))
                                                             for x in RIPARIAN]},
                                pname="evt_riparian", filename="kzoo.evt_riparian.evt")
    if soil:  # groundwater up into SWAT+'s soil: a drain at land surface on each cell of the two-model map
        flopy.mf6.ModflowGwfdrn(gwf, stress_period_data={0: [(c((0, r, k)), u(land), u(10.0, *COND))
                                                             for r in range(1, 9) for k in range(1, 9)]},
                                pname="soil", filename="kzoo.soil.drn")
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="kzoo.cbc", head_filerecord="kzoo.hds",
                           saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")], printrecord=[("BUDGET", "LAST")])
    if maw:  # a farm's multi-aquifer well, screened through both layers, idle in its own files
        r, col = 7, 4
        flopy.mf6.ModflowGwfmaw(gwf, pname="maw_farm", nmawwells=1, packagedata=[(0, u(0.15), u(BOTTOM), u(STAGE),
                                                                                  "thiem", 2)],
                                connectiondata=[(0, 0, c((0, r, col)), u(TOP), u(MIDDLE), 0.0, 0.0),
                                                (0, 1, c((1, r, col)), u(MIDDLE), u(BOTTOM), 0.0, 0.0)],
                                perioddata={0: [(0, "rate", 0.0)]}, filename="kzoo.maw_farm.maw")
    if lake:  # a pond in the top layer's inactive ring, as a lake's cells are, its bed connected down to the cells below
        conn = [(0, k, c(x), "vertical", u(0.01, 0, -1), 0.0, 0.0, 0.0, 0.0) for k, x in enumerate(POND)]
        flopy.mf6.ModflowGwflak(gwf, pname="lak_pond", nlakes=1, packagedata=[(0, u(MIDDLE + 0.5), len(conn), "pond")],
                                connectiondata=conn, boundnames=True, noutlets=0,
                                perioddata={0: [(0, "rainfall", u(0.002, *LT)), (0, "evaporation", u(0.003, *LT))]},
                                filename="kzoo.lak_pond.lak")
    if neighbor:
        _neighbor(sim, ims, u, lu, neighbor_recharge, neighbor_river, neighbor_drains, neighbor_lake, soil, land)
    sim.write_simulation(silent=True)
    return ws


def reservoir_project(dest: str | Path, refdata: str | Path) -> Path:
    """Osu_1hru with a reservoir, res01, below its channel (cha01 flows into it, and it out of the watershed), for
    the lake (``lake=True``): built from the start, its own seepage off (k = 0), a principal spillway at 1 ha and
    1 ha-m, an emergency one at 1.5 ha and 2 ha-m, released by res_rel.dtl's drawdown_days. The formats follow SWAT+
    62.0.1's readers, as swatplus-bmi's tests/reservoir.py does."""
    import shutil
    dest = Path(dest)
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(Path(refdata) / "Osu_1hru", dest)
    lines = (dest / "object.cnt").read_text().splitlines()
    f = lines[2].split()
    f[3], f[10] = str(int(f[3]) + 1), "1"                      # obj, res
    lines[2] = "  ".join(f)
    (dest / "object.cnt").write_text("\n".join(lines) + "\n")
    lines = (dest / "chandeg.con").read_text().splitlines()
    f = lines[2].split()
    f[12] = "1"                                                # out_tot: all of it to res 1
    lines[2] = "  ".join(f + ["res", "1", "tot", "1.00000"])
    (dest / "chandeg.con").write_text("\n".join(lines) + "\n")
    (dest / "reservoir.con").write_text("reservoir.con: a reservoir below cha01 (existing.py)\n"
                                        "id name gis_id area lat lon elev res wst cst ovfl rule out_tot\n"
                                        "1 res01 1 10.0 35.52 127.33 0.0 1 s35610n127290e 0 0 0 0\n")
    (dest / "reservoir.res").write_text("reservoir.res: existing.py\nid name init hyd rel sed nut\n"
                                        "1 res01 initwet1 hydres01 drawdown_days sedwet1 nutwet1\n")
    (dest / "hydrology.res").write_text("hydrology.res: existing.py\n"
                                        "name yr_op mon_op area_ps vol_ps area_es vol_es k evap_co shp_co1 shp_co2\n"
                                        "hydres01 0 1 1.0 1.0 1.5 2.0 0.0 0.6 0.0 0.0\n")
    lines = (dest / "file.cio").read_text().splitlines()
    for i, line in enumerate(lines):
        f = line.split()
        if f and f[0] == "connect":
            f[8] = "reservoir.con"
        elif f and f[0] == "reservoir":
            f[2], f[3] = "reservoir.res", "hydrology.res"
        else:
            continue
        lines[i] = "  ".join(f)
    (dest / "file.cio").write_text("\n".join(lines) + "\n")
    return dest


EAST_COLS = 4


def _neighbor(sim, ims, u, lu, recharge=False, river=False, drains=False, lake=False, soil=False, land=TOP):
    """``kzoo_east``: the aquifer east of the grid, its far column at the regional head, joined to kzoo's east
    column (bottom layer, where kzoo's top layer is inactive) by a GWF-GWF exchange. With ``recharge``, its own
    recharge package, ``rch_swb`` as kzoo's, for coupling both models together; with ``river``, a river above the
    regional head down its second column, ``riv_east`` (or the name ``river`` gives), losing to the aquifer; with ``drains``, tile drains of its
    own named as kzoo's (``drn_tile``) down its third column; with ``lake``, a lake at a fixed stage above the head
    at the river's foot (``lakes_east``, a RIV package of its own), gaining the aquifer."""
    LT = (1, -1)
    east = flopy.mf6.ModflowGwf(sim, modelname="kzoo_east", save_flows=True)
    cell = u(CELL)
    flopy.mf6.ModflowGwfdis(east, length_units=lu, nlay=2, nrow=N, ncol=EAST_COLS, delr=cell, delc=cell, top=u(land),
                            botm=[u(MIDDLE), u(BOTTOM)], pname="dis_east")
    flopy.mf6.ModflowGwfnpf(east, icelltype=[1, 0], k=[u(8.0, *LT), u(2.0, *LT)], k33=[u(0.8, *LT), u(0.2, *LT)])
    flopy.mf6.ModflowGwfsto(east, iconvert=[1, 0], sy=0.2, ss=[u(1e-5, -1), u(5e-6, -1)], transient={0: True})
    flopy.mf6.ModflowGwfic(east, strt=u(STAGE + 0.5))
    flopy.mf6.ModflowGwfchd(east, stress_period_data={0: [((k, r, EAST_COLS - 1), u(STAGE + 1.0)) for k in range(2)
                                                          for r in range(N)]}, pname="chd_regional")
    if recharge:
        flopy.mf6.ModflowGwfrcha(east, recharge=u(0.0008, *LT), pname="rch_swb", filename="kzoo_east.rch_swb.rcha")
    if river:
        flopy.mf6.ModflowGwfriv(east, stress_period_data={0: [((0, r, 1), u(STAGE + 2.0), u(10.0, 2, -1),
                                                                u(STAGE - 1.0)) for r in range(2, 10)]},
                                pname=river if isinstance(river, str) else "riv_east",
                                filename=f"kzoo_east.{river if isinstance(river, str) else 'riv_east'}.riv")
    if drains:
        flopy.mf6.ModflowGwfdrn(east, stress_period_data={0: [((0, r, 2), u(STAGE + 0.3), u(4.0, 2, -1))
                                                                for r in range(2, 10)]},
                                pname="drn_tile", filename="kzoo_east.drn_tile.drn")
    if soil:  # a drain at land surface on each top cell off the regional head, as kzoo's ``soil``
        flopy.mf6.ModflowGwfdrn(east, stress_period_data={0: [((0, r, k), u(land), u(10.0, 2, -1)) for r in range(N)
                                                                for k in range(EAST_COLS - 1)]},
                                pname="soil", filename="kzoo_east.soil.drn")
    if lake:
        flopy.mf6.ModflowGwfriv(east, stress_period_data={0: [((0, 10, c), u(STAGE + 1.5), u(10.0, 2, -1),
                                                                u(STAGE - 1.0)) for c in (1, 2)]},
                                pname="lakes_east", filename="kzoo_east.lakes_east.riv")
    flopy.mf6.ModflowGwfoc(east, budget_filerecord="kzoo_east.cbc", saverecord=[("BUDGET", "LAST")])
    flopy.mf6.ModflowGwfgwf(sim, exgtype="GWF6-GWF6", exgmnamea="kzoo", exgmnameb="kzoo_east", nexg=N,
                            exchangedata=[((1, r, N - 1), (1, r, 0), 1, cell / 2, cell / 2, cell * u(MIDDLE - BOTTOM))
                                          for r in range(N)], filename="kzoo_kzoo_east.gwfgwf")
    sim.register_solution_package(ims, ["kzoo", "kzoo_east"])


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
