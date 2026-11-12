"""Building SFR reaches for the exchange from channel lines (coupling/sfr.py), on structured and vertex grids.

- lines cut into reaches, in order downstream, the line's whole length kept; a line back into a cell makes a
  second reach there; on a refined vertex grid, the reaches are its cells of each size;
- the package: reaches chained within each channel only, Manning's constant for meters and days, and what the
  exchange takes without a problem; its map read back by the exchange;
- refused: a streambed rising downstream, a flat one without min_slope, a one-reach channel without its slope,
  another Manning's constant, other units, values that don't match the reaches (tier 1).
The synthetic SFR example (sfr.toml) is built with it, so tests/test_sfr.py runs it coupled.
"""

from __future__ import annotations

import datetime as dt
import importlib.util

import flopy
import numpy as np
import pytest

shapely = pytest.importorskip("shapely")

from basinledger.coupling.exchange import ChannelMap, StreamExchange, sfr_reach_links  # noqa: E402
from basinledger.coupling.modflow6 import Modflow6  # noqa: E402
from basinledger.coupling.sfr import add_sfr, channel_reaches, slopes, write_reach_map  # noqa: E402
from conftest import ROOT, RUNS  # noqa: E402

SYNTHETIC = ROOT / "examples" / "synthetic"


def synthetic_model():
    spec = importlib.util.spec_from_file_location("synthetic_model", SYNTHETIC / "model.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


M = synthetic_model()
SIDE = M.SMALL_SIDE
WEST = shapely.LineString([(M.SMALL_CELL / 2, SIDE), (M.SMALL_CELL / 2, 0.0)])       # down the west column
DIAGONAL = shapely.LineString([(SIDE * 0.95, SIDE * 0.95), (SIDE * 0.35, SIDE * 0.05)])


def small(name, disv=False):
    ws = RUNS / name
    (M.build_small_disv if disv else M.build_small)(ws, 2, dt.date(2010, 1, 1))
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(ws), verbosity_level=0)
    gwf = sim.get_model("gwf")
    gwf.remove_package("riv")
    return sim, gwf


def test_lines_into_reaches_in_order():
    _, gwf = small("sfrb-order")
    reaches = channel_reaches(gwf.modelgrid, {"cha01": WEST, "cha02": DIAGONAL})
    # the grid read back from its files has its cell size to 8 digits, so it ends ~2e-8 m inside the west line's ends
    inside = WEST.intersection(shapely.box(*gwf.modelgrid.extent[::2], *gwf.modelgrid.extent[1::2]))
    for ch, line in (("cha01", inside), ("cha02", DIAGONAL)):
        mine = [r for r in reaches if r.channel == ch]
        assert sum(r.length for r in mine) == pytest.approx(line.length, rel=1e-12)
        d = [r.distance for r in mine]
        assert d == sorted(d) and 0 < d[0] < d[-1] < line.length
    west = [r for r in reaches if r.channel == "cha01"]
    assert [r.cell for r in west] == [(row, 0) for row in range(10)]
    assert all(r.length == pytest.approx(gwf.modelgrid.delc[0], rel=1e-12) for r in west)


def test_a_line_back_into_a_cell_makes_a_second_reach():
    _, gwf = small("sfrb-back")
    c = M.SMALL_CELL
    zig = shapely.LineString([(0.5 * c, 5.5 * c), (1.5 * c, 5.5 * c), (1.5 * c, 4.5 * c), (0.6 * c, 4.5 * c),
                              (0.6 * c, 5.2 * c)])  # into the cell east of it and back into the first cell's row
    reaches = channel_reaches(gwf.modelgrid, {"cha01": zig})
    cells = [r.cell for r in reaches]
    assert sum(r.length for r in reaches) == pytest.approx(zig.length, rel=1e-12)
    assert len(cells) > len(set(cells))  # a cell visited twice has two reaches


def test_reaches_on_a_refined_vertex_grid():
    _, gwf = small("sfrb-disv", disv=True)
    reaches = channel_reaches(gwf.modelgrid, {"cha01": WEST})
    grid = {c["icell2d"]: c for c in M.small_disv_grid()}
    assert all(len(r.cell) == 1 and grid[r.cell[0]]["west"] for r in reaches)
    # the grid read back has its vertices to ~10 digits, so its edge is ~1e-8 m inside the line's ends
    assert sum(r.length for r in reaches) == pytest.approx(WEST.length, rel=1e-9)
    # down the west column of the refined grid: 79 m cells all the way (the refined block is in cols 2-3)
    assert [r.cell[0] for r in reaches] == [0, 4, 14, 24]


def build(name, lines, disv=False, **kw):
    sim, gwf = small(name, disv)
    reaches = channel_reaches(gwf.modelgrid, lines)
    tops = [88.0 - 0.001 * r.distance for r in reaches]
    pkg = add_sfr(gwf, reaches, **{**dict(top=tops, width=2.0, thickness=1.0, k=0.3, manning=0.035), **kw})
    sim.write_simulation(silent=True)
    return sim, reaches, pkg


@pytest.mark.parametrize("disv", [False, True])
def test_the_package_is_what_the_exchange_takes(libmf6, tmp_path, disv):
    sim, reaches, pkg = build(f"sfrb-ok-{disv}", {"cha01": WEST, "cha02": DIAGONAL}, disv=disv)
    assert float(pkg.unit_conversion.get_data()) == 86400.0
    write_reach_map(tmp_path / "sfr.csv", reaches)
    mf = Modflow6(sim.sim_path)
    try:
        mf.prepare_day()
        cmap = ChannelMap.from_spec("sfr.csv", tmp_path, mf, "sfr", ["cha01", "cha02"])
        assert [cmap.channels[i] for i in cmap.entry_channel] == [r.channel for r in reaches]
        ex = StreamExchange(mf, "sfr", cmap)  # no problem refused
        first = {ch: next(i for i, r in enumerate(reaches) if r.channel == ch) for ch in ("cha01", "cha02")}
        assert ex._top == first
        links = sfr_reach_links(mf, "sfr")
        assert all(reaches[a].channel == reaches[b].channel and b == a + 1 for a, b in links)
        assert len(links) == len(reaches) - 2
        ex.set_inflow({"cha01": 500.0, "cha02": 300.0})
        mf.solve_day()
    finally:
        mf.finalize()


def test_slopes_from_the_tops():
    _, gwf = small("sfrb-slope")
    reaches = channel_reaches(gwf.modelgrid, {"cha01": WEST})
    top = np.array([88.0 - 0.002 * r.distance for r in reaches])
    np.testing.assert_allclose(slopes(reaches, top), 0.002, rtol=1e-9)
    top[4] = top[3]  # a flat reach
    with pytest.raises(ValueError, match="reaches \\[3\\] are flat"):
        slopes(reaches, top)
    assert slopes(reaches, top, min_slope=1e-4)[3] == 1e-4
    top[6] = top[5] + 0.05
    with pytest.raises(ValueError, match="rises 0.0[0-9]+ m downstream from reach 5"):
        slopes(reaches, top, min_slope=1e-4)


@pytest.mark.parametrize("kw, refused", [
    (dict(unit_conversion=1.0), "unit_conversion is set here"),
    (dict(width={"cha01": 2.0}), "width: no value for channels \\['cha02'\\]"),
    (dict(k=[0.3, 0.3]), "k: 2 values for"),
    (dict(manning=0.0), "manning: every reach's must be above 0"),
])
def test_what_the_helper_refuses(kw, refused):
    with pytest.raises(ValueError, match=refused):
        build("sfrb-refused", {"cha01": WEST, "cha02": DIAGONAL}, **kw)


def test_one_reach_needs_its_slope():
    c = M.SMALL_CELL
    short = shapely.LineString([(0.5 * c, 5.8 * c), (0.5 * c, 5.2 * c)])
    with pytest.raises(ValueError, match="one reach, so its slope can't come from the tops"):
        build("sfrb-one", {"cha01": short})
    build("sfrb-one", {"cha01": short}, slope=0.001)


def test_other_units_are_refused():
    sim, gwf = small("sfrb-units")
    gwf.dis.length_units = "feet"
    reaches = channel_reaches(gwf.modelgrid, {"cha01": WEST})
    with pytest.raises(ValueError, match="meters and days"):
        add_sfr(gwf, reaches, top=88.0, width=2.0, thickness=1.0, k=0.3, manning=0.035, slope=0.001)
