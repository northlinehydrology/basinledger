"""The water-table link (``[groundwater_to_soil] rule = "water_table"``; docs/design.md, A water-table link between
MODFLOW 6 and SWAT+'s soils): SWAT+'s layers below MODFLOW 6's water table held, filled from MODFLOW 6's drain at the
profile's bottom (MODFLOW 6's top), what no layer takes returned as recharge the same day.

The first part needs no models: the flag's rule, the table's height, the top check, the return's placement, and the
link's day against stand-ins for swatplus-bmi and MODFLOW 6. The second runs Osu_1hru's one HRU on one MODFLOW 6 cell
whose head a general head boundary raises and lowers; it needs swatplus-bmi's held soil layers and skips without
them."""

from __future__ import annotations

import numpy as np
import pytest

from basinledger.coupling import gwsoil as G
from basinledger.mapping import AreaMap

# one HRU of six layers (Osu_1hru's: bottoms 0.01 ... 0.75 m), one of two
BOTTOMS = np.array([0.01, 0.17, 0.30, 0.45, 0.60, 0.75, 0.2, 1.0])
LAYER_HRU = np.array([0, 0, 0, 0, 0, 0, 1, 1])


def test_each_layers_top_is_the_bottom_of_the_one_above():
    tops = G.layer_tops(LAYER_HRU, BOTTOMS)
    np.testing.assert_allclose(tops, [0.0, 0.01, 0.17, 0.30, 0.45, 0.60, 0.0, 0.2])


def test_the_table_is_weighted_by_the_hrus_ground_in_each_cell():
    # HRU a: 25 m2 on cell 0, 75 m2 on cell 1; HRU b: cell 2, dry; cell tops 100, heads 100.4 and 100.8
    m = AreaMap(["a", "b"], np.array([0, 0, 1]), np.array([0, 1, 2]), np.array([25.0, 75.0, 50.0]))
    h = G.hru_table_height(m, np.array([100.4, 100.8, 1e30]), np.array([100.0, 100.0, 100.0]), 2)
    assert h[0] == pytest.approx(0.25 * 0.4 + 0.75 * 0.8)
    assert h[1] < -1e3  # a dry cell holds nothing
    assert np.isnan(G.hru_table_height(m, np.full(3, 100.0), np.full(3, 100.0), 3)[2])  # an HRU the map misses


def test_a_layer_is_held_above_its_middle_and_released_below_its_bottom():
    depth = np.array([0.75, 1.0])
    tops = G.layer_tops(LAYER_HRU, BOTTOMS)
    held = lambda table_a, before=None: G.held_layers(  # noqa: E731
        np.array([0.75 - table_a, -5.0]), depth, LAYER_HRU, tops, BOTTOMS, before)
    # the table 0.40 m down: layer 4 (0.30-0.45) has it below its middle, so only 5 and 6 are held from nothing
    np.testing.assert_array_equal(held(0.40)[:6], [0, 0, 0, 0, 1, 1])
    # risen to 0.36 m, above layer 4's middle: held
    rose = held(0.36)
    np.testing.assert_array_equal(rose[:6], [0, 0, 0, 1, 1, 1])
    # falling back to 0.40, then 0.44 (still above its bottom): layer 4 stays held, the gap
    np.testing.assert_array_equal(held(0.40, rose), rose)
    np.testing.assert_array_equal(held(0.44, rose), rose)
    # below its bottom: released; below the profile: nothing held
    np.testing.assert_array_equal(held(0.46, rose)[:6], [0, 0, 0, 0, 1, 1])
    np.testing.assert_array_equal(held(0.80, rose)[:6], [0, 0, 0, 0, 0, 0])
    # above the surface: everything
    np.testing.assert_array_equal(held(-0.1)[:6], [1, 1, 1, 1, 1, 1])
    np.testing.assert_array_equal(held(-0.1)[6:], [0, 0])  # the other HRU's table far below


def test_held_layers_run_up_from_the_bottom():
    depth = np.array([0.75, 1.0])
    tops = G.layer_tops(LAYER_HRU, BOTTOMS)
    # a stale flag on layer 2 with nothing held below it is dropped
    before = np.array([0, 1, 0, 0, 0, 0, 0, 0], float)
    out = G.held_layers(np.array([0.75 - 0.05, -5.0]), depth, LAYER_HRU, tops, BOTTOMS, before)
    # the table at 0.05 m: layers 2 (0.01-0.17, middle 0.09) and below held; layer 1 not (its middle 0.005)
    np.testing.assert_array_equal(out[:6], [0, 1, 1, 1, 1, 1])
    out = G.held_layers(np.array([0.75 - 0.9, -5.0]), depth, LAYER_HRU, tops, BOTTOMS, before)
    np.testing.assert_array_equal(out[:6], [0, 0, 0, 0, 0, 0])
    for f in (out, G.held_layers(np.array([0.3, 0.5]), depth, LAYER_HRU, tops, BOTTOMS)):
        for h in (0, 1):
            x = f[LAYER_HRU == h]
            assert np.all(np.diff(x) >= 0), x  # held never above an unheld layer


def test_the_top_check_names_the_hrus_off_the_profiles_bottom():
    m = AreaMap(["a", "b"], np.array([0, 0, 1]), np.array([0, 1, 2]), np.array([50.0, 50.0, 10.0]))
    top = np.array([99.0, 99.5, 90.0])
    land, depth = np.array([100.2, 95.0]), np.array([0.75, 1.0])  # a 0.2 m off, b 4 m
    assert G.top_problems(m, top, land, depth, ["a", "b"], 1.0) == [
        "b: MODFLOW 6's top 90.00 m, its land surface less its profile 94.00 m (-4.00 m)"]
    assert len(G.top_problems(m, top, land, depth, ["a", "b"], 0.1)) == 2


def test_what_no_layer_took_goes_back_to_the_cells_it_came_from():
    m = AreaMap(["a", "b"], np.array([0, 0, 1]), np.array([0, 1, 1]), np.array([30.0, 70.0, 40.0]))
    rows = np.array([1.0, 3.0, 0.0])  # a's water came 1:3 from cells 0 and 1; b's none
    out = G.split_back(rows, m, np.array([8.0, 4.0]), 3)
    np.testing.assert_allclose(out, [2.0, 6.0 + 4.0, 0.0])  # b's spread by its ground: all on cell 1
    assert out.sum() == pytest.approx(12.0)


def test_a_land_surface_per_cell_gives_each_cells_profile_bottom(tmp_path):
    # HRU a (0.75 m profile) on cells 0-2, HRU b (1.0 m) shares cell 2; cell 3 unmapped
    m = AreaMap(["a", "b"], np.array([0, 0, 0, 1]), np.array([0, 1, 2, 2]), np.array([10.0, 10.0, 5.0, 5.0]))
    depth = np.array([0.75, 1.0])
    land = np.array([100.0, 101.0, 102.0, 50.0])
    # cell 0 at its profile's bottom, cell 1 a stream's cell keeping the land surface, cell 2 at its bottom
    top = np.array([99.25, 101.0, 102.0 - 0.875, 40.0])
    bottom, kept, bad = G.profile_bottom_cells(m, top, land, depth, 4, 0.1)
    np.testing.assert_allclose(bottom[:3], [99.25, 100.25, 101.125])
    assert np.isnan(bottom[3])
    np.testing.assert_array_equal(kept, [False, True, False, False])
    assert bad == []
    top[2] = 101.6  # neither its bottom nor its land surface
    land[0] = np.nan
    _, _, bad = G.profile_bottom_cells(m, top, land, depth, 4, 0.1)
    assert bad == ["node 0: no land surface given",
                   "node 2: top 101.60 m, land surface 102.00 m, profile's bottom 101.12 m"]
    # the file: a map table's cell columns and land_m
    f = tmp_path / "land.csv"
    f.write_text("layer,row,col,land_m\n0,0,1,101.0\n")
    g = tmp_path / "nodes.csv"
    g.write_text("node,land_m\n2,102.0\n")
    node_of = lambda *cell, node=None: node if node is not None else cell[2]  # noqa: E731
    np.testing.assert_array_equal(G.read_land_surface(f, node_of, 3), [np.nan, 101.0, np.nan])
    np.testing.assert_array_equal(G.read_land_surface(g, node_of, 3), [np.nan, np.nan, 102.0])


def test_a_model_built_for_water_table_is_told_apart_from_one_at_the_land_surface():
    m = AreaMap(["a", "b"], np.array([0, 1]), np.array([0, 1]), np.array([10.0, 30.0]))
    land, depth = np.array([100.0, 200.0]), np.array([0.75, 1.5])
    assert G.looks_like_profile_bottom(m, np.array([99.25, 198.5]), land, depth) == 1.0
    assert G.looks_like_profile_bottom(m, np.array([100.0, 198.6]), land, depth) == 0.75
    assert G.looks_like_profile_bottom(m, np.array([100.0, 200.0]), land, depth) == 0.0
    assert G.looks_like_profile_bottom(m, np.array([88.5, 88.5]), land, depth) == 0.0  # far from both


# --- the link's day, against stand-ins ---------------------------------------------------------------------------

class StandInSwat:
    """swatplus-bmi as the link uses it: two HRUs, the layers above; Hook 1 takes what the held layers have room
    for and the rest goes to the lowest unheld layer, as swatplus-bmi's hook will."""

    def __init__(self, vars_=None):
        self.vars = set(vars_ if vars_ is not None else (G.INFLOW, G.FLAG, G.NOT_TAKEN,
                                                          G.SATURATION_EXCESS, G.BOTTOM))
        self.values = {G.BOTTOM: BOTTOMS.copy(), G.INFLOW: np.zeros(8), G.FLAG: np.zeros(8),
                       G.NOT_TAKEN: np.zeros(2), G.SATURATION_EXCESS: np.zeros(2)}
        self.sets = []
        self.room = np.full(8, 0.01)  # m of room per layer

    def has(self, v):
        return v in self.vars

    def names(self, v):
        return [f"{h}:{i}" for h, i in zip(np.array(["a", "b"])[LAYER_HRU], range(8))]

    def get(self, v):
        return self.values[v].copy()

    def set(self, v, x):
        self.sets.append(v)
        self.values[v] = np.asarray(x, float).copy()

    def day(self, area):
        """A day: the inflow set (m s-1 per layer) taken by held layers up to their room, the rest not taken."""
        q = self.values[G.INFLOW] * G.DAY * area[LAYER_HRU]          # m3 per layer
        held = self.values[G.FLAG] > 0.5
        took = np.where(held, np.minimum(q, self.room * area[LAYER_HRU]), q)
        nt = np.bincount(LAYER_HRU, weights=q - took, minlength=2)
        self.values[G.INFLOW] = took / area[LAYER_HRU] / G.DAY
        self.values[G.NOT_TAKEN] = nt / area / G.DAY


class StandInModflow:
    def __init__(self, top):
        self.top = np.asarray(top, float)
        self.nodes = len(top)
        self.h = self.top - 1.0
        self.elev = None
        self.out = np.zeros(self.nodes)

    def node_values(self, name, pkg):
        return self.top.copy()

    def package_nodes(self, pkg):
        return np.arange(self.nodes)

    def head(self):
        return self.h.copy()

    def set_entries(self, name, pkg, v):
        self.elev = np.asarray(v, float)

    def simvals(self, pkg):
        return -self.out


def stand_ins():
    m = AreaMap(["a", "b"], np.array([0, 1]), np.array([0, 1]), np.array([100.0, 100.0]))
    return m, StandInSwat(), StandInModflow([99.25, 99.0])


def test_water_table_refuses_an_older_swatplus_bmi_and_a_cap():
    m, sw, _ = stand_ins()
    with pytest.raises(ValueError, match="newer swatplus-bmi"):
        G.GroundwaterToSoil(StandInSwat({G.INFLOW, G.BOTTOM}), m, np.array([100.0, 100.0]), "water_table", "none",
                            "soil")
    with pytest.raises(ValueError, match='cap = "none"'):
        G.GroundwaterToSoil(sw, m, np.array([100.0, 100.0]), "water_table", "room", "soil")
    with pytest.raises(ValueError, match="held_lateral"):
        G.GroundwaterToSoil(sw, m, np.array([100.0, 100.0]), "water_table", "none", "soil", held_lateral="all")


def test_the_link_holds_fills_and_returns_what_isnt_taken():
    m, sw, mf = stand_ins()
    area = np.array([100.0, 100.0])
    land = np.array([100.0, 100.0])  # a's top 99.25 = 100 - 0.75; b's 99.0 = 100 - 1.0
    gs = G.GroundwaterToSoil(sw, m, area, "water_table", "none", "soil", held_lateral="swatplus", land=land,
                             n_cells=2)
    gs.attach(mf)
    np.testing.assert_allclose(mf.top, gs.elev0)  # the drain at the top: the profile's bottom
    # day 0: nothing due, nothing held
    assert gs.before_swat() == 0.0
    sw.day(area)
    assert gs.received() == 0.0 and gs.not_taken == 0.0 and gs.returned is None
    # MODFLOW 6's day: a's head 0.4 m above its top (the table 0.35 m down), b's below; the drain gives 2 m3 on a
    gs.before_solve()
    np.testing.assert_allclose(mf.elev, mf.top)
    mf.h = np.array([99.25 + 0.4, 98.0])
    mf.out = np.array([2.0, 0.0])
    sent, at_bottom = gs.after_solve()
    assert sent == 2.0 and at_bottom == 0.0
    np.testing.assert_array_equal(gs.flags[:6], [0, 0, 0, 1, 1, 1])  # layer 4's middle (0.375) is below 0.35
    np.testing.assert_array_equal(gs.flags[6:], [0, 0])
    # its water went into a's layers below the table (0.35 m down), by thickness
    due = gs.due * G.DAY * area[LAYER_HRU]
    assert due.sum() == pytest.approx(2.0)
    np.testing.assert_allclose(due[:6] / 2.0, [0, 0, 0, 0.10 / 0.40, 0.15 / 0.40, 0.15 / 0.40])
    # day 1: flags and the inflow set; the held layers take 3 x 0.01 m x 100 m2 = 3 m3 of room, all of it
    assert gs.before_swat() == pytest.approx(2.0)
    np.testing.assert_array_equal(sw.values[G.FLAG], gs.flags)
    assert gs.held_count() == 3
    sw.room[:] = 0.001  # only 0.3 m3 of room in a's held layers: 1.7 m3 isn't taken
    sw.day(area)
    got = gs.received()
    assert got + gs.not_taken == pytest.approx(2.0)
    assert gs.not_taken == pytest.approx(1.7)
    np.testing.assert_allclose(gs.returned, [1.7, 0.0])  # back on cell 0, where a's water came from
    # a restart carries the flags and the rows the due came from
    st = gs.state()
    again = G.GroundwaterToSoil(sw, m, area, "water_table", "none", "soil", held_lateral="swatplus", land=land,
                                saved=st, due=gs.due, n_cells=2)
    np.testing.assert_array_equal(again.flags, gs.flags)
    np.testing.assert_array_equal(again.rows, gs.rows)


def test_a_stream_cell_keeping_the_land_surface_has_its_drain_a_profile_down():
    m, sw, mf = stand_ins()
    mf.top = np.array([99.25, 100.0])  # cell 1, under b (1.0 m), keeps its land surface
    gs = G.GroundwaterToSoil(sw, m, np.array([100.0, 100.0]), "water_table", "none", "soil",
                             land_cells=np.array([100.0, 100.0]), n_cells=2)
    gs.attach(mf)
    np.testing.assert_allclose(gs.elev0, [99.25, 99.0])
    assert "1 cells keep the land surface" in gs.notes[0]
    gs.before_solve()
    mf.h = np.array([98.0, 99.5])  # b's table 0.5 m above its profile's bottom, 0.5 m down
    mf.out = np.array([0.0, 1.0])
    gs.after_solve()
    np.testing.assert_array_equal(gs.flags[6:], [0, 1])  # b's layers: 0-0.2 and 0.2-1.0 (middle 0.6)


def test_drain_warns_on_a_model_built_for_water_table():
    m, sw, mf = stand_ins()
    gs = G.GroundwaterToSoil(sw, m, np.array([100.0, 100.0]), "drain", "none", "soil", land=np.array([100.0, 100.0]))
    gs.attach(mf)
    assert "looks like the soil profile's bottom under 100%" in gs.notes[0]


def test_the_link_refuses_a_top_off_the_profiles_bottom():
    m, sw, mf = stand_ins()
    gs = G.GroundwaterToSoil(sw, m, np.array([100.0, 100.0]), "water_table", "none", "soil",
                             land=np.array([100.0, 103.0]), n_cells=2)
    with pytest.raises(ValueError, match="1 HRUs are not: b"):
        gs.attach(mf)


def test_the_config_takes_water_table_and_its_settings(tmp_path):
    from basinledger import config as C
    base = ('name = "x"\n[swatplus]\nproject = "p"\n[modflow6]\nbuilder = "m:b"\n'
            '[recharge]\nsource = "hru"\nmap = "all:hru0001"\n[exchange]\nmap = "all:cha01"\n')

    def load(extra):
        p = tmp_path / "c.toml"
        p.write_text(base + "\n[groundwater_to_soil]\n" + extra)
        return C.load(p)

    gs = load('rule = "water_table"\ncap = "none"\npackage = "soil"\nheld_lateral = "swatplus"\ntop_tolerance = 0.5\n')
    assert gs.groundwater_to_soil == {"rule": "water_table", "cap": "none", "package": "soil",
                                      "held_lateral": "swatplus", "top_tolerance": 0.5}
    gs = load('rule = "water_table"\ncap = "none"\npackage = "soil"\nland_surface = "land.csv"\n')
    assert gs.groundwater_to_soil["land_surface"] == tmp_path / "land.csv"
    for extra, match in [
        ('rule = "water_table"\ncap = "room"\npackage = "soil"\n', 'takes cap = "none"'),
        ('rule = "water_table"\ncap = "none"\npackage = "soil"\nheld_lateral = "some"\n', "held_lateral is one of"),
        ('rule = "water_table"\ncap = "none"\npackage = "soil"\ntop_tolerance = 0\n', "top_tolerance"),
        ('rule = "drain"\ncap = "none"\npackage = "soil"\nheld_lateral = "none"\n', "water_table's"),
    ]:
        with pytest.raises(ValueError, match=match):
            load(extra)
    p = tmp_path / "t.toml"
    p.write_text(base + '\n[transport.nitrate]\n\n[groundwater_to_soil]\nrule = "water_table"\ncap = "none"\n'
                 'package = "soil"\n')
    with pytest.raises(ValueError, match=r"not built with \[transport\]"):
        C.load(p)


# --- the land surface per map row (a land_surface file per row; docs/design.md, The land surface per map row) ---

def two_pieces():
    """HRUs a (profile 0.75 m) and b (1.0 m) on one 200 m2 cell, 100 m2 each: a on a rise (land 100.5 m), b in a
    hollow (99.5 m). The cell's land is 100.0 m, its profile's bottom 100 - (0.75 + 1.0) / 2 = 99.125 m, MODFLOW 6's
    top; a's piece bottoms at 99.75 m, b's at 98.5 m."""
    m = AreaMap(["a", "b"], np.array([0, 1]), np.array([0, 0]), np.array([100.0, 100.0]))
    return m, np.array([100.5, 99.5]), StandInSwat(), StandInModflow([99.125])


def test_each_piece_has_its_own_table_and_each_hru_its_pieces():
    m = AreaMap(["a", "b"], np.array([0, 0, 1]), np.array([0, 1, 1]), np.array([25.0, 75.0, 50.0]))
    bottom = np.array([99.0, 100.0, 98.0])  # a's two pieces, and b's
    z = G.row_table_height(m, np.array([99.5, 100.5]), bottom)
    np.testing.assert_allclose(z, [0.5, 0.5, 2.5])
    h = G.hru_table_height(m, np.array([99.5, 100.5]), np.zeros(2), 2, row_bottom=bottom)
    np.testing.assert_allclose(h, [0.5, 2.5])
    assert G.row_table_height(m, np.array([1e30, 100.5]), bottom)[0] < -1e3  # a dry cell holds nothing
    # with every piece at its cell's top, the per-cell table
    top = np.array([99.0, 100.0])
    np.testing.assert_array_equal(G.hru_table_height(m, np.array([99.5, 100.7]), top, 2, row_bottom=top[m.node]),
                                  G.hru_table_height(m, np.array([99.5, 100.7]), top, 2))


def test_the_drain_water_goes_to_the_pieces_with_soil_below_the_table():
    m = AreaMap(["a", "b", "c"], np.array([0, 1, 2]), np.array([0, 0, 1]), np.array([100.0, 300.0, 50.0]))
    # cell 0: a's piece above its table, b's 0.5 m into its soil; cell 1: c's piece above it, nowhere wet
    s = G.table_shares(m, np.array([-0.2, 0.5, -1.0]), 2)
    np.testing.assert_allclose(s, [0.0, 1.0, 1.0])  # cell 1 falls back to ground (it gives nothing anyway)
    s = G.table_shares(m, np.array([0.25, 0.5, 0.3]), 2)
    np.testing.assert_allclose(s, [25.0 / 175.0, 150.0 / 175.0, 1.0])  # by ground x height in the soil
    # at one height, by ground: cell_shares
    np.testing.assert_allclose(G.table_shares(m, np.array([0.4, 0.4, 0.4]), 2), G.cell_shares(m, 2))


def test_the_cells_land_is_its_pieces_and_a_file_per_row_carries_it(tmp_path):
    m = AreaMap(["a", "b"], np.array([0, 1, 1]), np.array([0, 0, 5]), np.array([100.0, 300.0, 50.0]))
    land = G.cell_land(m, np.array([101.0, 99.0, 50.0]), 6)
    np.testing.assert_allclose(land[[0, 5]], [99.5, 50.0])
    assert np.isnan(land[1])
    node_of = lambda layer, row, col: 3 * row + col  # noqa: E731
    f = tmp_path / "land.csv"  # the map's rows, in its order, by cell (the study's form) or by node
    f.write_text("object,layer,row,col,land_m\na,0,0,0,101.0\nb,0,0,0,99.0\nb,0,1,2,50.0\n")
    assert G.is_per_row(f)
    np.testing.assert_array_equal(G.read_row_land(f, m, node_of), [101.0, 99.0, 50.0])
    g = tmp_path / "map.csv"  # the recharge map's own table, with a land_m column
    g.write_text("object,node,area_m2,land_m\na,0,100,101.0\nb,0,300,99.0\nb,5,50,50.0\n")
    np.testing.assert_array_equal(G.read_row_land(g, m, node_of), [101.0, 99.0, 50.0])
    cells = tmp_path / "cells.csv"
    cells.write_text("layer,row,col,land_m\n0,0,0,100.0\n")
    assert not G.is_per_row(cells)
    for text, match in [
        ("object,node,land_m\na,0,101.0\nb,0,99.0\n", "2 rows, the recharge map has 3"),
        ("object,node,land_m\na,0,101.0\nb,0,99.0\nb,5,50.0\na,0,1.0\n", "more rows than"),
        ("object,node,land_m\na,0,101.0\nb,5,50.0\nb,0,99.0\n", "row 2 .* isn't the recharge map's"),
        ("object,node,land_m\na,0,101.0\nb,0,\nb,5,50.0\n", "row 2 has no land_m"),
    ]:
        f.write_text(text)
        with pytest.raises(ValueError, match=match):
            G.read_row_land(f, m, node_of)


def test_a_hollow_is_held_while_its_cell_is_below_the_top_and_takes_the_drains_water():
    m, land, sw, mf = two_pieces()
    area = np.array([100.0, 100.0])
    gs = G.GroundwaterToSoil(sw, m, area, "water_table", "none", "soil", n_cells=1, land_rows=land)
    gs.attach(mf)
    np.testing.assert_allclose(gs.elev0, [99.125])  # the drain per cell, at its top
    np.testing.assert_allclose(gs.row_bottom, [99.75, 98.5])
    gs.before_swat()
    sw.day(area)
    gs.received()
    # the head 99.0 m, below the cell's top: the drain gives nothing, but b's table is 0.5 m into its soil (0.5 m
    # down, above its lower layer's middle at 0.6 m): that layer is held, and fills from rain alone
    gs.before_solve()
    mf.h, mf.out = np.array([99.0]), np.array([0.0])
    assert gs.after_solve() == (0.0, 0.0)
    np.testing.assert_array_equal(gs.flags, [0, 0, 0, 0, 0, 0, 0, 1])
    assert gs.due.sum() == 0.0
    # the per-cell table at the same head holds nothing: the cell's table is below its top
    cell = G.GroundwaterToSoil(StandInSwat(), m, area, "water_table", "none", "soil", n_cells=1,
                               land_cells=np.array([100.0]))
    cell.attach(StandInModflow([99.125]))
    cell.before_solve()
    cell.mf6.h, cell.mf6.out = np.array([99.0]), np.array([0.0])
    cell.after_solve()
    assert cell.flags.sum() == 0
    # the head 99.5 m, above the top: the drain's 2 m3 all to b, whose table is at its surface; none to a, whose
    # piece is above its table (Andrew's A)
    gs.before_solve()
    mf.h, mf.out = np.array([99.5]), np.array([2.0])
    assert gs.after_solve() == (2.0, 0.0)
    np.testing.assert_allclose(gs.rows, [0.0, 2.0])
    due = gs.due * G.DAY * area[LAYER_HRU]
    assert due[:6].sum() == 0.0 and due[6:].sum() == pytest.approx(2.0)
    np.testing.assert_allclose(due[6:] / 2.0, [0.2, 0.8])  # b's layers by thickness below its table
    np.testing.assert_array_equal(gs.flags, [0, 0, 0, 0, 0, 0, 1, 1])
    # and it closes: what b's held layers can't take goes back to the cell
    gs.before_swat()
    sw.room[:] = 0.005
    sw.day(area)
    got = gs.received()
    assert got + gs.not_taken == pytest.approx(2.0)
    np.testing.assert_allclose(gs.returned, [gs.not_taken])


def test_pieces_at_their_cells_land_are_the_per_cell_link():
    def run(**land):
        m, sw, mf = stand_ins()
        area = np.array([100.0, 100.0])
        gs = G.GroundwaterToSoil(sw, m, area, "water_table", "none", "soil", n_cells=2, **land)
        gs.attach(mf)
        out = []
        for h, q in [([99.25 + 0.4, 98.0], [2.0, 0.0]), ([99.25 + 0.1, 99.6], [0.5, 1.5]), ([98.0, 98.5], [0, 0])]:
            gs.before_swat()
            sw.day(area)
            gs.received()
            gs.before_solve()
            mf.h, mf.out = np.array(h), np.array(q, float)
            gs.after_solve()
            out.append((gs.flags.copy(), gs.due.copy(), gs.rows.copy()))
        return out
    for (f1, d1, r1), (f2, d2, r2) in zip(run(land_rows=np.array([100.0, 100.0])),
                                          run(land_cells=np.array([100.0, 100.0]))):
        np.testing.assert_array_equal(f1, f2)
        np.testing.assert_array_equal(d1, d2)
        np.testing.assert_array_equal(r1, r2)


def test_a_land_surface_per_row_is_water_tables_in_place_of_one_per_cell():
    m, land, sw, _ = two_pieces()
    with pytest.raises(ValueError, match="in place of one per cell"):
        G.GroundwaterToSoil(sw, m, np.array([100.0, 100.0]), "water_table", "none", "soil", land_rows=land,
                            land_cells=np.array([100.0]))
    with pytest.raises(ValueError, match="water_table's"):
        G.GroundwaterToSoil(sw, m, np.array([100.0, 100.0]), "drain", "none", "soil", land_rows=land)


# --- on the models: Osu_1hru's one HRU on one MODFLOW 6 cell, its head raised and lowered ---------------------

import warnings  # noqa: E402

LAND = 350.96629   # m: Osu_1hru's hru.con elev
PROFILE = 0.75     # m: its profile's depth
TOP = LAND - PROFILE
DAYS = 150


def need_water_table():
    import shutil
    import tempfile
    from swatplus_bmi import SwatPlusBmi
    from basinledger.config import swatplus_refdata
    with tempfile.TemporaryDirectory() as d:
        shutil.copytree(swatplus_refdata() / "Osu_1hru", f"{d}/p")
        m = SwatPlusBmi()
        m.initialize(f"{d}/p")
        names = set(m.get_input_var_names()) | set(m.get_output_var_names())
        m.finalize()
    missing = [v for v in (G.FLAG, G.NOT_TAKEN, G.SATURATION_EXCESS) if v not in names]
    if missing:
        pytest.skip(f"swatplus-bmi has no held soil layers ({', '.join(missing)})")


def stages(kind):
    """The general head boundary's stage, m, day by day: the table below the profile, risen to the surface and
    above it (every layer held), down to mid-profile, and below the profile again; or held mid-layer, or below."""
    if kind == "below":
        return [TOP - 1.0] * (DAYS + 1)
    if kind == "mid_layer":  # 0.40 m down, inside layer 4 (0.30-0.45 m) below its middle
        return [TOP + PROFILE - 0.40] * (DAYS + 1)
    if kind == "mid_profile":  # the stage 0.30 m down; the stiff drain holds the head about 0.34 m down: layers 4-6
        return [TOP + PROFILE - 0.30] * (DAYS + 1)  # held (4's middle 0.375 m), 1-3 not (3's bottom 0.30 m)
    up = np.interp(np.arange(DAYS + 1), [0, 15, 35, 70, 85, 110, 125, DAYS],
                   [TOP - 0.5, TOP - 0.5, TOP + 0.9, TOP + 0.9, TOP + 0.4, TOP + 0.4, TOP - 0.6, TOP - 0.6])
    return [float(x) for x in up]


def not_ponded():
    """Osu_1hru with its HRU a forest with no management (hru-data.hru's surf_stor null, lu_mgt frst_tems_lum), not
    the irrigated rice paddy it ships as: a paddy's saturation excess and seepage go through its ponded water, which
    the column's land balance would have to carry (swatplus-bmi's note, Nov 27, 2026), and its irrigation would run
    off a held profile. Written once under tests/runs."""
    import shutil
    from basinledger.config import swatplus_refdata
    from conftest import RUNS
    p = RUNS / "osu-forest"
    if not (p / "hru-data.hru").exists():
        shutil.rmtree(p, ignore_errors=True)
        shutil.copytree(swatplus_refdata() / "Osu_1hru", p)
        f = p / "hru-data.hru"
        lines = f.read_text().splitlines()
        head = lines[1].split()
        i, lu = head.index("surf_stor"), head.index("lu_mgt")
        for k in range(2, len(lines)):
            cols = lines[k].split()
            if len(cols) > i:
                cols[i], cols[lu] = "null", "frst_tems_lum"
                lines[k] = "  ".join(cols)
        f.write_text("\n".join(lines) + "\n")
    return p


def column_config(kind="cycle", link=True, held_lateral="none", ponded=False):
    from basinledger import config as C
    from conftest import ROOT, RUNS
    cfg = C.load(ROOT / "examples" / "synthetic" / "two_way.toml")
    cfg.run_folder = RUNS
    if not ponded:
        cfg.swat_project = not_ponded()
    cfg.recharge_source, cfg.recharge_map = "hru", "all:hru0001"
    cfg.mf6_builder = cfg.mf6_builder.replace("build_small", "build_column")
    from basinledger.coupling import swatplus as S
    first = S._period(cfg.swat_project)[0]  # the run's first day: a restart's model takes the stages up from it
    cfg.mf6_builder_args = {"top": TOP, "stages": stages(kind), "stages_from": first.isoformat()}
    cfg.wells, cfg.users = [], []
    cfg.outputs = {"balance": "cha01"}
    cfg.groundwater_to_soil = ({"rule": "water_table", "cap": "none", "package": "soil",
                                "held_lateral": held_lateral} if link else None)
    return cfg


def column_run(cfg, name, days=DAYS, restart=None, save_at=None):
    """The run day by day, with SWAT+'s layers and the flags read at each pause."""
    from basinledger.driver import Basin
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run=name, users=False, days=days, restart=restart)
    rec = {k: [] for k in ("water", "ul", "flags", "perc")}
    stop = save_at if save_at else days
    while b.day < stop:
        b.run(1)
        rec["water"].append(b.sw.get(G.WATER).astype(float))
        rec["ul"].append(b.sw.get(G.SATURATION).astype(float))
        rec["flags"].append(b.sw.get(G.FLAG).astype(float) if b.sw.has(G.FLAG) else np.zeros(6))
        rec["perc"].append(float(b.sw.get("soil_water__percolation_rate")[0]))
    if save_at:
        b.save(b.folder / "state")
    b.write_outputs()
    out = {"table": b.ledger.table(), "folder": b.folder, "series": dict(b.series), "heads": b.mf6.head(),
           **{k: np.array(v) for k, v in rec.items()}}
    b.finalize()
    return out


@pytest.fixture(scope="module")
def column(swatplus):
    need_water_table()
    return column_run(column_config(), "water-table-column")


@pytest.mark.slow
def test_the_held_layers_fill_from_the_drain_and_hold(column):
    held = column["flags"] > 0.5
    assert held.any(axis=0).all()                 # every layer held at some point (the table above the surface)
    # a layer held two pauses running is full at the second: filled from the drain, not drained
    both = held[1:] & held[:-1]
    full = column["water"][1:] >= column["ul"][1:] * (1 - 1e-5)
    assert both.sum() > 50
    assert full[both].all()
    t = column["table"]
    assert t["gw_soil:sent"].sum() > 0
    np.testing.assert_allclose(t["gw_soil:due"], t["gw_soil:received"] + t["gw_soil:not_taken"], rtol=1e-6,
                               atol=1e-6)


@pytest.mark.slow
def test_held_layers_run_up_from_the_bottom_in_swatplus(column):
    assert np.all(np.diff(column["flags"], axis=1) >= 0)


@pytest.mark.slow
def test_a_full_profile_spills_and_what_no_layer_takes_returns_as_recharge(column):
    t = column["table"]
    every = (column["flags"] > 0.5).all(axis=1)
    assert every.sum() > 10
    assert t["soil:saturation_excess"].sum() > 0  # rain on a held profile ran off as saturation excess
    np.testing.assert_allclose(t["gw_soil:returned"], t["gw_soil:not_taken"], rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(t["recharge:received"], t["recharge:sent"] + t["gw_soil:returned"], rtol=1e-9,
                               atol=1e-6)


@pytest.mark.slow
def test_a_table_mid_profile_puts_no_groundwater_above_it(swatplus):
    """Andrew's choice (A), Nov 27, 2026: with layers held, drain water beyond their room isn't put above the table.
    It goes back to MODFLOW 6 as recharge, so a dry day spills nothing, however stiff the drain."""
    need_water_table()
    from basinledger.balance import _table
    r = column_run(column_config("mid_profile"), "water-table-mid-profile", days=90)
    t = r["table"]
    np.testing.assert_array_equal(r["flags"][5:], np.tile([0, 0, 0, 1, 1, 1], (85, 1)))
    head, rows = _table(r["folder"] / "swatplus" / "basin_wb_day.txt")
    dry = np.array([float(x[head.index("precip")]) == 0.0 for x in rows])
    late = np.arange(len(dry)) >= 15  # the held layers full
    assert (dry & late).sum() > 20
    np.testing.assert_array_equal(t["soil:saturation_excess"][dry & late], 0.0)
    assert t["gw_soil:not_taken"][late].min() > 0  # the drain's excess, every day
    np.testing.assert_allclose(t["gw_soil:returned"], t["gw_soil:not_taken"], rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(t["recharge:received"], t["recharge:sent"] + t["gw_soil:returned"], rtol=1e-9,
                               atol=1e-6)
    # the layers above the table never take groundwater: what went in is the held layers' room
    assert t["gw_soil:received"][late].max() < 0.05 * t["gw_soil:due"][late].min()


@pytest.mark.slow
def test_on_a_paddy_the_saturation_excess_goes_to_its_pond(swatplus):
    """Osu_1hru as it ships, a paddy: SWAT+ puts what a full profile spills into the ponded water, not the runoff,
    and the ledger books it there."""
    need_water_table()
    r = column_run(column_config(ponded=True), "water-table-paddy", days=90)
    t = r["table"]
    assert t["soil:saturation_excess_to_ponds"].sum() > 0
    np.testing.assert_array_equal(t["soil:saturation_excess"], 0.0)


@pytest.mark.slow
def test_a_released_layer_drains_as_recharge(column):
    held = column["flags"] > 0.5
    last = np.flatnonzero(held.any(axis=1))[-1]  # the last pause with any layer held
    after = column["water"][last + 1:last + 15]
    assert after.sum(axis=1)[-1] < after.sum(axis=1)[0]  # the profile drains
    assert column["perc"][last + 1:last + 15].sum() > 0  # through the bottom: recharge


@pytest.mark.slow
def test_the_column_closes(column):
    """MODFLOW 6's side and the link to the bit's tolerance; SWAT+'s land day by day from its own records, with the
    runoff before its surface lag (surq_cont), so the whole run's residual is what the lag still holds at the end."""
    from basinledger.balance import _table, basin_balance
    r = basin_balance(column["folder"], "cha01")
    assert abs(r["groundwater_residual"]) <= 1e-6 * (r["recharge_received"] + r["groundwater_to_soil_sent"])
    head, rows = _table(column["folder"] / "swatplus" / "basin_wb_day.txt")
    c = {k: np.array([float(x[head.index(k)]) for x in rows])
         for k in ("precip", "irr", "gwsoilq", "et", "surq_cont", "latq", "perc", "qtile", "sw_init", "sw_final",
                   "sno_init", "sno_final", "lagsurf")}
    day = (c["precip"] + c["irr"] + c["gwsoilq"] - c["et"] - c["surq_cont"] - c["latq"] - c["perc"] - c["qtile"]
           - (c["sw_final"] - c["sw_init"]) - (c["sno_final"] - c["sno_init"]))  # mm, to SWAT+'s printed digits
    assert c["gwsoilq"].sum() > 100 and np.abs(np.cumsum(day)).max() < 0.5
    assert r["land_residual"] / r["land_m2"] * 1000 == pytest.approx(c["lagsurf"][-1], abs=0.5)
    # SWAT+'s groundwater into the soil is the ledger's, and the not-taken and saturation excess are reported
    assert r["groundwater_into_soil"] == pytest.approx(column["table"]["gw_soil:received"].sum(), rel=1e-3)
    assert r["groundwater_to_soil_not_taken"] == pytest.approx(column["table"]["gw_soil:not_taken"].sum())
    assert 0 < r["runoff_saturation_excess"] <= r["runoff_to_channels"] + c["lagsurf"][-1] / 1000 * r["land_m2"]


@pytest.mark.slow
def test_a_table_held_mid_layer_doesnt_flap(swatplus):
    need_water_table()
    r = column_run(column_config("mid_layer"), "water-table-mid", days=60)
    f = r["flags"][1:]
    assert (f == f[0]).all()
    np.testing.assert_array_equal(f[0], [0, 0, 0, 0, 1, 1])


@pytest.mark.slow
def test_a_table_below_the_profile_is_the_run_without_the_link(swatplus):
    need_water_table()
    linked = column_run(column_config("below"), "water-table-below", days=90)
    alone = column_run(column_config("below", link=False), "water-table-below-alone", days=90)
    assert np.all(linked["table"]["gw_soil:sent"] == 0.0)
    assert np.all(linked["flags"] == 0.0)
    np.testing.assert_array_equal(linked["series"]["flow:cha01"], alone["series"]["flow:cha01"])
    np.testing.assert_array_equal(linked["heads"], alone["heads"])


@pytest.mark.slow
def test_a_restart_carries_the_flags_and_what_is_due(swatplus):
    need_water_table()
    cfg = column_config()
    whole = column_run(cfg, "water-table-whole", days=90)
    first = column_run(cfg, "water-table-first", days=90, save_at=45)
    rest = column_run(cfg, "water-table-rest", days=90, restart=first["folder"] / "state")
    assert (first["flags"][-1] > 0).any() and rest["table"]["gw_soil:due"][45] > 0
    np.testing.assert_array_equal(rest["table"]["gw_soil:sent"], whole["table"]["gw_soil:sent"])
    np.testing.assert_array_equal(rest["series"]["flow:cha01"], whole["series"]["flow:cha01"])
    np.testing.assert_array_equal(rest["heads"], whole["heads"])
    np.testing.assert_array_equal(rest["flags"], whole["flags"][45:])


@pytest.mark.slow
def test_lateral_flow_from_held_layers_is_an_option(swatplus, column):
    need_water_table()
    from basinledger.balance import _table
    own = column_run(column_config(held_lateral="swatplus"), "water-table-lateral")

    def latq(run):
        head, rows = _table(run["folder"] / "swatplus" / "basin_wb_day.txt")
        return sum(float(r[head.index("latq")]) for r in rows)
    assert latq(own) > latq(column)


@pytest.mark.slow
def test_a_land_surface_per_row_at_the_cells_land_is_the_run_without_it(swatplus, column, tmp_path):
    """A land surface file per map row read through the driver: the one piece at the cell's land surface is the
    column's run, to round-off (the cell's land comes back as its piece's mean)."""
    need_water_table()
    COLUMN_CELL = 100_000.0 ** 0.5  # examples/synthetic/model.py's, the column's one cell of 10 ha
    f, g = tmp_path / "map.csv", tmp_path / "land.csv"
    f.write_text(f"object,node,area_m2\nhru0001,0,{COLUMN_CELL ** 2!r}\n")
    g.write_text(f"object,node,land_m\nhru0001,0,{LAND!r}\n")
    cfg = column_config()
    cfg.recharge_map = str(f)
    cfg.groundwater_to_soil = dict(cfg.groundwater_to_soil, land_surface=g)
    r = column_run(cfg, "water-table-rows")
    import json
    assert json.loads((r["folder"] / "outputs" / "summary.json").read_text())["water_table"]["land_surface"] == "rows"
    np.testing.assert_array_equal(r["flags"], column["flags"])
    for k in ("gw_soil:sent", "gw_soil:not_taken", "soil:saturation_excess"):
        np.testing.assert_allclose(r["table"][k], column["table"][k], rtol=1e-9, atol=1e-6)
    np.testing.assert_allclose(r["series"]["flow:cha01"], column["series"]["flow:cha01"], rtol=1e-9)
