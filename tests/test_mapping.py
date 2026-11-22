"""Area maps move depth rates between SWAT+ objects and MODFLOW 6 cells without gaining or losing volume."""

from __future__ import annotations

import numpy as np
import pytest

from basinledger.mapping import AreaMap


def two_aquifers():
    # 4 cells of 100 m2; aquifer a covers cells 0, 1 and half of 2; aquifer b the other half of 2 and cell 3
    cell_area = np.full(4, 100.0)
    m = AreaMap(["a", "b"], np.array([0, 0, 0, 1, 1]), np.array([0, 1, 2, 2, 3]),
                np.array([100.0, 100.0, 50.0, 50.0, 100.0]))
    return m, cell_area


def test_volume_is_kept():
    m, cell_area = two_aquifers()
    rate = np.array([0.002, 0.005])  # m d-1
    cells = m.to_cells(rate, cell_area)
    assert np.sum(cells * cell_area) == pytest.approx(m.volume(rate), rel=1e-15)
    assert cells[2] == pytest.approx((0.002 * 50 + 0.005 * 50) / 100)
    assert m.mapped_area().tolist() == [250.0, 150.0]


def test_over_covered_cell_is_refused():
    m, cell_area = two_aquifers()
    m.area[2] = 80.0  # cell 2 now 130 m2 of 100
    with pytest.raises(ValueError, match="cell 2"):
        m.check(cell_area)


def test_area_must_match_the_source_model_when_asked():
    m, cell_area = two_aquifers()
    m.check(cell_area, np.array([250.0, 150.0]), rtol=1e-3)
    with pytest.raises(ValueError, match="b: mapped area"):
        m.check(cell_area, np.array([250.0, 200.0]), rtol=1e-3)


def test_from_csv(tmp_path):
    f = tmp_path / "map.csv"
    f.write_text("object,layer,row,col,area_m2\na,0,0,1,100\nb,0,1,0,40\n")
    m = AreaMap.from_csv(f, ["a", "b"], lambda l, r, c: r * 2 + c)
    assert m.node.tolist() == [1, 2] and m.obj.tolist() == [0, 1]
    with pytest.raises(KeyError):
        AreaMap.from_csv(f, ["a"], lambda l, r, c: 0)


def test_from_csv_finds_each_cell_once(tmp_path):
    """A map from pixels repeats its cells: each distinct cell's node is looked up once, and every row keeps its own
    object, node and area; a cell the grid refuses is still refused."""
    f = tmp_path / "map.csv"
    f.write_text("object,layer,row,col,area_m2\na,0,0,1,100\nb,0,0,1,40\na,0,1,0,60\nb,0,0,1,5\n")
    calls = []

    def node_of(l, r, c):
        calls.append((l, r, c))
        return r * 2 + c

    m = AreaMap.from_csv(f, ["a", "b"], node_of)
    assert sorted(calls) == [(0, 0, 1), (0, 1, 0)]
    assert m.obj.tolist() == [0, 1, 0, 1] and m.node.tolist() == [1, 1, 2, 1]
    assert m.area.tolist() == [100.0, 40.0, 60.0, 5.0]

    def inactive(l, r, c):
        raise ValueError(f"cell {(l, r, c)} is inactive")

    with pytest.raises(ValueError, match="inactive"):
        AreaMap.from_csv(f, ["a", "b"], inactive)


def test_cells_back_to_objects():
    m, cell_area = two_aquifers()
    drawdown = np.array([1.0, 2.0, 4.0, 8.0])
    # a: (100*1 + 100*2 + 50*4) / 250; b: (50*4 + 100*8) / 150
    np.testing.assert_allclose(m.to_objects(drawdown), [2.0, 1000 / 150])
    assert m.outside(6).tolist() == [4, 5]
