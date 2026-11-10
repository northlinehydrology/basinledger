"""Vertex grids (DISV): cells named as layer, icell2d through the coupling (coupling/modflow6.py, mapping.py).

- cells: a DISV cell to its reduced node, and a cell in the other grid's form, outside the grid or inactive refused;
  map rows and configured cells in either form (tier 1);
- the synthetic aquifer as a uniform 10 x 10 DISV grid, which is build_small's grid cell for cell, gives the DIS
  run's numbers through the whole coupling (tier 2);
- the coupled refined run (examples/synthetic/disv.toml): recharge and the river through icell2d maps, the well,
  the head output and a drawdown rule on a 40 m cell, chloride carried; MODFLOW 6's budget closes every day and
  matches the ledger, and the transport's mass closes (tier 1 and 2).
"""

from __future__ import annotations

import datetime as dt
import importlib.util

import flopy
import numpy as np
import pytest

from basinledger import config as C
from basinledger.coupling.modflow6 import Modflow6
from basinledger.mapping import AreaMap, cell_node
from conftest import ROOT, RUNS

SYNTHETIC = ROOT / "examples" / "synthetic"
DAYS = 365


def synthetic_model():
    spec = importlib.util.spec_from_file_location("synthetic_model", SYNTHETIC / "model.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ---------- cells ----------

def test_the_grid_is_a_quadtree_with_its_hanging_vertices():
    m = synthetic_model()
    cells = m.small_disv_grid()
    assert len(cells) == 28
    assert sum(c["side"] ** 2 for c in cells) == pytest.approx(1e5, rel=1e-12)  # Osu_1hru's 10 ha aquifer
    # the large cells beside the refined block list the vertex mid-edge: six of them
    assert sum(len(c["vertices"]) == 5 for c in cells) == 6
    assert m.small_disv_well() == (0, 8) and cells[8]["side"] < cells[0]["side"]


def test_disv_cells_to_nodes_and_what_is_refused(libmf6):
    m = synthetic_model()
    ws = RUNS / "disv-cells"
    m.build_small_disv(ws, 2, dt.date(2010, 1, 1))
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws, verbosity_level=0)
    idomain = np.ones((1, 28), int)
    idomain[0, 27] = 0  # the south-east cell inactive
    sim.get_model("gwf").disv.idomain = idomain
    sim.write_simulation(silent=True)
    mf = Modflow6(ws)
    try:
        assert mf.grid == "disv" and mf.shape == (1, 28) and mf.nodes == 27
        assert mf.node(0, 8) == 8 and mf.node(0, 26) == 26
        with pytest.raises(ValueError, match="DISV, whose cells are layer, icell2d"):
            mf.node(0, 5, 6)
        with pytest.raises(ValueError, match="outside the grid"):
            mf.node(0, 28)
        with pytest.raises(ValueError, match="inactive"):
            mf.node(0, 27)
        mf.prepare_day()
        assert mf.package_nodes("wel").tolist() == [8]
    finally:
        mf.finalize()


def test_a_dis_model_refuses_disv_cells(libmf6):
    m = synthetic_model()
    ws = RUNS / "dis-cells"
    m.build_small(ws, 2, dt.date(2010, 1, 1))
    mf = Modflow6(ws)
    try:
        assert mf.grid == "dis" and mf.node(*m.SMALL_WELL) == 56
        with pytest.raises(ValueError, match="DIS, whose cells are layer, row, col"):
            mf.node(0, 56)
        with pytest.raises(ValueError, match="outside the grid"):
            mf.node(0, 10, 0)  # would have wrapped to (0, 0, 0) of the next layer's numbering before
    finally:
        mf.finalize()


def test_map_rows_in_either_form():
    seen = []

    def node_of(*cell):
        seen.append(cell)
        if len(cell) != 2:
            raise ValueError("this model's grid is DISV")
        return 100 + cell[1]

    assert cell_node({"layer": "0", "icell2d": "7"}, node_of) == 107
    assert cell_node({"node": "12"}, node_of) == 12
    assert cell_node({"node": "12", "row": "", "col": ""}, node_of) == 12  # empty columns are not given
    with pytest.raises(ValueError, match="m.csv: this model's grid is DISV"):
        cell_node({"layer": "0", "row": "5", "col": "6"}, node_of, "m.csv")
    with pytest.raises(ValueError, match="gives \\['icell2d', 'layer', 'node'\\]"):
        cell_node({"layer": "0", "icell2d": "7", "node": "3"}, node_of)
    with pytest.raises(ValueError, match="no cell"):
        cell_node({"object": "aqu001"}, node_of)
    assert seen == [(0, 7), (0, 5, 6)]


def test_an_icell2d_map_read_against_the_model(libmf6, tmp_path):
    m = synthetic_model()
    ws = RUNS / "disv-map"
    m.build_small_disv(ws, 2, dt.date(2010, 1, 1))
    m.write_small_disv_maps(tmp_path)
    mf = Modflow6(ws)
    try:
        amap = AreaMap.from_csv(tmp_path / "disv_recharge.csv", ["aqu001"], mf.node)
        assert amap.node.tolist() == list(range(28))
        np.testing.assert_allclose(amap.area, mf.cell_area, rtol=1e-9)  # MODFLOW 6 reads vertices to ~10 digits
        amap.check(mf.cell_area, np.array([1e5]), 1e-9)
        (tmp_path / "rowcol.csv").write_text("object,layer,row,col,area_m2\naqu001,0,5,6,100.0\n")
        with pytest.raises(ValueError, match="rowcol.csv: cell \\(0, 5, 6\\): this model's grid is DISV"):
            AreaMap.from_csv(tmp_path / "rowcol.csv", ["aqu001"], mf.node)
    finally:
        mf.finalize()


def test_the_committed_maps_are_the_grids():
    m = synthetic_model()
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as d:
        m.write_small_disv_maps(d)
        for name in ("disv_recharge.csv", "disv_river.csv"):
            assert (Path(d) / name).read_text() == (SYNTHETIC / "maps" / name).read_text(), name


@pytest.mark.parametrize("cell", ["[0]", "[0, 1, 2, 3]", "[0, 1.5]", '"0,5"'])
def test_configured_cells_are_two_or_three_integers(tmp_path, cell):
    text = (SYNTHETIC / "disv.toml").read_text()
    q = tmp_path / "disv.toml"
    q.write_text(text.replace("cell = [0, 8]                    #", f"cell = {cell}  #", 1))
    with pytest.raises(ValueError, match=r"well 'pw1' cell must be \[layer, row, col\] \(DIS\) or \[layer, icell2d\]"):
        C.load(q)


def test_the_example_config_reads_its_cells():
    cfg = C.load(SYNTHETIC / "disv.toml")
    assert cfg.wells[0].cell == (0, 8)
    assert cfg.rules[0]["cell"] == [0, 8] and cfg.outputs["heads"][0]["cell"] == [0, 8]


# ---------- coupled runs ----------

def run(cfg, name, days=DAYS):
    from basinledger.driver import Basin
    cfg.run_folder = RUNS
    b = Basin(cfg, run=name, days=days)
    b.run()
    out = dict(table=b.ledger.table(), series={k: np.array(v) for k, v in b.series.items()}, ws=b.mf6.ws)
    b.finalize()
    return out


@pytest.fixture(scope="module")
def refined(swatplus):
    return run(C.load(SYNTHETIC / "disv.toml"), "disv")


def budget(ws, texts=("STO-SS", "STO-SY", "RCH", "WEL", "RIV")):
    cbc = flopy.utils.CellBudgetFile(ws / "gwf.cbc", precision="double")
    return {t: np.array([float(d["q"].sum()) if d.dtype.names else float(np.sum(d)) for d in cbc.get_data(text=t)])
            for t in texts}


@pytest.mark.slow
def test_refined_budget_closes_and_matches_the_ledger(refined):
    t, terms = refined["table"], budget(refined["ws"])
    net = sum(terms.values())
    gross = sum(np.abs(v) for v in terms.values())
    assert len(net) == DAYS and np.max(np.abs(net) / gross) < 1e-7
    np.testing.assert_allclose(terms["RCH"], t["recharge:sent"], rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(-terms["RIV"], t["exchange:cha01:to_stream"], rtol=1e-12)
    np.testing.assert_allclose(-terms["WEL"], t["route:well:pw1:delivered"], rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(t["recharge:sent"], t["recharge:received"], rtol=1e-12, atol=1e-9)
    assert (terms["RIV"] > 0).sum() > 10 and (terms["RIV"] < 0).sum() > 10  # the river lost and gained
    assert (t["route:well:pw1:delivered"] > 0).all()


@pytest.mark.slow
def test_refined_heads_well_and_rule(refined):
    t, s = refined["table"], refined["series"]
    hds = flopy.utils.HeadFile(refined["ws"] / "gwf.hds").get_alldata()  # (days, 1, 1, 28)
    np.testing.assert_allclose(s["head:pw1"], hds[:, 0, 0, 8], rtol=0, atol=1e-9)  # the output is icell2d 8's head
    active = t["rule:plant well limit:active"].astype(bool)
    assert 10 < active.sum() < DAYS - 10                         # it acted on some days, not all
    np.testing.assert_array_equal(active[1:], s["head:pw1"][:-1] < 87.9)   # the day after the head was below
    limited = t["route:well:pw1:delivered"][1:][active[1:]]
    assert (limited <= 0.0002 * 86400 * (1 + 1e-9)).all()


@pytest.mark.slow
def test_refined_transport_mass_closes(refined):
    t = refined["table"]
    k = "transport:chloride:"
    assert (np.abs(t[k + "closure_share"]) < 1e-9).all()
    np.testing.assert_allclose(t[k + "recharge:received"], t[k + "recharge:sent"], rtol=1e-9, atol=1e-12)
    total = (t[k + "recharge:received"] + t[k + "stream:cha01:leakage_received"] - t[k + "stream:cha01:to_stream"]
             - t[k + "well:pw1:withdrawn"] + t[k + "storage_release"]).sum()
    assert abs(total) < 1e-9 * t[k + "recharge:received"].sum()
    assert not any(a.startswith(k + "boundary:") for a in t)


@pytest.mark.slow
def test_uniform_disv_is_the_dis_run(swatplus):
    """build_small's grid as DISV, with its cells as icell2d, gives the DIS run through the whole coupling. MODFLOW 6
    computes each grid's geometry its own way, so the two agree closely but not bit for bit."""
    dis = C.load(SYNTHETIC / "two_way.toml")
    disv = C.load(SYNTHETIC / "two_way.toml")
    disv.mf6_builder = "model:build_small_disv"
    disv.mf6_builder_args = {"refine": False}
    disv.wells[0].cell = (0, 56)
    disv.outputs = {"heads": [{"name": "pw1", "cell": [0, 56]}]}
    a, b = run(dis, "disv-dis", 180), run(disv, "disv-uniform", 180)
    assert set(a["table"]) == set(b["table"])
    # MODFLOW 6 reads the vertices to ~10 digits, so cell areas differ by ~4e-10: ~1e-8 m3 a day in the map's terms
    off = [k for k in a["table"] if not np.allclose(b["table"][k], a["table"][k], rtol=1e-7, atol=1e-6)]
    assert not off, {k: float(np.max(np.abs(b["table"][k] - a["table"][k]))) for k in off}
    np.testing.assert_allclose(b["series"]["head:pw1"], a["series"]["head:pw1"], rtol=0, atol=1e-7)
