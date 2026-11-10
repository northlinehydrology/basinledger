"""The Plaster Creek groundwater model's building blocks, and (where its data are cached) a short coupled run.

The unit tests need the ``gis`` dependency group (shapely, rasterio, scipy). The coupled run needs the cached NextGen
rerun and the model arrays (examples/plaster-creek/groundwater/fetch.py, prepare.py, rerun_nextgen.py) and is slow.
"""

from __future__ import annotations

import json
import struct
import sys

import numpy as np
import pytest

shapely = pytest.importorskip("shapely")
pytest.importorskip("rasterio")
pytest.importorskip("scipy")

from shapely.geometry import LineString, box  # noqa: E402

from conftest import ROOT, RUNS  # noqa: E402

GW = ROOT / "examples" / "plaster-creek" / "groundwater"
sys.path.insert(0, str(GW))
import geometry as G  # noqa: E402
import prepare as P  # noqa: E402


def gpkg_blob(geom, srs=5070, envelope=True):
    flags = 0b00000011 if envelope else 0b00000001  # little endian; envelope xy (code 1) or none
    head = b"GP" + bytes([0, flags]) + struct.pack("<i", srs)
    env = struct.pack("<4d", *geom.bounds) if envelope else b""
    return head + env + geom.wkb


def test_geopackage_geometry():
    p = box(0, 0, 100, 50)
    for env in (True, False):
        g = G._gpkg_geometry(gpkg_blob(p, envelope=env))
        assert g.equals(p)
    with pytest.raises(ValueError):
        G._gpkg_geometry(b"XX" + p.wkb)


def test_maps_keep_area_and_length():
    import flopy
    mg = flopy.discretization.StructuredGrid(delr=np.full(8, 100.0), delc=np.full(6, 100.0), xoff=0, yoff=0, nlay=1,
                                             top=np.ones((6, 8)), botm=np.zeros((1, 6, 8)))
    gi = G.grid_intersect(mg)
    polys = {"cat-a": box(30, 20, 420, 330), "cat-b": box(420, 20, 790, 590)}
    cc = G.catchment_cells(gi, polys)
    for k, p in polys.items():
        assert sum(r[3] for r in cc if r[0] == k) == pytest.approx(p.area, rel=1e-12)
    line = {"wb-1": LineString([(10, 10), (790, 590)])}
    fc = G.flowpath_cells(gi, line)
    assert sum(r[3] for r in fc) == pytest.approx(line["wb-1"].length, rel=1e-12)
    x0, y0, nrow, ncol = G.extent(box(120, 130, 380, 290), 100, 50)
    assert (x0, y0, nrow, ncol) == (0.0, 0.0, 8, 10)  # x 20-480 and y 30-390, snapped to 50 m



def test_maps_on_a_vertex_grid(tmp_path):
    """GridIntersect's icell2d on a vertex grid gives rows of (object, icell2d, value), written as layer, icell2d."""
    import csv
    import flopy
    verts = [[i * 4 + j, j * 100.0, 300.0 - i * 100.0] for i in range(4) for j in range(4)]
    cell2d = [[r * 3 + c, c * 100 + 50.0, 250.0 - r * 100, 4, r * 4 + c, r * 4 + c + 1, (r + 1) * 4 + c + 1,
               (r + 1) * 4 + c] for r in range(3) for c in range(3)]
    mg = flopy.discretization.VertexGrid(vertices=verts, cell2d=cell2d, nlay=1, top=np.ones(9), botm=np.zeros((1, 9)))
    gi = G.grid_intersect(mg)
    p = box(30, 20, 270, 280)
    cc = G.catchment_cells(gi, {"cat-a": p})
    assert all(len(r) == 3 for r in cc) and sum(r[2] for r in cc) == pytest.approx(p.area, rel=1e-12)
    fc = G.flowpath_cells(gi, {"wb-1": LineString([(10, 10), (290, 290)])})
    assert {r[1] for r in fc} == {2, 4, 6}  # south-west to north-east, icell2d 0 at the north-west
    G.write_map(tmp_path / "m.csv", cc, "area_m2")
    rows = list(csv.DictReader(open(tmp_path / "m.csv")))
    assert list(rows[0]) == ["object", "layer", "icell2d", "area_m2"]
    assert sum(float(r["area_m2"]) for r in rows) == pytest.approx(p.area, rel=1e-12)


def test_wellogic_filter(tmp_path):
    def well(i, **p):
        base = dict(WELLID=str(i), AQ_CODE="D", SWL=10, SWL_FLAG=None, ELEVATION=700, ELEV_DEM=702,
                    CONST_DATE="2005-06-01 00:00:00", SCREEN_FRM=40, SCREEN_TO=50, WELL_DEPTH=50)
        base.update(p)
        return {"type": "Feature", "geometry": {"type": "Point", "coordinates": [1.0, 2.0]}, "properties": base}
    feats = [well(1), well(2, AQ_CODE="R"), well(3, SWL=0), well(4, SWL_FLAG="X"), well(5, ELEV_DEM=720)]
    (tmp_path / "wellogic.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    w = P.wells(tmp_path)
    assert [v["id"] for v in w] == ["1"]
    assert w[0]["head"] == pytest.approx(690 * 0.3048) and w[0]["screen_mid"] == pytest.approx(45 * 0.3048)


def tiny_arrays(path):
    nrow, ncol = 5, 6
    top = np.full((nrow, ncol), 250.0)
    botm = np.stack([top - 20, top - 40, top - 70])
    active = np.ones((nrow, ncol), bool)
    active[0, 0] = False
    gh = np.full((nrow, ncol), np.nan)
    gh[:, -1] = 245.0
    np.savez_compressed(path, x0=0.0, y0=0.0, nrow=nrow, ncol=ncol, cell=250.0, top=top, botm=botm,
                        k=np.full((3, nrow, ncol), 3.0), k33=np.full((3, nrow, ncol), 0.1),
                        riv_cells=np.array([[r, 0] for r in range(1, nrow)]), riv_bed=np.full(nrow - 1, 244.0),
                        riv_cond=np.full(nrow - 1, 100.0), riv_obj=np.array(["wb-1"] * (nrow - 1)), ghb_head=gh,
                        active=active, campus_well=np.array([1, 2, 3]), campus_well_distance=750.0)


def test_model_writer_runs_steady_and_transient(tmp_path, libmf6):
    import subprocess
    import flopy
    import model as M
    npz = tmp_path / "a.npz"
    tiny_arrays(npz)
    ws = tmp_path / "steady"
    M.build(ws, 1, None, steady=True, recharge=np.full((5, 6), 5e-4), depth={"wb-1": 0.3}, npz=npz)
    r = subprocess.run([str(libmf6.with_name("mf6"))], cwd=ws, capture_output=True, text=True)
    assert "Normal termination" in r.stdout
    h = flopy.utils.HeadFile(ws / "gwf.hds").get_data()
    assert np.all(h[:, 1:, :] < 250.1) and np.all(h[:, 1:, :] > 243.0)  # between the river and the land surface
    lst = flopy.utils.Mf6ListBudget(ws / "gwf.lst").get_incremental()
    assert abs(lst["PERCENT_DISCREPANCY"][-1]) < 0.01


CACHE_OK = None


def cached():
    from basinledger.config import cache_dir
    run = cache_dir() / "plaster-creek" / "nextgen-recharge" / "outputs" / "ngen"
    npz = cache_dir() / "plaster-creek" / "groundwater.npz"
    if not (run.is_dir() and npz.exists()):
        pytest.skip("Plaster Creek groundwater data not cached (groundwater/fetch.py, rerun_nextgen.py, prepare.py)")


@pytest.mark.slow
def test_plaster_creek_coupled_month(libmf6):
    cached()
    pytest.importorskip("waterledger_models.datacenter.bmi")
    from basinledger import config as C
    from basinledger.driver import Basin
    import model as M

    cfg = C.load(ROOT / "examples" / "plaster-creek" / "groundwater-wells.toml")
    cfg.run_folder = RUNS
    assert tuple(cfg.wells[0].cell) == M.campus_well()  # the config's well is the one prepare.py placed
    b = Basin(cfg, run="plaster-gw-month", days=30)
    assert b.status == "screening"  # groundwater isn't calibrated
    b.run()
    t = b.ledger.table()
    s = {k: np.array(v) for k, v in b.series.items()}
    b.finalize()
    np.testing.assert_allclose(t["recharge:received"], t["recharge:sent"], rtol=1e-12)
    np.testing.assert_allclose(t["route:well:campus:delivered"], t["user:campus:withdrawal"], rtol=1e-9)
    assert (s["head:campus-well"][1:] <= s["head:campus-well"][0] + 1e-6).sum() > 20  # pumping draws it down
    import flopy
    lst = flopy.utils.Mf6ListBudget(RUNS / "plaster-gw-month" / "mf6" / "gwf.lst").get_incremental()
    assert np.max(np.abs(lst["PERCENT_DISCREPANCY"])) < 0.1
