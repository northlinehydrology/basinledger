"""Recharge HRU by HRU (``[recharge] source = "hru"``): each HRU's percolation placed on the cells under its ground.

Units: the pixel overlay splits pixels exactly among cells, and the conserving map keeps each HRU's volume. On the
models (needs swatplus-bmi 0.6.0's ``soil_water__percolation_rate``, skipped otherwise):

- Osu_1hru (one HRU recharging one aquifer): the volume sent equals SWAT+'s own aquifer recharge every day (the
  driver holds ``recharge:made`` to ``recharge:swat_aquifers``) and the aquifer path's volume, and MODFLOW 6 takes it
  whole, its budget closing;
- Ames_sub1 (12 HRUs, one year, each given its own known percolation): each HRU's percolation times its area reaches
  its own block of cells and no other, every day.
"""

from __future__ import annotations

import csv
import shutil
import warnings

import flopy
import numpy as np
import pytest

from basinledger.mapping import AreaMap, pixel_rows


# --- units ---------------------------------------------------------------------------------------------------------

def test_a_pixel_across_four_cells_is_split_exactly():
    # cells of 50 m from x = 0 and y = 100 down; a 30 m pixel with its lower-left corner at (40, 30)
    o, r, c, a = pixel_rows(np.array([0]), np.array([1.0]), np.array([40.0]), np.array([30.0]), 30.0,
                            np.array([0.0, 50.0, 100.0]), np.array([100.0, 50.0, 0.0]))
    got = {(int(ri), int(ci)): float(ai) for ri, ci, ai in zip(r, c, a)}
    assert got == {(0, 0): 10.0 * 10.0, (0, 1): 20.0 * 10.0, (1, 0): 10.0 * 20.0, (1, 1): 20.0 * 20.0}


def test_pixels_keep_their_area_and_never_overfill_a_cell():
    rng = np.random.default_rng(1)
    edges_x, edges_y = np.arange(0.0, 501.0, 50.0), np.arange(500.0, -1.0, -50.0)
    gx, gy = np.meshgrid(np.arange(0.0, 480.0, 30.0), np.arange(0.0, 480.0, 30.0))  # 16 x 16 pixels inside the grid
    x0, y0 = gx.ravel() + 5.0, gy.ravel() + 7.0
    # each pixel held by one of three objects, and some pixels folded into two (weights adding to 1)
    obj = rng.integers(0, 3, len(x0))
    folded = rng.random(len(x0)) < 0.2
    w = np.where(folded, 0.4, 1.0)
    o2, x2, y2, w2 = (obj + 1) % 3, x0[folded], y0[folded], np.full(folded.sum(), 0.6)
    o, r, c, a = pixel_rows(np.r_[obj, o2[folded]], np.r_[w, w2], np.r_[x0, x2], np.r_[y0, y2], 30.0, edges_x, edges_y)
    want = np.bincount(np.r_[obj, o2[folded]], weights=np.r_[w, w2] * 900.0, minlength=3)
    np.testing.assert_allclose(np.bincount(o, weights=a, minlength=3), want, rtol=1e-12)
    cover = np.zeros((10, 10))
    np.add.at(cover, (r, c), a)
    assert cover.max() <= 2500.0 * (1 + 1e-12)
    assert abs(cover.sum() - len(x0) * 900.0) < 1e-6


def test_conserving_map_places_each_objects_own_volume():
    m = AreaMap(["h1", "h2"], np.array([0, 0, 1]), np.array([0, 1, 1]), np.array([300.0, 100.0, 500.0]))
    source = np.array([1000.0, 250.0])  # each HRU's area in SWAT+, m2
    c = m.conserving(source)
    np.testing.assert_allclose(c.mapped_area(), source, rtol=1e-15)
    rate = np.array([0.002, 0.010])  # m d-1
    cells = c.to_cells(rate, np.array([2500.0, 2500.0]))
    assert np.sum(cells * 2500.0) == pytest.approx(np.sum(rate * source), rel=1e-15)
    # h1 takes its share by its ground: three quarters to cell 0
    assert cells[0] * 2500.0 == pytest.approx(0.002 * 1000.0 * 0.75, rel=1e-15)


def test_conserving_map_refuses_an_hru_with_no_cells():
    m = AreaMap(["h1", "h2"], np.array([0]), np.array([0]), np.array([300.0]))
    with pytest.raises(ValueError, match="no cells"):
        m.conserving(np.array([1000.0, 250.0]))


# --- on the models -------------------------------------------------------------------------------------------------

def need_hru_percolation():
    from swatplus_bmi import SwatPlusBmi
    from basinledger.config import swatplus_refdata
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        shutil.copytree(swatplus_refdata() / "Osu_1hru", f"{d}/p")
        m = SwatPlusBmi()
        m.initialize(f"{d}/p")
        has = "soil_water__percolation_rate" in m.get_output_var_names()
        m.finalize()
    if not has:
        pytest.skip("swatplus-bmi has no soil_water__percolation_rate (0.6.0)")


def osu_one_year(source, recharge_map, name):
    from basinledger.driver import Basin
    from conftest import synthetic_config
    cfg = synthetic_config(recharge_source=source, recharge_map=recharge_map)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run=name, users=False, days=365)
    b.run()
    table, ws = b.ledger.table(), b.mf6.ws
    b.finalize()
    return table, ws


@pytest.fixture(scope="module")
def osu(swatplus):
    need_hru_percolation()
    return osu_one_year("hru", "all:hru0001", "hru-osu"), osu_one_year("aquifer", "all:aqu001", "hru-osu-aquifer")


@pytest.mark.slow
def test_osu_hru_recharge_is_swats_aquifer_recharge(osu):
    (hru, ws), (aqu, _) = osu
    assert hru["recharge:made"].sum() > 0
    # what the HRU percolated, times its area, is what SWAT+ put in its aquifer, every day
    np.testing.assert_allclose(hru["recharge:made"], hru["recharge:swat_aquifers"], rtol=1e-5, atol=1e-9)
    np.testing.assert_allclose(hru["recharge:made"], aqu["recharge:made"], rtol=1e-5, atol=1e-9)
    # the map placed it and added nothing; MODFLOW 6 took it whole
    np.testing.assert_allclose(hru["recharge:by_the_map"], 0.0, atol=1e-9)
    np.testing.assert_allclose(hru["recharge:received"], hru["recharge:sent"], rtol=1e-12, atol=1e-9)


@pytest.mark.slow
def test_osu_hru_recharge_budget_closes_every_day(osu):
    from test_one_way import budget
    (hru, ws), _ = osu
    b = budget(ws)
    np.testing.assert_allclose(b["RCH"], hru["recharge:sent"], rtol=1e-12, atol=1e-9)
    net = sum(b.values())
    gross = sum(np.abs(v) for v in b.values())
    assert np.max(np.abs(net) / gross) < 1e-7


AMES_DAYS = 365
AMES_BAND = 3  # rows of the synthetic grid per HRU: Ames_sub1's 12 HRUs over rows 0-35, rows 36-40 left out


class KnownPercolation:
    """Ames_sub1 percolates nothing (an irrigated field whose evapotranspiration takes all its water; stock SWAT+
    prints zero for every HRU), so its percolation is replaced with a known value per HRU and day, different for
    each: the test is of Basin Ledger's placement. That the value read is SWAT+'s is swatplus-bmi's test, and Osu_1hru's
    above."""

    def __init__(self, m):
        self._m, self.given = m, []

    def get_value(self, name, dest):
        if name != "soil_water__percolation_rate":
            return self._m.get_value(name, dest)
        t = self._m.get_current_time() - self._m.get_start_time()
        mm_d = (np.arange(len(dest)) + 1.0) * (1.0 + 0.5 * np.sin(2 * np.pi * t / 365.0))
        dest[:] = mm_d / 1000.0 / 86400.0
        self.given.append(np.asarray(dest, np.float64) * 86400.0)  # m d-1, as handed over
        return dest

    def __getattr__(self, name):
        return getattr(self._m, name)


@pytest.fixture(scope="module")
def ames(swatplus, tmp_path_factory):
    """Ames_sub1's first year, each of its 12 HRUs recharging its own band of three rows."""
    from basinledger.driver import Basin
    from conftest import RUNS, synthetic_config
    need_hru_percolation()
    proj = RUNS / "Ames_sub1-hru-recharge"
    if proj.exists():
        shutil.rmtree(proj)
    shutil.copytree(swatplus / "Ames_sub1", proj)
    names = [ln.split()[1] for ln in (proj / "hru.con").read_text().splitlines()[2:] if ln.split()]
    path = tmp_path_factory.mktemp("ames") / "hru_cells.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["object", "layer", "row", "col", "area_m2"])
        for i, h in enumerate(names):
            for r in range(AMES_BAND * i, AMES_BAND * (i + 1)):
                for c in range(41):
                    w.writerow([h, 0, r, c, 2500.0])
    cfg = synthetic_config(swat_project=proj, recharge_source="hru", recharge_map=str(path), area_rtol=None)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run="hru-ames", users=False, days=AMES_DAYS, wrap_surface=KnownPercolation)
    b.run()
    table, ws, swat_dir, given = b.ledger.table(), b.mf6.ws, b.swat.run_dir, np.array(b.swat.given)
    b.finalize()
    return names, table, ws, swat_dir, given


@pytest.mark.slow
def test_ames_each_hru_recharges_its_own_cells(ames):
    from basinledger.coupling.swatplus import hru_areas
    names, table, ws, swat_dir, given = ames
    area = np.array([hru_areas(swat_dir)[h] for h in names])
    assert given.shape == (AMES_DAYS, len(names))
    cbc = flopy.utils.CellBudgetFile(ws / "gwf.cbc", precision="double")
    rch = cbc.get_data(text="RCH")
    assert len(rch) == AMES_DAYS
    by_band = np.zeros((AMES_DAYS, len(names) + 1))
    for d, x in enumerate(rch):
        band = np.minimum((x["node"] - 1) // 41 // AMES_BAND, len(names))
        by_band[d] = np.bincount(band, weights=x["q"], minlength=len(names) + 1)
    assert np.all(by_band[:, -1] == 0)  # rows no HRU covers take nothing
    # each HRU's percolation times its area in SWAT+ reaches its own cells, every day, and the total is the ledger's
    np.testing.assert_allclose(by_band[:, :-1], given * area, rtol=1e-12)
    np.testing.assert_allclose(table["recharge:sent"], (given * area).sum(axis=1), rtol=1e-12)
    np.testing.assert_allclose(table["recharge:received"], table["recharge:sent"], rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(table["recharge:by_the_map"], 0.0, atol=1e-6)
    assert "recharge:swat_aquifers" not in table  # Ames_sub1 has no aquifers: nothing to hold it to
