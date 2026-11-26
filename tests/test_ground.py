"""Returns to the ground (route ``ground:<name>``, ``[[ground]]``): a septic drainfield or a system's leaks, water into
an aquifer cell through a WEL entry Basin Ledger adds, injecting it.

- the entries read, and what the start-up refuses (tier 1; their mass with transport and salt is tests/test_ground_mass.py, in a replay tests/test_replay_routes.py);
- coupled on the synthetic example, a house pumping its own well and returning 80% of it to a drainfield three cells
  away, against the same house returning it outside the basin (tier 2): the drainfield gets exactly the return every
  day, MODFLOW 6's budget closes with it, the capture analysis puts the return in the aquifer's water
  (``from_returns``), and the water table under the drainfield is higher.
"""

from __future__ import annotations

import datetime as dt

import numpy as np
import pytest
from conftest import ROOT, RUNS

from basinledger import config as C
from basinledger.analysis import capture
from basinledger.config import UserConfig
from basinledger.coupling import swatplus as S
from basinledger.driver import Basin

DAY = 86400.0
DAYS = 200
PUMP = 0.002              # m3 s-1, the house (a large one, so the aquifer notices)
RETURNED = 0.8            # to the drainfield; the rest consumed
FIELD = (0, 20, 31)       # three cells east of the well pw1 at (0, 20, 28)


def test_ground_entries_read(tmp_path):
    base = (ROOT / "examples" / "synthetic" / "surface.toml").read_text().replace(
        '"model:build"', f'"{(ROOT / "examples" / "synthetic").as_posix()}/model.py:build"')
    (tmp_path / "x.toml").write_text(base + '\n[[ground]]\nname = "field"\ncell = [0, 20, 31]\n')
    (g,) = C.load(tmp_path / "x.toml").ground
    assert (g.name, g.package, g.cell, g.add) == ("field", "wel_ground", (0, 20, 31), True)
    (tmp_path / "x.toml").write_text(base + '\n[[ground]]\nname = "field"\ncell = [0, 20, 31]\nrate = 3\n')
    with pytest.raises(ValueError, match="takes name and cell"):
        C.load(tmp_path / "x.toml")
    (tmp_path / "x.toml").write_text(base + '\n[[ground]]\nname = "field"\n')
    with pytest.raises(ValueError, match="cell must be"):
        C.load(tmp_path / "x.toml")


def _config(tmp_path, returns):
    cfg = C.load(ROOT / "examples" / "synthetic" / "surface.toml")
    cfg.run_folder = RUNS
    cfg.places_of_use, cfg.point_sources, cfg.utilities, cfg.rules = [], [], [], []
    cfg.ground = [C.WellConfig("field", "wel_ground", FIELD, True)]
    cfg.outputs = dict(cfg.outputs, heads=[{"name": "field", "cell": list(FIELD)}])
    start = S._period(cfg.swat_project)[0]
    w = PUMP * DAY
    rows = ["date,withdrawal_m3,consumption_m3,return_m3"]
    rows += [f"{start + dt.timedelta(days=d)},{w},{(1 - RETURNED) * w},{RETURNED * w}" for d in range(DAYS)]
    (tmp_path / "house.csv").write_text("\n".join(rows) + "\n")
    (tmp_path / "house.toml").write_text(f'series = "house.csv"\nstart = "{start}"\n')
    cfg.users = [UserConfig(name="house", component="basinledger.generic:ReportedUser",
                            config=tmp_path / "house.toml", withdrawal="well:pw1", returns=returns)]
    return cfg


def _run(tmp_path, name, returns):
    b = Basin(_config(tmp_path, returns), run=name, days=DAYS)
    b.run()
    out = b.ledger.table(), np.array(b.series["head:field"])
    b.finalize()
    return out


@pytest.mark.slow
def test_a_drainfield_returns_the_house_water_to_the_aquifer(swatplus, tmp_path):
    t, h = _run(tmp_path, "ground-drainfield", "ground:field")
    t0, h0 = _run(tmp_path, "ground-offsite", "offsite:house")
    r = "route:ground:field:"
    np.testing.assert_allclose(t[r + "requested"], RETURNED * PUMP * DAY, rtol=1e-12)
    np.testing.assert_allclose(t[r + "delivered"], t[r + "requested"], rtol=1e-9)
    np.testing.assert_allclose(t[r + "delivered"], t["user:house:return"], rtol=1e-9)
    np.testing.assert_allclose(t["route:well:pw1:delivered"], PUMP * DAY, rtol=1e-9)
    # MODFLOW 6's budget closes every day with the drainfield in it
    assert np.abs(t["groundwater:closure_share"]).max() < 1e-6
    # against the house returning outside the basin, the aquifer holds the drainfield's water
    c = capture(t0, t)
    assert c["total"]["pumping"] == pytest.approx(0.0, abs=1e-6)
    assert c["total"]["from_returns"] == pytest.approx(RETURNED * PUMP * DAY * DAYS, rel=1e-9)
    assert abs(c["total"]["residual"]) < 1e-6 * c["total"]["from_returns"]
    assert (h[20:] > h0[20:]).all()  # the water table under the field is higher
    assert "groundwater:boundary:wel_ground" not in t or not t["groundwater:boundary:wel_ground"].any()


def test_cells_read_inline_and_from_a_file(tmp_path):
    from basinledger.config import _cells
    got = _cells(tmp_path / "x.toml", "w", [[0, 20, 30, 1], [0, 20, 31, 3]])
    assert got == (((0, 20, 30), 0.25), ((0, 20, 31), 0.75))
    (tmp_path / "c.csv").write_text("layer,icell2d,weight\n0,412,2\n0,413,2\n")
    assert _cells(tmp_path / "x.toml", "w", "c.csv") == (((0, 412), 0.5), ((0, 413), 0.5))
    for bad in ([], [[0, 20]], [[0, 20, 30, 0]], [[0, 20, 30, 1], [0, 20, 30, 2]], [[0, 20, 30, True]]):
        with pytest.raises(ValueError):
            _cells(tmp_path / "x.toml", "w", bad)
    base = (ROOT / "examples" / "synthetic" / "surface.toml").read_text().replace(
        '"model:build"', f'"{(ROOT / "examples" / "synthetic").as_posix()}/model.py:build"')
    (tmp_path / "x.toml").write_text(base + '\n[[ground]]\nname = "f"\ncells = [[0, 20, 31, 1], [0, 21, 31, 1]]\n')
    (g,) = C.load(tmp_path / "x.toml").ground
    assert g.cell is None and g.cells == (((0, 20, 31), 0.5), ((0, 21, 31), 0.5))
    (tmp_path / "x.toml").write_text(base + '\n[[ground]]\nname = "f"\ncell = [0, 20, 31]\ncells = [[0, 21, 31, 1]]\n')
    with pytest.raises(ValueError, match="cell or cells, not both"):
        C.load(tmp_path / "x.toml")


@pytest.mark.slow
def test_a_well_and_a_drainfield_over_several_cells(swatplus, tmp_path):
    """A house's pumping spread over two cells (3:1) and its drainfield over three (1:2:1), as a group of domestic
    wells and their septic water would be over a HUC-12: each cell's WEL entry gets its weight's share of the day's
    request, the route's delivered is their sum, and MODFLOW 6's budget closes."""
    cfg = _config(tmp_path, "ground:field")
    cfg.wells = [C.WellConfig("pw1", "wel_house", None, True, cells=(((0, 20, 28), 0.75), ((0, 22, 28), 0.25)))]
    cfg.ground = [C.WellConfig("field", "wel_ground", None, True,
                               cells=(((0, 20, 31), 0.25), ((0, 21, 31), 0.5), ((0, 22, 31), 0.25)))]
    cfg.outputs = {}
    b = Basin(cfg, run="ground-spread", days=60)  # the pump's cells in a WEL package of Basin Ledger's own
    cells = {}
    for _ in range(60):
        b.step()
        for route in ("well:pw1", "ground:field"):
            for sub, w in b.router.groups[route]:
                kind_, name = sub.split(":", 1)
                pkg, i = (b.router._entries if kind_ == "well" else b.router._ground_entries)[name]
                cells.setdefault(sub, []).append((abs(float(b.mf6.simvals(pkg)[i])), w))
    t = b.ledger.table()
    b.finalize()
    for route, total in (("route:well:pw1:", PUMP * DAY), ("route:ground:field:", RETURNED * PUMP * DAY)):
        np.testing.assert_allclose(t[route + "requested"], total, rtol=1e-12)
        np.testing.assert_allclose(t[route + "delivered"], total, rtol=1e-9)
    for sub, rows in cells.items():
        got, w = np.array([r[0] for r in rows]), rows[0][1]
        total = PUMP * DAY if sub.startswith("well:") else RETURNED * PUMP * DAY
        np.testing.assert_allclose(got, w * total, rtol=1e-9)
    assert np.abs(t["groundwater:closure_share"]).max() < 1e-6
