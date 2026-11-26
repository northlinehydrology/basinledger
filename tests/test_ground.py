"""Returns to the ground (route ``ground:<name>``, ``[[ground]]``): a septic drainfield or a system's leaks, water into
an aquifer cell through a WEL entry Basin Ledger adds, injecting it.

- the entries read, and what the start-up refuses (tier 1);
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


def test_a_return_to_the_ground_needs_modflow_and_no_transport(swatplus, tmp_path):
    cfg = _config(tmp_path, "ground:field")
    cfg.salt = {"cl": object()}
    with pytest.raises(NotImplementedError, match=r"\[\[ground\]\] with transport or salt isn't built"):
        Basin(cfg, run="ground-refused", days=2)
