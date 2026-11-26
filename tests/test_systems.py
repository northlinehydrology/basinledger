"""A water system with customers, coupled (chain.py; docs/design.md, "Water systems with customers").

On the synthetic example, a town draws from an intake on Osu_1hru's channel first (up to 0.003 m3 s-1, often more
than the channel gives) and then from the well pw1:
- its homes' demand changes by month, 70% to the sewer;
- its edge is on septic, 90% to a drainfield;
- half its losses leak to the ground;
- a data center on the town returns 20% of what it gets to the sewer;
- a second system cools with the town's plant effluent.

Each day the system's day scales by the share of production its sources delivered, and its returns leave the next
day. Septic water and leaks go into the aquifer through ``ground:`` routes, and the plant sends yesterday's sewer
less what the reuse took. Every pair closes, and MODFLOW 6's budget closes with the ground returns.
"""

from __future__ import annotations

import datetime as dt

import numpy as np
import pytest
from conftest import ROOT, RUNS

from basinledger import config as C
from basinledger.config import UserConfig
from basinledger.coupling import swatplus as S
from basinledger.driver import Basin

pytestmark = pytest.mark.slow

DAY = 86400.0
DAYS = 240


@pytest.fixture(scope="module")
def run(swatplus, tmp_path_factory):
    tmp = tmp_path_factory.mktemp("systems")
    cfg = C.load(ROOT / "examples" / "synthetic" / "surface.toml")
    cfg.run_folder = RUNS
    cfg.places_of_use = [("town_intake", "cha01")]
    cfg.point_sources = [("wwtp", "cha01")]
    cfg.rules = []
    cfg.ground = [C.WellConfig("edge", "wel_ground", (0, 20, 31), True),
                  C.WellConfig("mains", "wel_ground", (0, 12, 20), True)]
    cfg.utilities = [
        {"name": "town", "outfall": "swat:point_source:wwtp", "losses": 0.1,
         "leaks": {"share": 0.5, "to": "ground:mains"},
         "sources": [{"route": "swat:place_of_use:town_intake", "priority": 1, "capacity": 0.003},
                     {"route": "well:pw1", "priority": 2}],
         "customers": [{"name": "homes", "demand": [0.0015] * 4 + [0.0025] * 5 + [0.0015] * 3,
                        "sewered": 0.7, "consumed": 0.3},
                       {"name": "edge", "demand": 0.0004, "septic": 0.9, "septic_to": "ground:edge",
                        "consumed": 0.1}]},
        {"name": "cooling", "outfall": "outfall:cooling", "sources": [{"route": "effluent:town", "share": 1.0}],
         "customers": [{"name": "towers", "demand": 0.0003, "consumed": 1.0}]}]
    start = S._period(cfg.swat_project)[0]
    w = 0.001 * DAY
    rows = ["date,withdrawal_m3,consumption_m3,return_m3"]
    rows += [f"{start + dt.timedelta(days=d)},{w},{0.8 * w},{0.2 * w}" for d in range(DAYS)]
    (tmp / "dc.csv").write_text("\n".join(rows) + "\n")
    (tmp / "dc.toml").write_text(f'series = "dc.csv"\nstart = "{start}"\n')
    cfg.users = [UserConfig(name="dc", component="basinledger.generic:ReportedUser", config=tmp / "dc.toml",
                            withdrawal="utility:town", returns="utility:town")]
    b = Basin(cfg, run="systems-town", days=DAYS)
    b.run()
    t = b.ledger.table()
    b.finalize()
    return t


def lag(x):
    return np.r_[0.0, x[:-1]]


def test_the_sources_in_priority_and_the_share(run):
    t = run
    k = "utility:town:"
    intake, well = (k + "source:swat:place_of_use:town_intake:", k + "source:well:pw1:")
    np.testing.assert_allclose(t[intake + "requested"], np.minimum(t[k + "asked"], 0.003 * DAY), rtol=1e-12)
    np.testing.assert_allclose(t[well + "requested"], t[k + "asked"] - t[intake + "requested"], rtol=1e-9, atol=1e-9)
    share = t[k + "share"]
    assert (share < 0.99).sum() > 20 and (share > 0.999).sum() > 5  # days the channel fell short, and not
    np.testing.assert_allclose(t[k + "produced"], t[intake + "delivered"] + t[well + "delivered"], rtol=1e-9)
    np.testing.assert_allclose(t[k + "delivered"], share * t[k + "demand"], rtol=1e-9)
    # the homes' demand by month: April's and May's
    month = np.array([(dt.date(2010, 1, 1) + dt.timedelta(days=d)).month for d in range(DAYS)])
    np.testing.assert_allclose(t[k + "demand:homes"][month == 4], 0.0015 * DAY, rtol=1e-12)
    np.testing.assert_allclose(t[k + "demand:homes"][month == 5], 0.0025 * DAY, rtol=1e-12)


def test_returns_leave_the_next_day(run):
    t = run
    k = "utility:town:"
    np.testing.assert_allclose(t["route:ground:edge:delivered"], lag(t[k + "septic:edge"]), rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(t["route:ground:mains:delivered"], lag(t[k + "leaks"]), rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(t[k + "plant_outflow"], lag(t[k + "sewer"]), rtol=1e-12)
    np.testing.assert_allclose(t["route:swat:point_source:wwtp:delivered"], t[k + "outfall"], rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(t[k + "sewer"], 0.7 * t[k + "delivered:homes"] + t[k + "users_return"], rtol=1e-9)
    assert np.abs(t["groundwater:closure_share"]).max() < 1e-6


def test_the_data_center_on_the_town(run):
    t = run
    share = t["utility:town:share"]
    np.testing.assert_allclose(t["user:dc:delivered"], share * 0.001 * DAY, rtol=1e-9)
    np.testing.assert_allclose(t["user:dc:return_due"], 0.2 * t["user:dc:delivered"], rtol=1e-9)
    np.testing.assert_allclose(t["utility:town:users_return"], t["user:dc:return_due"], rtol=1e-12)
    np.testing.assert_allclose(t["utility:town:users_delivered"], t["user:dc:delivered"], rtol=1e-9)


def test_the_reuse_takes_the_plants_effluent(run):
    t = run
    want = 0.0003 * DAY
    out = t["utility:town:plant_outflow"]
    np.testing.assert_allclose(t["utility:town:reuse:cooling"], np.minimum(want, out), rtol=1e-12)
    np.testing.assert_allclose(t["utility:cooling:produced"], np.minimum(want, out), rtol=1e-12)
    np.testing.assert_allclose(t["utility:town:outfall"], out - t["utility:town:reuse"], rtol=1e-12)
    # short only when the town's plant sends less than it asks (none on the first day)
    np.testing.assert_allclose(t["utility:cooling:share"], np.minimum(out / want, 1.0), rtol=1e-12)
    assert t["utility:cooling:share"][0] == 0.0 and (t["utility:cooling:share"] == 1.0).sum() > 100
