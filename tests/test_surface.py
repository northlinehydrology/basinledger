"""Surface routes, the utility chain and rules on the synthetic example (examples/synthetic/surface.toml).

A mill takes from Osu_1hru's channel through a SWAT+ water allocation place of use and returns through a point
source. The plant is on a utility that pumps a well and discharges to the channel through the treatment plant's
point source. Three rules act on the users. Two years, with the users and without (the baseline still has the
utility's existing customers).
"""

from __future__ import annotations

import numpy as np
import pytest

from basinledger import config as C
from basinledger.driver import Basin
from conftest import ROOT, RUNS

pytestmark = pytest.mark.slow

DAYS = 730
DAY = 86400.0


def run(name, users):
    cfg = C.load(ROOT / "examples" / "synthetic" / "surface.toml")
    cfg.run_folder = RUNS
    b = Basin(cfg, run=name, users=users, days=DAYS)
    b.run()
    out = b.ledger.table(), {k: np.array(v) for k, v in b.series.items()}
    b.finalize()
    return out


@pytest.fixture(scope="module")
def with_users(swatplus):
    return run("surface-users", True)


@pytest.fixture(scope="module")
def baseline(swatplus):
    return run("surface-base", False)


def test_routes_deliver_what_was_asked(with_users):
    t, _ = with_users
    r = "route:swat:place_of_use:mill_intake:"
    assert t[r + "delivered"].sum() > 0
    assert (t[r + "delivered"] <= t[r + "requested"] * (1 + 1e-6) + 1e-6).all()
    np.testing.assert_allclose(t[r + "requested"], t["user:mill:withdrawal"], rtol=0)
    # the mill's day scales by the share its intake delivered, and its return leaves the next day: no return of
    # water the channel didn't give (0.13.0)
    share = t["user:mill:delivered_share"]
    assert (share <= 1).all()
    np.testing.assert_allclose(t["user:mill:delivered"], t[r + "delivered"], rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(t["user:mill:unmet"], t[r + "shortfall"], rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(t["user:mill:return_due"], share * t["user:mill:return"], rtol=1e-12)
    sent = t["route:swat:point_source:mill_outfall:delivered"]
    assert sent[0] == 0
    np.testing.assert_allclose(sent[1:], t["user:mill:return_due"][:-1], rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(t["user:mill:return_sent"][1:], t["user:mill:return_due"][:-1], rtol=1e-12)
    assert (sent <= np.r_[0.0, t[r + "delivered"][:-1]] + 1e-6).all()  # never more than the intake gave


def test_the_channel_carries_the_users_water(with_users, baseline):
    """Same day: the channel's outflow changes by what the users' routes put in and took out."""
    t, s = with_users
    tb, sb = baseline
    put = (t["route:swat:point_source:mill_outfall:delivered"] + t["route:swat:point_source:wwtp:delivered"]
           - tb["route:swat:point_source:wwtp:delivered"] - t["route:swat:place_of_use:mill_intake:delivered"])
    dq = (s["flow:cha01"] - sb["flow:cha01"]) * DAY
    np.testing.assert_allclose(dq, put, rtol=1e-4, atol=0.05)


def test_the_utility_balances(with_users):
    t, _ = with_users
    np.testing.assert_allclose(t["route:well:pw1:delivered"], t["utility:city:produced"], rtol=1e-12)
    existing = 0.002 * DAY
    np.testing.assert_allclose(t["utility:city:produced"], (existing + t["user:plant:withdrawal"]) / 0.9, rtol=1e-12)
    # the well delivers all it's asked, so the system's share is 1, and the plant sends yesterday's sewer (0.15.0)
    np.testing.assert_allclose(t["utility:city:share"], 1.0, rtol=1e-12)
    sewer = existing * 0.8 + t["user:plant:return"] + 0.0001 * DAY
    np.testing.assert_allclose(t["utility:city:sewer"], sewer, rtol=1e-12)
    np.testing.assert_allclose(t["route:swat:point_source:wwtp:requested"], np.r_[0.0, sewer[:-1]], rtol=1e-12)


def test_rules_act_on_yesterday(with_users):
    t, s = with_users
    q, h = s["flow:cha01"], s["head:pw1"]
    yq, yh = np.r_[np.inf, q[:-1]], np.r_[np.inf, h[:-1]]
    cut, dry, lim = (t[f"rule:{n}:active"].astype(bool) for n in
                     ("mill low-flow cut", "plant dry on low flow", "plant well limit"))
    assert np.array_equal(cut, yq < 0.004) and 50 < cut.sum() < DAYS - 50
    assert np.array_equal(dry, yq < 0.0035) and dry.sum() > 20
    assert np.array_equal(lim, yh < 88.3) and lim.sum() > 20
    assert (t["user:mill:withdrawal"][cut] == 0).all() and (t["user:mill:withdrawal"][~cut] > 0).all()
    assert (t["user:plant:withdrawal"][lim] <= 0.002 * DAY * (1 + 1e-12)).all()
    w = t["user:plant:withdrawal"]
    assert np.median(w[dry & ~lim]) < 0.2 * np.median(w[~dry & ~lim])  # dry: a tenth of the water, roughly
    e = t["user:plant:MWh"]
    assert e[dry].min() > e[~dry].max()  # and more electricity


def test_the_baseline_has_no_users(baseline):
    t, _ = baseline
    assert "user:plant:withdrawal" not in t
    assert t["utility:city:users_withdrawal"].sum() == 0
