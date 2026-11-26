"""Water put on SWAT+'s fields (route ``land:<name>``, ``[[land]]``; docs/design.md, Mass through the human-system
routes, part 4), through swatplus-bmi 0.15.0's irrigation applied from outside.

On the synthetic salt and nitrate examples, Osu_1hru made upland, two reported irrigators on its HRU: a golf course
drawing from an intake on the channel (30% lost in the application, 20% of what lands running off) and a farm on a
well over two cells.

- the read and the refusals (a paddy, an HRU irrigated from wells for SWAT+'s decision tables too, the runoff's
  bounds, a source with no known concentration);
- each day's delivered water reaches the fields the next day, split into the soil and the runoff, and SWAT+ reports
  taking it all;
- its mass is the whole withdrawal's: the intake's at salt's mixed concentration (what the stream lost), the well's
  as the transport model took it out, or a stated one;
- salt's land takes it as a source, and every closure holds; nitrate and soluble P go into the top layer, SWAT+'s
  read-back equal to what was set;
- a restart and a replay carry it.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from conftest import ROOT, RUNS
from land_fixtures import DAY, land_config

from basinledger import config as C
from basinledger.config import LandConfig
from basinledger.driver import Basin

DAYS = 90


@pytest.fixture(scope="module")
def swatplus15(swatplus):
    """SWAT+ through swatplus-bmi 0.15.0 or later (its irrigation applied from outside)."""
    import swatplus_bmi
    if tuple(int(x) for x in swatplus_bmi.__version__.split(".")[:2]) < (0, 15):
        pytest.skip(f"the land route needs swatplus-bmi 0.15.0 (this is {swatplus_bmi.__version__})")
    return swatplus


def _run(cfg, name, days=DAYS, save_at=None, restart=None, replay=None):
    more = {"replay": replay} if replay is not None else {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run=name, days=days, restart=restart, **more)
    try:
        while b.day < (save_at or days):
            b.step()
        if save_at:
            b.save(b.folder / "state")
        return b.ledger.table()
    finally:
        b.finalize()


def lag(x):
    return np.r_[0.0, np.nan_to_num(x)[:-1]]


# ---------- read and refused ----------

def test_land_entries_read(tmp_path):
    from basinledger.config import _land
    g = _land(tmp_path / "x.toml", {"name": "golf", "hrus": [{"hru": "a", "weight": 1}, {"hru": "b", "weight": 3}],
                                    "runoff": 0.1, "concentration": {"chloride": 0.02}})
    assert g == LandConfig("golf", (("a", 0.25), ("b", 0.75)), 0.1, {"chloride": 0.02})
    (tmp_path / "h.csv").write_text("hru,weight\nhru0001,2\nhru0002,2\n")
    assert _land(tmp_path / "x.toml", {"name": "f", "hrus": "h.csv"}).hrus == (("hru0001", 0.5), ("hru0002", 0.5))
    for bad in ({"name": "f"}, {"name": "f", "hrus": []}, {"name": "f", "hrus": [{"hru": "a", "weight": 0}]},
                {"name": "f", "hrus": [{"hru": "a", "weight": 1}], "runoff": 1.0},
                {"name": "f", "hrus": [{"hru": "a", "weight": 1}, {"hru": "a", "weight": 1}]},
                {"name": "f", "hrus": [{"hru": "a", "weight": 1}], "rate": 2}):
        with pytest.raises(ValueError):
            _land(tmp_path / "x.toml", bad)


@pytest.mark.slow
def test_what_the_start_up_refuses(swatplus15, tmp_path):
    cfg = land_config("salt.toml", tmp_path, 3)
    cfg.swat_project = C.load(ROOT / "examples" / "synthetic" / "salt.toml").swat_project  # Osu_1hru as shipped
    for s in cfg.salt["chloride"].sources:  # (as shipped, the septic source's land use is the paddy's)
        if s.by == "landuse":
            s.rate = {"rice140_lum": v for v in s.rate.values()}
    with pytest.raises(ValueError, match="HRU hru0001 is a wetland or paddy"):
        Basin(cfg, run="land-refused", days=3)
    cfg = land_config("salt.toml", tmp_path, 3)
    cfg.irrigation = {"hru0001": [("farm", 1.0)]}
    with pytest.raises(ValueError, match="irrigated by SWAT.'s decision tables from wells"):
        Basin(cfg, run="land-refused", days=3)
    cfg = land_config("nitrate.toml", tmp_path, 3)
    cfg.users[-2].withdrawal = "offsite:golf_supply"
    with pytest.raises(ValueError, match="water onto the land comes from a well, an intake or a water system"):
        Basin(cfg, run="land-refused", days=3)


# ---------- on the models ----------

@pytest.fixture(scope="module")
def salted(swatplus15, tmp_path_factory):
    tmp = tmp_path_factory.mktemp("land-salt")
    cfg = land_config("salt.toml", tmp, DAYS)
    cfg.replay_record, cfg.replay_precision = RUNS / "record-land", "float64"
    return _run(cfg, "land-salt"), tmp


@pytest.mark.slow
def test_the_water_reaches_the_fields_the_next_day(salted):
    t, _ = salted
    for user, runoff in (("golf", 0.2), ("farm", 0.0)):
        r = f"route:land:{user}:"
        np.testing.assert_allclose(t[r + "requested"], lag(t[f"user:{user}:return_due"]), rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(t[r + "delivered"], t[r + "requested"], rtol=1e-6)  # SWAT+ took all of it
        np.testing.assert_allclose(t[r + "runoff"], runoff * t[r + "delivered"], rtol=1e-6, atol=1e-9)
        assert np.nansum(t[r + "delivered"]) > 0
    # the golf course's intake fell short on some days: what it returns is what was delivered
    assert np.nanmin(t["user:golf:delivered_share"]) < 1.0


@pytest.mark.slow
def test_the_mass_is_the_whole_withdrawals(salted):
    t, _ = salted
    for sp in ("chloride", "weathering"):
        k = f"user:golf:land:{sp}:"
        # the intake's water left the channel at salt's mixed concentration: the golf course is its only taker
        np.testing.assert_allclose(t[k + "kg_held"], t[f"salt:{sp}:streams:withdrawn"], rtol=1e-6, atol=1e-12)
        # the farm's left the aquifer as the transport model took it out of its two cells
        np.testing.assert_allclose(t[f"user:farm:land:{sp}:kg_held"], t[f"transport:{sp}:well:farm:withdrawn"],
                                   rtol=1e-9, atol=1e-12)
        # and reaches the land's stores the next day, as a source from outside
        np.testing.assert_allclose(t[f"salt:{sp}:land:source:irrigation"],
                                   lag(t[k + "kg_held"]) + lag(t[f"user:farm:land:{sp}:kg_held"]), rtol=1e-12,
                                   atol=1e-15)
        assert np.nanmax(np.abs(t[f"salt:{sp}:closure_share"])) < 1e-9
        assert np.nanmax(np.abs(t[f"salt:{sp}:land:closure_share"])) < 1e-12
        assert np.nanmax(np.abs(t[f"transport:{sp}:closure_share"])) < 1e-6


@pytest.mark.slow
def test_a_stated_concentration_and_the_nutrients_into_the_top_layer(swatplus15, tmp_path):
    cfg = land_config("nitrate.toml", tmp_path, 40, stated={"nitrate": 0.004})
    t = _run(cfg, "land-nitrate", 40)
    np.testing.assert_allclose(t["user:golf:land:nitrate:kg_held"], 0.004 * t["user:golf:delivered"], rtol=1e-12)
    for sp in ("nitrate", "phosphorus"):
        for user in ("golf", "farm"):
            k = f"transport:{sp}:land:{user}:"
            np.testing.assert_allclose(t[k + "sent"], lag(t[f"user:{user}:land:{sp}:kg_held"]), rtol=1e-12,
                                       atol=1e-15)
            np.testing.assert_allclose(t[k + "received"], t[k + "sent"], rtol=1e-6, atol=1e-12)  # SWAT+ took it
        assert np.nanmax(np.abs(t[f"transport:{sp}:closure_share"])) < 1e-6
    assert np.nansum(t["transport:nitrate:land:golf:sent"]) > 0


@pytest.mark.slow
def test_a_restart_carries_what_is_held(swatplus15, tmp_path):
    n, at = 40, 17
    a = _run(land_config("salt.toml", tmp_path, n), "land-unbroken", n)
    _run(land_config("salt.toml", tmp_path, n), "land-warm", n, save_at=at)
    b = _run(land_config("salt.toml", tmp_path, n), "land-restarted", n, restart=RUNS / "land-warm" / "state")
    for key in a:
        if key.startswith(("route:land:", "user:golf:land:", "user:farm:land:", "salt:chloride:land:")):
            np.testing.assert_array_equal(a[key], b[key], err_msg=key)


@pytest.mark.slow
def test_a_replay_puts_the_same_mass_on_the_land(salted):
    t, tmp = salted
    r = _run(land_config("salt.toml", tmp, DAYS, stated={}), "land-replay", replay=RUNS / "record-land")
    for key in r:
        if key.startswith(("salt:", "user:golf:land:", "user:farm:land:", "transport:chloride:")):
            np.testing.assert_array_equal(r[key], t[key], err_msg=key)
    # the stated concentration may change in a replay: it moves only mass
    s = _run(land_config("salt.toml", tmp, DAYS, stated={"chloride": 0.05}), "land-replay-stated",
             replay=RUNS / "record-land")
    np.testing.assert_allclose(s["user:golf:land:chloride:kg_held"], 0.05 * t["user:golf:delivered"], rtol=1e-12)
    assert np.nanmax(np.abs(s["salt:chloride:closure_share"])) < 1e-9
