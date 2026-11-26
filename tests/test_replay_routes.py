"""Replay with the human-system routes (replay.py; docs/design.md, Mass through the human-system routes).

The coupled town of tests/test_systems_mass.py, recorded and replayed: an intake and a well over two cells, homes on
the sewer, an edge on septic over three cells, leaks to the ground, a cooling user on the town's water (run again in
the replay, since it takes in the town's water), the plant at the point source:

- from a float64 record the replay is the full run bit for bit: the town's water and mass, the ground routes' mass,
  the salt's and the transport's accounts;
- from a float32 record, the same to float32's rounding of the flows;
- the stated concentrations may change in a replay (they move only mass), and a salt multiplier on septic scales the
  routes' septic mass; every account closes every day;
- a replay saved and restarted is the unbroken replay.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from conftest import RUNS
from test_systems_mass import EFFLUENT_CL, _config

from basinledger.calibration import Parameter
from basinledger.driver import Basin

pytestmark = pytest.mark.slow

DAYS = 60
EFFLUENT = {"chloride": EFFLUENT_CL, "weathering": "sewer"}
SEPTIC = {"chloride": 0.08, "weathering": "tap"}
MASS = ("utility:town:chloride:", "utility:town:weathering:", "salt:", "transport:", "user:cooling:")


def _cfg(effluent=EFFLUENT, septic=SEPTIC, record=None, precision="float64"):
    cfg = _config("salt.toml", effluent, septic)
    if record is not None:
        cfg.replay_record, cfg.replay_precision = RUNS / record, precision
    return cfg


def _run(cfg, name, replay=None, parameters=None, save_at=None, restart=None, days=DAYS):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        more = {"replay": replay} if replay is not None else {}
        b = Basin(cfg, run=name, days=days, parameters=parameters, restart=restart, **more)
    try:
        while b.day < (save_at or days):
            b.step()
        if save_at:
            b.save(b.folder / "state")
        return b.ledger.table()
    finally:
        b.finalize()


@pytest.fixture(scope="module")
def full64(swatplus):
    return _run(_cfg(record="record-routes64"), "replay-routes-full64")


@pytest.fixture(scope="module")
def replay64(full64):
    return _run(_cfg(), "replay-routes-64", replay=RUNS / "record-routes64")


def test_a_float64_record_replays_the_routes_mass_bit_for_bit(full64, replay64):
    keys = [k for k in replay64 if k.startswith(MASS)]
    assert any(k.startswith("utility:town:chloride:to_ground:") for k in keys)
    assert "transport:chloride:ground:edge:received" in keys and "transport:chloride:well:pw2:withdrawn" in keys
    for k in keys:
        np.testing.assert_array_equal(replay64[k], full64[k], err_msg=k)
    # the town's water, again from the routes' recorded lines
    for k in ("utility:town:produced", "utility:town:outfall", "utility:town:septic:edge", "utility:town:leaks"):
        np.testing.assert_array_equal(replay64[k], full64[k], err_msg=k)


def test_a_float32_record_replays_to_its_rounding(swatplus):
    full = _run(_cfg(record="record-routes32", precision="float32"), "replay-routes-full32")
    rep = _run(_cfg(), "replay-routes-32", replay=RUNS / "record-routes32")
    for k in rep:
        if k.startswith(MASS) and not k.endswith("closure_share"):
            a, b = np.nan_to_num(full[k]), np.nan_to_num(rep[k])
            np.testing.assert_allclose(b, a, rtol=1e-6, atol=1e-6 * max(np.abs(a).max(), 1e-3), err_msg=k)  # (1 nanogram on an account near zero)


def test_mass_settings_change_in_a_replay_and_close(full64):
    t = _run(_cfg(septic={"chloride": 0.04, "weathering": "tap"}), "replay-routes-septic", replay=RUNS / "record-routes64",
             parameters=[Parameter("septic", "salt", "chloride:septic", 3.0, 0.1, 10.0)])
    k = "utility:town:chloride:"
    np.testing.assert_allclose(t[k + "septic:edge:kg"], 3.0 * 0.04 * t["utility:town:septic:edge"], rtol=1e-12)
    np.testing.assert_array_equal(t["utility:town:septic:edge"], full64["utility:town:septic:edge"])
    assert not np.allclose(np.nan_to_num(t[k + "septic:edge:kg"]), np.nan_to_num(full64[k + "septic:edge:kg"]))
    for sp in ("chloride", "weathering"):
        assert np.nanmax(np.abs(t[f"utility:town:{sp}:closure_share"])) < 1e-12
        assert np.nanmax(np.abs(t[f"salt:{sp}:closure_share"])) < 1e-9
        assert np.nanmax(np.abs(t[f"transport:{sp}:closure_share"])) < 1e-6


def test_a_replay_restarted_is_the_unbroken_replay(replay64):
    _run(_cfg(), "replay-routes-warm", replay=RUNS / "record-routes64", save_at=25)
    b = _run(_cfg(), "replay-routes-restarted", replay=RUNS / "record-routes64",
             restart=RUNS / "replay-routes-warm" / "state")
    for k in replay64:
        if k.startswith(MASS):
            np.testing.assert_array_equal(b[k], replay64[k], err_msg=k)
