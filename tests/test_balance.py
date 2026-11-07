"""Tier 2 and 3: the whole basin closes over a coupled run (balance.py, testing G3), on the synthetic two-way example."""

from __future__ import annotations

import numpy as np
import pytest

from basinledger import config as C
from basinledger.balance import basin_balance
from basinledger.driver import Basin
from conftest import ROOT, RUNS

pytestmark = pytest.mark.slow


def test_the_whole_basin_closes(swatplus):
    cfg = C.load(ROOT / "examples" / "synthetic" / "two_way.toml")
    cfg.run_folder = RUNS
    cfg.outputs = dict(cfg.outputs, balance="cha01")
    b = Basin(cfg, run="two-way-balance", days=1096)
    b.run()
    b.write_outputs()
    b.finalize()
    r = basin_balance(RUNS / "two-way-balance", "cha01")
    p = r["precipitation_land"]
    assert p > 0
    # SWAT+'s own land budget: Osu_1hru's HRU has a wetland whose storage basin_wb doesn't break out (1.4% of
    # precipitation and irrigation when written, Nov 6, 2026)
    assert abs(r["land_residual"]) <= 0.02 * (p + r["irrigation"])
    assert abs(r["inflow_from_elsewhere"]) <= 0.001 * r["outflow"]   # Osu_1hru has none; SWAT+ prints rounded flows
    assert abs(r["groundwater_residual"]) <= 1e-6 * r["recharge_received"]
    assert r["recharge_made"] == pytest.approx(r["recharge_received"], rel=1e-6)
    assert abs(r["channels_residual"]) <= 0.001 * r["outflow"]          # channel storage isn't printed
    assert r["max_daily_closure_share"] <= 1e-6
