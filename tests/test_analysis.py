"""Capture of a well's pumping and a return's concentration change (analysis.py), on the synthetic examples."""

from __future__ import annotations

import numpy as np
import pytest

from basinledger import analysis as A
from basinledger import config as C
from basinledger.driver import Basin
from conftest import ROOT, RUNS

DAYS = 365


def run(example: str, run_name: str, users: bool):
    cfg = C.load(ROOT / "examples" / "synthetic" / example)
    cfg.run_folder = RUNS
    b = Basin(cfg, run=run_name, users=users, days=DAYS)
    b.run()
    t, series = b.ledger.table(), {k: np.array(v) for k, v in b.series.items()}
    b.finalize()
    return t, series


@pytest.mark.slow
@pytest.mark.parametrize("example", ["two_way.toml", "basin.toml"])
def test_pumping_is_all_captured(swatplus, example):
    """The extra pumping equals the extra release from storage, the water the streams lose and the extra boundary
    inflow, every day: nothing made or lost between the runs."""
    tag = example.split(".")[0]
    base, _ = run(example, f"{tag}-capture-base", users=False)
    case, _ = run(example, f"{tag}-capture-case", users=True)
    c = A.capture(base, case)
    assert c["total"]["pumping"] > 0
    assert abs(c["total"]["residual"]) <= 1e-6 * c["total"]["pumping"]
    assert np.max(np.abs(c["daily"]["residual"])) <= 1e-6 * np.max(c["daily"]["pumping"])
    assert c["total"]["from_recharge"] == pytest.approx(0.0, abs=1e-6)
    assert c["share"]["from_storage"] > 0
    if example == "two_way.toml":  # the streams supply part of it, through the coupled exchange
        assert c["share"]["from_streams"] > 0
    else:  # one-way: the river is a boundary of the groundwater model
        assert c["share"]["from_boundary"] > 0


def test_concentration_change_by_hand():
    out = A.concentration_change(np.array([100.0, 50.0, 0.0]), np.array([10.0, 10.0, 0.0]),
                                 np.array([5.0, 5.0, 0.0]), 0.2)
    np.testing.assert_allclose(out["concentration_kg_m3"][:2], [(0.2 * 90 + 5) / 100, (0.2 * 40 + 5) / 50])
    np.testing.assert_allclose(out["change_percent"][:2], [15.0, 30.0])
    assert np.isnan(out["concentration_kg_m3"][2])


@pytest.mark.slow
def test_a_return_raises_the_receiving_reach(swatplus):
    """The surface example's mill returns dissolved solids to Osu_1hru's channel through a point source: the
    receiving reach's concentration rises on the days it returns, and not on the others."""
    t, series = run("surface.toml", "surface-concentration", users=True)
    ret = t["route:swat:point_source:mill_outfall:delivered"]
    kg = t["user:mill:dissolved_solids:kg"]
    flow = series["flow:cha01"] * 86400.0
    out = A.concentration_change(flow, ret, kg, 0.2)
    on = ret > 0
    assert on.sum() > 100
    assert np.nanmin(out["change_percent"][on]) > 0
    assert np.allclose(out["change_percent"][~on & (flow > 0)], 0.0)
