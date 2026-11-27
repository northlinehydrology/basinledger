"""``[rivers] stage_floor`` (driver.py, docs/design.md, River stage floor): the day's stage from the surface model's
depth, no lower than a floor the rivers' map gives cell by cell. On the synthetic two-way example, its river along
the west edge (cells (0, r, 0), bed bottom 86 m, so the stage is 87 m + SWAT+'s centimeters of depth):
- a floor that never binds changes nothing: the same run to the bit;
- the floor holds where it is above the stage, cell by cell, the ledger counts it and outputs/stage_floor.csv lists it;
- the basin still closes;
- a bad setting is refused."""

from __future__ import annotations

import csv
import warnings

import numpy as np
import pytest
from conftest import ROOT, RUNS

from basinledger import config as C
from basinledger.balance import basin_balance
from basinledger.coupling import swatplus as S
from basinledger.driver import Basin

DAYS = 365
HIGH, LOW = 88.0, 80.0   # above every day's stage (87 m + centimeters), below every one


def _map(path, floors):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["object", "layer", "row", "col", "surface"])
        w.writerows(["cha01", 0, r, 0, "" if x is None else x] for r, x in enumerate(floors))
    return str(path)


def _cfg(rivers=None, balance=False):
    cfg = C.load(ROOT / "examples" / "synthetic" / "two_way.toml")
    cfg.run_folder = RUNS
    cfg.rivers = rivers
    if balance:
        cfg.outputs = dict(cfg.outputs, balance="cha01")
    return cfg


def _rivers(path, floor=True, **kw):
    return {"package": "riv", "map": path, "stage": "surface_depth", "bed_thickness": 1.0,
            **({"stage_floor": "surface"} if floor else {}), **kw}


def _run(cfg, name, days=DAYS, each_day=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run=name, days=days)
        try:
            for _ in range(days):
                b.step()
                if each_day:
                    each_day(b)
            b.write_outputs()
            return b.ledger.table()
        finally:
            b.finalize()


def _same(a: dict, b: dict) -> bool:
    return a.keys() == b.keys() and all(np.asarray(a[k]).tobytes() == np.asarray(b[k]).tobytes() for k in a)


@pytest.fixture
def depth(swatplus):
    if not S.has_groundwater_exchange():
        pytest.skip("needs swatplus-bmi 0.2.0 (channel depth)")
    return swatplus


def test_a_floor_that_never_binds_changes_nothing(depth, tmp_path):
    path = _map(tmp_path / "rivers.csv", [LOW] * 5 + [None] * 5)
    off = _run(_cfg(_rivers(path, floor=False)), "floor-off")
    low = _run(_cfg(_rivers(path)), "floor-low")
    assert (low.pop("rivers:stage_floor:cells") == 0).all() and (low.pop("rivers:stage_floor:raise_m") == 0).all()
    assert _same(off, low)
    assert not (RUNS / "floor-off" / "outputs" / "stage_floor.csv").exists()


def test_the_floor_holds_cell_by_cell(depth, tmp_path):
    """Rows 0-3 at 88 m (the input stage, above every day's), rows 4-6 with none, rows 7-9 at 80 m. With the
    exchange's supply limit off, the stage set is what MODFLOW 6 sees."""
    floors = [HIGH] * 4 + [None] * 3 + [LOW] * 3
    cfg = _cfg(_rivers(_map(tmp_path / "rivers.csv", floors)))
    cfg.exchange = dict(cfg.exchange, limit=False)
    days = 120
    seen = []

    def check(b):
        depth = float(b.sw.get(S.DEPTH)[0])  # the day's, which set the day's stage
        stage = b.mf6.get("STAGE", "GWF", "RIV")
        rbot = b.mf6.get("RBOT", "GWF", "RIV")
        nominal = rbot + 1.0 + depth
        np.testing.assert_array_equal(stage[:4], HIGH)
        np.testing.assert_allclose(stage[4:], nominal[4:], rtol=0, atol=1e-5)
        seen.append(float(np.sum(HIGH - nominal[:4])))

    t = _run(cfg, "floor-holds", days=days, each_day=check)
    assert (t["rivers:stage_floor:cells"] == 4).all()
    np.testing.assert_allclose(t["rivers:stage_floor:raise_m"], seen, rtol=0, atol=1e-4)
    with open(RUNS / "floor-holds" / "outputs" / "stage_floor.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    assert [int(r["entry"]) for r in rows] == [0, 1, 2, 3, 7, 8, 9]      # each cell with a floor; none for 4-6
    assert [int(r["days_at_floor"]) for r in rows] == [days] * 4 + [0] * 3
    assert [float(r["floor_m"]) for r in rows] == [HIGH] * 4 + [LOW] * 3
    assert sum(float(r["raise_m_days"]) for r in rows) == pytest.approx(sum(seen), rel=1e-9)


@pytest.mark.slow
def test_the_basin_closes_with_the_floor(depth, tmp_path):
    path = _map(tmp_path / "rivers.csv", [HIGH] * 5 + [None] * 5)
    t = _run(_cfg(_rivers(path), balance=True), "floor-balance", days=730)
    assert (t["rivers:stage_floor:cells"] > 0).all()
    r = basin_balance(RUNS / "floor-balance", "cha01")
    p = r["precipitation_land"]
    assert abs(r["land_residual"]) <= 0.02 * (p + r["irrigation"])      # as test_balance: Osu_1hru's wetland
    assert abs(r["inflow_from_elsewhere"]) <= 0.001 * r["outflow"]
    assert abs(r["groundwater_residual"]) <= 1e-6 * r["recharge_received"]
    assert abs(r["channels_residual"]) <= 0.001 * r["outflow"]
    assert r["max_daily_closure_share"] <= 1e-6


@pytest.mark.parametrize("rivers, match", [
    (lambda p: _rivers(p, stage_floor=3), "names a column"),
    (lambda p: _rivers(p, stage="fixed"), "fixed stage is MODFLOW 6's input"),
    (lambda p: _rivers("all:cha01"), "has none"),
    (lambda p: _rivers(p, stage_floor="dem"), "no column 'dem'"),
])
def test_a_bad_floor_is_refused(depth, tmp_path, rivers, match):
    path = _map(tmp_path / "rivers.csv", [HIGH] * 10)
    problems, _ = Basin.check(_cfg(rivers(path)))
    assert any(match in p for p in problems), problems


def test_a_floor_that_isnt_a_number_is_refused(depth, tmp_path):
    path = _map(tmp_path / "rivers.csv", ["high"] + [None] * 9)
    with pytest.raises(ValueError, match="'high' isn't a number"):
        _run(_cfg(_rivers(path)), "floor-bad", days=1)
