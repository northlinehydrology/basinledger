"""A model's time discretization rewritten to one-day steps in a run's copy (coupling/timesteps.py; docs/design.md,
"Extending an existing MODFLOW 6 model": stress periods).

- tier 1, on TDIS files alone: monthly periods of one step, and of five steps growing by 1.2, become daily steps,
  each period kept and its original row after it; a model already in daily steps is left untouched; periods that
  aren't whole days end on the nearest day, measured from the start, so nothing builds up; a model in hours;
  a period too short for a day of its own refused; units undefined, stated or in years; adaptive time stepping left
  out; the start date read; a one-period model cut for a restart once rewritten; Pleasant Lake's 85 periods;
- restarts across stress periods (``remaining_days``): period blocks renumbered, the one in effect carried as period
  1; time series and time-array series moved back; the synthetic foreign model cut at a day in its second month;
- tier 1 with libmf6: the synthetic foreign model with monthly growing steps, rewritten, runs the same days as its
  daily twin, to the last digit; and restarted in its second month, with a time series on a well, it continues the
  unbroken run to the last digit.
"""

from __future__ import annotations

import datetime as dt
import importlib.util
import shutil

import numpy as np
import pytest
from conftest import ROOT, RUNS

from basinledger.coupling.timesteps import (
    _rekey,
    _shift_times,
    daily_steps,
    remaining_days,
    start_date,
    tdis_file,
)

flopy = pytest.importorskip("flopy")
spec = importlib.util.spec_from_file_location("existing", ROOT / "examples" / "synthetic" / "existing.py")
E = importlib.util.module_from_spec(spec)
spec.loader.exec_module(E)
START = dt.date(2010, 1, 1)

NAM = """BEGIN options
END options

BEGIN timing
  TDIS6  model.tdis
END timing
"""


def tdis(tmp_path, rows, units="days", extra=""):
    """A simulation folder holding only mfsim.nam and a TDIS file with these (perlen, nstp, tsmult) rows."""
    (tmp_path / "mfsim.nam").write_text(NAM)
    opts = (f"  TIME_UNITS  {units}\n" if units else "") + "  START_DATE_TIME  2012-01-01\n" + extra
    body = "\n".join(f"  {p} {n} {s}  to somewhere" for p, n, s in rows)
    (tmp_path / "model.tdis").write_text(f"BEGIN options\n{opts}END options\n\nBEGIN dimensions\n  NPER  {len(rows)}\n"
                                         f"END dimensions\n\nBEGIN perioddata\n{body}\nEND perioddata\n")
    return tmp_path


def periods(ws):
    """The rewritten (perlen, nstp, tsmult) rows."""
    from basinledger.coupling.timesteps import _read
    return [(p, n, s) for _, p, n, s in _read(tdis_file(ws))[2]]


def test_monthly_periods_of_one_step_become_daily_steps(tmp_path):
    ws = tdis(tmp_path, [(1.0, 1, 1.0), (31.0, 1, 1.0), (29.0, 1, 1.0), (31.0, 1, 1.0)])
    notes = daily_steps(ws)
    assert periods(ws) == [(1.0, 1, 1.0), (31.0, 31, 1.0), (29.0, 29, 1.0), (31.0, 31, 1.0)]
    assert notes == [("MODFLOW 6's time steps rewritten to one day each in the run's copy: 4 stress periods, 92 days "
                      "(they had 4 steps); each period's input stays as given")]
    text = (ws / "model.tdis").read_text()
    assert "  31.0 31 1.0  was 31.0 1 1.0 to somewhere" in text and "NPER  4" in text
    assert text.startswith("# Basin Ledger: every time step one day")


def test_growing_steps_become_daily_steps(tmp_path):
    ws = tdis(tmp_path, [(1.0, 1, 1.0), (31.0, 5, 1.2), (30.0, 5, 1.2)])
    notes = daily_steps(ws)
    assert periods(ws) == [(1.0, 1, 1.0), (31.0, 31, 1.0), (30.0, 30, 1.0)]
    assert "62 days (they had 11 steps, growing in 2 periods)" in notes[0]


def test_flopy_and_modflow_read_the_rewritten_file(tmp_path):
    ws = RUNS / "timesteps-monthly-text"
    if ws.exists():
        shutil.rmtree(ws)
    E.build(ws, 59, START, daily=False, steps=5, tsmult=1.2)
    daily_steps(ws)
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(ws), verbosity_level=0, load_only=["tdis"])
    assert [tuple(r) for r in sim.tdis.perioddata.get_data()] == [(31.0, 31, 1.0), (28.0, 28, 1.0)]


def test_a_model_in_daily_steps_is_left_untouched(tmp_path):
    ws = tdis(tmp_path, [(1.0, 1, 1.5), (31.0, 31, 1.0)])  # one step: its multiplier means nothing
    before = (ws / "model.tdis").read_text()
    assert daily_steps(ws) == []
    assert (ws / "model.tdis").read_text() == before


def test_periods_not_whole_days_end_on_the_nearest_day_with_nothing_building_up(tmp_path):
    ws = tdis(tmp_path, [(30.4375, 1, 1.0)] * 12)  # a year of 365.25 days in equal months
    notes = daily_steps(ws)
    lengths = [n for _, n, _ in periods(ws)]
    ends = np.cumsum(lengths)
    assert all(abs(e - 30.4375 * (k + 1)) <= 0.5 for k, e in enumerate(ends))  # never more than 12 hours off
    assert ends[-1] == 365 and set(lengths) == {30, 31}
    assert "12 of 12 period ends moved to the nearest day, at most 12 hours" in notes[1]


def test_a_model_in_hours(tmp_path):
    ws = tdis(tmp_path, [(36.0, 1, 1.0), (36.0, 3, 1.0), (24.0, 1, 1.0)], units="hours")
    notes = daily_steps(ws)
    assert periods(ws) == [(48.0, 2, 1.0), (24.0, 1, 1.0), (24.0, 1, 1.0)]  # ends at 1.5 -> 2, 3, 4 days
    assert "1 of 3 period ends moved to the nearest day, at most 12 hours" in notes[1]


def test_a_period_too_short_for_a_day_of_its_own_is_refused(tmp_path):
    ws = tdis(tmp_path, [(31.0, 1, 1.0), (0.25, 1, 1.0), (30.0, 1, 1.0)])
    before = (ws / "model.tdis").read_text()
    with pytest.raises(ValueError, match="period 2 is 0.25 days, too short to keep a whole day"):
        daily_steps(ws)
    assert (ws / "model.tdis").read_text() == before


def test_time_units_undefined_stated_or_years(tmp_path):
    ws = tdis(tmp_path, [(31.0, 1, 1.0)], units=None)
    with pytest.raises(ValueError, match="time units are undefined .* \\[modflow6\\] time_units"):
        daily_steps(ws)
    assert daily_steps(ws, "days") and periods(ws) == [(31.0, 31, 1.0)]
    ws = tdis(tmp_path, [(86400.0 * 2, 1, 1.0)], units=None)
    daily_steps(ws, "seconds")
    assert periods(ws) == [(172800.0, 2, 1.0)]
    ws = tdis(tmp_path, [(1.0, 1, 1.0)], units="years")
    with pytest.raises(ValueError, match="in years .* a one-day step isn't a whole number of them"):
        daily_steps(ws)


def test_adaptive_time_stepping_is_left_out(tmp_path):
    ws = tdis(tmp_path, [(31.0, 31, 1.0)], extra="  ATS6  FILEIN  model.ats\n")
    notes = daily_steps(ws)
    assert notes == ["MODFLOW 6's adaptive time stepping (ATS6) left out: every step is a day"]
    assert "# ATS6  FILEIN  model.ats  (adaptive time stepping left out" in (ws / "model.tdis").read_text()


def test_the_start_date(tmp_path):
    assert start_date(tdis(tmp_path, [(1.0, 1, 1.0)])) == dt.date(2012, 1, 1)
    (tmp_path / "model.tdis").write_text((tmp_path / "model.tdis").read_text().replace("2012-01-01", "2012-01-01T06"))
    assert start_date(tmp_path) == dt.date(2012, 1, 1)
    (tmp_path / "model.tdis").write_text((tmp_path / "model.tdis").read_text().replace("START_DATE_TIME", "#"))
    assert start_date(tmp_path) is None


def test_a_one_period_model_restarts_once_rewritten():
    from basinledger.coupling.modflow6 import remaining_simulation
    ws = RUNS / "timesteps-one-period"
    if ws.exists():
        shutil.rmtree(ws)
    E.build(ws, 31, START, daily=False)  # January in one step
    with pytest.raises(ValueError, match="steps aren't one day each; a run rewrites them first"):
        remaining_simulation(ws, RUNS / "timesteps-one-period-restart", 10)
    daily_steps(ws)
    remaining_simulation(ws, ws, 10)  # in place, in a run's copy
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(ws), verbosity_level=0, load_only=["tdis"])
    assert [tuple(r) for r in sim.tdis.perioddata.get_data()] == [(21.0, 21, 1.0)]


def test_pleasant_lakes_monthly_periods():
    spec = importlib.util.spec_from_file_location("cs_fetch", ROOT / "examples" / "central-sands" / "fetch.py")
    F = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(F)
    src = F.workspace()
    if not (src / "mfsim.nam").exists():
        pytest.skip("Pleasant Lake not fetched: uv run python examples/central-sands/fetch.py")
    ws = RUNS / "timesteps-pleasant"
    if ws.exists():
        shutil.rmtree(ws)
    ws.mkdir(parents=True)
    for f in ("mfsim.nam", tdis_file(src).name):
        shutil.copy(src / f, ws / f)
    before = periods(ws)
    assert before[0] == (1.0, 1, 1.0) and before[1] == (31.0, 5, 1.2)  # a steady day, then months of 5 steps
    notes = daily_steps(ws)
    after = periods(ws)
    assert len(after) == 85 and all(n == p and s == 1.0 for p, n, s in after)
    assert [p for p, _, _ in after] == [p for p, _, _ in before]  # every period its own length
    assert start_date(ws) == dt.date(2012, 1, 1)
    assert notes == [(f"MODFLOW 6's time steps rewritten to one day each in the run's copy: 85 stress periods, "
                      f"{int(sum(p for p, _, _ in before)):,} days (they had 421 steps, growing in 84 periods); "
                      f"each period's input stays as given")]


def run_alone(ws, days):
    """Heads after each day, MODFLOW 6 as given (its own recharge and wells)."""
    from basinledger.coupling.modflow6 import Modflow6
    mf = Modflow6(ws, "kzoo")
    try:
        out = []
        for _ in range(days):
            mf.prepare_day()
            mf.solve_day()
            out.append(mf.head().copy())
        return np.array(out)
    finally:
        mf.finalize()


def test_growing_monthly_steps_rewritten_run_as_the_daily_twin(libmf6):
    daily, monthly = RUNS / "timesteps-daily", RUNS / "timesteps-monthly"
    for ws in (daily, monthly):
        if ws.exists():
            shutil.rmtree(ws)
    E.build(daily, 59, START)
    E.build(monthly, 59, START, daily=False, steps=5, tsmult=1.2)
    from basinledger.coupling.modflow6 import Modflow6
    mf = Modflow6(monthly, "kzoo")
    try:
        with pytest.raises(ValueError, match="time step is 4.16"):  # 31 days in 5 steps growing by 1.2: 4.17 to 8.6
            mf.prepare_day()
    finally:
        mf.finalize()
    assert daily_steps(daily) == []
    assert daily_steps(monthly)
    a, b = run_alone(daily, 59), run_alone(monthly, 59)
    assert np.array_equal(a, b)  # the same simulation: January and February, day by day


# ---------- restarts across stress periods ----------

PACKAGE = """BEGIN options
END options

BEGIN period  1
  a
END period  1

BEGIN PERIOD 3
  c
END PERIOD 3

begin period 5  # five
  e
end period
"""


@pytest.mark.parametrize("first, kept", [
    (1, None),
    (2, ["BEGIN period  1", "  a", "END period  1", "BEGIN PERIOD 2", "END PERIOD 2", "begin period 4  # five"]),
    (3, ["BEGIN PERIOD 1", "  c", "END PERIOD 1", "begin period 3  # five"]),
    (4, ["BEGIN PERIOD 1", "  c", "END PERIOD 1", "begin period 2  # five"]),
    (6, ["begin period 1  # five"]),  # the last block, in effect from period 5 on
])
def test_period_blocks_renumbered(tmp_path, first, kept):
    f = tmp_path / "x.wel"
    f.write_text(PACKAGE)
    assert _rekey(f, first) == (kept is not None)
    text = f.read_text()
    if kept is None:
        assert text == PACKAGE
        return
    for line in kept:
        assert line in text.splitlines()
    if first >= 3:
        assert "  a" not in text  # left out: it was over before the restart
    if first == 6:
        assert "  e" in text and "  c" not in text and text.count("BEGIN") + text.count("begin") == 2


def test_time_series_moved_back(tmp_path):
    ts = tmp_path / "q.ts"
    ts.write_text("BEGIN attributes\n  NAMES q\n  METHODS linear\nEND attributes\n\nBEGIN timeseries\n"
                  "  0.0 -25.0\n  40.5  -60.0\n  7.5d1 -10.0\nEND timeseries\n")
    assert _shift_times(ts, 40.0)
    rows = [line.split() for line in ts.read_text().splitlines() if line.startswith("  ") and "NAMES" not in line
            and "METHODS" not in line]
    assert rows == [["-40.0", "-25.0"], ["0.5", "-60.0"], ["35.0", "-10.0"]]
    tas = tmp_path / "r.tas"
    tas.write_text("BEGIN attributes\n  NAME r\n  METHOD linear\nEND attributes\n\nBEGIN TIME 0.0\n  CONSTANT 1.0\n"
                   "END TIME\n\nBEGIN TIME 60.0\n  CONSTANT 2.0\nEND TIME\n")
    _shift_times(tas, 40.0)
    assert "BEGIN TIME -40.0" in tas.read_text() and "BEGIN TIME 20.0" in tas.read_text()


def with_a_time_series(ws):
    """The synthetic foreign model's second city well at a rate from a time series rather than its own -25."""
    wel = ws / "kzoo.wel_city.wel"
    text = wel.read_text().replace("BEGIN options\n", "BEGIN options\n  TS6  FILEIN  kzoo.wel.ts\n")
    wel.write_text(text.replace("-2.50000000E+01", "qcity"))
    (ws / "kzoo.wel.ts").write_text("BEGIN attributes\n  NAMES qcity\n  METHODS linear\nEND attributes\n\n"
                                    "BEGIN timeseries\n  0.0 -25.0\n  40.5 -60.0\n  75.0 -10.0\nEND timeseries\n")
    return ws


def test_the_days_left_from_a_day_in_the_second_month():
    ws = RUNS / "timesteps-remaining"
    if ws.exists():
        shutil.rmtree(ws)
    with_a_time_series(E.build(ws, 75, START))  # January, February and half of March, daily
    changed = remaining_days(ws, 40)  # February 10th
    assert periods(ws) == [(19.0, 19, 1.0), (16.0, 16, 1.0)]
    assert "NPER  2" in (ws / "kzoo_regional.tdis").read_text()
    assert {"kzoo_regional.tdis", "kzoo.wel_city.wel", "kzoo.rch_swb.rcha", "kzoo.wel.ts"} <= set(changed)
    wel = (ws / "kzoo.wel_city.wel").read_text()
    assert "BEGIN period  1" in wel and "BEGIN period  2" in wel and "period  3" not in wel
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(ws), verbosity_level=0)
    assert sim.get_model("kzoo").rch_swb.recharge.get_data(0) is not None


def test_a_restart_in_the_second_month_continues_the_unbroken_run(libmf6):
    from basinledger.coupling.modflow6 import Modflow6, remaining_simulation
    src = RUNS / "timesteps-restart-src"
    if src.exists():
        shutil.rmtree(src)
    with_a_time_series(E.build(src, 75, START, daily=False, steps=5, tsmult=1.2))
    daily_steps(src)
    mf = Modflow6(src, "kzoo")
    try:
        unbroken = []
        for _ in range(75):
            mf.prepare_day()
            mf.solve_day()
            unbroken.append((mf.head().copy(), mf.simvals("wel_city")[1]))
    finally:
        mf.finalize()
    dst = remaining_simulation(src, RUNS / "timesteps-restart-dst", 40)
    mf = Modflow6(dst, "kzoo")
    try:
        mf.set_heads(unbroken[39][0])
        assert mf.end == 35.0
        for d in range(40, 75):
            mf.prepare_day()
            mf.solve_day()
            assert np.array_equal(mf.head(), unbroken[d][0]), d
            assert mf.simvals("wel_city")[1] == unbroken[d][1], d
    finally:
        mf.finalize()
    # day 41 (days 40 to 41): the series averaged over the step, through its point at 40.5
    a, b = -25.0 - 35.0 * 40.0 / 40.5, -60.0 + 50.0 * 0.5 / 34.5
    assert unbroken[40][1] == pytest.approx(0.25 * (a - 60.0) + 0.25 * (b - 60.0), rel=1e-9)
