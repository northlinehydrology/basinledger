"""Weather written into SWAT+'s gage records (``[weather]``, weather.py; swatplus-bmi 0.6's ``set_weather``).

Units: CSV series, a scenario's change (factors, shifts, monthly values, humidity capped), gridded values to gages
(nearest cell; the cells under each gage's HRUs, by area), the config's checks. On Osu_1hru (one gage per variable,
SWAT+ alone):

- a series equal to the project's own records gives the run the files give, bit for bit;
- a change reaches SWAT+: the gage's record read back each day is the changed value from the first day that can be
  set (the third for precipitation, the second for the rest), the files' before;
- a case from a saved warm-up takes its own weather from two days after its first (precipitation), the warm-up's
  before, and a case with the warm-up's weather matches the unbroken run;
- a user reading the station's weather gets the series, not the files;
- a series that ends before the run stops it, rather than leaving SWAT+ on its files.
"""

from __future__ import annotations

import datetime as dt
import json

import numpy as np
import pytest

from basinledger import config as C
from basinledger import weather as W
from basinledger.driver import Basin
from basinledger.scenarios import run_case
from conftest import RUNS

DAY = 86400.0


# --- units ---------------------------------------------------------------------------------------------------------

def test_a_series_reads_back_from_csv(tmp_path):
    s = W.Series(dt.date(2010, 1, 30), ["a.pcp", "b.pcp"], np.array([[1.0, 2.0], [0.0, 3.25], [7.5, 0.0]]))
    W.write_csv(tmp_path / "pcp.csv", s)
    t = W.read_csv(tmp_path / "pcp.csv")
    assert t.start == s.start and t.end == dt.date(2010, 2, 1) and t.gages == s.gages
    np.testing.assert_array_equal(t.values, s.values)
    np.testing.assert_array_equal(t.columns(["b.pcp"]).values[:, 0], [2.0, 3.25, 0.0])
    with pytest.raises(KeyError, match="c.pcp"):
        t.columns(["c.pcp"])
    with pytest.raises(ValueError, match="2010-01-29"):
        t.between(dt.date(2010, 1, 29), dt.date(2010, 1, 31))


def test_a_series_with_a_gap_is_refused(tmp_path):
    (tmp_path / "x.csv").write_text("date,a\n2010-01-01,1\n2010-01-03,2\n")
    with pytest.raises(ValueError, match="every day"):
        W.read_csv(tmp_path / "x.csv")


def test_a_change_scales_shifts_and_caps():
    days = 400
    s = W.Series(dt.date(2010, 1, 1), ["g"], np.full((days, 1), 0.9))
    monthly = [1.0 + m / 10 for m in range(12)]
    dates = s.dates()
    pc = W.changed("precipitation", s, {"precipitation": monthly}).values[:, 0]
    np.testing.assert_allclose(pc, [0.9 * monthly[d.month - 1] for d in dates])
    t = W.changed("temperature_min", s, {"temperature": 2.0}).values[:, 0]
    np.testing.assert_allclose(t, 2.9)
    h = W.changed("humidity", s, {"humidity": 1.5}).values[:, 0]
    np.testing.assert_array_equal(h, 1.0)                        # a fraction
    assert W.changed("wind", s, {"precipitation": 2.0}) is s      # nothing for wind
    with pytest.raises(ValueError, match="12 monthly"):
        W.changed("wind", s, {"wind": [1.0, 2.0]})


def test_gridded_values_go_to_the_nearest_cell_or_by_area():
    # four cells; gage A sits nearest cell 1, gage B nearest cell 3
    cell_lon, cell_lat = np.array([-85.0, -84.9, -85.0, -84.9]), np.array([44.0, 44.0, 43.9, 43.9])
    v = np.array([[1.0, 2.0, 3.0, 4.0], [10.0, 20.0, 30.0, 40.0]])
    s = W.from_grid(dt.date(2010, 1, 1), v, cell_lon, cell_lat, ["A", "B"], gage_lon=[-84.91, -84.89],
                    gage_lat=[44.01, 43.88])
    np.testing.assert_array_equal(s.values, [[2.0, 4.0], [20.0, 40.0]])
    # by area: A's HRUs are 3 ha near cell 0 and 1 ha near cell 2
    pts = {"A": ([-85.0, -85.0], [44.0, 43.9], [3e4, 1e4]), "B": ([-84.9], [43.9], [5e4])}
    s = W.from_grid(dt.date(2010, 1, 1), v, cell_lon, cell_lat, ["A", "B"], points=pts)
    np.testing.assert_allclose(s.values, [[0.75 * 1 + 0.25 * 3, 4.0], [0.75 * 10 + 0.25 * 30, 40.0]])


def test_the_config_checks_weather(tmp_path):
    base = '[swatplus]\nproject = "p"\n'
    (tmp_path / "a.toml").write_text('name = "a"\n' + base + '[weather]\nrain = "x.csv"\n')
    with pytest.raises(ValueError, match="rain"):
        C.load(tmp_path / "a.toml")
    (tmp_path / "b.toml").write_text('name = "b"\n' + base + '[weather.change]\nwind = 1.1\nsnow = 2\n')
    with pytest.raises(ValueError, match="snow"):
        C.load(tmp_path / "b.toml")
    (tmp_path / "c.toml").write_text('name = "c"\n[nextgen]\nrun = "r"\n[weather.change]\nprecipitation = 1.1\n')
    with pytest.raises(ValueError, match="NextGen"):
        C.load(tmp_path / "c.toml")
    (tmp_path / "d.toml").write_text('name = "d"\n' + base + '[weather]\nprecipitation = "w/pcp.csv"\n'
                                     '[weather.change]\ntemperature = [1,2,3,4,5,6,7,8,9,10,11,12]\n')
    cfg = C.load(tmp_path / "d.toml")
    assert cfg.weather["precipitation"] == tmp_path / "w" / "pcp.csv"
    assert cfg.weather["change"] == {"temperature": list(range(1, 13))}


# --- on Osu_1hru -----------------------------------------------------------------------------------------------------

PCP_FILE, TMP_FILE = "Imsilpcp.pcp", "Imsiltmp.tmp"


def osu_config(tmp_path, swatplus, name, weather=None, days=0) -> C.ExampleConfig:
    (tmp_path / f"{name}.toml").write_text(f'name = "{name}"\n[swatplus]\nproject = "{swatplus / "Osu_1hru"}"\n')
    cfg = C.load(tmp_path / f"{name}.toml")
    cfg.run_folder = RUNS
    cfg.weather = weather
    cfg.days = days
    return cfg


def run_days(cfg, days, read=(W.PCP,), **kw):
    """Run ``days`` and read the gages' records back each day (precipitation in mm/day)."""
    b = Basin(cfg, run=cfg.name, days=days, **kw)
    got = {v: [] for v in read}
    flow = []
    for _ in range(days - b.day):
        b.step()
        for v in read:
            got[v].append(float(b.sw.get(v)[0]) * (86_400_000.0 if v == W.PCP else 1.0))
        flow.append(float(b.sw.get("channel_water__volume_outflow_rate")[0]))
    return b, {v: np.array(x) for v, x in got.items()}, np.array(flow)


def file_values(swatplus, gage, start, days, column=0):
    return W.gage_records(swatplus / "Osu_1hru", gage, column).between(start, start + dt.timedelta(days - 1)).values[:, 0]


@pytest.mark.slow
def test_the_projects_own_records_give_the_files_run(swatplus, tmp_path):
    days = 120
    for kind in ("precipitation", "temperature_max", "temperature_min"):
        W.write_csv(tmp_path / f"{kind}.csv", W.project_series(swatplus / "Osu_1hru", kind,
                                                               [PCP_FILE if kind == "precipitation" else TMP_FILE]))
    spec = {k: tmp_path / f"{k}.csv" for k in ("precipitation", "temperature_max", "temperature_min")} | {"change": {}}
    b0, _, q0 = run_days(osu_config(tmp_path, swatplus, "wx-files"), days)
    b0.finalize()
    b1, _, q1 = run_days(osu_config(tmp_path, swatplus, "wx-series", spec), days)
    b1.write_outputs()
    b1.finalize()
    np.testing.assert_array_equal(q1, q0)
    s = json.loads((b1.folder / "outputs" / "summary.json").read_text())["weather"]
    assert s["precipitation"]["from"] == "2010-01-03" and s["temperature_max"]["from"] == "2010-01-02"
    assert s["precipitation"]["to"] == "2010-04-30" and s["precipitation"]["gages"] == [PCP_FILE]


@pytest.mark.slow
def test_a_change_reaches_swat(swatplus, tmp_path):
    days = 120
    start = dt.date(2010, 1, 1)
    cfg = osu_config(tmp_path, swatplus, "wx-wetter", {"change": {"precipitation": 1.5, "temperature": 2.0}})
    b, got, q = run_days(cfg, days, read=(W.PCP, W.TMAX))
    b.finalize()
    pcp, tmax = file_values(swatplus, PCP_FILE, start, days), file_values(swatplus, TMP_FILE, start, days)
    np.testing.assert_allclose(got[W.PCP][:2], pcp[:2], rtol=1e-6, atol=1e-6)       # read before any set
    np.testing.assert_allclose(got[W.PCP][2:], 1.5 * pcp[2:], rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(got[W.TMAX][:1], tmax[:1], atol=1e-5)
    np.testing.assert_allclose(got[W.TMAX][1:], tmax[1:] + 2.0, atol=1e-5)
    b0, _, q0 = run_days(osu_config(tmp_path, swatplus, "wx-base"), days)
    b0.finalize()
    assert pcp[2:].sum() > 0 and q.sum() > q0.sum()


@pytest.mark.slow
def test_a_case_from_a_warm_up_takes_its_own_weather(swatplus, tmp_path):
    days, warm = 90, 40
    start = dt.date(2010, 1, 1)
    from basinledger.scenarios import warm_up
    wet = {"change": {"precipitation": 1.5}}
    cfg = osu_config(tmp_path, swatplus, "wx-warm", wet, days=days)
    with pytest.warns(Warning):  # Osu_1hru's soils start empty and its warm-up is short
        state = warm_up(cfg, warm)
    pcp = file_values(swatplus, PCP_FILE, start, days)

    def dry(c):
        c.weather = {"change": {"precipitation": 0.0}}

    case = run_case(cfg, state, "wx-dry", users=False, change=dry)
    assert case.weather["precipitation"]["from"] == start + dt.timedelta(days=warm + 2)
    case.finalize()
    # the same case, stepped and read: its first two days are the warm-up's records (wetter), then dry
    b, got, _ = run_days(osu_config(tmp_path, swatplus, "wx-dry-read", {"change": {"precipitation": 0.0}},
                                    days=days), days, restart=state)
    b.finalize()
    np.testing.assert_allclose(got[W.PCP][:2], 1.5 * pcp[warm:warm + 2], rtol=1e-6, atol=1e-6)
    np.testing.assert_array_equal(got[W.PCP][2:], 0.0)
    # a case with the warm-up's weather is the unbroken run with it
    same = run_case(cfg, state, "wx-same", users=False)
    q_case = np.array(same.series["flow:cha01"])
    same.finalize()
    whole = Basin(cfg, run="wx-whole", users=False)
    whole.run()
    np.testing.assert_array_equal(q_case, np.array(whole.series["flow:cha01"]))
    whole.finalize()


def test_a_user_reads_the_series(swatplus, tmp_path):
    from basinledger.forcing import TEMPERATURE, swatplus_inputs
    start = dt.date(2010, 1, 1)
    proj = swatplus / "Osu_1hru"
    tmax = W.gage_records(proj, TMP_FILE, 0).between(start, start + dt.timedelta(9))
    tmin = W.gage_records(proj, TMP_FILE, 1).between(start, start + dt.timedelta(9))
    weather = {"temperature_max": {"series": W.changed("temperature_max", tmax, {"temperature": 3.0}).between(
                   start + dt.timedelta(1), start + dt.timedelta(9))},
               "temperature_min": {"series": W.changed("temperature_min", tmin, {"temperature": 3.0}).between(
                   start + dt.timedelta(1), start + dt.timedelta(9))}}
    t0 = dt.datetime(2010, 1, 1)
    base = swatplus_inputs(proj, "s35610n127290e", t0)[TEMPERATURE]
    warm = swatplus_inputs(proj, "s35610n127290e", t0, weather)[TEMPERATURE]
    assert warm(3600.0 * 15) == base(3600.0 * 15)                     # the first day is the file's
    for hour in range(24, 240, 7):
        assert warm(3600.0 * hour) == pytest.approx(base(3600.0 * hour) + 3.0)


def test_hru_points_name_each_gages_hrus(swatplus):
    pts = W.hru_points(swatplus / "Osu_1hru", "pcp")
    assert list(pts) == [PCP_FILE]
    lon, lat, area = pts[PCP_FILE]
    assert len(area) == 1 and area[0] > 0 and -180 <= lon[0] <= 180 and -90 <= lat[0] <= 90
    assert W.hru_points(swatplus / "Osu_1hru", "hmd") == {}           # humidity is generated


@pytest.mark.slow
def test_a_series_that_ends_early_stops_the_run(swatplus, tmp_path):
    s = W.project_series(swatplus / "Osu_1hru", "precipitation", [PCP_FILE])
    W.write_csv(tmp_path / "short.csv", s.between(dt.date(2010, 1, 1), dt.date(2010, 1, 20)))
    cfg = osu_config(tmp_path, swatplus, "wx-short", {"precipitation": tmp_path / "short.csv", "change": {}})
    b = Basin(cfg, run=cfg.name, days=40)
    try:
        with pytest.raises(ValueError, match="ends 2010-01-20"):
            b.run()
        assert b.day == 20
    finally:
        b.finalize()
