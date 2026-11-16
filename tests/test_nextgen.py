"""NextGen as the surface model: a finished NGIAB run replayed day by day, forcing for users, and the Plaster Creek
example where its run folder exists.

The synthetic run here has three flowpaths (wb-1 and wb-2 drain to nex-3, wb-3 to nex-4), 10 days of hourly t-route
output with known values, and hourly forcing for two catchments.
"""

from __future__ import annotations

import datetime as dt
import json
import sqlite3
from pathlib import Path

import numpy as np
import pytest

netCDF4 = pytest.importorskip("netCDF4")

from basinledger import config as C  # noqa: E402
from basinledger.coupling.nextgen import DEPTH, QOUT, NgenReplay  # noqa: E402
from basinledger.driver import Basin  # noqa: E402
from basinledger.forcing import HUMIDITY, TEMPERATURE, nextgen_inputs, relative_humidity  # noqa: E402
from conftest import ROOT, RUNS  # noqa: E402

DAYS = 10
START = dt.datetime(2010, 1, 1)


def flow(i, h):
    """Flowpath i's flow at the end of hour h (1-based), m3 s-1."""
    return (i + 1) * (1.0 + 0.5 * np.sin(2 * np.pi * h / (24 * 4)))


def make_run(folder):
    (folder / "config").mkdir(parents=True)
    (folder / "outputs" / "troute").mkdir(parents=True)
    (folder / "forcings").mkdir()
    (folder / "config" / "realization.json").write_text(json.dumps({"time": {
        "start_time": "2010-01-01 00:00:00", "end_time": "2010-01-11 00:00:00", "output_interval": 3600}}))
    with sqlite3.connect(folder / "config" / "test_subset.gpkg") as db:
        db.execute("create table flowpaths (id text, toid text)")
        db.executemany("insert into flowpaths values (?, ?)", [("wb-1", "nex-3"), ("wb-2", "nex-3"), ("wb-3", "nex-4")])
    hours = np.arange(1, DAYS * 24 + 1)
    with netCDF4.Dataset(folder / "outputs" / "troute" / "troute_output_201001010000.nc", "w") as d:
        d.createDimension("feature_id", 3)
        d.createDimension("time", len(hours))
        d.createVariable("feature_id", "i8", ("feature_id",))[:] = [1, 2, 3]
        t = d.createVariable("time", "f8", ("time",))
        t.units = "seconds since 2010-01-01 00:00:00"
        t[:] = hours * 3600.0
        d.createVariable("flow", "f4", ("feature_id", "time"))[:] = np.array([flow(i, hours) for i in range(3)])
        d.createVariable("depth", "f4", ("feature_id", "time"))[:] = np.array([0.1 * flow(i, hours) for i in range(3)])
    epoch = START.replace(tzinfo=dt.timezone.utc).timestamp()
    th = epoch + 3600.0 * np.arange(DAYS * 24 + 1)
    with netCDF4.Dataset(folder / "forcings" / "forcings.nc", "w") as d:
        d.createDimension("catchment-id", 2)
        d.createDimension("time", len(th))
        ids = d.createVariable("ids", str, ("catchment-id",))
        ids[0], ids[1] = "cat-1", "cat-2"
        d.createVariable("Time", "f8", ("catchment-id", "time"))[:] = np.array([th, th])
        temp = 273.15 + 20 + 10 * np.sin(2 * np.pi * (np.arange(len(th)) - 9) / 24)
        d.createVariable("TMP_2maboveground", "f4", ("catchment-id", "time"))[:] = np.array([temp, temp - 5])
        d.createVariable("SPFH_2maboveground", "f4", ("catchment-id", "time"))[:] = np.full((2, len(th)), 0.008)
        d.createVariable("PRES_surface", "f4", ("catchment-id", "time"))[:] = np.full((2, len(th)), 98000.0)
    return folder


@pytest.fixture
def run(tmp_path):
    return make_run(tmp_path / "ngiab")


def test_replay_gives_daily_means_keyed_to_flowpaths(run):
    m = NgenReplay()
    m.initialize(str(run))
    assert m.get_grid_object_names(0) == ["wb-1", "wb-2", "wb-3"]
    assert m.get_time_units() == "days since 2010-01-01 00:00:00" and m.get_end_time() == DAYS
    assert m.nexus == {"wb-1": "nex-3", "wb-2": "nex-3", "wb-3": "nex-4"}
    q, h = np.empty(3), np.empty(3)
    for d in range(DAYS):
        m.update()
        m.get_value(QOUT, q)
        m.get_value(DEPTH, h)
        hours = np.arange(d * 24 + 1, d * 24 + 25)  # hours ending within day d
        want = [np.float32(flow(i, hours)).astype(float).mean() for i in range(3)]
        np.testing.assert_allclose(q, want, rtol=1e-6)
        np.testing.assert_allclose(h, 0.1 * np.array(want), rtol=1e-5)
    with pytest.raises(RuntimeError):
        m.update()


def test_replay_saves_and_loads_its_day(run, tmp_path):
    m = NgenReplay()
    m.initialize(str(run))
    m.update_until(4)
    m.save_state(str(tmp_path / "s"))
    n = NgenReplay()
    n.initialize(str(run))
    n.load_state(str(tmp_path / "s"))
    a, b = np.empty(3), np.empty(3)
    m.update(), n.update()
    assert np.array_equal(m.get_value(QOUT, a), n.get_value(QOUT, b))


def test_relative_humidity():
    t = np.array([0.0, 20.0, 30.0])
    es = 611.2 * np.exp(17.67 * t / (t + 243.5))
    q_sat = 0.622 * es / (101325 - 0.378 * es)  # specific humidity at saturation
    np.testing.assert_allclose(relative_humidity(q_sat, np.full(3, 101325.0), t), 100.0, rtol=1e-9)
    assert relative_humidity(np.array([0.0]), np.array([1e5]), np.array([10.0]))[0] == 0.0


def test_forcing_for_a_catchment(run):
    f = nextgen_inputs(run, "cat-2", START)
    assert f[TEMPERATURE](0.0) == pytest.approx(20 + 10 * np.sin(2 * np.pi * -9 / 24) - 5, abs=1e-4)
    assert f[TEMPERATURE](3600 * 5 + 1800) == pytest.approx(20 + 10 * np.sin(2 * np.pi * -4 / 24) - 5, abs=1e-4)
    assert 0 < f[HUMIDITY](0.0) <= 100
    with pytest.raises(KeyError):
        nextgen_inputs(run, "cat-9", START)


def config(run, tmp_path, extra=""):
    (tmp_path / "basin.toml").write_text(f'''
name = "ngen-test"
[run]
folder = "{(RUNS / 'ngen-test').as_posix()}"
[nextgen]
run = "{Path(run).as_posix()}"
[[utilities]]
name = "city"
sources = [{{route = "offsite:lake", share = 1.0}}]
outfall = "offsite:plant"
[[users]]
name = "plant"
component = "basinledger.testing:SyntheticUser"
withdrawal = "utility:city"
return = "utility:city"
{extra}
''')
    return C.load(tmp_path / "basin.toml")


def test_driver_on_a_replay(run, tmp_path):
    cfg = config(run, tmp_path, '''
[[rules]]
name = "low flow"
kind = "curtailment"
user = "plant"
channel = "wb-1"
threshold = 1.0
''')
    assert cfg.surface == "nextgen" and cfg.models == ["nextgen"] and cfg.calibration_status == "screening"
    b = Basin(cfg, run="replay")
    b.run()
    t, s = b.ledger.table(), {k: np.array(v) for k, v in b.series.items()}
    out = b.write_outputs()
    b.finalize()
    assert len(s["flow:wb-1"]) == DAYS
    np.testing.assert_allclose(s["flow:nex-3"], s["flow:wb-1"] + s["flow:wb-2"], rtol=1e-12)
    np.testing.assert_allclose(s["flow:nex-4"], s["flow:wb-3"], rtol=1e-12)
    cut = t["rule:low flow:active"].astype(bool)
    assert np.array_equal(cut, np.r_[False, s["flow:wb-1"][:-1] < 1.0]) and 0 < cut.sum() < DAYS
    assert (t["user:plant:withdrawal"][cut] == 0).all()
    np.testing.assert_allclose(t["route:offsite:lake:delivered"], t["utility:city:produced"], rtol=0)
    assert json.loads((out / "summary.json").read_text())["models"] == {"nextgen": "screening"}


def test_cooling_user_on_nextgen_forcing(run, tmp_path):
    pytest.importorskip("waterledger_models.datacenter.bmi")
    (tmp_path / "campus.yaml").write_text("design: towers\nit_capacity_mw: 100\n")
    cfg = config(run, tmp_path, '''
[[users]]
name = "campus"
component = "waterledger_models.datacenter.bmi:CoolingWaterBmi"
config = "campus.yaml"
withdrawal = "utility:city"
return = "utility:city"
forcing = {source = "nextgen", catchment = "cat-1"}
variables = {withdrawal = "data_center_cooling_water__withdrawal_volume_flux", consumption = "data_center_cooling_water__consumption_volume_flux", return = "data_center_cooling_water__return_volume_flux", gain = "data_center_cooling_water__condensate_volume_flux", power = "data_center__electric_power"}
constituents = {}
''')
    b = Basin(cfg, run="campus", days=3)
    campus = [wu for u, wu in b.users if u.name == "campus"][0]
    assert set(campus.inputs) == {TEMPERATURE, HUMIDITY}
    b.run()
    w = b.ledger.table()["user:campus:withdrawal"]
    b.finalize()
    assert (w > 0).all()


def test_what_a_replay_refuses(run, tmp_path):
    cfg = config(run, tmp_path)
    cfg.places_of_use = [("x", "wb-1")]
    with pytest.raises(ValueError, match="no surface routes"):
        Basin(cfg, run="refuse")
    cfg = config(run, tmp_path)
    cfg.mf6_builder = "model:build"
    with pytest.raises(NotImplementedError, match="percolation"):
        Basin(cfg, run="refuse")


PLASTER = ROOT / "examples" / "plaster-creek"


@pytest.fixture(scope="module")
def plaster():
    cfg = C.load(PLASTER / "basin.toml")
    if not (cfg.nextgen_run / "outputs" / "troute").is_dir():
        pytest.skip(f"no NextGen run at {cfg.nextgen_run} (nextgen-michigan's steps 3 to 6; see the example's README)")
    cfg.run_folder = RUNS
    return cfg


def test_plaster_creek_fit_matches_its_record(plaster):
    import importlib.util
    import tomllib

    spec = importlib.util.spec_from_file_location("check_calibration", PLASTER / "check_calibration.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    obs = C.nextgen_michigan() / "runs" / "plaster-cal2" / "calibration" / "obs_hourly_discharge_cms.csv"
    if not obs.exists():
        pytest.skip(f"no gage record at {obs}")
    fit = mod.fit(plaster.nextgen_run, obs)
    rec = tomllib.loads((PLASTER / "calibration" / "record.toml").read_text())["models"]["nextgen"]
    for period in ("calibration", "validation"):
        want = rec[period]["recomputed"]
        assert fit[period]["days"] == want["days"]
        for k in ("kge", "nse", "log_nse"):
            assert fit[period][k] == pytest.approx(want[k], abs=1e-3), (period, k)
        assert fit[period]["pbias"] == pytest.approx(want["pbias"], abs=0.05)


def test_plaster_creek_runs_calibrated(plaster):
    pytest.importorskip("waterledger_models.datacenter.bmi")
    b = Basin(plaster, run="plaster-creek", days=60)
    assert b.status == "calibrated"
    b.run()
    out = b.write_outputs()
    b.finalize()
    s = json.loads((out / "summary.json").read_text())
    assert s["models"] == {"nextgen": "calibrated"} and "USGS 04119055" in s["calibrated_to"]["nextgen"]
    assert "not calibrated" in s["status_note"]  # the users aren't


def add_cfe_and_areas(run, cell_area):
    """CFE outputs with SOIL_TO_GW_FLUX for cat-1 and cat-2, and their areas: 50 cells each of the small grid."""
    hours = DAYS * 24
    (run / "outputs" / "ngen").mkdir(parents=True, exist_ok=True)
    for k, cat in enumerate(("cat-1", "cat-2")):
        rows = ["Time Step,Time,Q_OUT,SOIL_TO_GW_FLUX"]
        for h in range(hours):
            t = START + dt.timedelta(hours=h)
            perc = (k + 1) * 1e-5 * (1 + np.sin(2 * np.pi * h / 48))  # m per hour
            rows.append(f"{h},{t:%Y-%m-%d %H:%M:%S},0,{float(perc)!r}")
        (run / "outputs" / "ngen" / f"{cat}.csv").write_text("\n".join(rows) + "\n")
    with sqlite3.connect(run / "config" / "test_subset.gpkg") as db:
        db.execute("create table divides (divide_id text, areasqkm real)")
        db.executemany("insert into divides values (?, ?)", [("cat-1", 50 * cell_area / 1e6), ("cat-2", 50 * cell_area / 1e6)])


def test_nextgen_with_modflow6(run, tmp_path, libmf6):
    import importlib.util

    spec = importlib.util.spec_from_file_location("m", ROOT / "examples" / "synthetic" / "model.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    cell = m.SMALL_CELL ** 2
    add_cfe_and_areas(run, cell)
    # cat-1 covers columns 0-4, cat-2 columns 5-9; the river (column 0) belongs to wb-1
    (tmp_path / "rch.csv").write_text("object,layer,row,col,area_m2\n" + "".join(
        f"{'cat-1' if c < 5 else 'cat-2'},0,{r},{c},{cell}\n" for r in range(10) for c in range(10)))
    (tmp_path / "riv.csv").write_text("object,layer,row,col\n" + "".join(f"wb-1,0,{r},0\n" for r in range(10)))
    shutil_copy = __import__("shutil").copy
    shutil_copy(ROOT / "examples" / "synthetic" / "model.py", tmp_path / "model.py")
    cfg = config(run, tmp_path, '''
[modflow6]
builder = "model:build_small"
[recharge]
map = "rch.csv"
area_rtol = 1e-9
[rivers]
package = "riv"
map = "riv.csv"
stage = "nextgen_depth"
bed_thickness = 1.0
''')
    b = Basin(cfg, run="ngen-mf6")
    stages = []
    for _ in range(DAYS):
        b.step()
        stages.append(b.mf6.get("STAGE", "GWF", "RIV").copy())
    t, s = b.ledger.table(), {k: np.array(v) for k, v in b.series.items()}
    rbot = b.mf6.get("RBOT", "GWF", "RIV")
    b.finalize()
    # recharge: sent (CFE's percolation over each catchment's area) equals received, every day
    np.testing.assert_allclose(t["recharge:received"], t["recharge:sent"], rtol=1e-12)
    want = sum((k + 1) * 1e-5 * (1 + np.sin(2 * np.pi * np.arange(24) / 48)).sum() * 50 * cell for k in range(2))
    assert t["recharge:sent"][0] == pytest.approx(want, rel=1e-6)  # read through a 4-byte buffer, as for SWAT+
    # stage = RBOT + bed thickness + the day's t-route depth on wb-1
    depth = 0.1 * s["flow:wb-1"]  # the synthetic run's depth is a tenth of its flow
    np.testing.assert_allclose(np.array(stages)[:, 0], rbot[0] + 1.0 + depth, rtol=1e-6)
    assert np.isfinite(s["baseflow:wb-1"]).all() and len(s["baseflow:wb-1"]) == DAYS
