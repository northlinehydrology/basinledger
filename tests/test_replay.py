"""Replay (replay.py): the salt run again from a record of a coupled run's flow side.

- the record's storage: chunks of days round-trip exactly, arrays of changing length included;
- a replay from a float64 record is the full coupled run bit for bit, every day: both species' concentrations, every
  ledger line the replay books (the salt's, the transport's, the users'), the series;
- from a float32 record (the default), the same to float32's rounding of the flows: concentrations within a relative
  1e-6 (2^-24 is 6e-8 per value; seen: 3e-9), masses within 1e-6 of each account's scale;
- weathering in a top layer that dries and wets again replays bit for bit (the record keeps the saturation);
- a salt multiplier changes only the salt: the other species, the flows' series and the users' water stay as they
  were; the transport model's porosity and weathering rate may change; every account closes every day;
- fewer mantissa bits round each value to nearest; a record without the specific discharge replays without
  dispersion and refuses it with;
- anything that changes the flow is refused: a SWAT+ or MODFLOW 6 parameter, the flow model's files, the users, more
  days than the record holds, starting heads;
- restarts: a replay saved and restarted is the unbroken replay bit for bit, and a replay restarted from a full run's
  saved state is the unbroken full run; a recording run restarted carries its record on, day for day the same;
- the forward model (PEST++) replays.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from basinledger import replay as RP
from basinledger.calibration import Parameter
from conftest import RUNS
from test_salt import salt_config

DAYS = 365


def full_run(name, days, record=None, precision="float64", users=True, save_at=None, restart=None, **changes):
    from basinledger.driver import Basin
    cfg = salt_config(**changes)
    if record is not None:
        cfg.replay_record, cfg.replay_precision = RUNS / record, precision
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run=name, days=days, users=users, restart=restart)
    xs = []
    try:
        while b.day < (save_at or days):
            b.step()
            xs.append({sp: t.concentration().copy() for sp, t in b.transport.items()})
        if save_at:
            b.save(b.folder / "state")
        return {"table": b.ledger.table(), "series": {k: list(v) for k, v in b.series.items()}, "x": xs,
                "folder": b.folder}
    finally:
        b.finalize()


def replay_run(name, record, days=None, parameters=None, save_at=None, restart=None, **changes):
    from basinledger.driver import Basin
    cfg = salt_config(**changes)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run=name, days=days, parameters=parameters, restart=restart, replay=RUNS / record)
    xs = []
    try:
        stop = save_at or b.days
        while b.day < stop:
            b.step()
            xs.append({sp: t.concentration().copy() for sp, t in b.transport.items()})
        if save_at:
            b.save(b.folder / "state")
        return {"table": b.ledger.table(), "series": {k: list(v) for k, v in b.series.items()}, "x": xs,
                "times": dict(b.times), "users": b.users, "salt": b.salt}
    finally:
        b.finalize()


@pytest.fixture(scope="module")
def recorded64(swatplus):
    return full_run("replay-full64", DAYS, "record64")


@pytest.fixture(scope="module")
def recorded32(swatplus):
    return full_run("replay-full32", DAYS, "record32", precision="float32")


@pytest.fixture(scope="module")
def replayed64(recorded64):
    return replay_run("replay-64", "record64")


# ---------- storage ----------

def test_chunks_round_trip(tmp_path):
    rng = np.random.default_rng(0)
    days = [{"a": rng.random(50).astype(np.float32), "b": rng.random(7), "n": np.arange(i % 3, dtype=np.int32),
             "json": {"x": float(i), "s": {"p": [i, 2]}}} for i in range(9)]
    RP.write_chunk(tmp_path / "c.npz", days)
    back = RP.read_chunk(tmp_path / "c.npz")
    for d, e in zip(days, back):
        for k in ("a", "b", "n"):
            assert e[k].dtype == d[k].dtype
            np.testing.assert_array_equal(e[k], d[k])
        assert e["json"] == d["json"]


def test_fewer_mantissa_bits_round_to_nearest():
    x = np.random.default_rng(1).lognormal(0.0, 5.0, 10000) * np.where(np.arange(10000) % 2, 1.0, -1.0)
    y = RP.Precision("float32", 16)(x)
    assert y.dtype == np.float32
    assert np.max(np.abs(y - x) / np.abs(x)) <= 2.0 ** -17 * (1 + 1e-6)
    assert np.all((y.view(np.uint32) & 0x7F) == 0)
    np.testing.assert_array_equal(RP.Precision("float32")(x), x.astype(np.float32))


# ---------- exactness ----------

def same(a, b, key):
    np.testing.assert_array_equal(np.asarray(a, float), np.asarray(b, float), err_msg=key)


def test_float64_replay_is_the_full_run_bit_for_bit(recorded64, replayed64):
    full, rep = recorded64, replayed64
    assert len(rep["x"]) == DAYS
    for d in range(DAYS):
        for sp in ("chloride", "weathering"):
            same(full["x"][d][sp], rep["x"][d][sp], f"{sp} day {d}")
    booked = [k for k in rep["table"] if k.split(":")[0] in ("salt", "transport", "user")]
    assert len(booked) > 60
    for k in booked:
        same(full["table"][k], rep["table"][k], k)
    for k, v in rep["series"].items():
        same(full["series"][k], v, k)


def test_float32_replay_is_the_full_run_to_float32s_rounding(recorded32):
    full, rep = recorded32, replay_run("replay-32", "record32")
    for d in range(DAYS):
        for sp in ("chloride", "weathering"):
            np.testing.assert_allclose(rep["x"][d][sp], full["x"][d][sp], rtol=1e-6, atol=1e-15)
    for k, v in rep["table"].items():
        if k.split(":")[0] not in ("salt", "transport", "user") or k.endswith("closure_share"):
            continue  # a closure is a residual, checked on its own tolerance in the replay as in the run
        a, b = np.asarray(full["table"][k], float), np.asarray(v, float)
        assert np.nanmax(np.abs(a - b)) <= 1e-6 * np.nanmax(np.abs(a)) + 1e-9, k  # kg (a store that's ~1e-16)
    for k in ("chloride_mg_l:cha01", "weathering_mg_l:cha01", "conductance_us_cm:cha01"):
        np.testing.assert_allclose(rep["series"][k], full["series"][k], rtol=1e-6)


def test_every_account_closes_every_day_in_a_replay(replayed64):
    t = replayed64["table"]
    for sp in ("chloride", "weathering"):
        k = f"salt:{sp}:"
        assert np.nanmax(np.abs(t[k + "land:closure_share"])) < 1e-12
        assert np.nanmax(np.abs(t[k + "streams:closure_share"])) < 1e-12
        assert np.nanmax(np.abs(t[k + "closure_share"])) < 1e-9
        assert np.nanmax(np.abs(t[f"transport:{sp}:closure_share"])) < 1e-6


def test_weathering_in_cells_that_dry_replays_bit_for_bit(swatplus):
    dry = {"mf6_builder_args": dict(salt_config().mf6_builder_args, top=90.0, upper=88.2, newton=True, strt=88.6)}
    full = full_run("replay-dry-full", 160, "record-dry", users=False, **dry)
    rep = replay_run("replay-dry", "record-dry", **dry)
    for d in range(160):
        same(full["x"][d]["weathering"], rep["x"][d]["weathering"], f"day {d}")
    same(full["table"]["transport:weathering:boundary:weathering"], rep["table"]["transport:weathering:boundary:weathering"],
         "weathering")


# ---------- what may change, and what may not ----------

def test_a_salt_multiplier_changes_only_the_salt(recorded64, replayed64):
    base = replayed64
    more = replay_run("replay-road-x2", "record64",
                      parameters=[Parameter("road", "salt", "chloride:road_salt", 2.0, 0.5, 4.0)])
    t0, t1 = base["table"], more["table"]
    np.testing.assert_allclose(t1["salt:chloride:land:road_salt"], 2.0 * np.asarray(t0["salt:chloride:land:road_salt"]),
                               rtol=1e-12)
    assert np.nansum(t1["salt:chloride:land:road_salt"]) > 0
    assert not np.array_equal(more["x"][-1]["chloride"], base["x"][-1]["chloride"])
    for d in range(DAYS):  # the other species is untouched
        same(more["x"][d]["weathering"], base["x"][d]["weathering"], f"weathering day {d}")
    for k in t0:
        if k.startswith(("salt:weathering:", "transport:weathering:")):
            same(t0[k], t1[k], k)
    for k in ("withdrawal", "return", "consumption"):  # the users' water, rerun with their new intake
        same(t0[f"user:plant:{k}"], t1[f"user:plant:{k}"], k)
    assert not np.array_equal(t0["user:plant:chloride:kg"], t1["user:plant:chloride:kg"])
    for k, v in base["series"].items():  # the flow's series are the record's
        if k.startswith("flow:"):
            same(v, more["series"][k], k)
    assert np.nanmax(np.abs(t1["salt:chloride:closure_share"])) < 1e-9


def test_porosity_and_weathering_may_change(recorded64, replayed64):
    args = dict(salt_config().mf6_builder_args)
    args["transport"] = {sp: dict(v, porosity=0.15) for sp, v in args["transport"].items()}
    args["weathering"] = dict(args["weathering"], weathering=0.004)
    r = replay_run("replay-porosity", "record64", days=60, mf6_builder_args=args)
    assert np.nansum(r["table"]["transport:weathering:boundary:weathering"]) == pytest.approx(
        2.0 * np.nansum(replayed64["table"]["transport:weathering:boundary:weathering"][:60]), rel=1e-12)
    assert not np.array_equal(r["x"][-1]["chloride"], replayed64["x"][59]["chloride"])


@pytest.mark.parametrize("change, match", [
    ({"parameters": [Parameter("cn2", "swatplus", "cn2:pctchg", 5.0, -10.0, 10.0)]}, "parameters"),
    ({"parameters": [Parameter("k", "modflow6", "npf.k:all", 2.0, 0.5, 4.0)]}, "parameters"),
    ({"mf6_builder_args": dict(salt_config().mf6_builder_args, top=89.0)}, "flow model's input files"),
    ({"exchange": {"package": "riv", "map": "all:cha01", "limit": False}}, "config"),
    ({"users": False}, "users"),
    ({"days": DAYS + 1}, "record holds"),
    ({"heads": np.full(100, 88.0)}, "heads"),
])
def test_a_replay_refuses_what_would_change_the_flow(recorded64, change, match):
    from basinledger.driver import Basin
    change = dict(change)
    kw = {k: change.pop(k) for k in ("parameters", "users", "days", "heads") if k in change}
    with warnings.catch_warnings(), pytest.raises(RP.ReplayRefused, match=match):
        warnings.simplefilter("ignore")
        Basin(salt_config(**change), run="replay-refused", replay=RUNS / "record64", **kw)


def test_a_record_without_the_specific_discharge_replays_without_dispersion(swatplus):
    """Specific discharge is dispersion's only use and a third of a record: left out, a replay whose transport models
    disperse is refused, and one without dispersion runs."""
    from basinledger.driver import Basin
    cfg = salt_config(replay_record=RUNS / "record-nodsp", replay_specific_discharge=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        Basin(cfg, run="record-nodsp-run", days=30).run().finalize()
    with warnings.catch_warnings(), pytest.raises(RP.ReplayRefused, match="disperse"):
        warnings.simplefilter("ignore")
        Basin(salt_config(), run="replay-nodsp", replay=RUNS / "record-nodsp")
    args = dict(salt_config().mf6_builder_args)
    args["transport"] = {sp: dict(v, dispersivity=None) for sp, v in args["transport"].items()}
    r = replay_run("replay-nodsp", "record-nodsp", mf6_builder_args=args)
    assert len(r["x"]) == 30 and np.nanmax(np.abs(r["table"]["transport:chloride:closure_share"])) < 1e-6


# ---------- restarts ----------

def test_a_replay_restarts_bit_for_bit(recorded64, replayed64):
    replay_run("replay-warm", "record64", save_at=200)
    b = replay_run("replay-restarted", "record64", restart=RUNS / "replay-warm" / "state")
    for i, d in enumerate(range(200, DAYS)):
        for sp in ("chloride", "weathering"):
            same(replayed64["x"][d][sp], b["x"][i][sp], f"{sp} day {d}")
    for k, v in b["table"].items():
        if k.startswith(("salt:", "transport:")):
            same(replayed64["table"][k], v, k)


def test_a_replay_restarts_from_a_full_runs_state(recorded64):
    """The full run's saved state (its concentrations, salt and ledger) starts a replay: the replay is then the
    unbroken full run, as MODFLOW 6's flow before the restart day comes from the record."""
    full_run("replay-full-warm", DAYS, save_at=150)
    b = replay_run("replay-from-full", "record64", restart=RUNS / "replay-full-warm" / "state")
    for i, d in enumerate(range(150, DAYS)):
        same(recorded64["x"][d]["chloride"], b["x"][i]["chloride"], f"day {d}")


def test_a_restarted_recording_carries_its_record_on(swatplus):
    full_run("record-unbroken-run", 90, "record-unbroken")
    full_run("record-warm-run", 90, "record-pieces", save_at=50)
    full_run("record-rest-run", 90, "record-pieces", restart=RUNS / "record-warm-run" / "state")
    one, two = RP.Record(RUNS / "record-unbroken"), RP.Record(RUNS / "record-pieces")
    assert one.days == two.days == 90
    for d in range(90):
        x, y = one.day(d), two.day(d)
        assert sorted(x) == sorted(y)
        for k in x:
            if k != "json":
                np.testing.assert_array_equal(x[k], y[k], err_msg=f"{k} day {d}")
        assert x["json"]["series"] == y["json"]["series"]


# ---------- the forward model ----------

def test_the_forward_model_replays(recorded64, tmp_path):
    from basinledger.calibration import run_forward_model
    from test_salt import SALT
    obs = tmp_path / "obs.csv"
    obs.write_text("name,kind,object,date,value,weight,source\n"
                   "cl1,chloride_mg_l,cha01,2010-06-01,5.0,1,test\nq1,flow,cha01,2010-06-01,0.1,1,test\n")
    params = tmp_path / "par.csv"
    params.write_text("name,model,target,value,lower,upper,transform\nroad,salt,chloride:road_salt,1.0,0.5,2.0,none\n")
    out = run_forward_model(SALT, params, obs, tmp_path / "sim.csv", run_folder=tmp_path / "runs",
                            replay=RUNS / "record64")
    rows = dict(line.split(",") for line in out.read_text().splitlines()[1:])
    d = 151  # 2010-06-01
    assert float(rows["cl1"]) == recorded64["series"]["chloride_mg_l:cha01"][d]
    assert float(rows["q1"]) == recorded64["series"]["flow:cha01"][d]
