"""Scenarios of water users (scenarios.py, phases.py, attribution.py; docs/design.md, "Scenarios of water users").

The set is examples/synthetic/scenarios/data_centers.toml, on the synthetic town (examples/synthetic/town.toml): a
phased data center on the town's system, a second with its own well, the two together, and the town growing onto a
new well without the plant. The phased user matches the component built at each phase's size to the bit; the users'
own water adds exactly in the combination; every case's ledger closes and the capture's residual is zero; the
runner skips what is done and resumes what was stopped, to the same result.
"""

from __future__ import annotations

import datetime as dt
import json
import math

import numpy as np
import pytest
from conftest import ROOT, RUNS

from basinledger import attribution as A
from basinledger import scenarios as SC
from basinledger.phases import Phase, PhasedUser, combine, parse, scale_schedule
from basinledger.testing import SyntheticUser
from basinledger.users import WaterUser, load_component

SET = ROOT / "examples" / "synthetic" / "scenarios" / "data_centers.toml"


# ---------- phases ----------

def _user(tmp, name, text):
    (tmp / f"{name}.toml").write_text(text)
    return str(tmp / f"{name}.toml")


def test_a_phased_user_is_the_component_at_each_size(tmp_path):
    """Off before the first phase, then exactly the component built with each phase's config, times its
    multiplier."""
    start = dt.date(2020, 1, 1)
    phases = parse("t", [{"from": "2020-01-11", "config": {"withdrawal": 0.002}},
                         {"from": "2020-01-21", "config": {"withdrawal": 0.006}, "multiplier": 0.5},
                         {"from": "2020-01-26", "multiplier": 3.0}])
    p = WaterUser("p", PhasedUser(SyntheticUser, phases, start, tmp_path / "phases"),
                  _user(tmp_path, "u", "withdrawal = 0.004\nconsumed_fraction = 0.6\n"))
    a = WaterUser("a", SyntheticUser(), _user(tmp_path, "a", "withdrawal = 0.002\nconsumed_fraction = 0.6\n"))
    b = WaterUser("b", SyntheticUser(), _user(tmp_path, "b", "withdrawal = 0.006\nconsumed_fraction = 0.6\n"))
    for day in range(30):
        x, ya, yb = p.advance_day(day), a.advance_day(day), b.advance_day(day)
        m, y = (0.0, ya) if day < 10 else (1.0, ya) if day < 20 else (0.5, yb) if day < 25 else (3.0, yb)
        for k in ("withdrawal", "consumption", "returned", "gain", "energy"):
            if m in (0.0, 1.0, 0.5):  # the component itself, or a power of two: to the bit
                assert getattr(x, k) == m * getattr(y, k), (day, k)
            else:  # the hours scaled, then summed: to the sum's rounding
                assert getattr(x, k) == pytest.approx(m * getattr(y, k), rel=1e-14, abs=0), (day, k)
        assert x.constituents == pytest.approx({c: m * v for c, v in y.constituents.items()}, rel=1e-14)
        if m:
            assert x.return_temperature == pytest.approx(y.return_temperature, rel=1e-14)  # not scaled
    # each phase's config was written beside the run, in the original's format
    assert sorted(f.name for f in (tmp_path / "phases").iterdir()) == ["phase1.toml", "phase2.toml"]


def test_a_withdrawal_limit_binds_on_what_the_user_takes(tmp_path):
    phases = parse("t", [{"from": "2020-01-01", "multiplier": 4.0}])
    p = WaterUser("p", PhasedUser(SyntheticUser, phases, dt.date(2020, 1, 1), tmp_path), None)
    p.set_controls(limit=0.003)
    use = p.advance_day(0)
    assert use.withdrawal == pytest.approx(0.003 * 86400, rel=1e-12)  # 4 x 0.006 mean, held to 0.003
    assert use.imbalance == pytest.approx(0.0, abs=1e-9)


def test_cooling_water_in_phases(tmp_path):
    """CoolingWaterBmi at 60 MW, then 300 MW (its it_capacity_mw, set at initialize): each phase is the component
    built at that capacity, under a limit too, and the hours with heat not removed count the active phase's."""
    pytest.importorskip("waterledger_models")
    pytest.importorskip("yaml")
    spec = "waterledger_models.datacenter.bmi:CoolingWaterBmi"
    var = {"withdrawal": "data_center_cooling_water__withdrawal_volume_flux",
           "gain": "data_center_cooling_water__condensate_volume_flux",
           "consumption": "data_center_cooling_water__consumption_volume_flux",
           "return": "data_center_cooling_water__return_volume_flux", "power": "data_center__electric_power"}
    weather = {"atmosphere_bottom_air__temperature": lambda t: 22 + 8 * math.sin(2 * math.pi * t / 86400),
               "atmosphere_bottom_air_water_vapor__relative_saturation": lambda t: 60.0}
    report = {"lost_h": "data_center__setpoint_lost_time"}

    def user(name, comp, mw=100):
        (tmp_path / f"{name}.yaml").write_text(f"design: towers\nit_capacity_mw: {mw}\nend_time: 315360000\n")
        w = WaterUser(name, comp, str(tmp_path / f"{name}.yaml"), variables=var, report=report)
        w.inputs = weather
        return w

    phases = parse("t", [{"from": "2027-01-03", "mw": 60}, {"from": "2027-01-06", "mw": 300}])
    p = user("p", PhasedUser(lambda: load_component(spec), phases, dt.date(2027, 1, 1), tmp_path / "ph"))
    a, b = user("a", load_component(spec), 60), user("b", load_component(spec), 300)
    for day in range(9):
        if day == 7:
            for w in (p, b):
                w.set_controls(limit=0.002)
        x, ya, yb = p.advance_day(day), a.advance_day(day), b.advance_day(day)
        y = None if day < 2 else ya if day < 5 else yb
        assert x.withdrawal == (y.withdrawal if y else 0.0)
        assert x.energy == (y.energy if y else 0.0)
        if day >= 5:
            assert x.reported == yb.reported
    assert x.reported["lost_h"] == 24.0  # the limit binds at 300 MW
    assert json.loads((tmp_path / "ph" / "phase2.yaml").read_text())["it_capacity_mw"] == 300.0


def test_phases_refused():
    for bad, msg in (([{"from": "2027-01-01"}, {"from": "2026-01-01"}], "increase"),
                     ([{"from": "2027-01-01", "mw": 60, "config": {"it_capacity_mw": 60}}], "not both"),
                     ([{"from": "2027-01-01", "multiplier": -1}], "multiplier"),
                     ([{"from": "soon"}], "not a date"),
                     ([{"from": "2027-01-01", "size": 3}], "may have"),
                     ([], "list")):
        with pytest.raises(ValueError, match=msg):
            parse("t", bad)
    with pytest.raises(ValueError, match="multipliers only"):
        scale_schedule("t", phases=[{"from": "2027-01-01", "mw": 3}])


def test_a_scale_on_a_phased_user_multiplies_the_schedules():
    own = parse("t", [{"from": "2027-01-01", "mw": 60}, {"from": "2030-01-01", "mw": 300}])
    c = combine(own, scale_schedule("s", 0.5, dt.date(2028, 6, 1)))
    assert [(p.start, p.overrides, p.multiplier, p.off) for p in c] == [
        (None, {}, 0.0, True), (dt.date(2027, 1, 1), {"it_capacity_mw": 60.0}, 1.0, False),
        (dt.date(2028, 6, 1), {"it_capacity_mw": 60.0}, 0.5, False),
        (dt.date(2030, 1, 1), {"it_capacity_mw": 300.0}, 0.5, False)]
    assert combine(None, scale_schedule("s", 0.0)) == [Phase(None, (), 0.0)]


# ---------- the set and its cases ----------

def test_the_cases_configs():
    s = SC.load_set(SET)
    assert s.cases == ["baseline", "dc_a", "dc_b", "town_grows", "both"]
    configs = {c: SC.case_config(s, c) for c in s.cases}
    for c, (cfg, _) in configs.items():  # every case has every scenario's wells, and the heads compared
        assert [w.name for w in cfg.wells] == ["pw1", "pw2", "wb", "pw3"], c
        assert [h["name"] for h in cfg.outputs["heads"]] == ["pw1", "dc_a", "dc_b"]
    assert [u.name for u in configs["baseline"][0].users] == ["plant"]
    assert [u.name for u in configs["both"][0].users] == ["plant", "dc_a", "dc_b"]
    cfg, scales = configs["town_grows"]
    assert cfg.utilities[0]["sources"] == [{"route": "well:pw2", "share": 0.5}, {"route": "well:pw3", "share": 0.5}]
    assert [p.multiplier for p in scales["town"][0][0]] == [1.0, 1.5]
    assert cfg.users[0].phases == [Phase(None, (), 0.0)]  # the plant, removed
    # a case's changes don't leak into another case's
    b, _ = SC.case_config(s, "baseline")
    assert b.users[0].phases is None and b.utilities[0]["sources"] == [{"route": "well:pw2", "share": 1.0}]


def _set(tmp_path, scenarios: str, combos: str = "") -> SC.ScenarioSet:
    p = tmp_path / "set.toml"
    p.write_text(f'name = "t"\nexample = "{ROOT / "examples" / "synthetic" / "town.toml"}"\n'
                 f'[warm_up]\ndays = 10\n{scenarios}\n{combos}')
    return SC.load_set(p)


def test_changes_refused(tmp_path):
    town = '[[scenario]]\nname = "{n}"\n[[scenario.utilities]]\nname = "town"\nlosses = {x}\n'
    s = _set(tmp_path, town.format(n="a", x=0.2) + town.format(n="b", x=0.3),
             '[[combination]]\nname = "ab"\nscenarios = ["a", "b"]\n')
    SC.case_config(s, "a")
    with pytest.raises(ValueError, match="both change utility:town"):
        SC.case_config(s, "ab")
    s = _set(tmp_path, '[[scenario]]\nname = "a"\nremove = ["utility:town"]\n'
                       '[[scenario.users]]\nname = "dc"\ncomponent = "basinledger.testing:SyntheticUser"\n'
                       'withdrawal = "utility:town"\nreturn = "utility:town"\n')
    with pytest.raises(ValueError, match="users \\['dc'\\] aren't removed"):
        SC.case_config(s, "a")
    s = _set(tmp_path, '[[scenario]]\nname = "a"\n[[scenario.scale]]\nuser = "nobody"\nmultiplier = 2\n')
    with pytest.raises(ValueError, match="isn't there"):
        SC.case_config(s, "a")
    s = _set(tmp_path, '[[scenario]]\nname = "a"\n[[scenario.wells]]\nname = "pw1"\ncell = [0, 1, 1]\n')
    with pytest.raises(ValueError, match="defined differently"):
        SC.case_config(s, "a")
    s = _set(tmp_path, '[[scenario]]\nname = "a"\n[[scenario.users]]\nname = "plant"\n'
                       'component = "basinledger.testing:SyntheticUser"\nwithdrawal = "well:pw1"\nreturn = "outfall:x"\n')
    with pytest.raises(ValueError, match="already there"):
        SC.case_config(s, "a")
    for text, msg in (('[[scenario]]\nname = "baseline"\n', "the example itself"),
                      ('[[scenario]]\nname = "a"\nsize = 3\n', "unknown keys"),
                      ('[[scenario]]\nname = "a"\nremove = ["plant"]\n', "remove lists"),
                      ('[[scenario]]\nname = "a"\n[[combination]]\nname = "c"\nscenarios = ["a"]\n',
                       "two or more")):
        with pytest.raises(ValueError, match=msg):
            _set(tmp_path, text)


# ---------- the statistics ----------

def test_low_flow_statistics():
    dates = [dt.date(2000, 4, 1) + dt.timedelta(days=i) for i in range(366 + 365 * 4)]
    rng = np.random.default_rng(1)
    q = np.exp(rng.normal(0, 0.5, len(dates)))
    r = A.rolling_mean(q, 7)
    assert np.isnan(r[:6]).all() and r[6] == pytest.approx(q[:7].mean())
    years = A.climatic_years(dates)
    assert list(years) == [2000, 2001, 2002, 2003, 2004] and len(years[2000]) == 365 and len(years[2003]) == 366
    lows = A.low_flows(q, dates)
    assert lows["years"] == 5 and lows["method"] == "log-Pearson III"
    assert lows["lowest"] == pytest.approx(np.nanmin(r)) and lows["nQT"] < lows["mean_annual_minimum"]
    # with no skew, log-Pearson III is the log-normal quantile
    x = np.array([1.0, 10.0, 100.0])  # logs 0, 1, 2: skew 0
    assert A.lp3_low_flow(x, 10) == pytest.approx(10 ** (1 - 1.2815515655446004), rel=1e-12)
    # a zero minimum: the empirical quantile; fewer than 3 years: none
    assert A.low_flows(np.r_[np.zeros(400), q[400:]], dates)["method"].startswith("empirical")
    assert math.isnan(A.low_flows(q[:700], dates[:700])["nQT"])
    assert A.flow_duration(np.arange(101.0), [0.1, 0.9]).tolist() == pytest.approx([90.0, 10.0])


def test_interaction_arithmetic():
    i = A.interaction(-3.0, [-1.0, -1.5])
    assert i == {"sum_alone": -2.5, "together": -3.0, "interaction": -0.5, "interaction_share": 0.5 / 3}
    d = A.daily_interaction(np.array([-2.0, -2.0]), [np.array([-1.0, -1.0]), np.array([-1.0, -0.5])])
    assert d == {"interaction_share": 0.5 / 4, "largest": 0.5}


# ---------- the set, run ----------

@pytest.mark.slow
def test_the_set_runs_compares_and_resumes(swatplus, tmp_path):
    text = SET.read_text().replace('example = "../town.toml"',
                                   f'example = "{SET.parent.parent / "town.toml"}"\nfolder = "{RUNS / "user-scenarios"}"')
    for n in ("dc_a", "dc_b"):
        text = text.replace(f'config = "{n}.toml"', f'config = "{SET.parent / (n + ".toml")}"')
    path = tmp_path / "set.toml"
    path.write_text(text)
    log = []
    tables = SC.run_set(path, workers=2, fresh=True, log=log.append)
    s = SC.load_set(path)
    cases = {c: A.read_case(s.folder / "cases" / c) for c in s.cases}
    warm = json.loads((s.folder / "warm_up" / "done.json").read_text())["day"]
    dates = cases["baseline"]["dates"]
    t = {c: r["ledger"] for c, r in cases.items()}

    # the phased user: absent in the warm-up, off until its first phase, then on the town's system
    w = t["dc_a"]["user:dc_a:withdrawal"]
    first = dates.index(dt.date(2010, 11, 15))
    assert np.isnan(w[:warm]).all() and (w[warm:first] == 0).all() and (w[first:] > 0).all()
    k = "utility:town:"
    np.testing.assert_allclose(t["dc_a"][k + "demand"][warm:], t["dc_a"][k + "demand:homes"][warm:] + w[warm:],
                               rtol=1e-12)
    # its second phase is twice its first's config: the same days' withdrawal doubles from March 1
    second = dates.index(dt.date(2011, 3, 1))
    assert w[second] > 1.5 * w[second - 1]
    # dc_b doubles by its multiplier from April 1
    b = t["dc_b"]["user:dc_b:withdrawal"]
    april = dates.index(dt.date(2011, 4, 1))
    assert (b[warm:dates.index(dt.date(2010, 12, 1))] == 0).all() and b[april] > 1.8 * b[april - 1]

    # each case is the baseline up to the end of the warm-up
    for c in s.cases:
        for acct in ("recharge:received", "exchange:cha01:to_stream"):
            assert np.array_equal(t[c][acct][:warm], t["baseline"][acct][:warm]), (c, acct)
    # the users' own water adds exactly together: nothing feeds back on them
    for u in ("dc_a", "dc_b"):
        for acct in [a for a in t[u] if a.startswith(f"user:{u}:")]:
            assert np.array_equal(t["both"][acct], t[u][acct], equal_nan=True), acct
    # the town grows: its homes by half, half of its water from the new well; the plant closes
    tg = t["town_grows"]
    np.testing.assert_allclose(tg[k + "demand:homes"][warm:], 1.5 * t["baseline"][k + "demand:homes"][warm:],
                               rtol=1e-12)
    np.testing.assert_allclose(tg[k + "source:well:pw3:requested"][warm:], 0.5 * tg[k + "asked"][warm:], rtol=1e-12)
    assert (tg["user:plant:withdrawal"][warm:] == 0).all() and tg["user:plant:withdrawal"][warm - 1] > 0

    # closure: every case ran with its ledger checked each day; the capture's residual is zero
    p = slice(warm, None)
    for c in s.cases[1:]:
        cap = A.capture({a: v[p] for a, v in t["baseline"].items()}, {a: v[p] for a, v in t[c].items()})
        assert abs(cap["total"]["residual"]) <= 1e-8 * abs(cap["total"]["pumping"]), c
        assert cap["total"]["pumping"] != 0
    # the interaction is the effect together less the effects alone, from the flow table
    import csv
    flow = {(r["case"], r["statistic"]): float(r["change"]) for r in csv.DictReader(open(tables / "flow.csv"))}
    inter = {(r["kind"], r["point"], r["statistic"]): r for r in csv.DictReader(open(tables / "interaction.csv"))}
    row = inter[("flow", "cha01", "mean")]
    assert float(row["together"]) == pytest.approx(flow[("both", "mean")], rel=1e-12)
    assert float(row["interaction"]) == pytest.approx(
        flow[("both", "mean")] - flow[("dc_a", "mean")] - flow[("dc_b", "mean")], rel=1e-9, abs=1e-15)
    assert float(inter[("use", "dc_a", "consumption")]["interaction"]) == 0.0
    assert flow[("both", "mean")] < flow[("dc_a", "mean")] < 0
    md = (tables / "tables.md").read_text()
    for h in ("## Flow at cha01", "## Heads at dc_b", "## Capture", "## Do the effects add up?", "Screening run"):
        assert h in md

    # restartable: a second run does nothing; a case stopped part way resumes from its checkpoint, to the same days
    log.clear()
    SC.run_set(path, log=log.append)
    assert log[:6] == ["warm-up: done (300 days)"] + [f"{c}: done" for c in s.cases]
    before = (s.folder / "cases" / "dc_b" / "outputs" / "ledger.csv").read_text()
    (s.folder / "cases" / "dc_b" / "done.json").unlink()
    from basinledger.driver import Basin
    cfg, _ = SC.case_config(s, "dc_b")
    cfg.run_folder = tmp_path
    stopped = Basin(cfg, run="stopped", restart=s.state)
    stopped.run(100)
    ckpt = s.folder / "checkpoints" / "dc_b"
    stopped.save(ckpt)
    stopped.finalize()
    (ckpt / "checkpoint.json").write_text(json.dumps({"fingerprint": SC.case_fingerprint(s, "dc_b"), "day": 400}))
    log.clear()
    SC.run_set(path, cases=["dc_b"], log=log.append)
    assert "dc_b: finished (resumed from day 400)" in log
    assert (s.folder / "cases" / "dc_b" / "outputs" / "ledger.csv").read_text() == before
    assert not ckpt.exists()


# ---------- operating rules for systems and reservoirs ----------

def test_operating_rules_by_construction():
    from basinledger import rules as R
    st = R.State(day=3, flows={"river": 5.0, "town": 0.4}, head_at=None, volumes={"res": 6e5}, pools={"res": 1e6},
                 capacities={"fill": 2.0, "back": 1.0})
    # pumped storage: the least of the capacity, the flow above the threshold and the room over a day
    assert R.pumped_storage(st, transfer="fill", channel="river", above=4.0, reservoir="res") == \
        {"transfer:fill": {"rate": 1.0}}
    assert R.pumped_storage(st, transfer="fill", channel="river", above=1.0, reservoir="res")["transfer:fill"][
        "rate"] == 2.0
    assert R.pumped_storage(st, transfer="fill", channel="river", above=1.0, reservoir="res", below=0.6864)[
        "transfer:fill"]["rate"] == pytest.approx(1.0)  # 86,400 m3 of room
    assert R.pumped_storage(st, transfer="fill", channel="river", above=6.0, reservoir="res") == \
        {"transfer:fill": {"rate": 0.0}}  # the river below the threshold
    assert R.pumped_storage(st, transfer="fill", channel="river", above=1.0, reservoir="res", pool=5e5) == \
        {"transfer:fill": {"rate": 0.0}}  # the reservoir full
    # a release: the channel's deficit, up to the capacity, keeping a share of the pool
    assert R.release(st, transfer="back", channel="town", below=1.0) == {"transfer:back": {"rate": 0.6}}
    assert R.release(st, transfer="back", channel="town", below=0.3) == {"transfer:back": {"rate": 0.0}}
    assert R.release(st, transfer="back", channel="town", below=1.0, reservoir="res", keep=0.6) == \
        {"transfer:back": {"rate": 0.0}}
    # restrictions: the last stage whose trigger holds
    stages = [{"when": {"channel": "town", "below": 0.5}, "cut": {"outdoor": 0.5}},
              {"when": {"reservoir": "res", "below": 0.5}, "cut": {"outdoor": 1.0, "indoor": 0.1}}]
    assert R.restriction(st, utility="u", stages=stages) == {"utility:u": {"cut:outdoor": 0.5}}
    assert R.restriction(R.State(**{**st.__dict__, "volumes": {"res": 4e5}}), utility="u", stages=stages) == \
        {"utility:u": {"cut:outdoor": 1.0, "cut:indoor": 0.1}}
    assert R.restriction(R.State(**{**st.__dict__, "day": 0}), utility="u", stages=stages) == {}
    # the smallest rate, the largest cut
    assert R.combine([{"transfer:t": {"rate": 2.0}}, {"transfer:t": {"rate": 0.5}},
                      {"utility:u": {"cut:a": 0.1}}, {"utility:u": {"cut:a": 0.4}}]) == \
        {"transfer:t": {"rate": 0.5}, "utility:u": {"cut:a": 0.4}}
    for bad in ([], [{"when": {"channel": "x"}, "cut": {"a": 0.5}}], [{"when": {"channel": "x", "below": 1},
                                                                       "cut": {"a": 1.5}}]):
        with pytest.raises(ValueError):
            R.check("restriction", {"utility": "u", "stages": bad}, "t")


def test_transfers_and_point_sources_refused(tmp_path):
    from basinledger import config as C
    p = tmp_path / "x.toml"
    for x, msg in (({"name": "t", "from": "well:pw1", "to": "swat:point_source:a", "capacity": 1}, "SWAT\\+ intake"),
                   ({"name": "t", "from": "swat:place_of_use:a", "to": "ground:g", "capacity": 1}, "point source"),
                   ({"name": "t", "from": "swat:place_of_use:a", "to": "outfall:x", "capacity": -1}, "capacity"),
                   ({"name": "t", "from": "swat:place_of_use:a", "to": "outfall:x"}, "each has")):
        with pytest.raises(ValueError, match=msg):
            C._transfer(p, x)
    assert C._point_source(p, {"name": "in", "reservoir": "res01"}) == ("in", "res:res01")
    with pytest.raises(ValueError, match="one of channel or reservoir"):
        C._point_source(p, {"name": "in", "reservoir": "res01", "channel": "cha01"})


OPS_DAYS, OPS_SAVE = 240, 120


def _ops_config(tmp, swatplus):
    """Osu_1hru with res01 below its channel: a pump from the channel into the reservoir at high flow, a release
    from the reservoir back to the channel at low flow, and a town on a well whose outdoor use is cut by half when
    the channel is low and wholly (with a tenth of indoor use) when the reservoir is low."""
    import importlib.util
    from basinledger import config as C
    spec = importlib.util.spec_from_file_location("existing", ROOT / "examples" / "synthetic" / "existing.py")
    ex = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ex)
    cfg = C.load(ROOT / "examples" / "synthetic" / "surface.toml")
    cfg.run_folder = RUNS
    cfg.swat_project = ex.reservoir_project(tmp / "osu_res", swatplus)
    cfg.places_of_use = [("pump", "cha01"), ("out", "res:res01")]
    cfg.point_sources = [("fill", "res:res01"), ("back", "cha01")]
    cfg.users = []
    cfg.utilities = [{"name": "town", "sources": [{"route": "well:pw1", "share": 1.0}], "outfall": "outfall:town",
                      "customers": [{"name": "outdoor", "demand": 0.0004, "consumed": 1.0},
                                    {"name": "indoor", "demand": 0.0006, "sewered": 1.0}]}]
    cfg.transfers = [{"name": "fill", "from": "swat:place_of_use:pump", "to": "swat:point_source:fill",
                      "capacity": 0.002},
                     {"name": "release", "from": "swat:place_of_use:out", "to": "swat:point_source:back",
                      "capacity": 0.001}]
    cfg.rules = [{"name": "fill", "kind": "pumped_storage", "transfer": "fill", "channel": "cha01", "above": 0.001,
                  "reservoir": "res01", "below": 1.2},
                 {"name": "release", "kind": "release", "transfer": "release", "channel": "cha01", "below": 0.01,
                  "reservoir": "res01", "keep": 0.3},
                 {"name": "drought", "kind": "restriction", "utility": "town",
                  "stages": [{"when": {"channel": "cha01", "below": 0.003}, "cut": {"outdoor": 0.5}},
                             {"when": {"reservoir": "res01", "below": 0.5}, "cut": {"outdoor": 1.0, "indoor": 0.1}}]}]
    return cfg


@pytest.mark.slow
def test_pumped_storage_a_release_and_restrictions(swatplus, tmp_path):
    from basinledger.coupling import lakes as LK
    from basinledger.coupling import swatplus as S
    from basinledger.driver import Basin
    cfg = _ops_config(tmp_path, swatplus)
    b = Basin(cfg, run="ops", days=OPS_DAYS)
    assert b.pools == {"res01": 1e4}  # 1 ha-m, its principal spillway volume
    q, inflow, vol = [], [], []
    for d in range(OPS_DAYS):
        b.step()
        q.append(float(b.sw.get(S.QOUT)[0]))
        inflow.append(float(b.sw.get(LK.INFLOW)[0]) * 86400)
        vol.append(float(b.sw.get(LK.VOLUME)[0]))
        if d == OPS_SAVE - 1:
            b.save(tmp_path / "state")
    t = b.ledger.table()
    b.finalize()
    q, inflow, vol = np.array(q), np.array(inflow), np.array(vol)
    f, r = "transfer:fill:", "transfer:release:"

    # the pump: only on days after the river was above the threshold and the reservoir below its target, at the
    # least of capacity, flow above the threshold and room; what it got goes into the reservoir the next day
    rate = t[f + "rate"]
    assert (rate > 0).sum() > 20 and (rate == 0).sum() > 20
    expect = np.r_[0.0, np.clip(np.minimum(np.minimum(0.002, q[:-1] - 0.001), (1.2e4 - vol[:-1]) / 86400), 0, None)]
    np.testing.assert_allclose(rate, expect, rtol=1e-6, atol=1e-12)
    for k in (f, r):
        assert np.array_equal(t[k + "sent"][1:], t[k + "delivered"][:-1])
        assert t[k + "sent"][0] == 0
    # the reservoir takes it with the channel's water: its inflow is the channel's outflow and the pump's
    np.testing.assert_allclose(inflow, q * 86400 + t[f + "sent"], rtol=1e-5, atol=1e-3)
    # the release: only on days after low flow, never taking the reservoir below its keep
    rel = t[r + "rate"]
    assert (rel > 0).sum() > 5
    assert (rel[1:][q[:-1] >= 0.01] == 0).all()
    assert (rel[1:] * 86400 <= np.clip(vol[:-1] - 0.3e4, 0, None) * (1 + 1e-6)).all()
    assert t["route:swat:point_source:back:delivered"].sum() == pytest.approx(t[r + "sent"].sum(), rel=1e-6)

    # restrictions: booked apart from a shortfall, staged
    k = "utility:town:"
    cut_o = t[k + "restricted:outdoor"] / t[k + "wanted:outdoor"]
    cut_i = t[k + "restricted:indoor"] / t[k + "wanted:indoor"]
    stage2 = np.r_[False, vol[:-1] < 0.5e4]
    stage1 = np.r_[False, q[:-1] < 0.003] & ~stage2
    assert stage1.sum() > 5 and stage2.sum() > 5, (stage1.sum(), stage2.sum())
    np.testing.assert_allclose(cut_o, np.where(stage2, 1.0, np.where(stage1, 0.5, 0.0)), atol=1e-12)
    np.testing.assert_allclose(cut_i, np.where(stage2, 0.1, 0.0), atol=1e-12)
    np.testing.assert_allclose(t[k + "demand"], t[k + "wanted:outdoor"] + t[k + "wanted:indoor"] - t[k + "restricted"],
                               rtol=1e-12)
    assert np.allclose(t[k + "unmet"], 0.0, atol=1e-6)  # the well delivers all: nothing is a shortfall
    assert (t["rule:drought:active"] == (stage1 | stage2)).all()

    # a restart carries the water held overnight: the same days, bit for bit
    c = Basin(_ops_config(tmp_path / "again", swatplus), run="ops-restart", days=OPS_DAYS, restart=tmp_path / "state")
    c.run()
    tc = c.ledger.table()
    c.finalize()
    for a in [a for a in t if a.startswith(("transfer:", "utility:town:", "route:"))]:
        assert np.array_equal(t[a], tc[a], equal_nan=True), a
