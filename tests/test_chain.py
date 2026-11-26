"""Water systems (chain.py) balance by construction: production, delivery and losses; sources by share or priority;
customers by class; returns the next day; a shortfall scaling the day; reuse of another system's effluent."""

from __future__ import annotations

import datetime as dt

import numpy as np
import pytest

from basinledger.chain import Source, UtilityChain
from basinledger.ledger import FLOAT32, Ledger

DAY = 86400.0
JAN, JUL = dt.date(2010, 1, 15), dt.date(2010, 7, 15)


def city(**kw):
    c = dict(name="city", sources=[Source("well:a", 0.6), Source("swat:place_of_use:intake", 0.4)],
             outfall="swat:point_source:wwtp", losses=0.1, demand=0.002, sewered=0.8, infiltration=0.0001)
    c.update(kw)
    return UtilityChain(**c)


def town(**kw):
    c = dict(name="town", sources=[Source("swat:place_of_use:intake", None, 1, 0.003), Source("well:a", None, 2)],
             outfall="swat:point_source:wwtp", losses=0.1, leaks=0.5, leaks_to="ground:mains")
    c.update(kw)
    t = UtilityChain.from_config({
        "name": c["name"], "outfall": c["outfall"], "losses": c["losses"],
        "sources": [{"route": "swat:place_of_use:intake", "priority": 1, "capacity": 0.003},
                    {"route": "well:a", "priority": 2}],
        "leaks": {"share": c["leaks"], "to": c["leaks_to"]},
        "customers": [{"name": "homes", "demand": [0.002] * 6 + [0.004] * 6, "sewered": 0.7, "consumed": 0.3},
                      {"name": "edge", "demand": 0.001, "septic": 0.9, "septic_to": "ground:edge",
                       "consumed": 0.1}]})
    return t


def day(sys_, L, date, delivered: dict | None = None, systems=None):
    """One day of systems alone: plan, release, the routes delivering what they're asked unless ``delivered`` says
    otherwise (a share of it, by route), share, settle. Returns the requests."""
    systems = systems or {sys_.name: sys_}
    reqs = {}
    for s in systems.values():
        reqs[s.name] = s.plan(L, date, systems)
    for s in systems.values():
        for r, v in s.release(L).items():
            reqs[s.name][r] = reqs[s.name].get(r, 0.0) + v
    for s in systems.values():
        for r, v in s._day["requests"].items():
            if r.startswith("effluent:"):
                continue
            L.enter(f"route:{r}:requested", v)
            L.enter(f"route:{r}:delivered", v * (delivered or {}).get(r, 1.0))
    for s in systems.values():
        s.share(L, systems)
    for s in systems.values():
        s.settle(L)
    L.close_day()
    return reqs[sys_.name]


def test_a_day_on_the_chain():
    L = Ledger()
    u = city()
    u.declare(L)
    u.add_user(500.0, 0.0)   # m3 withdrawn from the tap
    u.add_user(0.0, 150.0)   # m3 returned to the sewer
    req = u.close_day(L)
    L.close_day()            # the chain's pairs hold
    existing = 0.002 * DAY
    produced = (existing + 500.0) / 0.9
    assert req["well:a"] == pytest.approx(0.6 * produced)
    assert req["swat:place_of_use:intake"] == pytest.approx(0.4 * produced)
    assert req["swat:point_source:wwtp"] == 0.0  # the plant sends yesterday's sewer: none on the first day
    sewer = existing * 0.8 + 150.0 + 0.0001 * DAY
    assert u.held["sewer"] == pytest.approx(sewer)
    assert L.totals()["utility:city:losses"] == pytest.approx(produced * 0.1)
    req = u.close_day(L)
    L.close_day()
    assert req["swat:point_source:wwtp"] == pytest.approx(sewer)


def test_users_start_fresh_each_day():
    L = Ledger()
    u = city(losses=0.0, demand=0.0, infiltration=0.0)
    u.add_user(100.0, 40.0)
    u.close_day(L)
    req = u.close_day(L)
    assert req["well:a"] == 0.0 and req["swat:point_source:wwtp"] == pytest.approx(40.0)
    assert u.close_day(L)["swat:point_source:wwtp"] == 0.0


def test_customers_by_month_sources_in_priority_and_returns_the_next_day():
    L = Ledger()
    t = town()
    t.declare(L, rtol=FLOAT32)
    req = day(t, L, JAN)
    demand = (0.002 + 0.001) * DAY
    asked = demand / 0.9
    assert req["swat:place_of_use:intake"] == pytest.approx(min(asked, 0.003 * DAY))
    assert req["well:a"] == pytest.approx(asked - 0.003 * DAY)
    assert "ground:edge" not in req or req["ground:edge"] == 0.0  # nothing held yet
    homes, edge = 0.002 * DAY, 0.001 * DAY
    assert t.held["sewer"] == pytest.approx(0.7 * homes)
    leaks = 0.5 * (asked - demand)
    assert t.held["ground"] == pytest.approx({"ground:edge": 0.9 * edge, "ground:mains": leaks})
    req = day(t, L, JUL)  # July's homes demand doubles; yesterday's returns leave today
    assert req["swat:point_source:wwtp"] == pytest.approx(0.7 * homes)
    assert req["ground:edge"] == pytest.approx(0.9 * edge) and req["ground:mains"] == pytest.approx(leaks)
    assert t.held["sewer"] == pytest.approx(0.7 * 0.004 * DAY)
    tot = L.totals()
    assert tot["utility:town:consumed"] == pytest.approx(0.3 * (0.002 + 0.004) * DAY + 0.1 * 2 * edge)


def test_a_short_source_scales_the_day():
    L = Ledger()
    t = town()
    t.declare(L, rtol=FLOAT32)
    day(t, L, JAN, delivered={"swat:place_of_use:intake": 0.25})  # the intake gives a quarter of its 0.003 m3/s
    demand = 0.003 * DAY
    asked = demand / 0.9
    produced = 0.25 * 0.003 * DAY + (asked - 0.003 * DAY)
    s = produced / asked
    tot = L.totals()
    assert tot["utility:town:share"] == pytest.approx(s)
    assert tot["utility:town:delivered"] == pytest.approx(s * demand)
    assert tot["utility:town:unmet"] == pytest.approx((1 - s) * demand)
    assert t.held["sewer"] == pytest.approx(0.7 * s * 0.002 * DAY)
    assert t.held["ground"]["ground:edge"] == pytest.approx(0.9 * s * 0.001 * DAY)


def test_daily_demand_from_a_file_and_monthly_shares(tmp_path):
    (tmp_path / "homes.csv").write_text("date,m3\n2010-01-15,300\n2010-07-15,700\n")
    t = UtilityChain.from_config({
        "name": "t", "outfall": "outfall:plant", "losses": 0.0,
        "sources": [{"route": "well:a", "share": [1.0] * 6 + [0.5] * 6},
                    {"route": "swat:place_of_use:river", "share": [0.0] * 6 + [0.5] * 6}],
        "customers": [{"name": "homes", "demand": "homes.csv", "sewered": 1.0}]}, tmp_path)
    L = Ledger()
    t.declare(L, rtol=FLOAT32)
    assert day(t, L, JAN) == {"well:a": pytest.approx(300.0), "swat:place_of_use:river": 0.0, "outfall:plant": 0.0}
    req = day(t, L, JUL)
    assert req["well:a"] == pytest.approx(350.0) and req["swat:place_of_use:river"] == pytest.approx(350.0)
    assert req["outfall:plant"] == pytest.approx(300.0)


def test_reuse_takes_another_plants_effluent_before_its_outfall():
    L = Ledger()
    t = town()
    dc = UtilityChain.from_config({"name": "dc", "outfall": "outfall:dc", "losses": 0.0,
                                   "sources": [{"route": "effluent:town", "share": 1.0}],
                                   "customers": [{"name": "cooling", "demand": 0.0005, "consumed": 1.0}]})
    systems = {"town": t, "dc": dc}
    for s in systems.values():
        s.declare(L, rtol=FLOAT32)
    day(t, L, JAN, systems=systems)          # nothing held yet: the data center gets none
    assert L.totals()["utility:dc:share"] == 0.0
    sewer = t.held["sewer"]
    req = day(t, L, JAN, systems=systems)    # yesterday's town sewer: the data center takes its 43.2 m3 first
    want = 0.0005 * DAY
    assert want < sewer
    assert req["swat:point_source:wwtp"] == pytest.approx(sewer - want)
    tot = L.totals()
    assert tot["utility:town:reuse:dc"] == pytest.approx(want)
    assert tot["utility:dc:source:effluent:town:delivered"] == pytest.approx(want)


@pytest.mark.parametrize("bad", [dict(sources=[Source("well:a", 0.5)]), dict(losses=1.0), dict(sewered=1.2),
                                 dict(sources=[Source("utility:other", 1.0)]),
                                 dict(sources=[Source("well:a", 0.5), Source("well:b", None, 1)]),
                                 dict(sources=[Source("effluent:city", 1.0)]),
                                 dict(leaks=0.5)])
def test_bad_systems_are_refused(bad):
    with pytest.raises(ValueError):
        city(**bad)


def test_bad_customers_are_refused():
    base = {"name": "t", "outfall": "outfall:p", "sources": [{"route": "well:a", "share": 1.0}]}
    for c in ({"name": "x", "demand": 1.0, "septic": 0.5, "sewered": 0.5},             # septic with nowhere to go
              {"name": "x", "demand": 1.0, "sewered": 0.7, "consumed": 0.4},           # shares beyond 1
              {"name": "x", "demand": [1.0] * 11},                                       # 11 months
              {"name": "x", "demand": 1.0, "outdoor": 0.2}):                             # an unknown key
        with pytest.raises(ValueError):
            UtilityChain.from_config(dict(base, customers=[c]))


def test_dated_shares_held_until_the_next_date(tmp_path):
    """Culpeper's wells start in 2015, before which all of it is Lake Pelham's: a table of dated shares, each row
    held until the next date."""
    (tmp_path / "shares.csv").write_text("date,swat:place_of_use:lake,well:a\n"
                                         "2010-01-01,1.0,0.0\n2015-03-01,0.6,0.4\n2016-01-01,0.5,0.5\n")
    t = UtilityChain.from_config({
        "name": "t", "outfall": "outfall:plant", "shares": "shares.csv",
        "sources": [{"route": "swat:place_of_use:lake"}, {"route": "well:a"}],
        "customers": [{"name": "homes", "demand": 0.001, "sewered": 1.0}]}, tmp_path)
    L = Ledger()
    t.declare(L, rtol=FLOAT32)
    asked = 0.001 * DAY
    for date, lake in ((dt.date(2014, 12, 31), 1.0), (dt.date(2015, 2, 28), 1.0), (dt.date(2015, 3, 1), 0.6),
                       (dt.date(2015, 12, 31), 0.6), (dt.date(2020, 6, 1), 0.5)):
        req = day(t, L, date)
        assert req["swat:place_of_use:lake"] == pytest.approx(lake * asked)
        assert req["well:a"] == pytest.approx((1 - lake) * asked, abs=1e-12)
    with pytest.raises(ValueError, match="no source shares on or before 2009-12-31"):
        day(t, L, dt.date(2009, 12, 31))


def test_bad_dated_shares_are_refused(tmp_path):
    base = {"name": "t", "outfall": "outfall:p", "shares": "s.csv",
            "sources": [{"route": "well:a"}, {"route": "well:b"}]}
    for text, sources, match in (
            ("date,well:a,well:b\n2010-01-01,0.5,0.4\n", None, "sum to 0.9"),
            ("date,well:a\n2010-01-01,1.0\n", None, "aren't its sources"),
            ("date,well:a,well:b\n2010-01-01,0.5,0.5\n", [{"route": "well:a", "share": 0.5},
                                                          {"route": "well:b", "share": 0.5}], "neither a share"),
            ("date,well:a,well:b\n", None, "has no shares")):
        (tmp_path / "s.csv").write_text(text)
        with pytest.raises(ValueError, match=match):
            UtilityChain.from_config(dict(base, sources=sources or base["sources"]), tmp_path)


# ---------- mass through a system (docs/design.md, Mass through the human-system routes, part 2) ----------

SP = ["chloride", "nitrate"]
C_INTAKE = {"chloride": 0.020, "nitrate": 0.001}    # kg m-3 of each source's water that day [A]
C_WELL = {"chloride": 0.010, "nitrate": 0.004}


def quality_town(tmp_path, **kw):
    (tmp_path / "dmr.csv").write_text("date,nitrate\n2010-01-01,0.006\n2010-01-03,0.008\n")
    c = {"name": "town", "outfall": "swat:point_source:wwtp", "losses": 0.1,
         "sources": [{"route": "swat:place_of_use:intake", "priority": 1, "capacity": 0.003},
                     {"route": "well:a", "priority": 2}],
         "leaks": {"share": 0.5, "to": "ground:mains"},
         "customers": [{"name": "homes", "demand": 0.004, "sewered": 0.7, "consumed": 0.3},
                       {"name": "edge", "demand": 0.001, "septic": 0.9, "septic_to": "ground:edge", "consumed": 0.1,
                        "septic_effluent": {"nitrate": 0.03}}],
         "effluent": {"chloride": "sewer", "nitrate": "dmr.csv"},
         "septic_effluent": {"chloride": 0.08, "nitrate": 0.05}}
    c.update(kw)
    return UtilityChain.from_config(c, tmp_path)


def mass_day(systems, L, date, delivered=None, users=(), scale=None):
    """``day`` with each system's mass: released after the water, settled after it. ``users``: (system, m3 from the
    tap, m3 returned to its sewer, its reported kg)."""
    for name, w, _, _ in users:
        systems[name].add_user(w, 0.0)
    for s in systems.values():
        s.plan(L, date, systems)
    out = {}
    for s in systems.values():
        s.release(L)
    for s in systems.values():
        out[s.name] = s.mass_release(L, date, SP, scale)
    for s in systems.values():
        for r, v in s._day["requests"].items():
            if not r.startswith("effluent:"):
                L.enter(f"route:{r}:requested", v)
                L.enter(f"route:{r}:delivered", v * (delivered or {}).get(r, 1.0))
    for s in systems.values():
        s.share(L, systems)
    for name, _, back, kg in users:
        share = systems[name]._day["share"]
        systems[name].add_return(share * back, {c: share * m for c, m in kg.items()})
    for s in systems.values():
        s.settle(L)
    for s in systems.values():
        cs = {}
        for src in s.sources:
            base = C_INTAKE if src.route.endswith("intake") else C_WELL
            cs[src.route] = {sp: (src.concentration.get(sp) if sp in src.concentration else
                                  systems[src.effluent].effluent_concentration(sp) if src.effluent else base[sp])
                             for sp in SP}
        s.mass_settle(L, date, SP, cs, scale)
    L.close_day()
    return out


def test_a_systems_mass_mixes_its_sources_and_closes_every_day(tmp_path):
    L = Ledger()
    t = quality_town(tmp_path)
    t.declare(L, rtol=FLOAT32)
    t.declare_mass(L, SP)
    systems = {"town": t}
    dates = [dt.date(2010, 1, 1) + dt.timedelta(days=i) for i in range(5)]
    released = [mass_day(systems, L, d, delivered={"swat:place_of_use:intake": 0.8},
                         users=[("town", 50.0, 20.0, {"chloride": 2.0})]) for d in dates]
    T = L.table()
    k = "utility:town:"
    for sp in SP:
        ks = k + f"{sp}:"
        # mixed by what each source delivered
        mix = (T[k + "source:swat:place_of_use:intake:delivered"] * C_INTAKE[sp]
               + T[k + "source:well:a:delivered"] * C_WELL[sp]) / T[k + "produced"]
        np.testing.assert_allclose(T[ks + "concentration"], mix, rtol=1e-12)
        # leaks at the mix; septic at its stated concentration (the class's own for nitrate)
        np.testing.assert_allclose(T[ks + "leaks:kg"], mix * T[k + "leaks"], rtol=1e-12)
        sep = {"chloride": 0.08, "nitrate": 0.03}[sp]
        np.testing.assert_allclose(T[ks + "septic:edge:kg"], sep * T[k + "septic:edge"], rtol=1e-12)
        # what goes to the ground today is what was held yesterday
        np.testing.assert_allclose(T[ks + "to_ground:ground:edge:kg"][1:], T[ks + "septic:edge:kg"][:-1], rtol=1e-12)
        np.testing.assert_allclose(T[ks + "to_ground:ground:mains:kg"][1:], T[ks + "leaks:kg"][:-1], rtol=1e-12)
        assert np.nanmax(np.abs(T[ks + "closure_share"])) < 1e-12
    # the sewer: the tap water's in the homes' sewered water, and the user's return (its chloride reported, its
    # nitrate at the tap's)
    c = "utility:town:chloride:"
    np.testing.assert_allclose(T[c + "users_return:kg"], 2.0 * T[k + "share"], rtol=1e-12)
    n = "utility:town:nitrate:"
    np.testing.assert_allclose(T[n + "users_return:kg"], T[n + "concentration"] * 20.0 * T[k + "share"], rtol=1e-12)
    # chloride passes the plant ("sewer"); nitrate leaves at the DMR's, held from its date, the day it leaves
    np.testing.assert_allclose(T[c + "plant_outflow:kg"][1:], T[c + "sewer:kg"][:-1], rtol=1e-12)
    np.testing.assert_allclose(T[c + "plant:added"], 0.0, atol=1e-12)
    np.testing.assert_allclose(T[n + "plant_outflow:kg"], [0.006, 0.006, 0.008, 0.008, 0.008] * T[k + "plant_outflow"],
                               rtol=1e-12)
    assert released[3]["town"]["outfall"]["nitrate"] == pytest.approx(0.008 * T[k + "outfall"][3])


def test_reuse_takes_the_plants_effluent_concentration_and_salts_multipliers_scale_stated_ones(tmp_path):
    L = Ledger()
    t = quality_town(tmp_path, effluent={"chloride": 0.12, "nitrate": 0.006})
    dc = UtilityChain.from_config({"name": "dc", "outfall": "outfall:dc", "losses": 0.0,
                                   "sources": [{"route": "effluent:town", "share": 1.0}],
                                   "customers": [{"name": "cooling", "demand": 0.0005, "consumed": 1.0}]})
    systems = {"town": t, "dc": dc}
    for s in systems.values():
        s.declare(L, rtol=FLOAT32)
        s.declare_mass(L, SP)
    scale = {"chloride": {"effluent": 2.0, "septic": 0.5}}
    for i in range(3):
        mass_day(systems, L, dt.date(2010, 1, 1 + i), scale=scale)
    T = L.table()
    np.testing.assert_allclose(T["utility:dc:chloride:concentration"][1:], 0.24, rtol=1e-12)
    np.testing.assert_allclose(T["utility:dc:nitrate:concentration"][1:], 0.006, rtol=1e-12)
    np.testing.assert_allclose(T["utility:town:chloride:reuse:dc:kg"], 0.24 * T["utility:town:reuse:dc"], rtol=1e-12)
    np.testing.assert_allclose(T["utility:town:chloride:septic:edge:kg"], 0.04 * T["utility:town:septic:edge"],
                               rtol=1e-12)
    np.testing.assert_allclose(T["utility:town:nitrate:septic:edge:kg"], 0.03 * T["utility:town:septic:edge"],
                               rtol=1e-12)
    for s in systems:
        for sp in SP:
            assert np.nanmax(np.abs(T[f"utility:{s}:{sp}:closure_share"])) < 1e-12


def test_a_restart_carries_the_held_mass(tmp_path):
    import json
    L1, L2 = Ledger(), Ledger()
    a, b = quality_town(tmp_path), quality_town(tmp_path)
    for L, s in ((L1, a), (L2, b)):
        s.declare(L, rtol=FLOAT32)
        s.declare_mass(L, SP)
    for i in range(2):
        for L, s in ((L1, a), (L2, b)):
            mass_day({"town": s}, L, dt.date(2010, 1, 1 + i))
    b2 = quality_town(tmp_path)
    b2.held = json.loads(json.dumps(b.held))
    b2.declare(L2, rtol=FLOAT32)
    mass_day({"town": a}, L1, dt.date(2010, 1, 3))
    mass_day({"town": b2}, L2, dt.date(2010, 1, 3))
    for x in ("chloride:plant_outflow:kg", "nitrate:to_ground:ground:edge:kg", "chloride:closure_share"):
        assert L1.table()[f"utility:town:{x}"][-1] == L2.table()[f"utility:town:{x}"][-1]


def test_what_a_systems_mass_needs_stated(tmp_path):
    t = quality_town(tmp_path, effluent={}, septic_effluent={},
                     customers=[{"name": "homes", "demand": 0.004, "sewered": 0.7, "consumed": 0.3},
                                {"name": "edge", "demand": 0.001, "septic": 0.9, "septic_to": "ground:edge",
                                 "consumed": 0.1, "septic_effluent": {"nitrate": "tap"}}])
    miss = t.missing(SP, False, {"well:a": set(SP), "swat:place_of_use:intake": {"chloride"}})
    assert miss == ["effluent chloride", "septic_effluent chloride (class edge)",
                    "effluent nitrate", "source swat:place_of_use:intake concentration nitrate"]
    t = quality_town(tmp_path, effluent={"chloride": 0.1, "nitrate": 0.1, "phosphorus": 0.1})
    assert t.missing(SP, False, {"well:a": set(SP), "swat:place_of_use:intake": set(SP)}) == [
        "['phosphorus'], which aren't transported"]
    with pytest.raises(ValueError, match="no column 'phosphorus'"):
        quality_town(tmp_path, effluent={"phosphorus": "dmr.csv"})
    with pytest.raises(ValueError, match="no concentration on or before 2009-12-31"):
        L = Ledger()
        t = quality_town(tmp_path)
        t.held["sewer"] = 1.0
        t.plan(L, dt.date(2009, 12, 31), {"town": t})
        t.release(L)
        t.mass_release(L, dt.date(2009, 12, 31), SP)
