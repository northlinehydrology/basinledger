"""Water systems (chain.py) balance by construction: production, delivery and losses; sources by share or priority;
customers by class; returns the next day; a shortfall scaling the day; reuse of another system's effluent."""

from __future__ import annotations

import datetime as dt

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
