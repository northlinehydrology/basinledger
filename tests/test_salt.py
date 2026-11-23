"""Salt (salt.py): chloride through the land, the aquifer and the streams of a coupled run.

Units (no models):
- the land's store closes every day; at steady state with uniform inputs its outflow carries the precipitation's
  concentration, raised by evapotranspiration as the chloride mass balance has it; the surface store holds snowfall's
  chloride until melt;
- the water bodies close every day; uniform inflow concentration reaches every body; a pulse passes downstream whole;
- the recharge's mass is placed on the cells exactly as the water is;
- land sources' timing (uniform, by month, the first fertilizer day of each year); an annual road salt amount spread
  over each calendar year's snow and ice days from the weather records, whole; tables by calendar year;
- SWAT+'s routing read from a project's files; the builder's weathering and stream-loss sources.

On the models (examples/synthetic/salt.toml; SWAT+ and libmf6):
- every account of both species (chloride and weathering's ions) closes every day, the land's, the streams' and the
  whole system's, with the links sent = received; conductance from both species by hand; each year's road salt and
  potash applied whole; multipliers scale an input group's mass;
- the recharge's mass is the land's, and GWT takes it whole; the gains reach the stream the next day; lakes at a
  fixed stage carry it through their outlet channel;
- a pulse of road salt arrives in the stream, passes through, and the rest stays in the land and the aquifer, the
  difference from a run without it adding to the pulse every day;
- a restart carries the salt bit for bit.
"""

from __future__ import annotations

import datetime as dt
import warnings

import numpy as np
import pytest

from basinledger import salt as SA
from basinledger.mapping import AreaMap
from conftest import ROOT, RUNS

SALT = ROOT / "examples" / "synthetic" / "salt.toml"


# ---------- the land ----------

def water(n=1, **kw):
    w = {k: np.zeros(n) for k in SA.HRU_READS}
    w["tmn"] += 5.0
    for k, v in kw.items():
        w[k] = np.full(n, float(v))
    return w


def run_land(days, w, perc_mm, store_mm, cp, surface=False, m0=0.0):
    area = np.array([1e4])
    m, ms = np.array([m0]), np.zeros(1)
    out = []
    for _ in range(days):
        dep = cp * w["precip"] / 1000.0 * area
        r = SA.land_day(m, ms, area, w, np.array([perc_mm]), np.array([store_mm]), dep, np.zeros(1), np.zeros(1),
                        surface)
        gone = sum(v.sum() for v in r["out"].values())
        assert r["m_soil"].sum() + r["m_surface"].sum() - m.sum() - ms.sum() == pytest.approx(dep.sum() - gone,
                                                                                               abs=1e-12)
        m, ms = r["m_soil"], r["m_surface"]
        out.append(r)
    return out


def test_land_closes_and_reaches_the_precipitations_concentration_without_et():
    # 3 mm a day of rain, all leaving (1 mm runoff, 2 mm percolation), 300 mm held: the outflow tends to c_p
    cp = 0.0004
    r = run_land(4000, water(precip=3.0, surq=1.0), 2.0, 300.0, cp)
    assert r[-1]["concentration"][0] == pytest.approx(cp, rel=1e-9)


def test_evapotranspiration_concentrates_the_outflow():
    # 3 mm of rain, 2 mm evapotranspired (gone with no chloride), 1 mm percolating: the chloride mass balance's
    # recharge concentration c_p P / R, three times the rain's
    cp = 0.0004
    r = run_land(12000, water(precip=3.0), 1.0, 300.0, cp)  # 40 turnovers of the store
    assert r[-1]["concentration"][0] == pytest.approx(3 * cp, rel=1e-9)
    assert r[-1]["out"]["perc"][0] == pytest.approx(3 * cp * 1e-3 * 1e4, rel=1e-9)  # mass out = mass in


def test_the_surface_store_holds_snowfalls_chloride_until_melt():
    area = np.array([1e4])
    m, ms = np.zeros(1), np.zeros(1)
    # snow: 10 mm falls, nothing melts; the mass waits on the surface
    r = SA.land_day(m, ms, area, water(precip=10.0, snofall=10.0, sno=10.0), np.zeros(1), np.array([300.0]),
                    np.array([1.0]), np.zeros(1), np.zeros(1), surface_store=True)
    assert r["m_surface"][0] == pytest.approx(1.0) and r["m_soil"][0] == 0.0 and r["out"]["sur"][0] == 0.0
    # it all melts, and 4 of the 10 mm run off: 40% of the mass reaches the stream that day, the rest the soil
    r = SA.land_day(r["m_soil"], r["m_surface"], area, water(snomlt=10.0, surq_cont=4.0, sno=0.0), np.zeros(1),
                    np.array([300.0]), np.zeros(1), np.zeros(1), np.zeros(1), surface_store=True)
    assert r["out"]["sur"][0] == pytest.approx(0.4) and r["m_soil"][0] == pytest.approx(0.6)
    assert r["m_surface"][0] == pytest.approx(0.0)


def test_half_the_snowpack_melting_releases_half_its_mass():
    area = np.array([1e4])
    r = SA.land_day(np.zeros(1), np.array([2.0]), area, water(snomlt=5.0, sno=5.0), np.zeros(1), np.array([300.0]),
                    np.zeros(1), np.zeros(1), np.zeros(1), surface_store=True)
    assert r["m_surface"][0] == pytest.approx(1.0) and r["m_soil"][0] == pytest.approx(1.0)


def test_road_salt_falls_on_snow_and_ice_days():
    road = SA.RoadSalt(rate=0.01, share=0.1, snowfall_mm=1.0, precipitation_mm=1.0, tmin_c=0.0)
    w = water(4)
    w["snofall"][:] = [0.0, 2.0, 0.0, 0.0]
    w["precip"][:] = [5.0, 2.0, 5.0, 0.5]
    w["tmn"][:] = [3.0, -1.0, -2.0, -5.0]  # rain above freezing; snow; freezing rain; too little to salt
    np.testing.assert_array_equal(SA.application_days(w, road, dt.date(2020, 1, 1)), [0.0, 0.01, 0.01, 0.0])
    road.series = {dt.date(2020, 1, 2): 0.03}
    np.testing.assert_array_equal(SA.application_days(w, road, dt.date(2020, 1, 2)), [0.03] * 4)
    np.testing.assert_array_equal(SA.application_days(w, road, dt.date(2020, 1, 3)), [0.0] * 4)


def test_land_sources_timing():
    days = [dt.date(2021, 1, 1) + dt.timedelta(days=i) for i in range(365)]
    assert sum(SA.source_share("uniform", d) for d in days) == pytest.approx(1.0, rel=1e-12)
    months = [0, 0, 1, 1, 0, 0, 0, 0, 2, 0, 0, 0]  # spring and fall
    share = [SA.source_share(months, d) for d in days]
    assert sum(share) == pytest.approx(1.0, rel=1e-12)
    assert sum(v for d, v in zip(days, share) if d.month == 9) == pytest.approx(0.5)
    assert all(v == 0 for d, v in zip(days, share) if d.month == 1)


def test_fertilizer_timing_applies_the_year_once():
    b = object.__new__(SA.SaltBalance)
    b.hrus = ["a", "b"]
    src = SA.LandSource("potash", 0.0, timing="fertilizer")
    b.sources = [(src, np.array([5.0, 7.0]))]
    b.fertilized = {"potash": (0, np.zeros(2, bool))}
    fert = {d: np.zeros(2) for d in ("n", "p")}
    reads = lambda v: fert["n"] if v.endswith("fertn") else fert["p"]
    got = {}
    for i in range(731):
        d = dt.date(2021, 1, 1) + dt.timedelta(days=i)
        fert["n"][:] = [1.0 if d in (dt.date(2021, 4, 10), dt.date(2021, 5, 1)) else 0.0, 0.0]
        fert["p"][:] = [0.0, 1.0 if d == dt.date(2022, 6, 1) else 0.0]
        out = b.sources_today(d, reads)["potash"]
        for h in range(2):
            if out[h]:
                got.setdefault(h, []).append((d, out[h]))
    # a: its first fertilizer day in 2021 (not the second), and 31 December 2022, which had none; b: 31 December 2021
    # (none that year), then its fertilizer day in 2022
    assert got[0] == [(dt.date(2021, 4, 10), 5.0), (dt.date(2022, 12, 31), 5.0)]
    assert got[1] == [(dt.date(2021, 12, 31), 7.0), (dt.date(2022, 6, 1), 7.0)]


def test_annual_road_salt_is_spread_over_each_years_days_whole():
    start, last = dt.date(2010, 1, 1), dt.date(2011, 12, 31)
    days = [start + dt.timedelta(days=i) for i in range((last - start).days + 1)]
    rng = np.random.default_rng(3)
    weather = {"s1": {}, "s2": {}}
    for d in days:
        cold = d.month in (12, 1, 2)
        p = float(rng.gamma(0.5, 6.0))
        weather["s1"][d] = (p, 2.0 if cold else 25.0, -6.0 if cold else 12.0)
        weather["s2"][d] = (p, 20.0, 8.0)  # never cold: no snow, no ice, no frost
    road = SA.RoadSalt(snowfall_mm=1.0, precipitation_mm=1.0, tmin_c=0.0, spread="precipitation")
    group, weights, total, frost, unapplied = SA.road_calendar(road, ["s1", "s1", "s2"], np.array([1.0, 1.0, 1.0]),
                                                               weather, start, last)
    assert group.tolist() == [0, 0, 1] and set(unapplied) == {(1, 2010), (1, 2011)} and not frost
    for y in (2010, 2011):
        in_year = [d for d in weights if d.year == y]
        # a calendar year's days: its January and February, and December, the next winter's first half
        assert {d.month for d in in_year} == {1, 2, 12}
        assert sum(weights[d][0] for d in in_year) == pytest.approx(total[(0, y)], rel=1e-12)
        assert all(weights[d][0] == pytest.approx(weather["s1"][d][0]) for d in in_year)  # weighted by precipitation
    # a year of only dry frost days takes its amount equally over them
    dry = {"s1": {d: (0.0, 2.0, -3.0 if d.month == 1 else 5.0) for d in days}}
    _, w2, t2, frost2, _ = SA.road_calendar(road, ["s1"], np.array([1.0]), dry, start, last)
    assert set(frost2) == {(0, 2010), (0, 2011)} and t2[(0, 2010)] == 31.0


def test_annual_tables_by_calendar_year(tmp_path):
    (tmp_path / "t.csv").write_text("hru,year,kg_m2\nh1,2000,1.0\nh1,2001,2.0\nh1,2002,4.0\nh2,2001,8.0\n")
    a = SA.Annual(tmp_path / "t.csv", repeat_last=2, repeat_first=1)
    got = SA.read_annual(a, ["h1", "h2", "h3"], range(1999, 2005), "test")
    np.testing.assert_array_equal(got[2001], [2.0, 8.0, 0.0])        # an HRU the table doesn't name: zero
    np.testing.assert_array_equal(got[2004], [3.0, 4.0, 0.0])        # after the last: the mean of the last 2 years
    np.testing.assert_array_equal(got[1999], [1.0, 0.0, 0.0])        # before the first: the first year
    with pytest.raises(ValueError, match="no year 2003"):
        SA.read_annual(SA.Annual(tmp_path / "t.csv"), ["h1", "h2"], range(2002, 2004), "test")
    with pytest.raises(ValueError, match="doesn't have"):
        SA.read_annual(a, ["h1"], range(2000, 2001), "test")
    (tmp_path / "b.csv").write_text("hru,kg_m2\nh1,0.5\n")
    (tmp_path / "f.csv").write_text("year,factor\n2000,1.0\n2001,2.0\n")
    got = SA.read_annual(SA.Annual(tmp_path / "b.csv", tmp_path / "f.csv"), ["h1"], range(2000, 2002), "test")
    assert got[2001][0] == 1.0


def test_conductance_is_a_sum_over_species():
    assert SA.conductance(20.0, {"chloride": 3.0, "weathering": 1.5}, {"chloride": 10.0, "weathering": 100.0}) == 200.0


# ---------- the water bodies ----------

# a -> c, b -> c, c -> r (a reservoir) -> d -> outlet
ORDER = [0, 1, 2, 4, 3]
DOWN = [[(2, 1.0)], [(2, 1.0)], [(4, 1.0)], [], [(3, 1.0)]]


def bodies_day(m, inflow, v_out, v_end=None, v_seep=None, v_loss=None):
    n = len(m)
    z = np.zeros(n)
    return SA.water_bodies_day(ORDER, DOWN, m, inflow, z if v_end is None else v_end, v_out, z,
                               z if v_seep is None else v_seep, z if v_loss is None else v_loss)


def flows(qa=100.0, qb=50.0, q_seep=10.0):
    # each body's outflow: the reservoir (4) loses q_seep to seepage
    out = np.array([qa, qb, qa + qb, qa + qb - q_seep, qa + qb - q_seep])
    seep = np.zeros(5)
    seep[4] = q_seep
    return out, seep


def test_uniform_inflow_concentration_reaches_every_body():
    c0 = 0.02
    v_out, v_seep = flows()
    v_end = np.array([30.0, 20.0, 60.0, 50.0, 5000.0])
    inflow = np.array([100.0, 50.0, 0.0, 0.0, 0.0]) * c0
    m = np.zeros(5)
    for _ in range(3000):
        q = bodies_day(m, inflow, v_out, v_end, v_seep)
        gone = q["outlet"].sum() + q["seepage"].sum() + q["withdrawn"].sum() + q["loss"].sum()
        assert q["m"].sum() - m.sum() == pytest.approx(inflow.sum() - gone, abs=1e-12)
        m = q["m"]
    np.testing.assert_allclose(q["concentration"], c0, rtol=1e-9)


def test_a_pulse_passes_downstream_whole():
    v_out, v_seep = flows(q_seep=0.0)
    v_end = np.array([30.0, 20.0, 60.0, 50.0, 500.0])
    m = np.zeros(5)
    pulse = np.array([5.0, 0, 0, 0, 0])
    outlet, peak_day, peak = 0.0, None, 0.0
    for day in range(400):
        q = bodies_day(m, pulse if day == 0 else np.zeros(5), v_out, v_end)
        m = q["m"]
        outlet += q["outlet"].sum()
        if q["concentration"][3] > peak:
            peak, peak_day = q["concentration"][3], day
        assert outlet + m.sum() == pytest.approx(5.0, rel=1e-12)
    assert peak_day > 0 and outlet == pytest.approx(5.0, rel=1e-9)  # the reservoir delays it; it all leaves


def test_a_dry_body_keeps_its_mass_and_losses_leave_at_the_mixed_concentration():
    v_out = np.array([0.0, 50.0, 50.0, 50.0, 50.0])
    q = bodies_day(np.array([2.0, 0, 0, 0, 0]), np.zeros(5), v_out, v_loss=np.array([0.0, 50.0, 0, 0, 0]))
    assert q["m"][0] == 2.0 and q["concentration"][0] == 0.0
    q = bodies_day(np.array([0, 1.0, 0, 0, 0]), np.zeros(5), v_out, v_loss=np.array([0.0, 50.0, 0, 0, 0]))
    assert q["loss"][1] == pytest.approx(0.5) and q["out"][1] == pytest.approx(0.5)


# ---------- the recharge's mass on the cells ----------

def test_recharge_mass_is_placed_like_the_water():
    # two HRUs of 2,000 and 1,000 m2; cell 1 holds 1,000 m2 of each
    amap = AreaMap(["h1", "h2"], np.array([0, 0, 1]), np.array([0, 1, 1]), np.array([1000.0, 1000.0, 1000.0]))
    area = np.array([2000.0, 1000.0])
    perc_m = np.array([0.002, 0.004])           # m today
    c_hru = np.array([0.01, 0.04])              # kg m-3
    perc_kg = c_hru * perc_m * area
    c = SA.cell_concentration(amap, perc_kg, perc_m, area)
    water_cell = np.bincount(amap.node, weights=perc_m[amap.obj] * amap.area)
    assert c[0] == pytest.approx(0.01)
    assert c[1] == pytest.approx((0.01 * 2.0 + 0.04 * 4.0) / 6.0)  # the flow-weighted mean
    assert float(np.sum(c * water_cell)) == pytest.approx(perc_kg.sum(), rel=1e-14)


# ---------- SWAT+'s routing ----------

def con(name, rows):
    head = "      id  name  gis_id  area  lat  lon  elev  obj  wst  cst  ovfl  rule  out_tot"
    lines = [f"{name}: test", head]
    for i, (n, area, outs) in enumerate(rows, 1):
        lines.append(f"{i} {n} {i} {area} 0 0 0 {i} w 0 0 0 {len(outs)} "
                     + " ".join(f"{t} {j} {h} {f}" for t, j, h, f in outs))
    return "\n".join(lines) + "\n"


def project(tmp_path, hru_outs=None, ru2_out=None):
    (tmp_path / "file.cio").write_text(
        "file.cio: test\nconnect hru.con null rout_unit.con null null null null reservoir.con null null null null "
        "chandeg.con\n")
    (tmp_path / "hru.con").write_text(con("hru.con", [("h1", 1.0, []), ("h2", 2.0, []), ("h3", 1.0, hru_outs or [])]))
    # h1 wholly in ru1; h2 split between ru1 (1.5 ha of it) and ru2 (0.5 ha); h3 routes itself
    (tmp_path / "rout_unit.def").write_text("rout_unit.def: test\nid name elem_tot elements\n1 ru1 2 1 -2\n2 ru2 1 3\n")
    (tmp_path / "rout_unit.ele").write_text("rout_unit.ele: test\nid name obj_typ obj_id frac dlr\n"
                                            "1 h1 hru 1 0.4 0\n2 h2 hru 2 0.6 0\n3 h2 hru 2 1.0 0\n")
    (tmp_path / "rout_unit.con").write_text(con("rout_unit.con", [
        ("ru1", 2.5, [("sdc", 1, "sur", 1.0), ("sdc", 2, "lat", 1.0), ("sdc", 2, "til", 1.0), ("aqu", 1, "rhg", 1.0)]),
        ("ru2", 0.5, ru2_out or [("ru", 1, "tot", 1.0)])]))
    (tmp_path / "chandeg.con").write_text(con("chandeg.con", [("c1", 0, [("sdc", 2, "tot", 1.0)]),
                                                              ("c2", 0, [("res", 1, "tot", 1.0)]), ("c3", 0, [])]))
    (tmp_path / "reservoir.con").write_text(con("reservoir.con", [("r1", 0, [("sdc", 3, "tot", 1.0)])]))
    return tmp_path


def test_routing_follows_units_hydrographs_and_bodies(tmp_path):
    p = project(tmp_path, hru_outs=[("sdc", 3, "tot", 1.0)])
    r = SA.routing(p, ["h1", "h2", "h3"], ["c1", "c2", "c3"], ["r1"])
    assert r.bodies == ["c1", "c2", "c3", "r1"]
    order = r.order
    assert order.index(0) < order.index(1) < order.index(3) < order.index(2)  # c1, c2, r1, c3
    assert r.down[1] == [(3, 1.0)] and r.down[2] == []

    def shares(comp):
        h, b, s = r.land[comp]
        out = {}
        for hi, bi, si in zip(h, b, s):
            out.setdefault(int(hi), {})[r.bodies[bi]] = float(si)
        return out
    # h2 is in ru1 (ef 0.6 * 2.5 / 2 = 0.75) and ru2 (ef 1, all of it: frac 1), so its water is split 0.75:1 between
    # them, and ru2 sends all its water through ru1: every path of h2 follows ru1's
    assert shares("sur") == {0: {"c1": 1.0}, 1: {"c1": 1.0}, 2: {"c3": 1.0}}
    assert shares("lat") == {0: {"c2": 1.0}, 1: {"c2": 1.0}, 2: {"c3": 1.0}}


def test_routing_refuses_what_isnt_built(tmp_path):
    p = project(tmp_path, ru2_out=[("hru", 1, "tot", 1.0)])
    with pytest.raises(ValueError, match="only channels, reservoirs and routing units"):
        SA.routing(p, ["h1", "h2", "h3"], ["c1", "c2", "c3"], ["r1"])
    p = project(tmp_path, ru2_out=[("sdc", 1, "sur", 1.0)])  # h2's lateral water half lost through ru2
    with pytest.raises(ValueError, match="adding to"):
        SA.routing(p, ["h1", "h2", "h3"], ["c1", "c2", "c3"], ["r1"])


def test_the_virginia_project_routes_to_one_outlet():
    swat = ROOT.parent / "virginia-basin-study" / "rappahannock" / "coupled" / "build" / "swatplus"
    if not (swat / "rout_unit.con").exists():
        pytest.skip("the Virginia study's SWAT+ build isn't on this machine")
    names = lambda f: [x["name"] for x in SA._con(swat / f)]
    r = SA.routing(swat, names("hru.con"), names("chandeg.con"), names("reservoir.con"))
    assert [r.bodies[b] for b in range(len(r.bodies)) if not r.down[b]] == ["cha1029"]  # Fredericksburg
    for comp in SA.COMPONENTS:
        h, _, s = r.land[comp]
        np.testing.assert_allclose(np.bincount(h, weights=s), 1.0, rtol=1e-12)


# ---------- the builder's sources ----------

def test_weathering_and_stream_sources_in_the_builder(tmp_path):
    import flopy

    from basinledger.coupling.transport import add_transport
    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=str(tmp_path), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, time_units="days", nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf")
    flopy.mf6.ModflowGwfdis(gwf, nlay=2, nrow=2, ncol=3, delr=10.0, delc=20.0, top=50.0, botm=[40.0, 20.0])
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=45.0)
    flopy.mf6.ModflowGwfriv(gwf, stress_period_data={0: [((0, 0, 0), 45.0, 1.0, 44.0), ((0, 1, 0), 45.0, 1.0, 44.0)]},
                            pname="riv")
    gwt = add_transport(sim, gwf, "chloride", 0.3, 1.0, sources=[("riv", "chloride")])
    zones = np.array([[[1, 1, 2], [1, 1, 2]], [[2, 2, 2], [2, 2, 2]]])
    src = SA.add_weathering(gwt, {1: 0.0365, 2: 0.0}, zones=zones)
    rec = src.stress_period_data.get_data(0)
    assert len(rec) == 4  # zone 1's cells in layer 1
    np.testing.assert_allclose([r[1] for r in rec], 0.0365 * 10 * 20 * 10 / 365.25)
    s2 = SA.add_stream_sources(gwt, gwf, ["riv"])
    assert [tuple(r[0]) for r in s2.stress_period_data.get_data(0)] == [(0, 0, 0), (0, 1, 0)]


# ---------- on the models ----------

def salt_config(**changes):
    from basinledger import config as C
    cfg = C.load(SALT)
    cfg.run_folder = RUNS
    for k, v in changes.items():
        setattr(cfg, k, v)
    return cfg


def salt_run(cfg, name, days, users=True, restart=None, save_at=None, parameters=None):
    from basinledger.driver import Basin
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run=name, days=days, users=users, restart=restart, parameters=parameters)
    if save_at:
        b.run(save_at - b.day)
        b.save(b.folder / "state")
    else:
        b.run()
    out = {"table": b.ledger.table(), "series": {k: list(v) for k, v in b.series.items()}, "folder": b.folder,
           "salt": b.salt["chloride"], "x": b.transport["chloride"].concentration()}
    b.finalize()
    return out


@pytest.fixture(scope="module")
def coupled(swatplus):
    return salt_run(salt_config(), "salt", 730)  # two whole calendar years


@pytest.mark.parametrize("sp", ["chloride", "weathering"])
def test_every_account_closes_every_day(coupled, sp):
    t = coupled["table"]
    k = f"salt:{sp}:"
    for a in ("land:closure_share", "streams:closure_share"):
        assert np.nanmax(np.abs(t[k + a])) < 1e-12
    assert np.nanmax(np.abs(t[k + "closure_share"])) < 1e-9
    assert np.nanmax(np.abs(t[f"transport:{sp}:closure_share"])) < 1e-6
    np.testing.assert_allclose(t[k + "land:to_streams"], t[k + "streams:from_land"], rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(t[k + "streams:losses"], t[k + "aquifer:from_streams"], rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(t[k + "soil:sent"], t[k + "soil:mapped"], rtol=1e-12, atol=1e-15)
    # yesterday's soil link arrives in the land today
    np.testing.assert_allclose(t[k + "land:from_groundwater"][1:], t[k + "soil:mapped"][:-1], rtol=1e-12, atol=1e-15)
    assert np.nansum(t[f"transport:{sp}:boundary:weathering"]) > 0
    assert np.nansum(t[k + "streams:outlet"]) > 0


def test_the_chloride_paths_carry_mass(coupled):
    t = coupled["table"]
    k = "salt:chloride:"
    # every path carried mass: road salt, rain, the soil link, recharge, runoff, the river's losses, the user's return
    for a in ("land:road_salt", "land:precipitation", "land:from_groundwater", "land:recharge", "land:to_streams",
              "land:source:septic", "land:source:potash", "streams:losses", "streams:point_sources", "streams:outlet"):
        assert np.nansum(t[k + a]) > 0, a
    # weathering's ions come only from the aquifer, the user's point source and nothing on the land
    w = "salt:weathering:"
    for a in ("land:precipitation", "land:road_salt", "streams:precipitation"):
        assert np.nansum(t[w + a]) == 0, a


def test_sources_and_road_salt_applied_whole(coupled):
    t = coupled["table"]
    k = "salt:chloride:"
    days = [dt.date(2010, 1, 1) + dt.timedelta(days=i) for i in range(len(t[k + "land:road_salt"]))]
    area = 1e5  # Osu_1hru's 10 ha
    # septic: 2e-5 kg m-2 a year by land use, evenly
    assert np.nansum(t[k + "land:source:septic"][:365]) == pytest.approx(2e-5 * area, rel=1e-12)
    road = np.nan_to_num(t[k + "land:road_salt"])
    potash = np.nan_to_num(t[k + "land:source:potash"])
    for y, kg_road, kg_potash in ((2010, 0.0039, 0.004), (2011, 0.0021, 0.003)):
        # potash: the year's amount once, on the HRU's first fertilizer day
        assert [p for d, p in zip(days, potash) if d.year == y and p > 0] == [pytest.approx(kg_potash * area)]
        # road salt: each calendar year's amount over that year's days, whole
        assert sum(r for d, r in zip(days, road) if d.year == y) == pytest.approx(kg_road * area, rel=1e-12)
    assert all(r == 0 or d.month in (11, 12, 1, 2, 3, 4) for d, r in zip(days, road))


def test_multipliers_scale_the_inputs_mass(swatplus):
    from basinledger.calibration import Parameter
    base = salt_run(salt_config(), "salt-x1", 400, users=False)["table"]
    params = [Parameter("road", "salt", "chloride:road_salt", 3.0, 0.1, 10.0),
              Parameter("septic", "salt", "chloride:septic", 2.0, 0.1, 10.0)]
    scaled = salt_run(salt_config(), "salt-x2", 400, users=False, parameters=params)["table"]
    k = "salt:chloride:"
    np.testing.assert_allclose(scaled[k + "land:road_salt"], 3.0 * base[k + "land:road_salt"], rtol=1e-12)
    np.testing.assert_allclose(scaled[k + "land:source:septic"], 2.0 * base[k + "land:source:septic"], rtol=1e-12)
    np.testing.assert_array_equal(scaled[k + "land:source:potash"], base[k + "land:source:potash"])
    assert np.nansum(base[k + "land:road_salt"]) > 0 and np.nanmax(np.abs(scaled[k + "closure_share"])) < 1e-9


def test_the_recharges_mass_is_the_lands(coupled):
    t = coupled["table"]
    np.testing.assert_allclose(t["transport:chloride:recharge:sent"], t["salt:chloride:land:recharge"], rtol=0)
    np.testing.assert_allclose(t["transport:chloride:recharge:received"], t["transport:chloride:recharge:sent"],
                               rtol=1e-9, atol=1e-12)
    # one HRU over every cell: each cell's recharge takes its store's concentration
    s = coupled["salt"]
    c = s.c_cell[s.c_cell > 0]
    assert c.size == 0 or np.ptp(c) <= 1e-15 * c.max()


def test_the_point_source_carries_the_users_chloride(coupled):
    t = coupled["table"]
    np.testing.assert_allclose(t["salt:chloride:streams:point_sources"], t["user:plant:chloride:kg"], rtol=1e-12,
                               atol=1e-15)


def test_conductance_from_both_species(coupled):
    s = coupled["series"]
    cl, wx = np.array(s["chloride_mg_l:cha01"]), np.array(s["weathering_mg_l:cha01"])
    assert wx.max() > 0 and cl.max() > 0
    np.testing.assert_allclose(s["conductance_us_cm:cha01"], 20.0 + 3.0 * cl + 1.5 * wx, rtol=1e-12)


def test_gains_reach_the_stream_the_next_day(swatplus):
    cfg = salt_config()
    cfg.mf6_builder_args = dict(cfg.mf6_builder_args, top=90.0)  # the water table above the river: it gains
    t = salt_run(cfg, "salt-gaining", 120, users=False)["table"]
    sent = t["transport:chloride:stream:cha01:to_stream"]
    assert np.nansum(sent) > 0
    np.testing.assert_allclose(t["salt:chloride:streams:gains"][1:], sent[:-1], rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(t["salt:chloride:transit:gains"], sent, rtol=1e-12, atol=1e-15)
    assert np.nanmax(np.abs(t["salt:chloride:closure_share"])) < 1e-9


@pytest.mark.slow
def test_a_pulse_of_road_salt_arrives_and_passes_through(swatplus):
    """Two runs, the same water, one with 50 kg of road salt on one day: chloride is conservative and its transport
    linear, so the difference between them is the pulse's own path."""
    days, on = 300, 40
    base = salt_config()
    base.salt["chloride"].road_salt = None
    pulse = salt_config()
    start = dt.date(2010, 1, 1)
    rs = pulse.salt["chloride"].road_salt = SA.RoadSalt(share=0.05)
    rs.series = {start + dt.timedelta(days=on): 50.0 / pulse_road_m2(rs)}
    a = salt_run(base, "salt-no-pulse", days, users=False)["table"]
    b = salt_run(pulse, "salt-pulse", days, users=False)["table"]
    k = "salt:chloride:"
    d = {x: np.nan_to_num(b[k + x]) - np.nan_to_num(a[k + x]) for x in
         ("land:road_salt", "land:storage", "streams:storage", "transit:gains", "transit:soil", "streams:outlet")}
    assert d["land:road_salt"].sum() == pytest.approx(50.0) and d["land:road_salt"][on] == pytest.approx(50.0)
    aquifer = np.cumsum(-(np.nan_to_num(b["transport:chloride:storage_release"])
                          - np.nan_to_num(a["transport:chloride:storage_release"])))
    held = d["land:storage"] + d["streams:storage"] + d["transit:gains"] + d["transit:soil"] + aquifer
    total = held + np.cumsum(d["streams:outlet"])
    np.testing.assert_allclose(total[on:], 50.0, rtol=1e-9)
    np.testing.assert_allclose(total[:on], 0.0, atol=1e-9)
    # it arrives in the stream, passes, and leaves a legacy in the land and the aquifer
    assert d["streams:outlet"][on:].max() > 0 and d["streams:outlet"][:on].max() == 0
    assert aquifer[-1] > 0 and 0 < d["land:storage"][-1] < 50.0


def pulse_road_m2(rs):
    return rs.share * 1e5  # Osu_1hru's one HRU is 10 ha


def test_restart_carries_the_salt(swatplus):
    n, at = 40, 20
    a = salt_run(salt_config(), "salt-unbroken", n)
    salt_run(salt_config(), "salt-warm", n, save_at=at)
    b = salt_run(salt_config(), "salt-restarted", n, restart=RUNS / "salt-warm" / "state")
    for key in a["table"]:
        if key.startswith("salt:"):
            np.testing.assert_array_equal(a["table"][key], b["table"][key], err_msg=key)
    np.testing.assert_array_equal(a["series"]["chloride_mg_l:cha01"], b["series"]["chloride_mg_l:cha01"])
    np.testing.assert_array_equal(a["x"], b["x"])


# ---------- the config ----------

def test_config_refuses_salt_it_cant_carry(tmp_path):
    from basinledger import config as C
    text = SALT.read_text()
    for change, match in (
            (lambda s: s.replace("[transport.chloride]", "[transport.nitrate]"), "needs \\[transport.chloride\\]"),
            (lambda s: s.replace('source = "hru"', 'source = "aquifer"'), "source = \"hru\""),
            (lambda s: s.replace("[transport.chloride]", "[transport.chloride]\nrecharge = 0.01"), "come from"),
            (lambda s: s.replace("precipitation = 0.0003", "precip = 0.0003"), "unknown settings"),
            (lambda s: s.replace('spread = "precipitation"', 'share = 0.05'), "annual amount is per m2"),
            (lambda s: s.replace('annual = "maps/salt_road_annual.csv", repeat_last = 3, ', "rate = 0.002, "),
             "needs share"),
            (lambda s: s.replace("septic = 1.0, point_sources", "sewage = 1.0, point_sources"), "no input group"),
            (lambda s: s.replace('timing = "fertilizer"', 'timing = "weekly"'), "timing is"),
            (lambda s: s.replace("slopes = {chloride", "slopes = {nitrate = 1.0, chloride"), "conductance slopes")):
        bad = tmp_path / "bad.toml"
        # forward slashes: a Windows path's backslashes are escapes in a TOML string
        bad.write_text(change(text).replace("intake_user.toml", (SALT.parent / "intake_user.toml").as_posix())
                       .replace('"maps/', f'"{SALT.parent.as_posix()}/maps/'))
        with pytest.raises(ValueError, match=match):
            C.load(bad)


def test_lakes_at_a_fixed_stage_carry_salt(swatplus):
    """A lake at a fixed stage (the Virginia study's form) couples through its outlet channel: its gains join the
    channel's the next day at the cells' concentration, and its losses take the channel's, by the stream-loss source's
    entries on the lake's cells."""
    cfg = salt_config()
    cfg.mf6_builder_args = dict(cfg.mf6_builder_args, lake_stage=88.5)
    cfg.exchange = dict(cfg.exchange, lakes={"package": "lakes", "map": "maps/lakes.csv"})
    r = salt_run(cfg, "salt-lakes", 120)
    t = r["table"]
    k = "salt:chloride:"
    assert np.nanmax(np.abs(t["exchange:cha01:lakes"])) > 0
    assert np.nanmax(np.abs(t[k + "closure_share"])) < 1e-9
    assert np.nansum(t[k + "aquifer:from_streams"]) > 0
    np.testing.assert_allclose(t[k + "streams:losses"], t[k + "aquifer:from_streams"], rtol=1e-9, atol=1e-12)
    assert len(r["salt"].loss_q) == 10 + 4  # the river's entries, then the lake's
