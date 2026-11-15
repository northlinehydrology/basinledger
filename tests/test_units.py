"""A MODFLOW 6 model in its own units, converted at every link (coupling/modflow6.py, Units; docs/design.md,
"Extending an existing MODFLOW 6 model", Units). The model stays as given; Basin Ledger reads and sets it in meters
and days and the ledger stays in cubic meters.

Each test runs twins: the same model written in meters and days and in feet (and seconds), every length, area, rate,
conductance and tolerance converted, so MODFLOW 6 solves the same problem. Whatever Basin Ledger sets or reads in
meters and days must then be the same in both:

- the units read from DIS and TDIS, or stated where the model leaves them undefined (refused if not stated, or
  stated otherwise than the model's own); a real-valued variable of unknown dimension refused (tier 1);
- MODFLOW 6 alone: recharge, a well, river stages set; heads, each package's flows and storage read (tier 1);
- the stream exchange on SFR, routed (A) and at SWAT+'s stage (B), with a cut between channels (tier 1);
- solute mass: concentrations set and read, the mass budget (tier 1);
- maps on the grid in square meters and meters, whatever its units (tier 1);
- a restart's remaining simulation in the model's time units (tier 1);
- the synthetic foreign model coupled to SWAT+ for 75 days in feet and days, and in feet and seconds left undefined
  and stated in the config: the same ledger as its meters twin, closing across each converted link (tier 2);
- the start-up check prints the units and their factors (tier 2).
"""

from __future__ import annotations

import datetime as dt
import importlib.util
import io
from contextlib import redirect_stdout

import flopy
import numpy as np
import pytest

from basinledger import config as C
from basinledger.coupling.exchange import ChannelMap, StreamExchange
from basinledger.coupling.modflow6 import DAYS, METERS, Modflow6, remaining_simulation
from basinledger.coupling.transport import Transport, add_transport
from conftest import ROOT, RUNS

SYNTHETIC = ROOT / "examples" / "synthetic"
START = dt.date(2010, 1, 1)
FT = 0.3048
SI = ("meters", "days")
TWINS = [("feet", "days"), ("feet", "seconds"), ("centimeters", "hours")]


def example():
    spec = importlib.util.spec_from_file_location("existing_model", SYNTHETIC / "existing.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


E = example()


def tag(units) -> str:
    return "-".join(units)


# ---------- the units ----------

def test_units_read_from_the_model_and_their_factors(libmf6):
    ws = E.build(RUNS / "units-read", 31, START, length_units="feet")
    mf = Modflow6(ws, "kzoo")
    try:
        assert (mf.length_units, mf.time_units, mf.si) == ("feet", "days", False)
        assert mf.factor((1, 0)) == FT and mf.factor((3, -1)) == pytest.approx(FT ** 3)
        assert mf.cell_area[0] == pytest.approx(1000.0, rel=1e-9)  # m2 (flopy writes 9 digits)
        assert mf.units_summary().startswith("feet and days; converted at every link to meters and days: lengths "
                                             "and heads x 0.3048, areas x 0.092903, flows x 0.0283168")
        with pytest.raises(ValueError, match="KZOO/DIS/XORIGIN: Basin Ledger doesn't know its dimension"):
            mf.get("XORIGIN", mf.model, "DIS")
        mf.prepare_day()
        assert mf.get("NBOUND", mf.model, "WEL_CITY").tolist() == [2]  # integers are never converted
    finally:
        mf.finalize()


def test_meters_and_days_convert_nothing(libmf6):
    ws = E.build(RUNS / "units-si", 31, START)
    mf = Modflow6(ws, "kzoo")
    try:
        assert mf.si and mf.units_summary() == "meters and days, the ledger's units: nothing converted"
        assert mf.get("XORIGIN", mf.model, "DIS")[0] == E.ORIGIN[0]  # passed through, as before
    finally:
        mf.finalize()


def test_undefined_units_are_stated_or_refused(libmf6):
    ws = E.build(RUNS / "units-undefined", 31, START, length_units="feet", time_units="seconds", state_units=False)
    with pytest.raises(ValueError, match="length units are undefined \\(DIS LENGTH_UNITS\\): state them in the "
                                         "config, \\[modflow6\\] length_units"):
        Modflow6(ws, "kzoo")
    with pytest.raises(ValueError, match="time units are undefined \\(TDIS TIME_UNITS\\)"):
        Modflow6(ws, "kzoo", length_units="feet")
    mf = Modflow6(ws, "kzoo", length_units="feet", time_units="seconds")
    try:
        assert mf.units_summary().startswith("feet and seconds (length and time stated in [modflow6], undefined in "
                                             "the model)")
        assert mf.factor((3, -1)) == pytest.approx(FT ** 3 * 86400.0)
        mf.prepare_day()  # a step of 86,400 seconds is a day
        assert mf.time == pytest.approx(1.0) and mf.end == pytest.approx(31.0)
    finally:
        mf.finalize()


def test_a_unit_stated_otherwise_than_the_models_is_refused(libmf6, tmp_path):
    ws = E.build(RUNS / "units-conflict", 31, START, length_units="feet")
    with pytest.raises(ValueError, match="\\[modflow6\\] length_units = 'meters', but the model's DIS LENGTH_UNITS "
                                         "is feet"):
        Modflow6(ws, "kzoo", length_units="meters")
    with pytest.raises(ValueError, match="\\[modflow6\\] length_units = 'yards'"):
        Modflow6(ws, "kzoo", length_units="yards")
    cfg = tmp_path / "x.toml"
    cfg.write_text('name = "x"\n[modflow6]\nworkspace = "mf6"\nlength_units = "foot"\n')
    with pytest.raises(ValueError, match="\\[modflow6\\] length_units = 'foot'; one of feet, meters, centimeters"):
        C.load(cfg)


def test_years_are_refused(libmf6, tmp_path):
    ws = tmp_path / "years"
    sim = flopy.mf6.MFSimulation(sim_ws=str(ws), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, time_units="years", nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf")
    flopy.mf6.ModflowGwfdis(gwf, length_units="meters", nlay=1, nrow=1, ncol=2, top=1.0, botm=0.0)
    flopy.mf6.ModflowGwfnpf(gwf)
    flopy.mf6.ModflowGwfic(gwf)
    sim.write_simulation(silent=True)
    with pytest.raises(ValueError, match="model in years \\(TDIS TIME_UNITS\\)"):
        Modflow6(ws)


# ---------- MODFLOW 6 alone ----------

def run_alone(units, days=40):
    """The foreign model, driven as the driver drives it, in meters and days: recharge, a well and the creek's stages
    set each day; heads, every package's flows and storage read."""
    ws = E.build(RUNS / f"units-alone-{tag(units)}", days, START, length_units=units[0], time_units=units[1],
                 complexity="simple")
    mf = Modflow6(ws, "kzoo")
    out = []
    try:
        stage0 = None
        for d in range(days):
            mf.prepare_day()
            mf.set_recharge("rch_swb", np.full(mf.nodes, 0.001 + 0.001 * np.sin(d / 4.0)))
            mf.set_well("wel_city", 0, -60.0 - d)
            stage0 = mf.get("STAGE", mf.model, "RIV_CREEK") if stage0 is None else stage0
            mf.set_stage("riv_creek", stage0 + 0.3 * np.sin(d / 3.0))
            mf.solve_day()
            flows = {p: mf.simvals(p) for p in ("rch_swb", "riv_creek", "drn_tile", "ghb_regional", "chd_lake",
                                                "wel_city", "evt_riparian")}
            sto = sum(float(mf.get(v, mf.model, "STO").sum()) for v in ("STRGSS", "STRGSY"))
            out.append((mf.head(), flows, sto))
    finally:
        mf.finalize()
    return out


@pytest.fixture(scope="module")
def alone_si(libmf6):
    return run_alone(SI)


@pytest.mark.parametrize("units", TWINS, ids=tag)
def test_modflow6_alone_gives_the_same_in_meters_and_days(alone_si, units):
    other = run_alone(units)
    for (h0, q0, s0), (h1, q1, s1) in zip(alone_si, other):
        np.testing.assert_allclose(h1, h0, rtol=0, atol=1e-6)  # m
        for p in q0:
            np.testing.assert_allclose(q1[p], q0[p], rtol=1e-6, atol=1e-6, err_msg=p)  # m3 d-1
        assert s1 == pytest.approx(s0, rel=1e-6, abs=1e-3)  # m3 d-1, to the solver's closure
    _, q, _ = other[-1]
    assert q["wel_city"][0] == pytest.approx(-99.0, rel=1e-12)  # the rate set, back as set
    assert q["wel_city"][1] == pytest.approx(-25.0, rel=1e-8)  # the model's own, read in m3 d-1 (9 digits written)


def test_a_restart_counts_days_in_the_models_time_units(libmf6):
    src = E.build(RUNS / "units-restart-src", 20, START, length_units="feet", time_units="hours")
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(src), verbosity_level=0)
    sim.tdis.nper = 1  # remaining_simulation takes one period of daily steps
    sim.tdis.perioddata = [(20 * 24.0, 20, 1.0)]
    sim.write_simulation(silent=True)
    dst = remaining_simulation(src, RUNS / "units-restart-dst", 5)
    pd = flopy.mf6.MFSimulation.load(sim_ws=str(dst), verbosity_level=0).tdis.perioddata.get_data()
    assert (float(pd[0][0]), int(pd[0][1])) == (15 * 24.0, 15)


# ---------- the stream exchange on SFR ----------

MAIN, TRIBUTARY, INTO = list(range(10)), [10, 11, 12], 5


def sfr_twin(ws, units, strt=88.5):
    """A 10 x 10 aquifer of 1,000 m2 cells with an SFR network: a main stem of ten reaches down the west column and
    a tributary of three along row 4 into reach 5 (tests/test_sfr_lineup.py's), in ``units``."""
    u = E.units(*units)
    cell = 1000.0 ** 0.5
    sim = flopy.mf6.MFSimulation(sim_ws=str(ws), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, time_units=units[1], nper=1, perioddata=[(u(5.0, 0, 1), 5, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="simple", outer_dvclose=u(1e-8), inner_dvclose=u(1e-10),
                         rcloserecord=[u(1e-8 if units[1] == "days" else 1e-7, 3, -1), "strict"])
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf")
    flopy.mf6.ModflowGwfdis(gwf, length_units=units[0], nlay=1, nrow=10, ncol=10, delr=u(cell), delc=u(cell),
                            top=u(100.0), botm=u(70.0))
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=u(2.0, 1, -1))
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, sy=0.15, ss=u(1e-5, -1), transient={0: True})
    flopy.mf6.ModflowGwfic(gwf, strt=u(strt))
    flopy.mf6.ModflowGwfrcha(gwf, recharge=0.0, pname="rch")
    cells = [(0, r, 0) for r in MAIN] + [(0, 4, 3), (0, 4, 2), (0, 4, 1)]
    down = {r: r + 1 for r in MAIN[:-1]} | {10: 11, 11: 12, 12: INTO}
    up: dict[int, list[int]] = {}
    for a, b in down.items():
        up.setdefault(b, []).append(a)
    conn = [[r] + up.get(r, []) + ([-down[r]] if r in down else []) for r in range(13)]
    k = 20.0 / (2.0 * cell)
    pd = [(r, cells[r], u(cell), u(2.0), 0.001, u(88.0 - 0.001 * cell * (r % 10)), u(1.0), u(k, 1, -1), 0.035,
           len(conn[r]) - 1, 1.0, 0) for r in range(13)]
    manning = (1.0 / METERS[units[0]]) ** (1.0 / 3.0) * 86400.0 * DAYS[units[1]]  # 1.486 in feet and seconds
    flopy.mf6.ModflowGwfsfr(gwf, pname="sfr", nreaches=13, packagedata=pd, connectiondata=conn,
                            unit_conversion=manning)
    sim.write_simulation(silent=True)
    return ws


def run_sfr(units, mode):
    mf = Modflow6(sfr_twin(RUNS / f"units-sfr-{mode}-{tag(units)}", units))
    out = []
    try:
        mf.prepare_day()
        cmap = ChannelMap(["up", "trib", "down"], np.array([0] * 5 + [2] * 5 + [1] * 3))
        ex = StreamExchange(mf, "sfr", cmap, sfr=mode)
        for d in range(5):
            if d:
                mf.prepare_day()
            mf.set_recharge("rch", np.full(mf.nodes, 0.002))
            ex.set_inflow({"up": 800.0 + 50 * d, "trib": 150.0, "down": 400.0},
                          depth={"up": 0.4, "trib": 0.2, "down": 0.5})
            mf.solve_day()
            out.append((ex.after_solve(), -mf.get("QOUTFLOW", mf.model, "SFR"), mf.get("STAGE", mf.model, "SFR"),
                        mf.head()))
    finally:
        mf.finalize()
    return out


@pytest.mark.parametrize("mode", ["route", "stage"])
@pytest.mark.parametrize("units", TWINS[:2], ids=tag)
def test_sfr_routed_or_at_swats_stage_gives_the_same_exchange(libmf6, units, mode):
    si, other = run_sfr(SI, mode), run_sfr(units, mode)
    for (n0, q0, s0, h0), (n1, q1, s1, h1) in zip(si, other):
        for c in n0:
            assert n1[c] == pytest.approx(n0[c], rel=1e-6, abs=1e-6), c  # m3 d-1, positive to the stream
        np.testing.assert_allclose(q1, q0, rtol=1e-6, atol=1e-6)        # m3 d-1 leaving each reach
        np.testing.assert_allclose(s1, s0, rtol=0, atol=1e-6)           # m
        np.testing.assert_allclose(h1, h0, rtol=0, atol=1e-6)
    if mode == "stage":  # SWAT+'s depth set over the streambed top, in meters
        _, _, s, _ = other[0]
        np.testing.assert_allclose(s[MAIN[:5]] - (88.0 - 0.001 * 1000.0 ** 0.5 * np.arange(5)), 0.4, atol=1e-9)


# ---------- solute mass ----------

def run_mass(units, days=15):
    u = E.units(*units)
    ws = RUNS / f"units-mass-{tag(units)}"
    sim = flopy.mf6.MFSimulation(sim_ws=str(ws), exe_name="mf6")
    cell = 1000.0 ** 0.5
    flopy.mf6.ModflowTdis(sim, time_units=units[1], nper=1, perioddata=[(u(days, 0, 1), days, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="simple", outer_dvclose=u(1e-8), inner_dvclose=u(1e-10),
                         rcloserecord=[u(1e-8 if units[1] == "days" else 1e-7, 3, -1), "strict"])
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf")
    flopy.mf6.ModflowGwfdis(gwf, length_units=units[0], nlay=1, nrow=10, ncol=10, delr=u(cell), delc=u(cell),
                            top=u(100.0), botm=u(70.0))
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=u(2.0, 1, -1))
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, sy=0.15, ss=u(1e-5, -1), transient={0: True})
    flopy.mf6.ModflowGwfic(gwf, strt=u(88.0))
    flopy.mf6.ModflowGwfrcha(gwf, recharge=0.0, pname="rch")
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data={0: [((0, 5, 6), 0.0)]}, pname="wel")
    flopy.mf6.ModflowGwfriv(gwf, stress_period_data={0: [((0, r, 0), u(88.0), u(20.0, 2, -1), u(86.0))
                                                         for r in range(10)]}, pname="riv")
    add_transport(sim, gwf, "chloride", 0.25, u(5.0), strt=u(0.01, -3),
                  sorption={"distcoef": u(1e-4, 3), "bulk_density": u(1800.0, -3)},
                  sources=[("rch", "chloride"), ("riv", "chloride")],
                  # the residual in kg d-1: in seconds the solver gets no closer than about 1e-11 kg s-1, 1e-6 kg d-1
                  ims=dict(outer_dvclose=u(1e-10, -3), inner_dvclose=u(1e-12, -3),
                           rcloserecord=[u(1e-10 if units[1] == "days" else 1e-6, 0, -1), "strict"]))
    sim.write_simulation(silent=True)
    mf = Modflow6(ws)
    t = Transport(mf, "chloride")
    out = []
    try:
        for d in range(days):
            mf.prepare_day()
            mf.set_recharge("rch", np.full(mf.nodes, 0.003))
            mf.set_well("wel", 0, -900.0 if d % 6 < 3 else 0.0)  # hard enough that the river loses
            t.set_inflow("rch", 0.03)
            t.set_inflow("riv", 0.05)
            mf.solve_day()
            out.append((t.budget("rch", {"pw1": ("wel", 0)}, [("riv", np.zeros(10, int), ["cha01"])]),
                        t.concentration(), float(mf.simvals("rch").sum())))
    finally:
        mf.finalize()
    return out


@pytest.mark.parametrize("units", TWINS[:2], ids=tag)
def test_solute_mass_is_the_same_in_kilograms(libmf6, units):
    si, other = run_mass(SI), run_mass(units)
    for (b0, c0, _), (b1, c1, r1) in zip(si, other):
        np.testing.assert_allclose(c1, c0, rtol=1e-6, atol=1e-12)  # kg m-3
        for k in ("recharge", "storage_release"):
            assert getattr(b1, k) == pytest.approx(getattr(b0, k), rel=1e-6, abs=1e-12), k
        for k in ("wells", "to_stream", "from_stream", "leakage_m3"):
            for name, v in getattr(b0, k).items():
                assert getattr(b1, k)[name] == pytest.approx(v, rel=1e-6, abs=1e-12), (k, name)
        assert abs(b1.closure_share) < 1e-8
        assert b1.recharge == pytest.approx(0.03 * r1, rel=1e-12)  # kg = kg m-3 x m3, across the converted link
    assert sum(b.leakage_m3["cha01"] > 0 for b, _, _ in other) > 3


# ---------- maps ----------

def placed(grid, points, meters):
    """Points given in meters from the grid's corner (E.world's, unrotated), in the grid's world coordinates and
    units, through its own origin and rotation (flopy writes a feet origin to nine figures, 0.03 ft off the exact)."""
    import shapely
    x, y = grid.get_coords(*(np.asarray(points, float).T / meters))
    return shapely.Polygon(np.c_[x, y]) if len(points) > 2 else shapely.LineString(np.c_[x, y])


@pytest.mark.parametrize("grid", ["dis", "disv"])
def test_maps_on_a_grid_in_feet_are_in_square_meters_and_meters(grid):
    pytest.importorskip("shapely")
    from basinledger import gridmaps as G
    maps = {}
    for units, state in ((SI, True), (("feet", "days"), True), (("feet", "days"), False)):
        ws = E.build(RUNS / f"units-maps-{grid}-{tag(units)}-{state}", 31, START, grid=grid, length_units=units[0],
                     state_units=state)
        _, gwf = G.load(ws)
        lu, g, m = None if state else "feet", gwf.modelgrid, METERS[units[0]]
        c = E.CELL
        aquifer = placed(g, [(c, 11 * c), (11 * c, 11 * c), (11 * c, c), (c, c)], m)  # E.aquifer_polygon
        creek = placed(g, [(1.5 * c, 11 * c), (1.5 * c, c)], m)                        # E.creek_line
        if not state:
            with pytest.raises(ValueError, match="length units are undefined: give length_units"):
                G.hru_map(g, {"aqu001": aquifer})
        hm = G.hru_map(g, {"aqu001": aquifer}, length_units=lu)
        sm = G.channel_map(g, G.stream_entries(gwf, "riv_creek"), {"cha01": creek}, length_units=lu)
        maps[(units, state)] = (hm, sm)
    hm0, sm0 = maps[(SI, True)]
    assert hm0.mapped()["aqu001"] == pytest.approx(1e5, rel=1e-9)  # 10 ha
    for key, (hm, sm) in maps.items():
        assert [r[:3] for r in hm.rows] == [r[:3] for r in hm0.rows], key
        np.testing.assert_allclose([r[3] for r in hm.rows], [r[3] for r in hm0.rows], rtol=1e-9)
        assert hm.problems({"aqu001": 1e5}, rtol=1e-6) == []
        np.testing.assert_allclose([r[2] for r in sm.rows], [r[2] for r in sm0.rows], rtol=1e-9)
        assert sm.in_grid["cha01"] == pytest.approx(sm0.in_grid["cha01"], rel=1e-9)
        assert sm.problems() == []
    with pytest.raises(ValueError, match="length_units 'meters', but the grid is in feet"):
        G.hru_map(G.load(RUNS / f"units-maps-{grid}-feet-days-True")[1].modelgrid, {}, length_units="meters")


def test_sfr_reach_lengths_in_meters(libmf6):
    pytest.importorskip("shapely")
    from basinledger import gridmaps as G
    _, gwf = G.load(sfr_twin(RUNS / "units-sfr-lengths", ("feet", "days")))
    assert [e.length for e in G.stream_entries(gwf, "sfr")] == pytest.approx([1000.0 ** 0.5] * 13, rel=1e-9)


# ---------- the foreign model coupled, in feet (tier 2) ----------

COUPLED = {"feet-days": dict(length_units="feet", complexity="simple"),
           "feet-seconds-undefined": dict(length_units="feet", time_units="seconds", state_units=False,
                                          complexity="simple")}


def config(**builder_args):
    cfg = C.load(SYNTHETIC / "existing.toml")
    cfg.run_folder = RUNS
    cfg.mf6_builder_args = dict(builder_args)
    if builder_args.get("state_units") is False:
        cfg.mf6_length_units = builder_args.get("length_units")
        cfg.mf6_time_units = builder_args.get("time_units", "days")
    return cfg


@pytest.fixture(scope="module")
def coupled(swatplus):
    from basinledger.driver import Basin
    out = {}
    for name, args in {"meters-days": dict(complexity="simple"), **COUPLED}.items():
        with pytest.warns(Warning):
            b = Basin(config(**args), run=f"units-coupled-{name}", days=75)
        try:
            b.run()
            out[name] = (b.ledger.totals(), b.mf6.head(), b.mf6.units_summary())
        finally:
            b.finalize()
    return out


@pytest.mark.parametrize("name", list(COUPLED))
def test_coupled_in_feet_the_ledger_is_the_meters_twins(coupled, name):
    t0, h0, _ = coupled["meters-days"]
    t1, h1, units = coupled[name]
    assert not units.startswith("meters")
    assert set(t1) == set(t0)
    for k, v in t0.items():
        assert t1[k] == pytest.approx(v, rel=1e-6, abs=1e-4), k  # m3 over 75 days
    np.testing.assert_allclose(h1, h0, rtol=0, atol=1e-6)        # m
    assert sum(abs(v) for k, v in t0.items() if k.startswith("exchange:")) > 0  # the creek exchanged water


@pytest.mark.parametrize("name", list(COUPLED))
def test_coupled_in_feet_the_ledger_closes_across_each_converted_link(coupled, name):
    t, _, _ = coupled[name]
    assert t["recharge:received"] == pytest.approx(t["recharge:made"], rel=1e-9)            # SWAT+ m3 -> ft/d -> m3
    assert t["route:well:city1:delivered"] == pytest.approx(t["route:well:city1:requested"], rel=1e-9)  # m3 d-1
    assert t["groundwater:boundary:wel_city"] == pytest.approx(-25.0 * 75, rel=1e-8)  # its own, 9 digits written
    assert abs(t["groundwater:closure_share"]) < 1e-6 * 75


def test_the_check_prints_the_units_and_their_factors(swatplus):
    from basinledger.driver import Basin
    units = []
    found, notes = Basin.check(config(**COUPLED["feet-seconds-undefined"]), run="units-check", units=units)
    assert found == []
    assert units == [next(n for n in notes if n.startswith("MODFLOW 6 in "))[len("MODFLOW 6 in "):]]
    assert units[0].startswith("feet and seconds (length and time stated in [modflow6], undefined in the model); "
                               "converted at every link")
    cfg = config(**COUPLED["feet-seconds-undefined"])
    cfg.mf6_length_units = None
    found, _ = Basin.check(cfg, run="units-check-unstated")
    assert len(found) == 1 and "length units are undefined" in found[0]


def test_the_command_line_check_prints_the_units(swatplus):
    from basinledger.__main__ import main
    out = io.StringIO()
    with redirect_stdout(out):
        main(["check", str(SYNTHETIC / "existing.toml")])
    assert "MODFLOW 6's units: meters and days, the ledger's units: nothing converted" in out.getvalue()
