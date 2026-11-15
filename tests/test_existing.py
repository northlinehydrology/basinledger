"""The second way in: a MODFLOW 6 model Basin Ledger didn't write, taken as given (docs/design.md, "Extending an
existing MODFLOW 6 model"; examples/synthetic/existing.py, a model in someone else's style: rotated, two layers with
inactive cells, monthly stress periods, its own package names and boundaries).

- the start-up check: each package listed with what Basin Ledger takes over and what it leaves as the model's;
  what would count water twice refused (coupling/existing.py; tier 1, and tier 2 on the example);
- maps on the model's own grid, DIS and DISV (gridmaps.py): HRUs to cells by area on the top active layer, with what
  falls outside or on inactive columns reported; channels to stream entries by length, entries no channel reaches
  reported and not guessed, and each channel's length checked (tier 1);
- the example coupled both ways across monthly stress periods, on DIS and DISV: the ledger closes every day, the
  model's own boundaries are booked as water from or to outside, and the two grids agree (tier 2).
"""

from __future__ import annotations

import copy
import csv
import datetime as dt
import importlib.util
import io
import shutil
import warnings
from contextlib import redirect_stdout

import numpy as np
import pytest

from basinledger import config as C
from basinledger.coupling.existing import PackageRole, problems, table
from conftest import ROOT, RUNS

shapely = pytest.importorskip("shapely")
from basinledger import gridmaps as G  # noqa: E402  (needs shapely)

SYNTHETIC = ROOT / "examples" / "synthetic"
CONFIG = SYNTHETIC / "existing.toml"
START = dt.date(2010, 1, 1)


def example():
    spec = importlib.util.spec_from_file_location("existing_model", SYNTHETIC / "existing.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


E = example()


@pytest.fixture(scope="module", params=["dis", "disv"])
def built(request):
    ws = RUNS / f"existing-{request.param}"
    E.build(ws, 75, START, grid=request.param)
    return request.param, *G.load(ws)


# ---------- the start-up check, on its own ----------

def test_a_package_is_basin_ledgers_the_models_or_partly_both():
    rows = [PackageRole("rch_swb", "rch", 144, 144, "recharge from SWAT+"),
            PackageRole("riv_creek", "riv", 12, 10, "stream exchange with SWAT+'s channels"),
            PackageRole("ghb_regional", "ghb", 10, 0, "")]
    assert [r.role for r in rows] == ["Basin Ledger's", "partly Basin Ledger's", "the model's"]
    assert "10 of 12 entries" in rows[1].describe() and "groundwater:boundary:riv_creek" in rows[1].describe()
    assert "groundwater:boundary:ghb_regional" in rows[2].describe()
    t = table(rows)
    assert t.splitlines()[0].startswith("package") and len(t.splitlines()) == 4


def test_what_would_count_water_twice_is_refused():
    left = lambda name, kind: PackageRole(name, kind, 5, 0, "")  # noqa: E731
    bad, notes = problems([left("rch_septic", "rch"), left("uzf", "uzf"), left("evt", "evt"),
                           left("sfr_model", "sfr"), left("lak_model", "lak"), left("ghb", "ghb")], False)
    assert [b.split(":")[0] for b in bad] == ["rch_septic", "uzf"]
    assert [n.split(":")[0] for n in notes] == ["evt", "sfr_model", "lak_model"]
    bad, _ = problems([left("evt", "evt")], groundwater_to_soil=True)  # SWAT+'s plants take it from the soil too
    assert bad and "EVT" in bad[0]
    # a recharge package Basin Ledger sets is no problem
    assert problems([PackageRole("rch", "rch", 5, 5, "recharge from SWAT+")], False) == ([], [])


# ---------- maps on the model's own grid ----------

def test_hrus_to_cells_on_the_top_active_layer_of_a_rotated_grid(built):
    grid_type, _, gwf = built
    hm = G.hru_map(gwf.modelgrid, E.hru_polygons())
    got = hm.mapped()
    assert got["hru01"] == pytest.approx(6e4, rel=1e-9) and got["hru02"] == pytest.approx(4e4, rel=1e-9)
    assert len(hm.rows) == 100 and {lay for _, lay, _, _ in hm.rows} == {0}  # 10 x 10 active top cells
    assert hm.problems() == []
    assert sum(hm.slivers.values()) < 1e-3  # rounding along the rotated cell edges, left out and reported


def test_ground_outside_the_grid_or_over_inactive_columns_is_reported(built):
    _, _, gwf = built
    grid = gwf.modelgrid
    # aqu001 two rows too far north: one row over the inactive ring (top layer inactive, bottom active), one off the grid
    poly = shapely.Polygon([E.world(-1, 1), E.world(-1, 11), E.world(9, 11), E.world(9, 1)])
    hm = G.hru_map(grid, {"aqu001": poly})
    assert hm.outside["aqu001"] == pytest.approx(1e4, rel=1e-6)
    assert {lay for _, lay, _, _ in hm.rows} == {0, 1}  # the ring's columns map to their top active layer, the bottom
    assert hm.inactive["aqu001"] == 0.0
    p = hm.problems({"aqu001": 1e5})
    assert len(p) == 1 and "outside the grid" in p[0]
    # a grid with a column inactive in every layer: its ground is reported, not placed
    holed = copy.deepcopy(grid)
    idomain = np.array(grid.idomain).reshape(2, -1)
    idomain[:, 5 * E.N + 5] = 0
    holed._idomain = idomain.reshape(np.shape(grid.idomain))
    hm = G.hru_map(holed, {"aqu001": E.aquifer_polygon()})
    assert hm.inactive["aqu001"] == pytest.approx(1e3, rel=1e-9) and len(hm.rows) == 99
    assert "over inactive columns" in hm.problems({"aqu001": 1e5})[0]


def test_channels_to_the_models_stream_entries(built):
    grid_type, _, gwf = built
    entries = G.stream_entries(gwf, "riv_creek")
    assert len(entries) == 12  # the creek's ten and the tributary's two, over all three periods
    sm = G.channel_map(gwf.modelgrid, entries, {"cha01": E.creek_line()})
    assert len(sm.rows) == 10 and sm.on_entries["cha01"] == pytest.approx(10 * E.CELL, rel=1e-9)
    assert sm.in_grid["cha01"] == pytest.approx(10 * E.CELL, rel=1e-9)
    tributary = [E.cellid(c, grid_type) for c in E.TRIBUTARY]
    assert [e.cell for e in sm.unreached] == tributary  # reported, not guessed
    assert sm.problems() == [] and "2 entries no channel reaches" in sm.report()
    # within reach_m, the tributary's first cell (half a cell from the creek's line) is reached; its second isn't
    sm = G.channel_map(gwf.modelgrid, entries, {"cha01": E.creek_line()}, reach_m=20.0)
    assert [e.cell for _, e, _ in sm.by_distance] == tributary[:1]
    assert [e.cell for e in sm.unreached] == tributary[1:]
    # tile drains: no channel near them
    assert len(G.channel_map(gwf.modelgrid, G.stream_entries(gwf, "drn_tile"), {"cha01": E.creek_line()}).rows) == 0


def test_a_channel_the_model_has_no_entries_along_fails_the_length_check(built):
    _, _, gwf = built
    longer = shapely.LineString([E.world(1.0, 1.5), E.world(10.5, 1.5), E.world(10.5, 9.0)])  # turns east on row 10
    sm = G.channel_map(gwf.modelgrid, G.stream_entries(gwf, "riv_creek"), {"cha01": longer})
    p = sm.problems()
    assert len(p) == 1 and p[0].startswith("cha01:") and "through cells with stream entries" in p[0]


def test_a_shared_cell_goes_to_the_longest_and_is_reported(built):
    _, _, gwf = built
    lines = {"cha01": E.creek_line(),
             "cha02": shapely.LineString([E.world(3.5, 4.0), E.world(3.5, 1.8)])}  # into the creek's row-3 cell
    sm = G.channel_map(gwf.modelgrid, G.stream_entries(gwf, "riv_creek"), lines)
    assert len(sm.shared) == 1 and list(sm.shared[0][1]) == ["cha01", "cha02"]
    assert {c for c, _, _ in sm.rows} == {"cha01"}


def test_sfr_reaches_to_channels_with_their_lengths(libmf6):
    from test_existing_api import _build
    ws = RUNS / "existing-sfr-map"
    _build(ws)
    _, gwf = G.load(ws)
    grid = gwf.modelgrid
    entries = G.stream_entries(gwf, "river")
    assert [e.reach for e in entries] == [0, 1, 2] and all(e.length == 100.0 for e in entries)
    line = shapely.LineString([(300.0, 750.0), (600.0, 750.0)])  # down row 2's middle, columns 3 to 5
    sm = G.channel_map(grid, entries, {"cha01": line})
    assert sm.sfr and sm.reach_length["cha01"] == pytest.approx(300.0) and sm.problems() == []
    out = RUNS / "existing-sfr-map.csv"
    sm.write(out)
    assert next(csv.reader(open(out))) == ["object", "reach", "length_m"]
    short = shapely.LineString([(300.0, 750.0), (560.0, 750.0)])  # stops 40 m into the last reach's cell
    p = G.channel_map(grid, entries, {"cha01": short}).problems()
    assert any("reaches are 300 m long" in x for x in p)


def test_the_examples_maps_are_what_the_helpers_write(built):
    grid_type, sim, _ = built
    folder = RUNS / f"existing-maps-{grid_type}"
    hm, sm = E.write_maps(folder, sim.simulation_data.mfpath.get_sim_path())
    tag = "" if grid_type == "dis" else "_disv"
    rows = list(csv.reader(open(folder / f"existing_recharge{tag}.csv")))
    assert rows[0] == (["object", "layer", "row", "col", "area_m2"] if grid_type == "dis"
                       else ["object", "layer", "icell2d", "area_m2"])
    assert len(rows) == 101
    if grid_type == "dis":  # the committed maps, which existing.toml reads
        for name in ("existing_recharge.csv", "existing_river.csv"):
            assert (folder / name).read_text() == (SYNTHETIC / "maps" / name).read_text()


def test_one_step_a_month_is_refused_by_modflow6_alone(libmf6):
    """A run rewrites the time discretization to daily steps first (coupling/timesteps.py, tests/test_timesteps.py);
    MODFLOW 6 alone, on the model as given, refuses a step longer than a day."""
    from basinledger.coupling.modflow6 import Modflow6
    ws = RUNS / "existing-monthly"
    E.build(ws, 59, START, daily=False)
    mf = Modflow6(ws, "kzoo")
    try:
        with pytest.raises(ValueError, match="time step is 31.0 days.*timesteps.daily_steps"):
            mf.prepare_day()
    finally:
        mf.finalize()


# ---------- the example, coupled (tier 2) ----------

def config(grid: str = "dis", **builder_args):
    cfg = C.load(CONFIG)
    cfg.run_folder = RUNS
    cfg.mf6_builder_args = {"grid": grid, **builder_args}
    if grid == "disv":
        ws = RUNS / "existing-disv-maps-model"
        E.build(ws, 31, START, grid="disv")
        E.write_maps(RUNS / "existing-disv-maps", ws)
        cfg.recharge_map = str(RUNS / "existing-disv-maps" / "existing_recharge_disv.csv")
        cfg.exchange = dict(cfg.exchange, map=str(RUNS / "existing-disv-maps" / "existing_river_disv.csv"))
        cfg.wells[0].cell = (1, 5 * E.N + 6)
    return cfg


EXPECTED = {"chd_lake": "the model's", "drn_tile": "the model's", "evt_riparian": "the model's",
            "ghb_regional": "the model's", "rch_swb": "Basin Ledger's", "riv_creek": "partly Basin Ledger's",
            "wel_city": "partly Basin Ledger's"}


@pytest.mark.parametrize("grid", ["dis", "disv"])
def test_the_start_up_check_lists_what_basin_ledger_takes_over(swatplus, grid):
    from basinledger.driver import Basin
    rows = []
    found, notes = Basin.check(config(grid), run=f"existing-check-{grid}", takeover=rows)
    assert found == []
    assert {r.package: r.role for r in rows} == EXPECTED
    riv = next(r for r in rows if r.package == "riv_creek")
    assert (riv.entries, riv.taken) == (12, 10)
    assert any(n.startswith("evt_riparian: groundwater ET stays the model's") for n in notes)


def test_a_second_recharge_package_stops_the_start_up(swatplus):
    from basinledger.driver import Basin
    found, _ = Basin.check(config(septic=True), run="existing-check-septic")
    assert len(found) == 1 and "rch_septic: a recharge package left as the model's recharges 1 cells" in found[0]


def test_a_second_recharge_package_off_swats_ground_is_the_models(swatplus):
    """A recharge package on ground SWAT+'s recharge doesn't cover (a cell of the bottom layer under the top's
    inactive ring, as a catchment outside SWAT+'s project) is the model's, booked as water from outside."""
    from basinledger.driver import Basin
    cfg = config(septic=(1, 0, 3))
    found, notes = Basin.check(cfg, run="existing-check-septic-outside")
    assert found == []
    assert "rch_septic: recharge left as the model's on 1 cells off SWAT+'s recharge, booked as water from outside" \
        in notes
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run="existing-run-septic-outside", days=40)
        try:
            b.run()
            t = b.ledger.totals()
        finally:
            b.finalize()
    assert t["groundwater:boundary:rch_septic"] == pytest.approx(0.002 * 1000.0 * 40, rel=1e-9)
    assert t["recharge:received"] == pytest.approx(t["recharge:made"], rel=1e-9)
    assert abs(t["groundwater:closure_share"]) < 1e-6 * 40


def test_where_a_recharge_package_lies_against_swats():
    from basinledger.coupling import existing as EX
    share = np.array([1.0, 0.6, 0.6, 0.0, 0.9995])
    area = np.full(5, 100.0)
    # node 0 wholly SWAT+'s; 1 with the rest (0.4) its own; 2 over by 0.1; 3 off SWAT+'s; 4 within the tolerance
    cover = EX.recharge_cover(share, np.array([0, 1, 2, 3, 4, -1]), np.array([1.0, 0.4, 0.5, 1.0, 0.001, 1.0]), area)
    assert cover[:3] == (6, 2, 2) and cover[3] == pytest.approx(100.0 + 10.0 + 0.05)
    # SWAT+'s map alone over a cell (its areas scaled to SWAT+'s): a package entry of no share there adds nothing
    assert EX.recharge_cover(np.array([1.013]), np.array([0]), np.array([0.0]), np.array([100.0])) == (1, 0, 1, 0.0)
    rows = [EX.PackageRole("rch_x", "rch", 3, 0, "")]
    bad, notes = EX.problems(rows, False, {"rch_x": (3, 0, 2, 0.0)})
    assert bad == [] and "2 of them shared with SWAT+'s recharge" in notes[0]
    bad, _ = EX.problems(rows, False, {"rch_x": (3, 2, 0, 110.0)})
    assert "rch_x: a recharge package left as the model's recharges 2 cells" in bad[0] and "(110 m2 recharged" in bad[0]
    bad, _ = EX.problems(rows, False)
    assert "rch_x: a recharge package left as the model's;" in bad[0]


def test_the_command_line_check_prints_the_packages(swatplus):
    from basinledger.__main__ import main
    out = io.StringIO()
    with redirect_stdout(out):
        code = main(["check", str(CONFIG)])
    text = out.getvalue()
    assert code == 0 and "MODFLOW 6's packages:" in text
    assert "10 of 12 entries Basin Ledger's" in text and "no problems found" in text


@pytest.fixture(scope="module")
def runs(swatplus):
    from basinledger.driver import Basin
    out = {}
    for grid in ("dis", "disv"):
        with pytest.warns(Warning):
            b = Basin(config(grid), run=f"existing-run-{grid}", days=75)
        try:
            b.run()
            out[grid] = (b.ledger.totals(), b.mf6.head(), b.takeover)
        finally:
            b.finalize()
    return out


def test_coupled_across_monthly_stress_periods_the_ledger_closes(runs):
    t, _, _ = runs["dis"]
    # the run went through January, February and half of March: three stress periods, the ledger checked every day
    assert t["recharge:received"] == pytest.approx(t["recharge:made"], rel=1e-9)
    assert t["route:well:city1:delivered"] == pytest.approx(t["route:well:city1:requested"])
    assert abs(t["groundwater:closure_share"]) < 1e-6 * 75


def test_monthly_steps_growing_give_the_daily_models_run(runs):
    """The same model with five steps a month growing by 1.2, as a study's model usually has: the run rewrites it to
    daily steps in its copy and says so, and the run is the daily model's."""
    from basinledger.driver import Basin
    with pytest.warns(Warning) as caught:
        b = Basin(config(daily=False, steps=5, tsmult=1.2), run="existing-run-monthly", days=75)
    try:
        b.run()
        t = b.ledger.totals()
    finally:
        b.finalize()
    notes = [str(w.message) for w in caught]
    assert ("MODFLOW 6's time steps rewritten to one day each in the run's copy: 3 stress periods, 75 days (they "
            "had 15 steps, growing in 3 periods); each period's input stays as given") in notes
    daily = runs["dis"][0]
    assert set(t) == set(daily)
    for k, v in daily.items():
        assert t[k] == pytest.approx(v, rel=1e-12, abs=1e-12), k


def test_a_models_own_folder_with_its_own_start_date(swatplus):
    """A model handed over as a folder (``workspace``), its periods monthly and its dates two years after the run's:
    the check rewrites its steps and notes that its stress periods fall before its dates."""
    from basinledger.driver import Basin
    ws = RUNS / "existing-handed-over"
    if ws.exists():
        shutil.rmtree(ws)
    E.build(ws, 75, dt.date(2012, 1, 1), daily=False)
    cfg = config()
    cfg.mf6_builder, cfg.mf6_workspace = None, ws
    found, notes = Basin.check(cfg, run="existing-check-handed-over")
    assert found == []
    assert any(n.startswith("MODFLOW 6's time steps rewritten to one day each") for n in notes)
    assert ("MODFLOW 6's time starts on 2012-01-01 (TDIS START_DATE_TIME) and the run on 2010-01-01: the model's "
            "first day is taken as the run's, so its stress periods fall 730 days earlier than its dates") in notes
    assert "Basin Ledger" not in (ws / "kzoo_regional.tdis").read_text()  # the model handed over is as it was


def test_a_restart_across_stress_periods_continues_the_unbroken_run(swatplus):
    """A model handed over as a folder, five growing steps a month: saved on February 10th and carried on in a new
    simulation of the days left (its January left out, February shortened, every package's periods renumbered), the
    run is the unbroken run's, to the last digit."""
    from basinledger.driver import Basin
    ws = RUNS / "existing-restart-model"
    if ws.exists():
        shutil.rmtree(ws)
    E.build(ws, 75, START, daily=False, steps=5, tsmult=1.2)
    cfg = config()
    cfg.mf6_builder, cfg.mf6_workspace = None, ws
    a = Basin(cfg, run="existing-restart-unbroken", days=75)
    a.run()
    ta, ha = a.ledger.table(), a.mf6.head()
    a.finalize()
    w = Basin(cfg, run="existing-restart-warm", days=75)
    w.run(40)
    w.save(RUNS / "existing-restart-state")
    w.finalize()
    b = Basin(cfg, run="existing-restart-restarted", days=75, restart=RUNS / "existing-restart-state")
    assert b.day == 40 and b.mf6.end == 35.0
    b.run()
    tb, hb = b.ledger.table(), b.mf6.head()
    b.finalize()
    assert np.array_equal(ha, hb)
    assert set(ta) == set(tb)
    for key in ta:
        assert np.array_equal(ta[key], tb[key], equal_nan=True), key


def test_another_groundwater_model_in_the_simulation_is_outside(swatplus):
    """A second groundwater model joined by a GWF-GWF exchange (``neighbor=True``), not coupled: what crosses to it is
    booked as water from or to outside, so the coupled model's ledger closes every day, and the check says so."""
    from basinledger.driver import Basin
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        b = Basin(config(neighbor=True), run="existing-run-neighbor", days=75)
        try:
            b.run()  # the takeover's notes come with the first day, once the packages have their entries
            t = b.ledger.totals()
        finally:
            b.finalize()
    notes = [str(w.message) for w in caught]
    assert ("gwf-gwf_1: groundwater to and from kzoo_east, another model of the simulation that Basin Ledger doesn't "
            "couple, booked as from or to outside") in notes
    assert t["groundwater:boundary:gwf-gwf_1"] > 0.0  # the regional head to the east feeds the aquifer
    assert abs(t["groundwater:closure_share"]) < 1e-6 * 75
    assert t["recharge:received"] == pytest.approx(t["recharge:made"], rel=1e-9)


def test_the_add_settings_read(tmp_path):
    base = CONFIG.read_text().replace('[recharge]\n', '[recharge]\nadd = true\n')
    base = base.replace('cell = [1, 5, 6]', 'cell = [1, 5, 6]\nadd = true')
    (tmp_path / "x.toml").write_text(base.replace('"maps/', f'"{SYNTHETIC}/maps/'))
    cfg = C.load(tmp_path / "x.toml")
    assert cfg.recharge_add and cfg.wells[0].add
    (tmp_path / "x.toml").write_text(base.replace("add = true", 'add = "yes"', 1))
    with pytest.raises(ValueError, match="\\[recharge\\] add must be true or false"):
        C.load(tmp_path / "x.toml")


def with_a_new_well(cfg, cell=(1, 3, 3), package="wel_users"):
    """A second plant drawing from a well the model doesn't have, which Basin Ledger adds."""
    cfg.wells.append(C.WellConfig("new1", package, cell, add=True))
    plant = copy.deepcopy(cfg.users[0])
    plant.name, plant.withdrawal = "plant2", "well:new1"
    cfg.users.append(plant)
    return cfg


def test_a_well_the_model_doesnt_have_is_added(swatplus):
    from basinledger.driver import Basin
    rows = []
    found, notes = Basin.check(with_a_new_well(config()), run="existing-check-new-well", takeover=rows)
    assert found == []
    assert ("wel_users: a well package Basin Ledger adds to the run's copy for wells the model doesn't have: new1 at "
            "[1, 3, 3]") in notes
    assert {r.package: (r.role, r.entries) for r in rows}["wel_users"] == ("Basin Ledger's", 1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(with_a_new_well(config()), run="existing-run-new-well", days=75)
        try:
            b.run()
            t = b.ledger.totals()
        finally:
            b.finalize()
    assert t["route:well:new1:requested"] > 0.0
    assert t["route:well:new1:delivered"] == pytest.approx(t["route:well:new1:requested"])
    assert t["route:well:city1:delivered"] == pytest.approx(t["route:well:city1:requested"])
    assert abs(t["groundwater:closure_share"]) < 1e-6 * 75
    assert "groundwater:boundary:wel_users" in t and t["groundwater:boundary:wel_users"] == 0.0  # all of it a user's


def test_an_added_well_refused(swatplus):
    from basinledger.driver import Basin
    found, _ = Basin.check(with_a_new_well(config(), package="wel_city"), run="existing-check-new-well-bad")
    assert len(found) == 1 and "the model has a package 'wel_city'" in found[0]
    found, _ = Basin.check(with_a_new_well(config(), cell=(1, 3)), run="existing-check-new-well-bad")
    assert len(found) == 1 and "cell [1, 3] isn't a DIS cell" in found[0]


def test_a_recharge_package_for_a_model_without_one(swatplus):
    """The model handed over without its recharge (rch_swb out of its name file): ``[recharge] add = true`` puts a
    package of Basin Ledger's own on the map's cells, and SWAT+'s recharge reaches the aquifer through it."""
    from basinledger.driver import Basin
    ws = RUNS / "existing-no-recharge"
    if ws.exists():
        shutil.rmtree(ws)
    E.build(ws, 75, START)
    nam = ws / "kzoo.nam"
    nam.write_text("\n".join(x for x in nam.read_text().splitlines() if "rch_swb" not in x) + "\n")
    cfg = config()
    cfg.mf6_builder, cfg.mf6_workspace = None, ws
    cfg.recharge_package, cfg.recharge_add = "rch_swat", True
    rows = []
    found, notes = Basin.check(cfg, run="existing-check-no-recharge", takeover=rows)
    assert found == []
    assert ("rch_swat: a recharge package Basin Ledger adds to the run's copy, the model having none: 100 cells of the "
            "recharge map") in notes
    assert {r.package: r.role for r in rows}["rch_swat"] == "Basin Ledger's"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run="existing-run-no-recharge", days=75)
        try:
            b.run()
            t = b.ledger.totals()
        finally:
            b.finalize()
    assert t["recharge:received"] > 0.0
    assert t["recharge:received"] == pytest.approx(t["recharge:made"], rel=1e-9)
    assert abs(t["groundwater:closure_share"]) < 1e-6 * 75
    cfg = config()  # the model with its own recharge: adding one would count it twice
    cfg.recharge_package, cfg.recharge_add = "rch_swat", True
    found, _ = Basin.check(cfg, run="existing-check-no-recharge-bad")
    assert len(found) == 1 and "rch_swb: a recharge package left as the model's recharges 100 cells" in found[0]


def test_the_models_own_boundaries_are_booked_as_from_or_to_outside(runs):
    t, _, _ = runs["dis"]
    assert t["groundwater:boundary:wel_city"] == pytest.approx(-25.0 * 75)  # the other city well, the model's rate
    assert t["groundwater:boundary:rch_swb"] == 0.0                         # all of it SWAT+'s
    for p in ("ghb_regional", "chd_lake", "evt_riparian", "riv_creek"):     # riv_creek: the tributary
        assert t[f"groundwater:boundary:{p}"] != 0.0
    assert t["groundwater:boundary:evt_riparian"] < 0.0


def test_dis_and_disv_give_the_same_run(runs):
    (td, hd, _), (tv, hv, _) = runs["dis"], runs["disv"]
    assert np.allclose(hd, hv, rtol=0, atol=1e-6)
    for k in ("recharge:received", "exchange:cha01:to_stream", "groundwater:boundary:ghb_regional",
              "groundwater:boundary:riv_creek", "groundwater:storage_release"):
        assert tv[k] == pytest.approx(td[k], rel=1e-6, abs=1e-6), k


def with_a_farm_well(cfg):
    """A second plant drawing from the farm's multi-aquifer well (``maw=True``), screened through both layers."""
    cfg.mf6_builder_args = dict(cfg.mf6_builder_args, maw=True)
    cfg.wells.append(C.WellConfig("farm1", "maw_farm", None, well=1))
    plant = copy.deepcopy(cfg.users[0])
    plant.name, plant.withdrawal = "plant2", "well:farm1"
    cfg.users.append(plant)
    return cfg


def test_a_multi_aquifer_well_as_a_users_well(swatplus):
    """MAW: the user's withdrawal is the well's rate each day; MODFLOW 6 splits it between the layers, and what it
    pumped is delivered, from the aquifer by its connections and from the wellbore's water (to MAW's own solve)."""
    from basinledger.driver import Basin
    rows = []
    found, _ = Basin.check(with_a_farm_well(config()), run="existing-check-maw", takeover=rows)
    assert found == []
    assert {r.package: (r.role, r.entries) for r in rows}["maw_farm"] == ("Basin Ledger's", 2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(with_a_farm_well(config()), run="existing-run-maw", days=75)
        try:
            b.run()
            t = b.ledger.table()
        finally:
            b.finalize()
    k = "route:well:farm1:"
    assert t[k + "requested"].sum() > 0
    np.testing.assert_allclose(t[k + "delivered"], t[k + "requested"], rtol=1e-12)
    np.testing.assert_allclose(t[k + "from_aquifer"] + t[k + "from_wellbore"], t[k + "delivered"], rtol=1e-4)
    assert np.abs(t[k + "from_wellbore"]).max() < 1e-3 * t[k + "delivered"].max()
    assert np.abs(t["groundwater:closure_share"]).max() < 1e-6
    assert np.nansum(t["groundwater:boundary:maw_farm"]) == 0.0  # all of it the user's


def test_a_multi_aquifer_well_refused(swatplus, tmp_path):
    from basinledger.driver import Basin
    cfg = config()
    cfg.wells.append(C.WellConfig("farm1", "wel_city", None, well=1))
    found, _ = Basin.check(cfg, run="existing-check-maw-bad")
    assert len(found) == 1 and "well = 1 names a multi-aquifer well, but wel_city is a WEL package" in found[0]
    base = CONFIG.read_text().replace('"maps/', f'"{SYNTHETIC}/maps/')
    (tmp_path / "x.toml").write_text(base.replace("cell = [1, 5, 6]", "cell = [1, 5, 6]\nwell = 1"))
    with pytest.raises(ValueError, match="a MAW well is named by its number .well., not a cell"):
        C.load(tmp_path / "x.toml")
    (tmp_path / "x.toml").write_text(base.replace("cell = [1, 5, 6]", "well = 0"))
    with pytest.raises(ValueError, match="well is the MAW well's number in its package, from 1"):
        C.load(tmp_path / "x.toml")


def test_storage_for_a_model_handed_over_steady(tmp_path):
    """``existing.add_storage``: a model with no storage package gets one of Basin Ledger's own, transient; a model
    whose own is steady has its periods made transient, its values kept; the refusals."""
    import flopy

    from basinledger.coupling import existing as X
    ws = tmp_path / "no_sto"
    E.build(ws, 31, START, steady="no_sto")
    assert X.steady_periods(ws, "kzoo") == "no STO"
    with pytest.raises(ValueError, match="kzoo has no storage package; give sy and ss"):
        X.add_storage(ws, "kzoo", None, None)
    with pytest.raises(ValueError, match="sy is a fraction"):
        X.add_storage(ws, "kzoo", 1.5, 1e-5)
    notes = X.add_storage(ws, "kzoo", 0.15, 1e-5)
    assert notes == ["sto_basinledger: a storage package Basin Ledger adds to the run's copy, the model having none: "
                     "ICONVERT 1, SY 0.15, SS 1e-05 per meter, transient from the first day"]
    assert X.steady_periods(ws, "kzoo") is None
    sto = flopy.mf6.MFSimulation.load(sim_ws=str(ws), verbosity_level=0).get_model("kzoo").get_package("sto")
    assert float(sto.sy.array.max()) == 0.15 and sto.transient.get_data() is True  # flopy reads it back
    ws = tmp_path / "feet"
    E.build(ws, 31, START, steady="no_sto", length_units="feet")
    X.add_storage(ws, "kzoo", 0.15, 1e-5)
    assert "CONSTANT 3.048e-06" in (ws / "kzoo.sto_basinledger.sto").read_text()  # per foot
    ws = tmp_path / "steady"
    E.build(ws, 31, START, steady="sto")
    assert X.steady_periods(ws, "kzoo") == "steady"
    with pytest.raises(ValueError, match="has its own storage package"):
        X.add_storage(ws, "kzoo", 0.15, 1e-5)
    assert X.add_storage(ws, "kzoo", None, None) == [
        "sto_glacial: the model's storage package made transient in the run's copy (1 steady period), its "
        "own SY and SS kept"]
    assert X.steady_periods(ws, "kzoo") is None


@pytest.mark.parametrize("storage", [None, {"sy": 0.2, "ss": 1e-5}])
def test_a_model_handed_over_steady_coupled(swatplus, storage):
    """With no storage package, each day is the day's steady state (noted) and nothing is released from storage;
    with [modflow6] storage, the model runs transient, storage releasing and taking water, and the ledger closes."""
    from basinledger.driver import Basin
    cfg = config(steady="no_sto")
    cfg.mf6_storage = storage
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        b = Basin(cfg, run=f"existing-run-steady-{bool(storage)}", days=75)
        try:
            b.run()
            t = b.ledger.table()
        finally:
            b.finalize()
    notes = [str(w.message) for w in caught]
    sto = t["groundwater:storage_release"]
    if storage is None:
        assert any(n.startswith("kzoo: no storage package, so MODFLOW 6 solves those days without storage") for n in notes)
        assert (sto == 0).all()
    else:
        assert any(n.startswith("sto_basinledger: a storage package Basin Ledger adds") for n in notes)
        assert (np.abs(sto) > 0).all()  # storage in every day's balance (heads falling from their start: released)
    assert np.abs(t["groundwater:closure_share"]).max() < 1e-6
