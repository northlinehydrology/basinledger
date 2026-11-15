"""The mover in the SFR line-up rule (coupling/mover.py, coupling/exchange.py ``sfr_problems``; docs/design.md,
"Extending an existing MODFLOW 6 model"). Tier 1, MODFLOW 6 alone.

The line-up test's network (tests/test_sfr_lineup.py), but the tributary joins the main stem through the mover, as a
model with an LGR inset or a diversion's return has it, rather than through SFR's connections:

- the movers read from their files (a model's own, and an exchange's on Pleasant Lake), each period's block in
  effect until the next;
- lined up, a move between channels, or from the model's own reaches into a channel's top reach, is cut each day
  (its value set to move nothing, for each of the four types): the tributary's water leaves SFR at its last reach,
  and the main stem carries only its own;
- refused by name: a move into the middle of a channel (A; B runs it, cut), from a coupled reach to a reach no
  channel has or to another package, and within a channel; a drain's move into a channel is cut and its water the
  channel's, and its move into the model's own reaches stays the model's.
"""

from __future__ import annotations

import datetime as dt
import importlib.util

import flopy
import numpy as np
import pytest
from conftest import ROOT, RUNS

from basinledger.coupling import mover as MV
from basinledger.coupling.exchange import ChannelMap, StreamExchange
from basinledger.coupling.modflow6 import Modflow6

spec = importlib.util.spec_from_file_location("synthetic_model", ROOT / "examples" / "synthetic" / "model.py")
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)
MAIN = list(range(10))
TRIBUTARY = [10, 11, 12]
INTO = 5
P = "SFR_RIVER"


def build(name: str, moves: dict[int, list], drain: bool = False):
    """The network with the tributary's last reach not connected to the main stem; ``moves``: the mover's period
    blocks (0-based period -> (provider, id, receiver, id, type, value) rows, ids 0-based). Three one-day periods."""
    ws = RUNS / name
    M.build_small(ws, 3, dt.date(2010, 1, 1))
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(ws), verbosity_level=0)
    sim.tdis.nper = 3
    sim.tdis.perioddata = [(1.0, 1, 1.0)] * 3
    gwf = sim.get_model("gwf")
    gwf.remove_package("riv")
    cells = [(0, r, 0) for r in MAIN] + [(0, 4, 3), (0, 4, 2), (0, 4, 1)]
    down = {r: r + 1 for r in MAIN[:-1]} | {10: 11, 11: 12}
    up: dict[int, list[int]] = {}
    for u, d in down.items():
        up.setdefault(d, []).append(u)
    conn = [[r] + up.get(r, []) + ([-down[r]] if r in down else []) for r in range(13)]
    pd = [(r, cells[r], M.SMALL_CELL, M.SFR_WIDTH, M.SFR_SLOPE, 88.0 - M.SFR_SLOPE * M.SMALL_CELL * (r % 10), 1.0,
           M.SFR_BED_K, 0.035, len(conn[r]) - 1, 1.0, 0) for r in range(13)]
    flopy.mf6.ModflowGwfsfr(gwf, pname="sfr_river", nreaches=13, packagedata=pd, connectiondata=conn,
                            perioddata={0: [(10, "inflow", 77.0)]}, unit_conversion=86400.0, mover=True,
                            filename="gwf.sfr_river.sfr")
    packages = [["sfr_river"]]
    if drain:  # tile drains at a cell above the water table's start
        flopy.mf6.ModflowGwfdrn(gwf, pname="drn_tile", stress_period_data={0: [((0, 6, 6), 85.0, 50.0)]},
                                mover=True, filename="gwf.drn_tile.drn")
        packages.append(["drn_tile"])
    flopy.mf6.ModflowGwfmvr(gwf, maxmvr=max(len(v) for v in moves.values()), maxpackages=len(packages),
                            packages=packages, perioddata=moves, filename="gwf.mvr")
    sim.write_simulation(silent=True)
    return ws


def channels(**spec) -> ChannelMap:
    e = np.full(13, -1)
    for i, reaches in enumerate(spec.values()):
        e[reaches] = i
    return ChannelMap(list(spec), e)


TRIB_INTO_MAIN = ("sfr_river", 12, "sfr_river", INTO, "FACTOR", 1.0)


@pytest.fixture
def run(libmf6, request):
    """A model built with the given moves, MODFLOW 6 started and its first day prepared."""
    made = []

    def start(moves, drain=False):
        mf = Modflow6(build(f"mover-{request.node.name}", moves, drain))
        made.append(mf)
        mf.prepare_day()
        return mf

    yield start
    for mf in made:
        mf.finalize()


def test_the_movers_read_from_their_files():
    ws = build("mover-read", {0: [TRIB_INTO_MAIN], 2: [("sfr_river", 12, "sfr_river", 6, "UPTO", 50.0)]})
    (mv,) = MV.movers(ws)
    assert mv.path == "GWF/MVR" and mv.file.name == "gwf.mvr" and sorted(mv.blocks) == [1, 3]
    (m,) = mv.active(1)
    assert (m.provider, m.receiver, m.kind, m.value) == (MV.Entry("gwf", "sfr_river", 12),
                                                        MV.Entry("gwf", "sfr_river", INTO), "FACTOR", 1.0)
    assert mv.active(2) == [m] and mv.active(3)[0].receiver.id == 6  # a block holds until the next


def test_the_mover_works_as_given(run):
    mf = run({0: [TRIB_INTO_MAIN]})
    mf.solve_day()
    assert mf.get("QFROMMVR", mf.model, P)[INTO] > 77.0  # the tributary's 77 m3 d-1 and what it gains


def test_lined_up_a_move_between_channels_is_cut(run):
    mf = run({0: [TRIB_INTO_MAIN]})
    ex = StreamExchange(mf, "sfr_river", channels(up=MAIN[:5], trib=TRIBUTARY, down=MAIN[5:]))
    assert ex._cut == [INTO] and [(m.provider.id, m.receiver.id) for m in ex._moves_cut] == [(12, INTO)]
    ex.set_inflow({"up": 1000.0, "trib": 200.0, "down": 500.0})
    assert mf.api.get_value("GWF/MVR/VALUE").tolist() == [0.0]
    mf.solve_day()
    q = -mf.get("QOUTFLOW", mf.model, P)
    leak = mf.simvals("sfr_river")
    assert mf.get("QFROMMVR", mf.model, P)[INTO] == 0.0
    assert q[INTO] == pytest.approx(500.0 - leak[INTO], rel=1e-9)        # only down's own water
    ext = -mf.get("QEXTOUTFLOW", mf.model, P)
    assert ext[12] == pytest.approx(200.0 - leak[10:13].sum(), rel=1e-9)  # the tributary's leaves SFR at its end


def test_the_models_own_tributary_by_the_mover_into_a_top_reach_is_cut(run):
    mf = run({0: [TRIB_INTO_MAIN]})
    ex = StreamExchange(mf, "sfr_river", channels(up=MAIN[:5], down=MAIN[5:]))
    ex.set_inflow({"up": 1000.0, "down": 500.0})
    mf.solve_day()
    assert mf.get("QFROMMVR", mf.model, P)[INTO] == 0.0
    assert mf.get("INFLOW", mf.model, P)[10] == 77.0  # the model's own tributary keeps its own inflow


@pytest.mark.parametrize("kind, value, nothing", [("EXCESS", 10.0, 1e30), ("THRESHOLD", 10.0, 1e30),
                                                  ("UPTO", 50.0, 0.0), ("FACTOR", 0.5, 0.0)])
def test_each_type_cut_each_day_in_each_period(run, kind, value, nothing):
    """Period 1 joins the tributary at the main stem's top reach by a factor; period 3's block by this type."""
    mf = run({0: [TRIB_INTO_MAIN], 2: [("sfr_river", 12, "sfr_river", INTO, kind, value), TRIB_INTO_MAIN]})
    ex = StreamExchange(mf, "sfr_river", channels(up=MAIN[:5], trib=TRIBUTARY, down=MAIN[5:]))
    assert len(ex._moves_cut) == 3
    for day in range(3):
        if day:
            mf.prepare_day()
        ex.set_inflow({"up": 1000.0, "trib": 200.0, "down": 500.0})
        n = int(mf.api.get_value("GWF/MVR/NMVR")[0])  # VALUE is as long as the most moves in any period
        assert mf.api.get_value("GWF/MVR/VALUE")[:n].tolist() == ([0.0] if day < 2 else [nothing, 0.0])
        mf.solve_day()
        assert mf.get("QFROMMVR", mf.model, P)[INTO] == 0.0, day


def test_a_move_into_the_middle_of_a_channel(run):
    mf = run({0: [("sfr_river", 12, "sfr_river", 7, "FACTOR", 1.0)]})
    with pytest.raises(ValueError, match="reach 12 \\(trib\\) gives water through the mover to reach 7 \\(down\\), "
                                         "inside channel down; every confluence must fall on a channel boundary"):
        StreamExchange(mf, "sfr_river", channels(up=MAIN[:5], trib=TRIBUTARY, down=MAIN[5:]))
    ex = StreamExchange(mf, "sfr_river", channels(up=MAIN[:5], trib=TRIBUTARY, down=MAIN[5:]), sfr="route_or_stage")
    assert ex.modes["down"][0] == "stage" and ex.modes["trib"] == ("route", [])  # B, the move cut
    ex.set_inflow({"up": 1000.0, "trib": 200.0, "down": 500.0}, depth=dict.fromkeys(["up", "trib", "down"], 0.2))
    mf.solve_day()
    assert mf.get("QFROMMVR", mf.model, P)[7] == 0.0


def test_moves_refused_by_name(run):
    mf = run({0: [TRIB_INTO_MAIN, ("sfr_river", 10, "sfr_river", 12, "FACTOR", 0.5),
                  ("drn_tile", 0, "sfr_river", 7, "FACTOR", 1.0)]}, drain=True)
    with pytest.raises(ValueError) as e:
        StreamExchange(mf, "sfr_river", channels(trib=TRIBUTARY, down=MAIN[7:]), sfr="route_or_stage")
    text = str(e.value)
    assert "reach 12 (trib) gives water to reach 5 (no channel) through the mover (GWF, period 1's move 1)" in text
    assert "reach 10 gives water to reach 12 in channel trib through the mover" in text
    assert "a mover within a channel isn't built" in text
    assert "drn_tile" not in text  # a drain's move into a channel is cut, its water the channel's (below)


def test_a_drains_move_into_a_channel_is_cut_and_its_water_the_channels(run):
    """A drain giving water through the mover to a coupled reach: the move is cut each day (SWAT+ carries the
    channel's water), and the drain's water is the channel's, as the exchange's drains' seepage is."""
    mf = run({0: [TRIB_INTO_MAIN, ("drn_tile", 0, "sfr_river", 7, "FACTOR", 1.0)]}, drain=True)
    ex = StreamExchange(mf, "sfr_river", channels(up=MAIN[:5], trib=TRIBUTARY, down=MAIN[5:]))
    assert {(m.provider.package, m.receiver.id) for m in ex._moves_cut} == {("sfr_river", INTO), ("drn_tile", 7)}
    assert ex._drains_in == {"down": [("GWF", "DRN_TILE", 0)]}
    ex.set_inflow({"up": 1000.0, "trib": 200.0, "down": 500.0})
    mf.solve_day()
    assert mf.get("QFROMMVR", mf.model, P)[7] == 0.0  # nothing moved into reach 7
    q = -float(mf.simvals("drn_tile")[0])
    due = ex.after_solve()
    assert q > 0.0 and ex.drains_moved == {"up": 0.0, "trib": 0.0, "down": pytest.approx(q, rel=1e-12)}
    assert due["down"] - q == pytest.approx(-float(mf.simvals("sfr_river")[MAIN[5:]].sum()), rel=1e-9)


def test_a_drains_move_into_the_models_own_reaches_stays_the_models(run):
    mf = run({0: [TRIB_INTO_MAIN, ("drn_tile", 0, "sfr_river", 11, "FACTOR", 1.0)]}, drain=True)
    ex = StreamExchange(mf, "sfr_river", channels(up=MAIN[:5], down=MAIN[5:]))
    assert [(m.provider.package, m.receiver.id) for m in ex._moves_cut] == [("sfr_river", INTO)]
    ex.set_inflow({"up": 1000.0, "down": 500.0})
    assert mf.api.get_value("GWF/MVR/VALUE").tolist() == [0.0, 1.0]
