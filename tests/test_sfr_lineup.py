"""The line-up rule for an SFR package, the same for a package Basin Ledger built from SWAT+'s channel lines and one a
study brought (coupling/exchange.py, ``sfr_problems``; docs/design.md, "Extending an existing MODFLOW 6 model", A
and B). Tier 1, MODFLOW 6 alone.

The network, built by hand the way someone else's model would have it: a main stem of ten reaches down the west
column of the small synthetic aquifer, and a tributary of three reaches along row 4 flowing west into reach 5.

- lined up (SWAT+ built on the network, an outlet at the confluence): the main stem as two channels, the tributary a
  third; connections between channels run from a last reach into a top reach, and are cut each day, so each channel
  carries only its own water and SWAT+ carries it on;
- a confluence inside a channel, from another channel or from the model's own reaches, is refused by name;
- the model's own tributary flowing into a channel's top reach is cut, and keeps its own inflow;
- a package built from channel lines has nothing to cut (the build-both path, unchanged);
- the confluences listed with their places, for building SWAT+'s outlets (gridmaps.sfr_confluences).
"""

from __future__ import annotations

import datetime as dt
import importlib.util

import flopy
import numpy as np
import pytest

from basinledger.coupling.exchange import ChannelMap, StreamExchange, sfr_reach_links
from basinledger.coupling.modflow6 import Modflow6
from conftest import ROOT, RUNS

SYNTHETIC = ROOT / "examples" / "synthetic"
MAIN = list(range(10))
TRIBUTARY = [10, 11, 12]
INTO = 5  # the main stem's reach the tributary joins


def synthetic_model():
    spec = importlib.util.spec_from_file_location("synthetic_model", SYNTHETIC / "model.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


M = synthetic_model()


def build(name: str):
    ws = RUNS / name
    M.build_small(ws, 3, dt.date(2010, 1, 1))
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(ws), verbosity_level=0)
    gwf = sim.get_model("gwf")
    gwf.remove_package("riv")
    cells = [(0, r, 0) for r in MAIN] + [(0, 4, 3), (0, 4, 2), (0, 4, 1)]
    down = {r: r + 1 for r in MAIN[:-1]} | {10: 11, 11: 12, 12: INTO}
    up: dict[int, list[int]] = {}
    for u, d in down.items():
        up.setdefault(d, []).append(u)
    conn = [[r] + up.get(r, []) + ([-down[r]] if r in down else []) for r in range(13)]
    pd = [(r, cells[r], M.SMALL_CELL, M.SFR_WIDTH, M.SFR_SLOPE, 88.0 - M.SFR_SLOPE * M.SMALL_CELL * (r % 10), 1.0,
           M.SFR_BED_K, 0.035, len(conn[r]) - 1, 1.0, 0) for r in range(13)]
    flopy.mf6.ModflowGwfsfr(gwf, pname="sfr_river", nreaches=13, packagedata=pd, connectiondata=conn,
                            perioddata={0: [(10, "inflow", 77.0)]}, unit_conversion=86400.0,
                            filename="gwf.sfr_river.sfr")
    sim.write_simulation(silent=True)
    return ws


@pytest.fixture
def model(libmf6, request):
    mf = Modflow6(build(f"sfr-lineup-{request.node.name}"))
    mf.prepare_day()
    yield mf
    mf.finalize()


def channels(**spec) -> ChannelMap:
    """ChannelMap from channel -> reaches; reaches not named have no channel."""
    e = np.full(13, -1)
    for i, reaches in enumerate(spec.values()):
        e[reaches] = i
    return ChannelMap(list(spec), e)


def test_the_network_as_built(model):
    links = sfr_reach_links(model, "sfr_river")
    assert (12, INTO) in links and (4, INTO) in links and len(links) == 12


def test_lined_up_each_channel_carries_its_own_water(model):
    ex = StreamExchange(model, "sfr_river", channels(up=MAIN[:5], trib=TRIBUTARY, down=MAIN[5:]))
    assert ex._top == {"up": 0, "trib": 10, "down": INTO} and ex._cut == [INTO]
    ex.set_inflow({"up": 1000.0, "trib": 200.0, "down": 500.0})
    model.solve_day()
    q = -model.get("QOUTFLOW", model.model, "SFR_RIVER")
    leak = model.simvals("sfr_river")  # + into the aquifer
    assert model.get("USTRF", model.model, "SFR_RIVER")[INTO] == 0.0
    assert q[INTO] == pytest.approx(500.0 - leak[INTO], rel=1e-9)        # only down's own water
    assert q[4] == pytest.approx(1000.0 - leak[:5].sum(), rel=1e-9)       # up's water leaves SFR at its last reach
    assert q[12] == pytest.approx(200.0 - leak[10:13].sum(), rel=1e-9)    # and the tributary's at its
    net = ex.after_solve()
    assert net["down"] == pytest.approx(-leak[5:10].sum(), rel=1e-12)


def test_a_confluence_inside_a_channel_is_refused(model):
    with pytest.raises(ValueError, match="reach 12 \\(trib\\) flows into reach 5 \\(main\\), inside channel main; "
                                         "every confluence must fall on a channel boundary"):
        StreamExchange(model, "sfr_river", channels(main=MAIN, trib=TRIBUTARY))


def test_the_models_own_tributary_into_a_top_reach_is_cut_and_keeps_its_inflow(model):
    ex = StreamExchange(model, "sfr_river", channels(up=MAIN[:5], down=MAIN[5:]))
    assert ex._cut == [INTO]
    ex.set_inflow({"up": 1000.0, "down": 500.0})
    assert model.get("INFLOW", model.model, "SFR_RIVER")[10] == 77.0  # the model's own
    model.solve_day()
    q = -model.get("QOUTFLOW", model.model, "SFR_RIVER")
    leak = model.simvals("sfr_river")
    assert q[INTO] == pytest.approx(500.0 - leak[INTO], rel=1e-9)


def test_the_models_own_tributary_into_the_middle_of_a_channel_is_refused(model):
    with pytest.raises(ValueError, match="reach 12 \\(no channel\\) flows into reach 5 \\(main\\), inside channel"):
        StreamExchange(model, "sfr_river", channels(main=MAIN))


def test_a_channel_that_isnt_one_chain_is_refused(model):
    with pytest.raises(ValueError, match="channel all has 2 top reaches \\(\\[0, 10\\]\\)"):
        StreamExchange(model, "sfr_river", channels(all=MAIN + TRIBUTARY))


def test_a_package_built_from_channel_lines_has_nothing_to_cut(libmf6):
    ws = RUNS / "sfr-lineup-built"
    M.build_small_sfr(ws, 3, dt.date(2010, 1, 1))
    mf = Modflow6(ws)
    try:
        ex = StreamExchange(mf, "sfr", ChannelMap(["cha01"], np.zeros(10, int)))
        assert ex._cut == []
    finally:
        mf.finalize()


def test_the_confluences_for_swat_outlets(libmf6):
    pytest.importorskip("shapely")
    from basinledger import gridmaps as G
    ws = build("sfr-lineup-confluences")
    _, gwf = G.load(ws)
    links = G.sfr_links(gwf, "sfr_river")
    mf = Modflow6(ws)
    try:
        mf.prepare_day()
        assert links == sorted(sfr_reach_links(mf, "sfr_river"))  # the file read as libmf6 reads it
    finally:
        mf.finalize()
    (c,) = G.sfr_confluences(gwf, "sfr_river")
    assert c.reach == INTO and c.upstream == [4, 12] and c.cell == (0, 5, 0)
    assert c.xy == pytest.approx((M.SMALL_CELL / 2, M.SMALL_CELL * 4.5))
