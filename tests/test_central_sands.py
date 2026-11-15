"""The second way in on a published model: Pleasant Lake, Wisconsin Central Sands (Fienen, Haserodt and Leaf, 2021,
doi:10.5066/P9BVFSGJ; examples/central-sands). A parent model and an LGR inset in one simulation, SFR built by
SFRmaker in both, a lake in the inset, the mover carrying the inset's streams into the parent's, a steady first day
and then monthly periods of five growing steps. Skipped until it is fetched (``uv run python
examples/central-sands/fetch.py``, 16 MB).

What Basin Ledger does with it today: reads it as given (meters, days, a georeferenced structured grid), lists both
models' packages, puts its SFR networks through the line-up rule with SFRmaker's lines as SWAT+'s channels (they
line up: SWAT+ built on the same lines would route every channel), lists their confluences, reads the mover from the
exchange's file and cuts it where it carries the inset's streams into the parent's channels, and refuses what isn't
built yet, by name: the inset's coupled streams given to the parent through the mover, and (on the model as given,
before a run rewrites them) time steps longer than a day. Tier 1: MODFLOW 6 reads it and solves its first, steady
day; nothing longer runs.
"""

from __future__ import annotations

import importlib.util

import numpy as np
import pytest

from basinledger.coupling import mover as MV
from basinledger.coupling.exchange import ChannelMap, StreamExchange, sfr_problems
from basinledger.coupling.existing import inventory
from basinledger.coupling.modflow6 import Modflow6
from conftest import ROOT

spec = importlib.util.spec_from_file_location("central_sands_fetch", ROOT / "examples" / "central-sands" / "fetch.py")
F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)
WS = F.workspace()
pytestmark = pytest.mark.skipif(not (WS / "mfsim.nam").exists(),
                                reason="Pleasant Lake not fetched: uv run python examples/central-sands/fetch.py")
INSET, PARENT = "plsnt_lgr_inset", "plsnt_lgr_parent"


@pytest.fixture(scope="module")
def mf(libmf6):
    m = Modflow6(WS, INSET)
    m.prepare_day()  # the steady first day: one step of one day
    yield m
    m.finalize()


@pytest.fixture(scope="module")
def sim():
    pytest.importorskip("shapely")
    from basinledger import gridmaps as G
    return G.load(WS, INSET)[0]


def lines(sim, model) -> ChannelMap:
    """SFRmaker's lines as SWAT+'s channels: each reach to the line it was cut from."""
    from basinledger import gridmaps as G
    ids = [str(int(x)) for x in G.sfr_aux(sim.get_model(model), "sfr_0", "line_id")]
    names = list(dict.fromkeys(ids))
    return ChannelMap(names, np.array([names.index(i) for i in ids]))


def test_read_as_given(mf):
    assert mf.grid == "dis" and mf.shape == (5, 100, 120)  # meters and days, or Modflow6 would have refused it


def test_both_models_packages_listed(mf):
    inset = {r.package: (r.type, r.entries) for r in inventory(mf)}
    parent = {r.package: (r.type, r.entries) for r in inventory(mf, PARENT)}
    assert inset == {"lak_0": ("lak", 1382), "rcha": ("rch", 12000), "sfr_0": ("sfr", 25)}
    assert parent == {"chd_0": ("chd", 1599), "rcha": ("rch", 10000), "sfr_0": ("sfr", 551), "wel_0": ("wel", 244)}
    assert all(r.role == "the model's" for r in inventory(mf))


def test_its_streams_line_up_with_their_own_lines(mf, sim):
    out = {}
    for model, reaches, n_lines, cuts in ((INSET, 25, 2, 0), (PARENT, 551, 48, 21)):
        mf.model = model.upper()
        try:
            cmap = lines(sim, model)
            out[model] = lay = sfr_problems(mf, "sfr_0", cmap)
        finally:
            mf.model = INSET.upper()
        assert len(cmap.entry_channel) == reaches and len(cmap.channels) == n_lines
        assert lay.channel_problems == {} and len(lay.top) == n_lines and len(lay.cut) == cuts
    # the mover carries the inset's two outlet reaches into two of the parent's top reaches: lined up on the parent's
    # side (cut each day), and on the inset's the one thing in the way, until both models are coupled
    parent = out[PARENT]
    assert parent.problems == []
    assert [(m.provider.id, m.receiver.id) for m in parent.moves_cut] == [(24, 218), (4, 213)]
    assert {218, 213} <= set(parent.top.values())
    assert [p.split(" through")[0] for p in out[INSET].problems] == [
        "sfr_0: reach 24 (104) gives water to plsnt_lgr_parent's sfr_0 entry 218",
        "sfr_0: reach 4 (102) gives water to plsnt_lgr_parent's sfr_0 entry 213"]


def test_the_mover_is_refused_by_name(mf, sim):
    with pytest.raises(ValueError, match="gives water to plsnt_lgr_parent's sfr_0 entry 218 through the mover "
                                         "\\(GWF-GWF_1, period 1's move 1\\); SWAT\\+ routes the channel's water"):
        StreamExchange(mf, "sfr_0", lines(sim, INSET), sfr="route_or_stage")


def renamed(cmap: ChannelMap, model: str) -> ChannelMap:
    """The lines as channels of their own model, as if SWAT+ had an outlet where the inset meets its parent."""
    return ChannelMap([f"{model}:{c}" for c in cmap.channels], cmap.entry_channel)


def test_both_models_streams_line_up_with_the_mover_between_them(mf, sim):
    """Both models coupled: each model's SFR to its own lines, and the mover from the inset's two outlet reaches into
    the parent's top reaches lined up on both sides and cut. Its lines as SFRmaker numbered them, the same line on
    both sides of the inset's edge is one channel, and the move a move within a channel, refused."""
    parent = mf.view(PARENT)
    inset_map, parent_map = renamed(lines(sim, INSET), "inset"), renamed(lines(sim, PARENT), "parent")
    li = sfr_problems(mf, "sfr_0", inset_map, {(PARENT, "sfr_0"): parent_map})
    lp = sfr_problems(parent, "sfr_0", parent_map, {(INSET, "sfr_0"): inset_map})
    assert li.problems == [] and li.channel_problems == {} and lp.problems == [] and lp.channel_problems == {}
    assert [(m.provider.id, m.receiver.id) for m in li.moves_cut] == [(24, 218), (4, 213)]
    assert li.moves_cut == lp.moves_cut
    assert {24, 4} <= set(li.last.values()) and {218, 213} <= set(lp.top.values())
    same = sfr_problems(mf, "sfr_0", lines(sim, INSET), {(PARENT, "sfr_0"): lines(sim, PARENT)})
    assert any("a mover within a channel isn't built" in p for p in same.problems)


def test_confluences_listed_for_swat_outlets(sim):
    from basinledger import gridmaps as G
    assert G.sfr_confluences(sim.get_model(INSET), "sfr_0") == []
    c = G.sfr_confluences(sim.get_model(PARENT), "sfr_0")
    assert len(c) == 10 and all(len(x.upstream) >= 2 and x.xy is not None for x in c)


def test_the_grid_is_georeferenced(sim):
    from basinledger import gridmaps as G
    grid = sim.get_model(INSET).modelgrid
    assert (grid.xoffset, grid.yoffset) == (554700.0, 389300.0)  # Wisconsin Transverse Mercator
    assert (G.top_active(grid) >= 0).sum() == 12000
    assert len(G.stream_entries(sim.get_model(INSET), "sfr_0")) == 25


def test_the_mover_cut_at_run_time_then_monthly_steps_refused(mf):
    """Last: it moves the shared model on (libmf6 holds one simulation per process). The mover here is the GWF-GWF
    exchange's, carrying the inset's two outlet reaches into the parent's SFR; cut as the parent's line-up has it
    (each move's factor set to 0), that water leaves the inset's SFR as its external outflow, as a cut upstream
    fraction does."""
    (mv,) = MV.movers(WS)
    assert mv.path == "GWF-GWF_1/MVR" and mv.file.name == "pleasant_lgr.mvr"
    assert mf.api.get_value("GWF-GWF_1/MVR/VALUE").tolist() == [1.0, 1.0]  # factors
    MV.cut(mf, [mv], mv.active(1))
    assert mf.api.get_value("GWF-GWF_1/MVR/VALUE").tolist() == [0.0, 0.0]
    mf.solve_day()
    ext = mf.get("QEXTOUTFLOW", mf.model, "SFR_0")
    assert ext[4] < -100.0 and ext[24] < -100.0  # m3 d-1 leaving at the two outlets
    with pytest.raises(ValueError, match="time step is 4.16"):
        mf.prepare_day()
