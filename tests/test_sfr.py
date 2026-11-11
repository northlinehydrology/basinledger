"""The two-way exchange through MODFLOW 6's streamflow routing package (coupling/exchange.py, SFR).

- MODFLOW 6 alone: a channel's water into its top reach, and a reach losing all it carries and no more; map rows by
  reach number; what stops an SFR package from being coupled (tier 1);
- the coupled run (examples/synthetic/sfr.toml, three years): MODFLOW 6's budget closes every day and its SFR term
  is the ledger's exchange, the exchange reaches SWAT+ the next day, no reach loses more than its channel put in, and
  nothing comes from outside the models (tier 1 and 2). Restarts: tests/test_restart.py.
"""

from __future__ import annotations

import datetime as dt
import importlib.util

import flopy
import numpy as np
import pytest

from basinledger import config as C
from basinledger.coupling.exchange import ChannelMap, StreamExchange, sfr_reach_links
from basinledger.coupling.modflow6 import Modflow6
from basinledger.driver import Basin
from conftest import ROOT, RUNS

SYNTHETIC = ROOT / "examples" / "synthetic"
DAYS = 1096


def synthetic_model():
    spec = importlib.util.spec_from_file_location("synthetic_model", SYNTHETIC / "model.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def sfr_model(name, strt=None, **kw):
    ws = RUNS / name
    synthetic_model().build_small_sfr(ws, 3, dt.date(2010, 1, 1), **kw)
    if strt is not None:
        sim = flopy.mf6.MFSimulation.load(sim_ws=ws, verbosity_level=0)
        sim.get_model("gwf").ic.strt = strt
        sim.write_simulation(silent=True)
    return Modflow6(ws)


# ---------- MODFLOW 6 alone ----------

def test_the_reaches_are_one_chain(libmf6):
    mf = sfr_model("sfr-chain")
    try:
        assert sfr_reach_links(mf, "sfr") == [(r, r + 1) for r in range(9)]
        ex = StreamExchange(mf, "sfr", ChannelMap(["cha01"], np.zeros(10, int)))
        assert ex.kind == "sfr" and ex._top == {"cha01": 0}
    finally:
        mf.finalize()


def test_a_reach_loses_what_it_carries_and_no_more(libmf6):
    """With the water table well below the streambed, the channel's water all leaks: the first reaches lose it and
    the ones below are dry. With more water than the bed can take, the loss is the bed's, less than the inflow."""
    mf = sfr_model("sfr-dry", strt=80.0)
    try:
        ex = StreamExchange(mf, "sfr", ChannelMap(["cha01"], np.zeros(10, int)))
        for inflow in (50.0, 1e5):
            mf.prepare_day()
            ex.set_inflow({"cha01": inflow})
            mf.solve_day()
            got = ex.after_solve()["cha01"]  # m3, positive to the stream
            depth = mf.get("DEPTH", mf.model, "SFR")
            if inflow == 50.0:
                assert got == pytest.approx(-50.0, rel=1e-9)
                assert depth[-1] == 0.0
            else:
                assert -inflow < got < 0 and depth[-1] > 0
    finally:
        mf.finalize()


def test_map_rows_by_reach(libmf6, tmp_path):
    mf = sfr_model("sfr-map")
    try:
        mf.prepare_day()  # a package's entries are there once MODFLOW 6 has read its first day
        (tmp_path / "m.csv").write_text("object,reach\n" + "".join(f"cha01,{r}\n" for r in range(10)))
        cmap = ChannelMap.from_spec("m.csv", tmp_path, mf, "sfr", ["cha01"])
        assert cmap.entry_channel.tolist() == [0] * 10
        (tmp_path / "bad.csv").write_text("object,reach\ncha01,10\n")
        with pytest.raises(ValueError, match="reach 10 isn't in 'sfr' \\(10 reaches"):
            ChannelMap.from_spec("bad.csv", tmp_path, mf, "sfr", ["cha01"])
        (tmp_path / "w.csv").write_text("object,reach\ncha01,0\n")
        with pytest.raises(ValueError, match="reach numbers are for an SFR package; 'wel' isn't one"):
            ChannelMap.from_spec("w.csv", tmp_path, mf, "wel", ["cha01"])
    finally:
        mf.finalize()


@pytest.mark.parametrize("kw, refused", [
    (dict(unit_conversion=None), "Manning's constant is 1, which isn't for meters and days"),
    (dict(storage=True), "STORAGE"),
    (dict(split=5), "channel cha01 has 2 top reaches \\(\\[0, 5\\]\\)"),
])
def test_what_stops_an_sfr_package(libmf6, kw, refused):
    mf = sfr_model("sfr-refused", **kw)
    try:
        with pytest.raises(ValueError, match=refused):
            StreamExchange(mf, "sfr", ChannelMap(["cha01"], np.zeros(10, int)))
    finally:
        mf.finalize()


def test_time_conversion_sets_mannings_constant(libmf6):
    mf = sfr_model("sfr-timeconv", unit_conversion=None, time_conversion=86400.0)
    try:
        StreamExchange(mf, "sfr", ChannelMap(["cha01"], np.zeros(10, int)))
    finally:
        mf.finalize()


def test_a_reach_flowing_out_of_its_channel_is_refused(libmf6):
    """Reaches 0-4 coupled, 5-9 not: reach 4 flows into reach 5, so SFR would route the channel's water on."""
    mf = sfr_model("sfr-cross")
    try:
        e = np.array([0] * 5 + [-1] * 5)
        with pytest.raises(ValueError, match="reach 4 \\(cha01\\) flows into reach 5 \\(no channel\\)"):
            StreamExchange(mf, "sfr", ChannelMap(["cha01"], e))
    finally:
        mf.finalize()


# ---------- settings refused in a coupled run ----------

def config(**exchange):
    cfg = C.load(SYNTHETIC / "sfr.toml")
    cfg.run_folder = RUNS
    cfg.exchange = dict(cfg.exchange, **exchange)
    return cfg


def test_limit_off_is_refused_with_sfr(swatplus):
    b = Basin(config(limit=False), run="sfr-limit-off", days=2)
    try:
        with pytest.raises(ValueError, match="SFR limits a reach's loss"):
            b.step()
    finally:
        b.finalize()


def test_river_stages_are_refused_on_sfr(swatplus):
    cfg = config()
    cfg.rivers = {"package": "sfr", "map": "all:cha01", "stage": "surface_depth"}
    b = Basin(cfg, run="sfr-rivers", days=2)
    try:
        with pytest.raises(ValueError, match="\\[rivers\\] sets a river \\(RIV\\) package's stage; 'sfr' is SFR"):
            b.step()
    finally:
        b.finalize()


# ---------- the coupled run ----------

@pytest.fixture(scope="module")
def sfr(swatplus):
    b = Basin(config(), run="sfr", days=DAYS)
    b.run()
    out = dict(table=b.ledger.table(), due=dict(b.due), ws=b.mf6.ws)
    b.finalize()
    return out


@pytest.mark.slow
def test_budget_closes_and_its_sfr_term_is_the_exchange(sfr):
    cbc = flopy.utils.CellBudgetFile(sfr["ws"] / "gwf.cbc", precision="double")
    terms = {t: np.array([float(d["q"].sum()) if d.dtype.names else float(np.sum(d)) for d in cbc.get_data(text=t)])
             for t in ("STO-SS", "STO-SY", "RCHA", "WEL", "SFR")}
    net = sum(terms.values())
    gross = sum(np.abs(v) for v in terms.values())
    assert len(net) == DAYS and np.max(np.abs(net) / gross) < 1e-7
    np.testing.assert_allclose(-terms["SFR"], sfr["table"]["exchange:cha01:to_stream"], rtol=1e-12, atol=1e-9)


@pytest.mark.slow
def test_no_reach_loses_more_than_its_channel_put_in(sfr):
    t = sfr["table"]
    x, inflow = t["exchange:cha01:to_stream"], t["exchange:cha01:sfr_inflow"]
    assert (x > 0).sum() > 100 and (x < 0).sum() > 10   # the reach gains and loses
    assert (-x <= inflow * (1 + 1e-9) + 1e-9).all()
    assert "exchange:cha01:stage_lowered_m" not in t   # RIV's limit isn't used


@pytest.mark.slow
def test_the_exchange_arrives_the_next_day(sfr):
    t = sfr["table"]
    k = "exchange:cha01:"
    np.testing.assert_allclose(t[k + "due"][1:], t[k + "to_stream"][:-1], rtol=0)
    np.testing.assert_allclose(t[k + "due"], t[k + "applied"] + t[k + "shortfall"], rtol=1e-6, atol=1e-6)
    # what remains short is the lag: today's loss is set against today's flow and taken from tomorrow's
    assert -t[k + "shortfall"].sum() < 1e-3 * np.abs(t[k + "to_stream"]).sum()


@pytest.mark.slow
def test_nothing_comes_from_outside(sfr):
    t = sfr["table"]
    for k in [k for k in t if k.startswith("groundwater:boundary:")]:
        assert np.abs(t[k]).max() < 1e-9, k
    assert np.abs(t["groundwater:closure_share"]).max() < 1e-6
