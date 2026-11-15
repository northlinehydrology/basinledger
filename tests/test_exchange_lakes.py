"""Lakes at a fixed stage coupled to channels (``[exchange] lakes``, coupling/exchange.py).

- MODFLOW 6 alone: what stops a lake package, and the lakes' map onto some of the exchange's channels (tier 0);
- the coupled runs (examples/synthetic/lakes.toml with a RIV exchange, lakes_sfr.toml with SFR; three years): the
  lake gains and loses; MODFLOW 6's budget closes every day and its lake term is the ledger's ``lakes``; the
  channel's exchange, the lake's part with it, reaches SWAT+ the next day; nothing comes from outside the models and
  the whole basin closes; a lake takes no more than its channel has, or under SFR the excess is a booked shortfall
  (tier 1 to 3);
- the refusals at start-up, listed by a check and stopping a run; a lake that exchanges nothing leaves the plain run
  as it was. Restarts: tests/test_restart.py.
"""

from __future__ import annotations

import datetime as dt
import importlib.util
import re

import flopy
import numpy as np
import pytest

from basinledger import config as C
from basinledger.balance import basin_balance
from basinledger.coupling.exchange import ChannelMap, StreamExchange, lake_map
from basinledger.coupling.modflow6 import Modflow6
from basinledger.driver import Basin
from conftest import ROOT, RUNS

SYNTHETIC = ROOT / "examples" / "synthetic"
DAYS = 1096
DAY = 86400.0


def synthetic_model():
    spec = importlib.util.spec_from_file_location("synthetic_model", SYNTHETIC / "model.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def small_model(name, **kw):
    ws = RUNS / name
    synthetic_model().build_small(ws, 3, dt.date(2010, 1, 1), lake_stage=88.5, **kw)
    mf = Modflow6(ws)
    mf.prepare_day()  # a package's entries are there once MODFLOW 6 has read its first day
    return mf


def config(example, **exchange):
    cfg = C.load(SYNTHETIC / example)
    cfg.run_folder = RUNS
    cfg.exchange = dict(cfg.exchange, **exchange)
    return cfg


# ---------- MODFLOW 6 alone ----------

def test_what_stops_a_lake_package(libmf6):
    mf = small_model("lakes-refused", drain_elevation=87.5)
    try:
        rivers = ChannelMap(["cha01"], np.zeros(10, int))
        lake = lambda p: (p, ChannelMap(["cha01"], np.zeros(mf.nbound(p), int)))  # noqa: E731
        drains = lake("drn")
        with pytest.raises(ValueError, match="'riv' is the exchange's own package"):
            StreamExchange(mf, "riv", rivers, drains, lakes=lake("riv"))
        with pytest.raises(ValueError, match="'drn' is the drains' package"):
            StreamExchange(mf, "riv", rivers, drains, lakes=lake("drn"))
        with pytest.raises(ValueError, match="'wel' is WEL; lakes at a fixed stage are a river \\(RIV\\) package"):
            StreamExchange(mf, "riv", rivers, drains, lakes=lake("wel"))
        ex = StreamExchange(mf, "riv", rivers, drains, lakes=lake("lakes"))
        assert ex.lakes[0] == "lakes"
    finally:
        mf.finalize()


def test_the_lakes_map_onto_some_of_the_channels(libmf6, tmp_path):
    """A lake maps to its outlet channel, in the exchange's order and indices; the lakes need not cover every coupled
    channel, but a channel the exchange doesn't couple is refused."""
    mf = small_model("lakes-map")
    try:
        swat = ["cha00", "cha01", "cha02"]
        coupled = ["cha00", "cha01"]
        (tmp_path / "lakes.csv").write_text("object,layer,row,col\n" + "".join(
            f"cha01,0,{r},{c}\n" for r, c in [(4, 8), (4, 9), (5, 8), (5, 9)]))
        m = lake_map("lakes.csv", tmp_path, mf, "lakes", swat, coupled)
        assert m.channels == coupled and m.entry_channel.tolist() == [1, 1, 1, 1]
        assert lake_map("all:cha01", tmp_path, mf, "lakes", swat, coupled).entry_channel.tolist() == [1, 1, 1, 1]
        with pytest.raises(ValueError, match="lakes: cha02 isn't a channel the exchange couples"):
            lake_map("all:cha02", tmp_path, mf, "lakes", swat, coupled)
        (tmp_path / "bad.csv").write_text("object,layer,row,col\ncha01,0,4,8\ncha02,0,4,9\n")
        with pytest.raises(ValueError, match="lakes: cha02 isn't a channel the exchange couples"):
            lake_map("bad.csv", tmp_path, mf, "lakes", swat, coupled)
        with pytest.raises(KeyError, match="'cha09' is not a SWAT\\+ channel"):
            lake_map("all:cha09", tmp_path, mf, "lakes", swat, coupled)
    finally:
        mf.finalize()


def test_the_config_names_a_package_and_a_map(tmp_path):
    text = (SYNTHETIC / "lakes.toml").read_text().replace(', map = "maps/lakes.csv"', "")
    (tmp_path / "lakes.toml").write_text(text)
    with pytest.raises(ValueError, match="\\[exchange\\] lakes is \\{package"):
        C.load(tmp_path / "lakes.toml")


# ---------- refused at start-up ----------

@pytest.mark.parametrize("change, refused", [
    (dict(lakes={"package": "riv", "map": "maps/lakes.csv"}), "'riv' is the exchange's own package"),
    (dict(lakes={"package": "wel", "map": "all:cha01"}), "'wel' is WEL"),
    (dict(lakes={"package": "drn", "map": "all:cha01"}, drains={"package": "drn", "map": "all:cha01"}),
     "'drn' is the drains' package"),
    (dict(lakes={"package": "lakes", "map": "all:cha09"}), "'cha09' is not a SWAT\\+ channel"),
    ("rivers", "\\[rivers\\] sets the stages of 'lakes'; a lake's stage is fixed"),
    ("transport", "transport with \\[exchange\\] lakes needs the lakes' own concentrations"),
])
def test_refused_at_start_up(swatplus, change, refused):
    """Each is listed by a check (``Basin.check`` builds the exchange), and stops a run on its first day."""
    if change == "transport":
        cfg = C.load(SYNTHETIC / "transport.toml")
        cfg.run_folder = RUNS
        cfg.mf6_builder_args = dict(cfg.mf6_builder_args, lake_stage=88.5)
        cfg.exchange = dict(cfg.exchange, lakes={"package": "lakes", "map": "maps/lakes.csv"})
    else:
        cfg = config("lakes.toml", **(change if isinstance(change, dict) else {}))
        if "drains" in cfg.exchange:
            cfg.mf6_builder_args = dict(cfg.mf6_builder_args, drain_elevation=87.5)
        if change == "rivers":
            cfg.rivers = {"package": "lakes", "map": "maps/lakes.csv", "stage": "surface_depth"}
    problems, _ = Basin.check(cfg)
    assert len(problems) == 1
    assert problems[0].startswith(("[exchange]", "[transport]")) and re.search(refused, problems[0])
    b = Basin(cfg, run="lakes-refused-run", days=2)
    try:
        with pytest.raises((ValueError, KeyError, NotImplementedError), match=refused):
            b.step()
    finally:
        b.finalize()


def test_a_good_lakes_example_has_no_problems(swatplus):
    for example in ("lakes.toml", "lakes_sfr.toml"):
        assert Basin.check(config(example))[0] == []


# ---------- the coupled runs ----------

def run(example, name, days=DAYS, **exchange):
    cfg = config(example, **exchange)
    cfg.outputs = dict(cfg.outputs, balance="cha01")
    b = Basin(cfg, run=name, days=days)
    b.run()
    b.write_outputs()
    out = dict(table=b.ledger.table(), due=dict(b.due), ws=b.mf6.ws, series=dict(b.series),
               heads=b.mf6.head(), run=b.folder)
    b.finalize()
    return out


@pytest.fixture(scope="module", params=["lakes.toml", "lakes_sfr.toml"])
def lakes(request, swatplus):
    out = run(request.param, request.param.split(".")[0])
    out["kind"] = "sfr" if "sfr" in request.param else "riv"
    return out


def budget(ws, text, paknam2=None):
    cbc = flopy.utils.CellBudgetFile(ws / "gwf.cbc", precision="double")
    kw = {"paknam2": paknam2} if paknam2 else {}
    return np.array([float(d["q"].sum()) if d.dtype.names else float(np.sum(d)) for d in cbc.get_data(text=text, **kw)])


@pytest.mark.slow
def test_the_lake_gains_and_loses(lakes):
    x = lakes["table"]["exchange:cha01:lakes"]
    assert (x > 0).sum() > 300 and (x < 0).sum() > 300


@pytest.mark.slow
def test_budget_closes_and_its_lake_term_is_the_ledgers(lakes):
    """MODFLOW 6's lake term is ``exchange:cha01:lakes`` every day, and with the river's (or SFR's) is the channel's
    ``to_stream``."""
    ws, t = lakes["ws"], lakes["table"]
    stream = "SFR" if lakes["kind"] == "sfr" else "RIV"
    terms = {"STO-SS": budget(ws, "STO-SS"), "STO-SY": budget(ws, "STO-SY"), "RCHA": budget(ws, "RCHA"),
             "WEL": budget(ws, "WEL"), "LAKES": budget(ws, "RIV", "LAKES"),
             stream: budget(ws, stream, "RIV") if stream == "RIV" else budget(ws, "SFR")}
    net = sum(terms.values())
    gross = sum(np.abs(v) for v in terms.values())
    assert len(net) == DAYS and np.max(np.abs(net) / gross) < 1e-7
    np.testing.assert_allclose(-terms["LAKES"], t["exchange:cha01:lakes"], rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(-terms["LAKES"] - terms[stream], t["exchange:cha01:to_stream"], rtol=1e-12, atol=1e-9)


@pytest.mark.slow
def test_the_lakes_water_reaches_swat_the_next_day(lakes):
    t = lakes["table"]
    k = "exchange:cha01:"
    np.testing.assert_allclose(t[k + "due"][1:], t[k + "to_stream"][:-1], rtol=0)
    assert t[k + "due"][0] == 0
    np.testing.assert_allclose(t[k + "due"], t[k + "applied"] + t[k + "shortfall"], rtol=1e-6, atol=1e-6)
    total = t[k + "to_stream"].sum()
    accounted = (t[k + "applied"] + t[k + "shortfall"]).sum() + lakes["due"]["cha01"]
    assert accounted == pytest.approx(total, rel=1e-6)


@pytest.mark.slow
def test_nothing_comes_from_outside_and_the_basin_closes(lakes):
    """The lake package is a coupled link, not a boundary flow; the aquifer closes with the links, and the whole basin
    closes as tightly as the two-way example's (tests/test_balance.py)."""
    t = lakes["table"]
    for k in [k for k in t if k.startswith("groundwater:boundary:")]:
        assert np.abs(t[k]).max() < 1e-9, k
    assert "groundwater:boundary:lakes" in t
    assert np.abs(t["groundwater:closure_share"]).max() < 1e-6
    links = (t["recharge:received"] - t["exchange:cha01:to_stream"] + t["groundwater:storage_release"]
             - t["route:well:pw1:delivered"])
    np.testing.assert_allclose(links, 0.0, atol=1e-6 * np.abs(t["recharge:received"]).max())
    r = basin_balance(lakes["run"], "cha01")
    p = r["precipitation_land"]
    assert abs(r["land_residual"]) <= 0.02 * (p + r["irrigation"])
    assert abs(r["inflow_from_elsewhere"]) <= 0.001 * r["outflow"]
    assert abs(r["groundwater_residual"]) <= 1e-6 * r["recharge_received"]
    assert abs(r["channels_residual"]) <= 0.001 * r["outflow"]
    assert r["max_daily_closure_share"] <= 1e-6


@pytest.mark.slow
def test_a_lake_takes_no_more_than_its_channel_has(lakes):
    """With RIV, the channel's rivers and lake together take no more than it has without today's exchange (its flow
    less what SWAT+ exchanged today), to the limit's stage tolerance; their stages are lowered by one depth. With SFR,
    the lake alone takes no more than the water SFR's top reach is given; SFR's losses from the same channel aren't
    counted, so a combined excess is a booked shortfall."""
    t, k = lakes["table"], "exchange:cha01:"
    available = np.maximum(np.array(lakes["series"]["flow:cha01"]) * DAY - t[k + "applied"], 0.0)
    lowered = t[k + "lake_stage_lowered_m"]
    assert (lowered > 0).sum() > 50
    tol = 1e-3 * (10 * 20.0 + 4 * 10.0)  # m3: the stage tolerance (1 mm) on every river and lake entry
    if lakes["kind"] == "riv":
        assert (-t[k + "to_stream"] <= available + tol).all()
        np.testing.assert_array_equal(lowered, t[k + "stage_lowered_m"])  # one depth for the channel's rivers and lake
    else:
        np.testing.assert_allclose(available, t[k + "sfr_inflow"], rtol=1e-12, atol=1e-9)
        assert (-t[k + "lakes"] <= t[k + "sfr_inflow"] + 1e-3 * 4 * 10.0).all()
        assert "exchange:cha01:stage_lowered_m" not in t
    # what remains short is booked: today's loss is set against today's flow and taken from tomorrow's
    assert -t[k + "shortfall"].sum() < 0.01 * np.abs(t[k + "to_stream"]).sum()


@pytest.mark.slow
@pytest.mark.parametrize("example, plain", [("lakes.toml", "two_way.toml"), ("lakes_sfr.toml", "sfr.toml")])
def test_a_lake_that_exchanges_nothing_leaves_the_plain_run(swatplus, monkeypatch, example, plain):
    """A lake of zero conductance, coupled, gives the run without lakes bit for bit: every account they share, the
    heads and the flows."""
    import basinledger.driver as D
    m = synthetic_model()
    monkeypatch.setattr(m, "LAKE_CONDUCTANCE", 0.0)
    monkeypatch.setattr(D, "_example_function", lambda cfg, spec, default: getattr(m, spec.partition(":")[2]))
    tag = plain.split(".")[0]
    a = run(plain, f"lakes-plain-{tag}", days=400)
    b = run(example, f"lakes-none-{tag}", days=400)
    assert np.abs(b["table"]["exchange:cha01:lakes"]).max() == 0
    assert np.array_equal(a["heads"], b["heads"])
    for key in a["table"]:
        assert np.array_equal(a["table"][key], b["table"][key], equal_nan=True), key
    for key in a["series"]:
        assert np.array_equal(a["series"][key], b["series"][key]), key
