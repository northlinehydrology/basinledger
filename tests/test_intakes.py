"""Intakes: SWAT+ water allocation places of use, several on one channel, on any channel (swatplus-bmi 0.14.0).

SWAT+ 62.0.1's wallo_withdraw.f90 sized a withdrawal from the wrong place of use's demand: one on any channel but the
first took nothing (or a stray value), and places sharing a channel took each other's demands. The Rappahannock
study's intakes got nothing past the first channel (Nov 24, 2026); swatplus-bmi 0.14.0 patches it.

- an older swatplus-bmi is refused places of use on more than one channel, with the reason (tier 1);
- several places of use on one channel share its point of diversion, listed in the config's order (tier 1);
- coupled on the synthetic example, two reported users on Osu_1hru's channel, their point of diversion the file's
  second (where 0.13.1 gave nothing): each day the senior (first listed) gets what it asks, or the day's outflow if
  less, and the junior what it asks of the rest; the route pairs close, nothing is overdrawn, and the channel's
  outflow is the outflow without the users less what they took, the same day (tier 2);
- a place of use on a reservoir read from the config and written as a point of diversion of type res (tier 1);
- a plant asking more than its channel carries: its day scales by the share delivered, and its return of that share
  leaves the next day, never more than its intake gave (tier 2);
- two users on a reservoir: senior first, and the reservoir's daily balance closing with what they took (tier 2).
"""

from __future__ import annotations

import datetime as dt
import shutil

import numpy as np
import pytest
from conftest import ROOT, RUNS

from basinledger import config as C
from basinledger.config import UserConfig
from basinledger.coupling import swatplus as S
from basinledger.driver import Basin

DAY = 86400.0
DAYS = 240
ASK = {"senior": 0.0015, "junior": 0.0008}  # m3 s-1; Osu_1hru's channel carries about 0.001-0.01


def test_older_swatplus_bmi_refused_intakes_on_several_channels(monkeypatch):
    two = [("a", "cha01"), ("b", "cha02")]
    one = [("a", "cha01"), ("b", "cha01")]
    monkeypatch.setattr(S, "_swatplus_bmi_version", lambda: (0, 13))
    assert not S.places_of_use_fixed()
    (msg,) = S.places_of_use_problems(two)
    assert "2 channels" in msg and "0.14.0" in msg and "wallo_withdraw" in msg
    assert S.places_of_use_problems(one) == []  # one channel: its places are the file's first, read right
    monkeypatch.setattr(S, "_swatplus_bmi_version", lambda: (0, 14))
    assert S.places_of_use_fixed() and S.places_of_use_problems(two) == []


def test_places_on_one_channel_share_its_point_of_diversion_in_order(swatplus, tmp_path):
    p = tmp_path / "Osu_1hru"
    shutil.copytree(swatplus / "Osu_1hru", p)
    S.add_places_of_use(p, [("senior", "cha01"), ("junior", "cha01")])
    pod = (p / "point_of_diver.wro").read_text().splitlines()
    assert pod[1] == "1"
    assert pod[3].split() == ["1", "cha01", "cha", "1", "2",
                              "1", "senior", "1", "out", "1", "senior",
                              "2", "junior", "1", "out", "2", "senior"]


def test_a_place_of_use_on_a_reservoir_read_and_written(swatplus, tmp_path):
    from basinledger.config import _place_of_use
    assert _place_of_use("x.toml", {"name": "a", "channel": "cha01"}) == ("a", "cha01")
    assert _place_of_use("x.toml", {"name": "b", "reservoir": "res30"}) == ("b", "res:res30")
    for bad in ({"name": "c"}, {"name": "c", "channel": "cha01", "reservoir": "res01"}, {"channel": "cha01"}):
        with pytest.raises(ValueError, match="one of channel or reservoir"):
            _place_of_use("x.toml", bad)
    assert S.place_source("res:res30") == ("res", "res30") and S.place_source("cha01") == ("cha", "cha01")
    p = _example("existing").reservoir_project(tmp_path / "osu_res", swatplus)
    S.add_places_of_use(p, [("river", "cha01"), ("lake_a", "res:res01"), ("lake_b", "res:res01")])
    pod = [ln.split() for ln in (p / "point_of_diver.wro").read_text().splitlines()]
    assert pod[1] == ["2"]
    assert pod[3][:5] == ["1", "cha01", "cha", "1", "1"]
    assert pod[4] == ["2", "res01", "res", "1", "2", "2", "lake_a", "1", "out", "2", "senior",
                      "3", "lake_b", "1", "out", "3", "senior"]
    pou = [ln.split() for ln in (p / "place_of_use.wro").read_text().splitlines()]
    assert pou[6][:5] == ["1", "2", "res01", "res", "1"]
    with pytest.raises(ValueError, match="no reservoir 'res09'"):
        S.add_places_of_use(p, [("lake", "res:res09")])


def _second_point_of_diversion(model, folder):
    """The run's point_of_diver.wro with a point of diversion before Basin Ledger's that serves nothing, so the
    places of use are read through the file's second, where SWAT+ 62.0.1 gave them nothing."""
    if model != "swatplus":
        return
    f = folder / "point_of_diver.wro"
    lines = f.read_text().splitlines()
    rows = [ln.split() for ln in lines[3:] if ln.split()]
    for r in rows:
        r[0] = str(int(r[0]) + 1)
    f.write_text("\n".join([lines[0], str(len(rows) + 1), lines[2], "1 nothing cha 1 0"]
                           + [" ".join(r) for r in rows]) + "\n")


def _config(tmp_path, users: bool):
    cfg = C.load(ROOT / "examples" / "synthetic" / "surface.toml")
    cfg.run_folder = RUNS
    cfg.places_of_use = [("senior_intake", "cha01"), ("junior_intake", "cha01")]
    cfg.point_sources, cfg.utilities, cfg.rules = [], [], []
    start = S._period(cfg.swat_project)[0]
    cfg.users = []
    if not users:
        return cfg
    for name, q in ASK.items():
        rows = ["date,withdrawal_m3,consumption_m3,return_m3"]
        rows += [f"{start + dt.timedelta(days=d)},{q * DAY},{q * DAY},0" for d in range(DAYS)]
        (tmp_path / f"{name}.csv").write_text("\n".join(rows) + "\n")
        (tmp_path / f"{name}.toml").write_text(f'series = "{name}.csv"\nstart = "{start}"\n')
        cfg.users.append(UserConfig(name=name, component="basinledger.generic:ReportedUser",
                                    config=tmp_path / f"{name}.toml",
                                    withdrawal=f"swat:place_of_use:{name}_intake", returns=f"offsite:{name}"))
    return cfg


def _run(tmp_path, name, users):
    b = Basin(_config(tmp_path, users), run=name, users=users, days=DAYS, prepare=_second_point_of_diversion)
    b.run()
    out = b.ledger.table(), np.array(b.series["flow:cha01"]) * DAY
    b.finalize()
    return out


@pytest.mark.slow
def test_intakes_on_a_shared_channel_get_their_own_demands_senior_first(swatplus, tmp_path):
    if not S.places_of_use_fixed():
        pytest.skip("needs swatplus-bmi 0.14.0 or later (SWAT+'s place-of-use fix)")
    t, q = _run(tmp_path, "intakes-users", True)
    _, q0 = _run(tmp_path, "intakes-base", False)
    r = {n: f"route:swat:place_of_use:{n}_intake:" for n in ASK}
    got = {n: t[r[n] + "delivered"] for n in ASK}
    for n, q_ask in ASK.items():
        np.testing.assert_allclose(t[r[n] + "requested"], q_ask * DAY, rtol=1e-12)
        np.testing.assert_allclose(t[r[n] + "requested"], t[r[n] + "delivered"] + t[r[n] + "shortfall"], rtol=1e-6)
        assert (t[r[n] + "overdraw"] == 0).all()

    # the outflow without the users is what the channel had; the senior takes first, the junior from the rest
    tol = 1e-5 * q0 + 1e-6
    senior = np.minimum(ASK["senior"] * DAY, q0)
    junior = np.minimum(ASK["junior"] * DAY, q0 - senior)
    assert (np.abs(got["senior"] - senior) <= tol).all()
    assert (np.abs(got["junior"] - junior) <= tol).all()
    # when the flow allows, each intake gets all it asks; days of both kinds are in the run
    plenty = q0 > 1.01 * sum(ASK.values()) * DAY
    assert plenty.sum() > 20 and (~plenty).sum() > 20
    for n, q_ask in ASK.items():
        np.testing.assert_allclose(got[n][plenty], q_ask * DAY, rtol=1e-5)
    # the channel lost exactly what the intakes took, the same day
    np.testing.assert_allclose(q, q0 - got["senior"] - got["junior"], rtol=1e-5, atol=1e-3)


def _example(name):
    """A module of the synthetic example's (examples/synthetic/<name>.py), loaded by path."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(f"synthetic_{name}", ROOT / "examples" / "synthetic" / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _reported(tmp_path, name, start, withdrawal, consumed):
    """A ReportedUser's files: ``withdrawal`` m3 s-1 every day, ``consumed`` of it consumed and the rest returned."""
    w = withdrawal * DAY
    rows = ["date,withdrawal_m3,consumption_m3,return_m3"]
    rows += [f"{start + dt.timedelta(days=d)},{w},{consumed * w},{w - consumed * w}" for d in range(DAYS)]
    (tmp_path / f"{name}.csv").write_text("\n".join(rows) + "\n")
    (tmp_path / f"{name}.toml").write_text(f'series = "{name}.csv"\nstart = "{start}"\n')
    return tmp_path / f"{name}.toml"


@pytest.mark.slow
def test_a_short_intake_returns_only_what_it_got_the_next_day(swatplus, tmp_path):
    """A plant asking more than the channel carries, consuming 5% and returning the rest to the same channel: each day
    its use scales by the share its intake delivered, and its return leaves the next day. So it never returns water
    the channel didn't give (the Virginia study's Lake Pelham intake returned 829,326 m3 of 8,697 delivered)."""
    if not S.places_of_use_fixed():
        pytest.skip("needs swatplus-bmi 0.14.0 or later (SWAT+'s place-of-use fix)")
    cfg = C.load(ROOT / "examples" / "synthetic" / "surface.toml")
    cfg.run_folder = RUNS
    cfg.places_of_use = [("plant_intake", "cha01")]
    cfg.point_sources, cfg.utilities, cfg.rules = [("plant_outfall", "cha01")], [], []
    start = S._period(cfg.swat_project)[0]
    cfg.users = [UserConfig(name="plant", component="basinledger.generic:ReportedUser",
                            config=_reported(tmp_path, "plant", start, 0.004, 0.05),
                            withdrawal="swat:place_of_use:plant_intake", returns="swat:point_source:plant_outfall")]
    b = Basin(cfg, run="intakes-short", days=DAYS)
    b.run()
    t = b.ledger.table()
    b.finalize()
    k, r = "user:plant:", "route:swat:place_of_use:plant_intake:"
    share = t[k + "delivered_share"]
    assert (share < 0.99).sum() > 20 and share.min() < 0.5 and (share > 0.999).sum() > 5  # short days, and not
    np.testing.assert_allclose(t[k + "delivered"], t[r + "delivered"], rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(t[k + "delivered"], t[k + "consumption_met"] + t[k + "return_due"], rtol=1e-9)
    np.testing.assert_allclose(t[k + "consumption_met"], 0.05 * t[k + "delivered"], rtol=1e-9)
    sent = t["route:swat:point_source:plant_outfall:delivered"]
    assert sent[0] == 0
    np.testing.assert_allclose(sent[1:], t[k + "return_due"][:-1], rtol=1e-6, atol=1e-6)
    assert (sent[1:] <= 0.95 * t[r + "delivered"][:-1] * (1 + 1e-6) + 1e-6).all()
    # over the run, what it returned and holds is what its intake gave less what it consumed
    np.testing.assert_allclose(sent.sum() + t[k + "return_held"][-1],
                               t[r + "delivered"].sum() - t[k + "consumption_met"].sum(), rtol=1e-6)  # float32


@pytest.mark.slow
def test_intakes_in_a_reservoir_take_from_its_storage_senior_first(swatplus, tmp_path):
    """Two users on one reservoir (Osu_1hru with res01 below its channel): SWAT+ takes their water out of the
    reservoir's storage after the day's inflow and before its release. Each day the reservoir's volume changes by its
    inflow, precipitation, release, evaporation and seepage, less exactly what the intakes delivered."""
    if not S.places_of_use_fixed():
        pytest.skip("needs swatplus-bmi 0.14.0 or later (SWAT+'s place-of-use fix)")
    from basinledger.coupling import lakes as LK
    reservoir_project = _example("existing").reservoir_project
    cfg = C.load(ROOT / "examples" / "synthetic" / "surface.toml")
    cfg.run_folder = RUNS
    cfg.swat_project = reservoir_project(tmp_path / "osu_res", swatplus)
    cfg.places_of_use = [("senior_intake", "res:res01"), ("junior_intake", "res:res01")]
    cfg.point_sources, cfg.utilities, cfg.rules = [], [], []
    start = S._period(cfg.swat_project)[0]
    cfg.users = [UserConfig(name=n, component="basinledger.generic:ReportedUser",
                            config=_reported(tmp_path, n, start, q, 1.0),
                            withdrawal=f"swat:place_of_use:{n}_intake", returns=f"offsite:{n}")
                 for n, q in ASK.items()]
    b = Basin(cfg, run="intakes-reservoir", days=DAYS)

    def get(var):  # the one reservoir
        return float(b.sw.get(var)[0])

    v = [get(LK.VOLUME)]
    flows = []
    for _ in range(DAYS):
        b.step()
        v.append(get(LK.VOLUME))
        flows.append([get(x) * DAY for x in (LK.INFLOW, LK.PRECIP, LK.OUTFLOW, LK.EVAP, LK.SEEP)])
    t = b.ledger.table()
    b.finalize()
    v, (q_in, p, q_out, e, seep) = np.array(v), np.array(flows).T
    r = {n: f"route:swat:place_of_use:{n}_intake:" for n in ASK}
    took = t[r["senior"] + "delivered"] + t[r["junior"] + "delivered"]
    assert took.sum() > 0.5 * sum(ASK.values()) * DAY * DAYS
    for n in ASK:
        assert (t[r[n] + "overdraw"] == 0).all()
        np.testing.assert_allclose(t[r[n] + "requested"], t[r[n] + "delivered"] + t[r[n] + "shortfall"], rtol=1e-6)
    # the junior gets water only on days the senior got all it asked
    junior_got = t[r["junior"] + "delivered"] > 0
    np.testing.assert_allclose(t[r["senior"] + "delivered"][junior_got], ASK["senior"] * DAY, rtol=1e-5)
    # the reservoir loses exactly what the intakes took, the same day (float32 volumes and rates)
    resid = np.diff(v) - (q_in + p - q_out - e - seep - took)
    scale = v[1:] + q_in + q_out + took
    assert (np.abs(resid) <= 1e-5 * scale + 1e-3).all(), np.abs(resid).max()
