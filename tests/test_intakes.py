"""Intakes: SWAT+ water allocation places of use, several on one channel, on any channel (swatplus-bmi 0.14.0).

SWAT+ 62.0.1's wallo_withdraw.f90 sized a withdrawal from the wrong place of use's demand: one on any channel but the
first took nothing (or a stray value), and places sharing a channel took each other's demands. The Rappahannock
study's intakes got nothing past the first channel (Nov 24, 2026); swatplus-bmi 0.14.0 patches it.

- an older swatplus-bmi is refused places of use on more than one channel, with the reason (tier 1);
- several places of use on one channel share its point of diversion, listed in the config's order (tier 1);
- coupled on the synthetic example, two reported users on Osu_1hru's channel, their point of diversion the file's
  second (where 0.13.1 gave nothing): each day the senior (first listed) gets what it asks, or the day's outflow if
  less, and the junior what it asks of the rest; the route pairs close, nothing is overdrawn, and the channel's
  outflow is the outflow without the users less what they took, the same day (tier 2).
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
