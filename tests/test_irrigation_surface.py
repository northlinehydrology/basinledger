"""Irrigation from SWAT+'s own channels and reservoirs, booked (coupling/irrigation.py, ``NativeIrrigation``;
docs/design.md, "Irrigation from surface water"; swatplus-bmi 0.14.0).

Osu_1hru's paddy table (``pad_irr_r``) switched to irrigate 1 mm (100 m3 on its 10 ha HRU) from its channel with
irr.ops' sprinkler_med (efficiency 0.85, no runoff), as swatplus-bmi's own tests do. Osu_1hru's channel routes all it
holds out each day, so its storage coefficient (``scoef``, parameters.bsn) is halved, so it holds some water some days.

- the HRUs that irrigate from a SWAT+ water body, and the water bodies, read from the project (tier 1);
- coupled both ways on the synthetic example for three years (tier 2): each day the channel gives what the HRU takes,
  which is what reached the soil plus the application's loss (15%); SWAT+'s own irrigation in its land budget is
  that water, not water from outside; and the whole-basin balance's channels close with the irrigation they gave.
"""

from __future__ import annotations

import shutil

import numpy as np
import pytest
from conftest import ROOT, RUNS

from basinledger import config as C
from basinledger.balance import basin_balance
from basinledger.coupling import irrigation as IR
from basinledger.driver import Basin

AMOUNT, AREA, EFF = 1.0, 10.0, 0.85  # mm an action, ha, sprinkler_med's efficiency
DEMAND = AMOUNT * AREA * 10.0         # m3 an action takes


def from_channel(folder, src="cha 1"):
    """The paddy table irrigates AMOUNT mm from ``src``; its channel keeps some water from day to day."""
    p = folder / "parameters.bsn"
    lines = p.read_text().splitlines()
    f = lines[2].split()
    f[lines[1].split().index("scoef")] = "0.50000"
    lines[2] = "  ".join(f)
    p.write_text("\n".join(lines) + "\n")
    p = folder / "lum.dtl"
    lines = p.read_text().splitlines()
    i = next(k for k, ln in enumerate(lines) if ln.split()[:1] == ["pad_irr_r"])
    j = next(k for k in range(i, i + 10) if lines[k].split()[:2] == ["irrigate", "null"])
    lines[j] = f"irrigate  {src}  irr_src  sprinkler_med  {AMOUNT}  50.  null  y  n"
    p.write_text("\n".join(lines) + "\n")


@pytest.fixture
def osu(swatplus, tmp_path):
    p = tmp_path / "Osu_1hru"
    shutil.copytree(swatplus / "Osu_1hru", p)
    return p


def test_the_hrus_and_water_bodies_irrigating(osu):
    assert IR.native(osu) == ([], set())  # the paddy irrigates from an outside source
    from_channel(osu)
    assert IR.native(osu) == (["hru0001"], {("cha", 1)})
    from_channel(osu, "res 2")  # the other irrigate row now from a reservoir
    hrus, objects = IR.native(osu)
    assert hrus == ["hru0001"] and ("res", 2) in objects


@pytest.mark.slow
def test_irrigation_from_a_channel_comes_out_of_it(swatplus):
    cfg = C.load(ROOT / "examples" / "synthetic" / "two_way.toml")
    cfg.run_folder = RUNS
    cfg.outputs = dict(cfg.outputs, balance="cha01")
    b = Basin(cfg, run="two-way-irrigation", days=1096,
              prepare=lambda model, folder: from_channel(folder) if model == "swatplus" else None)
    if b.native_irrigation is None:
        b.finalize()
        pytest.skip("needs swatplus-bmi 0.14.0 (its irrigation withdrawal reads)")
    b.run()
    b.write_outputs()
    t = b.ledger.table()
    b.finalize()

    k = "irrigation:native:"
    took = t[k + "withdrawn"]
    days = took > 0
    assert days.sum() >= 5                                              # actions that found water
    np.testing.assert_allclose(took[days] / DEMAND, np.round(took[days] / DEMAND), rtol=1e-5)  # 100 m3 each, or none
    np.testing.assert_allclose(t[k + "from_channels"], took, rtol=1e-6)
    np.testing.assert_allclose(t["irrigation:channel:cha01"], took, rtol=1e-6)
    assert (t[k + "from_reservoirs"] == 0).all() and (t[k + "excess"] == 0).all()
    np.testing.assert_allclose(t[k + "applied"][days], EFF * took[days], rtol=1e-5)  # sprinkler_med, no runoff
    np.testing.assert_allclose(t[k + "loss"][days], (1 - EFF) * took[days], rtol=1e-4)

    r = basin_balance(RUNS / "two-way-irrigation", "cha01")
    assert r["irrigation_from_streams"] == pytest.approx(took.sum(), rel=1e-9)
    assert r["irrigation_loss"] == pytest.approx((1 - EFF) * took.sum(), rel=1e-4)
    # SWAT+'s own irrigation in its land budget is the water the channel gave less the loss: none from outside
    assert r["irrigation"] == pytest.approx(EFF * took.sum(), rel=1e-3)
    # the channels close with what they gave; without it they'd be off by about all of it (1,300 m3 over the three
    # years, against -113 m3 with it, Nov 25, 2026)
    assert abs(r["channels_residual"]) <= 0.001 * r["outflow"]
    assert abs(r["channels_residual"]) < 0.2 * abs(r["channels_residual"] + r["irrigation_from_streams"])
