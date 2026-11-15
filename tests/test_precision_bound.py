"""A ledger pair's bound, and the precision of SWAT+'s HRU areas as written (Nov 13, 2026)."""

import pytest

from basinledger.coupling.swatplus import _half_ulp, hru_area_precision
from basinledger.ledger import Ledger


def test_half_ulp():
    assert _half_ulp("0.0000806") == pytest.approx(5e-8)
    assert _half_ulp("2.61659") == pytest.approx(5e-6)
    assert _half_ulp("1116") == pytest.approx(0.5)
    assert _half_ulp("1.5e-3") == pytest.approx(5e-5)


def test_hru_area_precision(tmp_path):
    (tmp_path / "hru.con").write_text(
        "hru.con: test\n      id  name                gis_id          area\n"
        "       1  hru0001         1       0.09000\n       2  hru0002         2       2.61659\n")
    (tmp_path / "rout_unit.ele").write_text(
        "rout_unit.ele: test\n      id  name                   obj_typ    obj_id          frac               dlr\n"
        "       1  hru0001                     hru         1     0.0000806                 0\n"
        "       2  hru0002                     hru         2     0.0023433                 0\n")
    p = hru_area_precision(tmp_path)
    # hru0001: half of 1e-5 ha in its area, plus half of 1e-7 of its unit's area (900 m2 / 0.0000806)
    assert p["hru0001"] == pytest.approx(0.05 + 5e-8 * 900 / 0.0000806)
    assert p["hru0001"] / 900 > 5e-4  # a small HRU in a large unit: more than 1e-5 of its area
    assert p["hru0002"] == pytest.approx(0.05 + 5e-8 * 26165.9 / 0.0023433)


def test_pair_bound():
    led = Ledger(rtol=1e-5, atol=1e-6)
    led.pair("a", "b", bound="why")
    led.enter("a", 4.758428)
    led.enter("b", 4.758341)
    with pytest.raises(AssertionError):
        led.close_day()  # 9e-5 apart, over 1e-6 + 1e-5 x 4.76
    led.enter("a", 4.758428)
    led.enter("b", 4.758341)
    led.enter("why", 2e-4)
    led.close_day()  # within the day's bound
    led.enter("a", 4.758428)
    led.enter("b", 4.0)
    led.enter("why", 2e-4)
    with pytest.raises(AssertionError):
        led.close_day()  # the bound doesn't excuse a real difference
