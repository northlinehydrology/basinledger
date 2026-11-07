"""The reach heat balance (temperature.py): energy closes, mixing is by volume, and the floor's heat is counted."""

from __future__ import annotations

import pytest

from basinledger.temperature import RHO_C, Coefficients, reach_heat_balance

GEO = {"wd": 10.0, "dp": 1.0, "slp": 0.001, "len": 3.0, "mann": 0.05}
SUMMER = dict(tave=22.0, rhum=0.7, wind=3.0, solar=24.0, solar_max=30.0)
WINTER = dict(tave=-12.0, rhum=0.8, wind=6.0, solar=4.0, solar_max=9.0)


def run(parts, weather, **coef):
    return reach_heat_balance(parts, **weather, geo=GEO, q_in=2.0, q_out=2.0, area_ha=20000.0,
                              coef=Coefficients(**coef))


@pytest.mark.parametrize("weather", [SUMMER, WINTER])
def test_energy_closes(weather):
    r = run([(100000.0, 14.0), (40000.0, 8.0), (5000.0, 20.0)], weather)
    assert r["in_J"] + r["exchanged_J"] + r["floor_J"] == pytest.approx(r["out_J"], rel=1e-12)
    assert r["out_J"] == pytest.approx(RHO_C * r["volume_m3"] * r["t_out"], rel=1e-12)


def test_mixing_is_by_volume():
    r = run([(100000.0, 14.0), (50000.0, 8.0)], SUMMER, heat_exchange=False)
    assert r["t_init"] == pytest.approx((100000 * 14 + 50000 * 8) / 150000)
    assert r["t_out"] == pytest.approx(r["t_init"]) and r["exchanged_J"] == 0.0


def test_more_groundwater_cools_a_summer_reach():
    warm = run([(100000.0, 18.0), (10000.0, 7.0)], SUMMER)
    cool = run([(100000.0, 18.0), (60000.0, 7.0)], SUMMER)
    assert cool["t_out"] < warm["t_out"]


def test_exchange_heads_toward_equilibrium():
    r = run([(100000.0, 10.0)], SUMMER)
    assert (r["t_eq"] - r["t_init"]) * (r["t_out"] - r["t_init"]) > 0


def test_the_floor_is_counted():
    r = run([(100000.0, 0.5)], WINTER)
    assert r["t_out"] == pytest.approx(0.001)
    assert r["floor_J"] > 0
    assert r["in_J"] + r["exchanged_J"] + r["floor_J"] == pytest.approx(r["out_J"], rel=1e-9, abs=1e-3)


@pytest.mark.slow
def test_on_a_coupled_run(swatplus):
    """The heat balance on the synthetic two-way run (testing G4): it closes every reach and day, follows the
    seasons, and colder groundwater makes a cooler summer stream."""
    import numpy as np
    from basinledger import config as C
    from basinledger.driver import Basin
    from basinledger.temperature import stream_temperature
    from conftest import ROOT, RUNS
    cfg = C.load(ROOT / "examples" / "synthetic" / "two_way.toml")
    cfg.run_folder = RUNS
    cfg.outputs = dict(cfg.outputs, temperature=True)
    b = Basin(cfg, run="two-way-temperature", days=730)
    b.run()
    b.write_outputs()
    b.finalize()
    run = RUNS / "two-way-temperature"
    r = stream_temperature(run)
    h = r["heat"]
    np.testing.assert_allclose(h["in_J"] + h["exchanged_J"] + h["floor_J"], h["out_J"], rtol=1e-9, atol=1.0)
    t = r["temperature"][:, 0]
    months = np.array([int(d[5:7]) for d in r["dates"]])
    assert np.isfinite(t).all()
    assert t[np.isin(months, (7, 8))].mean() > t[np.isin(months, (1, 2))].mean()
    cold = stream_temperature(run, Coefficients(gw_temperature=2.0))["temperature"][:, 0]
    summer = np.isin(months, (7, 8))
    assert cold[summer].mean() <= t[summer].mean()
