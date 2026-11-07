"""The water-user interface: hourly to daily, the user's balance, constituents, role mapping and units."""

from __future__ import annotations

import math

import numpy as np
import pytest

from basinledger.testing import SyntheticUser
from basinledger.users import BalanceError, WaterUser


def hourly_series(days, **p):
    """The synthetic user's hourly rates, stepped by hand."""
    u = SyntheticUser()
    u.initialize(None)
    u.p.update(p)
    rows = []
    for _ in range(24 * days):
        u.update()
        rows.append({k: float(v[0]) for k, v in u._v.items()})
    return rows


def test_daily_totals_are_hourly_rates_integrated():
    wu = WaterUser("u", SyntheticUser())
    days = [wu.advance_day(d) for d in range(3)]
    rows = hourly_series(3)
    for d, use in enumerate(days):
        hrs = rows[24 * d: 24 * (d + 1)]
        w = sum(r["water_user_water__withdrawal_volume_flux"] for r in hrs) * 3600
        c = sum(r["water_user_water__consumption_volume_flux"] for r in hrs) * 3600
        assert use.withdrawal == pytest.approx(w, rel=1e-12)
        assert use.consumption == pytest.approx(c, rel=1e-12)
        assert use.rate(use.withdrawal) == pytest.approx(w / 86400, rel=1e-12)
        assert use.energy == pytest.approx(50.0 * 24, rel=1e-12)  # MWh


def test_user_balance_closes_every_day():
    wu = WaterUser("u", SyntheticUser())
    wu.bmi.p["gain_fraction"] = 0.02
    for d in range(30):
        use = wu.advance_day(d)
        assert abs(use.imbalance) <= 1e-9 * use.withdrawal
        assert use.gain == pytest.approx(0.02 * use.withdrawal, rel=1e-12)


def test_constituent_mass_withdrawn_is_returned():
    """Evaporation carries no dissolved solids, so the mass withdrawn leaves in the return, concentrated."""
    wu = WaterUser("u", SyntheticUser())
    use = wu.advance_day(0)
    tds_in = wu.bmi.p["intake_tds"]
    assert use.constituents["dissolved_solids"] == pytest.approx(use.withdrawal * tds_in, rel=1e-12)
    assert use.constituents["dissolved_solids"] / use.returned > tds_in * 3  # about 1 / (1 - 0.7)


def test_return_temperature_is_flow_weighted():
    wu = WaterUser("u", SyntheticUser())
    use = wu.advance_day(200)  # mid-July: warm
    assert 25.0 < use.return_temperature < 30.0


class Broken(SyntheticUser):
    def update(self):
        super().update()
        self._v["water_user_water__return_volume_flux"][0] *= 1.01


def test_a_user_that_loses_water_is_refused():
    with pytest.raises(BalanceError):
        WaterUser("bad", Broken()).advance_day(0)


class Renamed(SyntheticUser):
    """The same user under its own names and units (m3 d-1, kW, days), as a model from another package would be."""

    MAP = {"water_user_water__withdrawal_volume_flux": "mine_water__intake_volume_flux",
           "water_user_water__consumption_volume_flux": "mine_water__loss_volume_flux",
           "water_user_water__return_volume_flux": "mine_water__discharge_volume_flux",
           "water_user_water__gain_volume_flux": "mine_water__dewatering_volume_flux",
           "water_user__electric_power": "mine__electric_power"}

    def get_output_var_names(self):
        return tuple(self.MAP.values())

    def _inner(self, name):
        return {v: k for k, v in self.MAP.items()}.get(name, name)

    def get_var_units(self, name):
        return {"mine__electric_power": "kW"}.get(name, "m3 d-1")

    def get_value(self, name, dest):
        k = self._inner(name)
        dest[:] = self._v[k] * (1000.0 if k == "water_user__electric_power" else 86400.0)
        return dest

    def get_var_grid(self, name):
        return 0

    def get_var_type(self, name):
        return "float64"

    def get_time_units(self):
        return "d"

    def get_current_time(self):
        return self._t / 86400.0


def test_any_component_plugs_in_by_role():
    variables = {"withdrawal": "mine_water__intake_volume_flux", "consumption": "mine_water__loss_volume_flux",
                 "return": "mine_water__discharge_volume_flux", "gain": "mine_water__dewatering_volume_flux",
                 "power": "mine__electric_power"}
    a = WaterUser("a", SyntheticUser()).advance_day(5)
    comp = Renamed()
    b = WaterUser("b", comp, variables=variables, constituents={}).advance_day(5)
    assert b.withdrawal == pytest.approx(a.withdrawal, rel=1e-12)
    assert b.returned == pytest.approx(a.returned, rel=1e-12)
    assert b.energy == pytest.approx(a.energy, rel=1e-12)


def test_cooling_water_bmi_plugs_in_without_engine_changes(tmp_path):
    """CoolingWaterBmi from waterledger-models, mapped by role; its condensate is a gain."""
    bmi = pytest.importorskip("waterledger_models.datacenter.bmi")
    cfg = tmp_path / "cooling.yaml"
    cfg.write_text("design: towers\nit_capacity_mw: 100\n")
    variables = {"withdrawal": "data_center_cooling_water__withdrawal_volume_flux",
                 "consumption": "data_center_cooling_water__consumption_volume_flux",
                 "return": "data_center_cooling_water__return_volume_flux",
                 "gain": "data_center_cooling_water__condensate_volume_flux",
                 "power": "data_center__electric_power"}
    # a warm, humid July day: hourly air temperature and humidity
    inputs = {"atmosphere_bottom_air__temperature": lambda t: 24 + 6 * math.sin(2 * math.pi * (t / 3600 - 9) / 24),
              "atmosphere_bottom_air_water_vapor__relative_saturation": lambda t: 65.0}
    wu = WaterUser("campus", bmi.CoolingWaterBmi(), str(cfg), variables=variables, constituents={}, inputs=inputs)
    use = wu.advance_day(0)
    assert use.withdrawal > 0 and use.consumption > 0
    assert abs(use.imbalance) <= 1e-9 * use.withdrawal
    assert use.energy > 0
    assert wu.bmi.get_current_time() == pytest.approx(86400.0)


def test_reported_user_replays_its_series(tmp_path):
    from basinledger.generic import ReportedUser, write_series
    write_series(tmp_path / "s.csv", ["2020-01-02", "2020-01-03"], [100.0, 250.0], [30.0, 50.0])
    (tmp_path / "u.toml").write_text('series = "s.csv"\nstart = "2020-01-01"\n')
    wu = WaterUser("ghost", ReportedUser(), str(tmp_path / "u.toml"), constituents={})
    days = [wu.advance_day(d) for d in range(4)]
    assert [round(x.withdrawal, 9) for x in days] == [0.0, 100.0, 250.0, 0.0]
    assert days[2].consumption == pytest.approx(50.0) and days[2].returned == pytest.approx(200.0)
    bad = tmp_path / "bad.csv"
    bad.write_text("date,withdrawal_m3,consumption_m3,return_m3\n2020-01-02,100,30,60\n")
    (tmp_path / "b.toml").write_text('series = "bad.csv"\nstart = "2020-01-01"\n')
    with pytest.raises(ValueError, match="withdrawal != consumption"):
        WaterUser("bad", ReportedUser(), str(tmp_path / "b.toml"), constituents={})


def test_synthetic_user_returns_the_mass_it_takes_in(tmp_path):
    """With intake_species, the synthetic user takes an intake concentration as an input and returns all of the mass
    (evaporation carries none), as the transport tests need."""
    (tmp_path / "u.toml").write_text('intake_species = ["chloride"]\n')
    var = "water_user_intake_water_chloride__mass_concentration"
    wu = WaterUser("u", SyntheticUser(), str(tmp_path / "u.toml"))
    assert var in wu.bmi.get_input_var_names() and "chloride" in wu.constituents
    wu.inputs = {var: lambda t: 0.012}
    use = wu.advance_day(0)
    assert use.constituents["chloride"] == pytest.approx(0.012 * use.withdrawal, rel=1e-12)
    assert var not in SyntheticUser().get_input_var_names()  # only when asked for
