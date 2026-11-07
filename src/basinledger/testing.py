"""A synthetic water user for the tests: a BMI component with the common water-user variables and nothing else.

It stands for any large user. Withdrawal follows a daily and a seasonal cycle; a fixed fraction is consumed (evaporated);
an optional gain enters with the return; all dissolved solids withdrawn leave in the return, since evaporation carries
none. Its numbers mean nothing physically; they only exercise the interface.

Config (TOML, all optional):
    withdrawal = 0.006          # m3 s-1, mean
    consumed_fraction = 0.7
    gain_fraction = 0.0         # of withdrawal
    intake_tds = 0.3            # kg m-3 of dissolved solids in the water withdrawn
    power = 50.0                # MW
    dry_withdrawal = 0.1        # in dry mode, the share of withdrawal still taken (staff and process water)
    dry_power = 1.15            # in dry mode, power x this
    time_step = 3600            # s
    intake_species = []         # e.g. ["chloride"]: an input water_user_intake_water_<c>__mass_concentration (kg m-3)
                                #   for each, all of it returned (water_user_return_water_<c>__mass_flux)
"""

from __future__ import annotations

import math
import tomllib

import numpy as np
from bmipy import Bmi

OUTPUTS = {
    "water_user_water__withdrawal_volume_flux": "m3 s-1",
    "water_user_water__gain_volume_flux": "m3 s-1",
    "water_user_water__consumption_volume_flux": "m3 s-1",
    "water_user_water__return_volume_flux": "m3 s-1",
    "water_user_return_water__temperature": "degC",
    "water_user_return_water_dissolved_solids__mass_flux": "kg s-1",
    "water_user__electric_power": "MW",
}
INPUTS = {
    "water_user_water__withdrawal_volume_flux_limit": "m3 s-1",  # a cap set by a feedback rule; inf = none
    "water_user__dry_mode_flag": "1",                            # 1 = dry: no evaporative cooling
}


class SyntheticUser(Bmi):
    _name = "Basin Ledger synthetic water user"

    def __init__(self):
        self._outputs, self._inputs = dict(OUTPUTS), dict(INPUTS)
        self._intake: list[str] = []
        self._v = {k: np.zeros(1) for k in [*OUTPUTS, *INPUTS]}
        self._v["water_user_water__withdrawal_volume_flux_limit"][:] = np.inf
        self._t = 0.0
        self._dt = 3600.0
        self._end = 100 * 365 * 86400.0
        self.p = dict(withdrawal=0.006, consumed_fraction=0.7, gain_fraction=0.0, intake_tds=0.3, power=50.0,
                      dry_withdrawal=0.1, dry_power=1.15)

    def initialize(self, config_file=None):
        cfg = {}
        if config_file:
            with open(config_file, "rb") as f:
                cfg = tomllib.load(f)
        self._dt = float(cfg.pop("time_step", 3600.0))
        self._intake = list(cfg.pop("intake_species", []))
        for c in self._intake:
            self._inputs[f"water_user_intake_water_{c}__mass_concentration"] = "kg m-3"
            self._outputs[f"water_user_return_water_{c}__mass_flux"] = "kg s-1"
        self._v.update({k: np.zeros(1) for k in [*self._outputs, *self._inputs] if k not in self._v})
        unknown = set(cfg) - set(self.p)
        if unknown:
            raise ValueError(f"unknown parameters: {sorted(unknown)}")
        self.p.update({k: float(v) for k, v in cfg.items()})
        self._t = 0.0

    def update(self):
        mid = self._t + self._dt / 2  # the step's midpoint sets its rates
        self._t += self._dt
        hour = (mid / 3600.0) % 24
        doy = (mid / 86400.0) % 365.25
        w = self.p["withdrawal"] * (1 + 0.2 * math.sin(2 * math.pi * (hour - 9) / 24)) \
            * (1 + 0.3 * math.sin(2 * math.pi * (doy - 110) / 365.25))
        dry = self._v["water_user__dry_mode_flag"][0] >= 0.5
        if dry:
            w *= self.p["dry_withdrawal"]
        w = min(w, float(self._v["water_user_water__withdrawal_volume_flux_limit"][0]))
        f = 0.2 if dry else self.p["consumed_fraction"]
        g = self.p["gain_fraction"] * w
        c = f * w
        self._v["water_user_water__withdrawal_volume_flux"][0] = w
        self._v["water_user_water__gain_volume_flux"][0] = g
        self._v["water_user_water__consumption_volume_flux"][0] = c
        self._v["water_user_water__return_volume_flux"][0] = (w - c) + g
        self._v["water_user_return_water_dissolved_solids__mass_flux"][0] = w * self.p["intake_tds"]
        self._v["water_user_return_water__temperature"][0] = 25.0 + 5.0 * math.sin(2 * math.pi * (doy - 110) / 365.25)
        self._v["water_user__electric_power"][0] = self.p["power"] * (self.p["dry_power"] if dry else 1.0)
        for c in self._intake:  # evaporation carries none: what is taken in is returned
            self._v[f"water_user_return_water_{c}__mass_flux"][0] = \
                w * self._v[f"water_user_intake_water_{c}__mass_concentration"][0]

    def update_until(self, time):
        while self._t + self._dt <= time + 1e-9:
            self.update()

    def finalize(self):
        pass

    def get_component_name(self):
        return self._name

    def get_input_item_count(self):
        return len(self._inputs)

    def get_output_item_count(self):
        return len(self._outputs)

    def get_input_var_names(self):
        return tuple(self._inputs)

    def get_output_var_names(self):
        return tuple(self._outputs)

    def get_var_grid(self, name):
        self._units(name)
        return 0

    def _units(self, name):
        return {**self._outputs, **self._inputs}[name]

    def get_var_type(self, name):
        return "float64"

    def get_var_units(self, name):
        return self._units(name)

    def get_var_itemsize(self, name):
        return 8

    def get_var_nbytes(self, name):
        return 8

    def get_var_location(self, name):
        return "node"

    def get_current_time(self):
        return self._t

    def get_start_time(self):
        return 0.0

    def get_end_time(self):
        return self._end

    def get_time_units(self):
        return "s"

    def get_time_step(self):
        return self._dt

    def get_value(self, name, dest):
        dest[:] = self._v[name]
        return dest

    def get_value_ptr(self, name):
        return self._v[name]

    def get_value_at_indices(self, name, dest, inds):
        dest[:] = self._v[name].take(inds)
        return dest

    def set_value(self, name, src):
        if name not in self._inputs:
            raise KeyError(f"{name} is not an input variable")
        self._v[name][:] = src

    def set_value_at_indices(self, name, inds, src):
        self.set_value(name, src)

    def get_grid_rank(self, grid):
        return 0

    def get_grid_size(self, grid):
        return 1

    def get_grid_type(self, grid):
        return "scalar"

    def get_grid_node_count(self, grid):
        return 1

    def get_grid_edge_count(self, grid):
        return 0

    def get_grid_face_count(self, grid):
        return 0

    def _none(self, *a):
        raise NotImplementedError("scalar grid")

    get_grid_shape = get_grid_spacing = get_grid_origin = _none
    get_grid_x = get_grid_y = get_grid_z = _none
    get_grid_edge_nodes = get_grid_face_edges = get_grid_face_nodes = get_grid_nodes_per_face = _none
