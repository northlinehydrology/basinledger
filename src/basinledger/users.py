"""The water-user interface: any BMI component that withdraws, consumes and returns water.

A water user (a data center's cooling, a mine, a battery plant, a fab) is a BMI component stepped at its own time step,
usually hourly. Basin Ledger reads it through roles, not through model-specific names:

    role          common variable name                                       units
    withdrawal    water_user_water__withdrawal_volume_flux                   m3 s-1   taken from a source
    gain          water_user_water__gain_volume_flux                         m3 s-1   water entering by other ways (condensate from air, ore moisture)
    consumption   water_user_water__consumption_volume_flux                  m3 s-1   evaporated, or leaving in product or waste
    return        water_user_water__return_volume_flux                       m3 s-1   discharged (to a sewer or an outfall)
    storage       water_user_water__volume                                   m3       water held on site (ponds, tailings), optional
    temperature   water_user_return_water__temperature                       degC     of the return, optional
    power         water_user__electric_power                                 MW       optional
    constituents  water_user_return_water_<name>__mass_flux                  kg s-1   in the return, one per constituent, optional

Inputs a feedback rule sets (optional; a rule asking for one the component lacks is refused):

    limit         water_user_water__withdrawal_volume_flux_limit             m3 s-1   a cap on withdrawal; inf = none
    dry_mode      water_user__dry_mode_flag                                  1        1 = run without evaporative cooling

The balance every user must close, every step:  withdrawal + gain = consumption + return + d(storage)/dt.

A component using other names is mapped role by role in the example's config (CoolingWaterBmi reports
``data_center_cooling_water__withdrawal_volume_flux`` and so on, and its condensate is a gain). Nothing here knows what
kind of user it is.

Values read after ``update()`` describe the step just ended. The daily hydrologic models get the day's total volumes:
``WaterUser.advance_day`` steps the component to the end of the day and integrates each rate over its steps.
"""

from __future__ import annotations

import importlib
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from . import units as U

COMMON = {
    "withdrawal": "water_user_water__withdrawal_volume_flux",
    "gain": "water_user_water__gain_volume_flux",
    "consumption": "water_user_water__consumption_volume_flux",
    "return": "water_user_water__return_volume_flux",
    "storage": "water_user_water__volume",
    "temperature": "water_user_return_water__temperature",
    "power": "water_user__electric_power",
}
CONTROLS = {
    "limit": "water_user_water__withdrawal_volume_flux_limit",
    "dry_mode": "water_user__dry_mode_flag",
}
VOLUME_ROLES = ("withdrawal", "gain", "consumption", "return")
CONSTITUENT = "water_user_return_water_{}__mass_flux"


class BalanceError(RuntimeError):
    pass


@dataclass
class DailyUse:
    """One day of a water user, as the daily models see it. Volumes in m3, masses in kg, energy in MWh."""

    day: int
    withdrawal: float = 0.0
    gain: float = 0.0
    consumption: float = 0.0
    returned: float = 0.0
    storage_change: float = 0.0
    constituents: dict[str, float] = field(default_factory=dict)  # kg returned
    return_temperature: float | None = None  # flow-weighted mean, degC
    energy: float | None = None  # MWh
    reported: dict[str, float] = field(default_factory=dict)  # other outputs the example asked for (``report``)

    def rate(self, volume: float) -> float:
        """A day's volume as the day's mean rate, m3 s-1."""
        return volume / U.DAY

    @property
    def imbalance(self) -> float:
        return self.withdrawal + self.gain - self.consumption - self.returned - self.storage_change


def load_component(spec: str):
    """'package.module:Class' -> a new instance."""
    mod, _, name = spec.partition(":")
    return getattr(importlib.import_module(mod), name)()


class WaterUser:
    """A BMI water-user component behind the common roles.

    variables: role -> the component's variable name, for components that don't use the common names.
    constituents: constituent name -> the component's variable (default: the common pattern, if the component has it).
    inputs: component input variable -> function of seconds since the run start, set before each step (forcing).
    report: name -> any other output to book each day, a variable name or {variable, how}: a power (MW or W) is
    integrated to MWh (``MW``); an amount of time (h or s) is taken as a running total, and the day's increase is
    booked (``change``: CoolingWaterBmi's ``data_center__setpoint_lost_time`` is summed since initialize); anything
    else is the day's time-weighted mean (``mean``). The campus's heat it couldn't reject under a limit
    (``data_center__unremoved_heat_power``) is the first use.
    """

    def __init__(self, name: str, component, config: str | None = None, *, variables: dict | None = None,
                 constituents: dict | None = None, inputs: dict[str, Callable[[float], float]] | None = None,
                 report: dict[str, str] | None = None, rtol: float = 1e-9):
        self.name = name
        self.bmi = component
        self.bmi.initialize(config)
        outs = set(self.bmi.get_output_var_names())
        ins = set(self.bmi.get_input_var_names())
        self.inputs = inputs or {}
        self.control_vars = {k: (variables or {}).get(k, v) for k, v in CONTROLS.items()
                             if (variables or {}).get(k, v) in ins}
        self.controls: dict[str, float] = {}
        self.rtol = rtol
        self.var: dict[str, str] = {}
        for role, common in COMMON.items():
            v = (variables or {}).get(role, common)
            if v in outs:
                self.var[role] = v
            elif variables and role in variables:
                raise KeyError(f"{name}: {v!r} (role {role}) is not an output of {self.bmi.get_component_name()}")
        if "withdrawal" not in self.var:
            raise KeyError(f"{name}: no withdrawal variable; map one under variables.withdrawal")
        self.constituents: dict[str, str] = {}
        for c, v in (constituents or {}).items():
            if v not in outs:
                raise KeyError(f"{name}: constituent variable {v!r} is not an output")
            self.constituents[c] = v
        if constituents is None:
            for v in outs:
                if v.startswith("water_user_return_water_") and v.endswith("__mass_flux"):
                    self.constituents[v[len("water_user_return_water_"):-len("__mass_flux")]] = v
        self.report: dict[str, tuple[str, str]] = {}
        for name, spec in (report or {}).items():
            v, how = (spec, None) if isinstance(spec, str) else (spec["variable"], spec.get("how"))
            if v not in outs:
                raise KeyError(f"{name}: reported variable {v!r} is not an output of {self.bmi.get_component_name()}")
            u = self.bmi.get_var_units(v)
            how = how or ("MW" if u in ("MW", "W") else "change" if u in ("h", "s") else "mean")
            if how not in ("MW", "change", "mean"):
                raise ValueError(f"{name}: how = {how!r} (MW, change or mean)")
            self.report[name] = (v, how)
        self.seconds = U.time_factor(self.bmi.get_time_units())
        self._factor = {}
        for role, v in self.var.items():
            u = self.bmi.get_var_units(v)
            if role in VOLUME_ROLES:
                self._factor[role] = U.factor(U.VOLUME_FLUX, u, v)
            elif role == "storage":
                self._factor[role] = U.factor(U.VOLUME, u, v)
            elif role == "power":
                self._factor[role] = U.factor(U.POWER, u, v)
            elif role == "temperature":
                self._factor[role] = U.factor(U.TEMPERATURE, u, v)
        for c, v in self.constituents.items():
            self._factor["c:" + c] = U.factor(U.MASS_FLUX, self.bmi.get_var_units(v), v)
        self._buf = np.zeros(1, dtype=np.float64)

    def _get(self, var: str) -> float:
        buf = np.zeros(self.bmi.get_grid_size(self.bmi.get_var_grid(var)),
                       dtype=self.bmi.get_var_type(var) or "float64")
        self.bmi.get_value(var, buf)
        return float(buf.sum())  # one user per component; a multi-point component is summed

    def _role(self, role: str) -> float:
        v = self.var.get(role)
        if v is None:
            return 0.0
        if role == "temperature":
            return self._get(v) + self._factor[role]
        return self._get(v) * self._factor[role]

    def set_controls(self, **controls: float):
        """Limits for the coming day (limit in m3 s-1, dry_mode 0 or 1). Controls not given go back to none."""
        for k in controls:
            if k not in self.control_vars:
                raise KeyError(f"{self.name}: the component has no {k!r} input ({CONTROLS.get(k, k)})")
        self.controls = dict(controls)

    def _apply_controls(self):
        for k, v in self.control_vars.items():
            default = np.inf if k == "limit" else 0.0
            self.bmi.set_value(v, np.array([self.controls.get(k, default)], dtype=np.float64))

    def fast_forward(self, day: int):
        """Step the component to the start of day ``day`` without recording, for a restart. Inputs are set as in a
        run, so a deterministic component reaches the same state."""
        t_end = day * U.DAY / self.seconds
        while self.bmi.get_current_time() < t_end - 1e-9 * max(1.0, abs(t_end)):
            t0 = self.bmi.get_current_time()
            for v, f in self.inputs.items():
                self.bmi.set_value(v, np.array([f(t0 * self.seconds)], dtype=np.float64))
            self.bmi.update()

    def advance_day(self, day: int) -> DailyUse:
        """Step the component to the end of day ``day`` (0-based, from the run start) and integrate the day."""
        t_end = (day + 1) * U.DAY / self.seconds
        out = DailyUse(day=day)
        energy = 0.0
        heat = 0.0
        rep = dict.fromkeys(self.report, 0.0)
        start = {n: float(self._get(v)) for n, (v, how) in self.report.items() if how == "change"}
        s0 = self._role("storage") if "storage" in self.var else 0.0
        eps = 1e-9 * max(1.0, abs(t_end))
        self._apply_controls()
        while self.bmi.get_current_time() < t_end - eps:
            t0 = self.bmi.get_current_time()
            for v, f in self.inputs.items():
                self.bmi.set_value(v, np.array([f(t0 * self.seconds)], dtype=np.float64))
            self.bmi.update()
            dt = (self.bmi.get_current_time() - t0) * self.seconds
            if dt <= 0:
                raise RuntimeError(f"{self.name}: update() did not advance time")
            w, g, c, r = (self._role(k) for k in VOLUME_ROLES)
            ds = 0.0
            if "storage" in self.var:
                s1 = self._role("storage")
                ds, s0 = s1 - s0, s1
            scale = max(w + g, c + r, 1e-12)
            if abs((w + g) * dt - (c + r) * dt - ds) > self.rtol * scale * dt:
                raise BalanceError(f"{self.name}: withdrawal + gain != consumption + return + storage change "
                                   f"at t={t0} ({w}+{g} vs {c}+{r}, ds={ds})")
            out.withdrawal += w * dt
            out.gain += g * dt
            out.consumption += c * dt
            out.returned += r * dt
            out.storage_change += ds
            for name, v in self.constituents.items():
                out.constituents[name] = out.constituents.get(name, 0.0) + self._get(v) * self._factor["c:" + name] * dt
            if "temperature" in self.var:
                heat += r * dt * self._role("temperature")
            if "power" in self.var:
                energy += self._role("power") * dt / 3600.0
            for name, (v, how) in self.report.items():
                x = float(self._get(v))
                if how == "MW":
                    rep[name] += x * (1e-6 if self.bmi.get_var_units(v) == "W" else 1.0) * dt / 3600.0
                elif how == "change":
                    rep[name] = x - start[name]
                else:
                    rep[name] += x * dt / U.DAY
        out.reported = rep
        if "temperature" in self.var and out.returned > 0:
            out.return_temperature = heat / out.returned
        if "power" in self.var:
            out.energy = energy
        return out

    def finalize(self):
        self.bmi.finalize()
