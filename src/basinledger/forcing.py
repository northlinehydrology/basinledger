"""Hourly weather for water users, from the same forcing the surface model uses.

A user component's weather inputs are set before each of its steps from a forcing source, by CSDMS standard name:

    atmosphere_bottom_air__temperature                        degC     from 2 m air temperature
    atmosphere_bottom_air_water_vapor__relative_saturation    percent  from 2 m specific humidity, temperature and
                                                                       surface pressure

Sources:
- ``nextgen``: an NGIAB run's ``forcings/forcings.nc`` (NWM retrospective or AORC, per catchment, hourly), the
  catchment named in the user's config.
- ``swatplus``: a SWAT+ weather station's daily files (``<station>.tmp``, ``.hmd``) in the run's project, made
  hourly [A]: air temperature on a cosine between the day's minimum at 03:00 and maximum at 15:00 (so the day's mean
  is SWAT+'s), relative humidity the day's value all day (SWAT+ keeps it as a fraction; here in percent).

Relative humidity is the ratio of vapor pressure, e = q p / (0.622 + 0.378 q), to saturation vapor pressure over
water, es = 611.2 exp(17.67 T / (T + 243.5)) Pa with T in degC (Bolton 1980), capped at 100%.
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import numpy as np

TEMPERATURE = "atmosphere_bottom_air__temperature"
HUMIDITY = "atmosphere_bottom_air_water_vapor__relative_saturation"


def relative_humidity(q: np.ndarray, p_pa: np.ndarray, t_c: np.ndarray) -> np.ndarray:
    e = q * p_pa / (0.622 + 0.378 * q)
    es = 611.2 * np.exp(17.67 * t_c / (t_c + 243.5))
    return np.minimum(100.0 * e / es, 100.0)


class HourlySeries:
    """Values at hourly times (seconds since the Unix epoch, UTC); value(t) is the hour containing t."""

    def __init__(self, epoch_s: np.ndarray, values: np.ndarray, run_start: dt.datetime):
        self.t0 = float(epoch_s[0])
        if np.any(np.diff(epoch_s) != 3600):
            raise ValueError("forcing must be hourly")
        self.v = values
        self.offset = run_start.replace(tzinfo=dt.timezone.utc).timestamp() - self.t0

    def __call__(self, t_run_s: float) -> float:
        i = int((t_run_s + self.offset) // 3600)
        if not 0 <= i < len(self.v):
            raise IndexError(f"forcing has no hour for {t_run_s} s after the run start")
        return float(self.v[i])


def nextgen_inputs(run_folder: Path, catchment: str, run_start: dt.datetime) -> dict[str, HourlySeries]:
    """Hourly temperature and humidity for one catchment of an NGIAB run, as user inputs."""
    import netCDF4

    with netCDF4.Dataset(Path(run_folder) / "forcings" / "forcings.nc") as d:
        ids = [str(i) for i in d["ids"][:]]
        if catchment not in ids:
            raise KeyError(f"no forcing for {catchment!r} (catchments: {', '.join(ids)})")
        k = ids.index(catchment)
        t = np.asarray(d["Time"][k, :], dtype=float)
        tk = np.asarray(d["TMP_2maboveground"][k, :], dtype=float)
        q = np.asarray(d["SPFH_2maboveground"][k, :], dtype=float)
        p = np.asarray(d["PRES_surface"][k, :], dtype=float)
    tc = tk - 273.15
    return {TEMPERATURE: HourlySeries(t, tc, run_start),
            HUMIDITY: HourlySeries(t, relative_humidity(q, p, tc), run_start)}


class DailyToHourly:
    """An hourly function of seconds since the run start from daily values keyed by date."""

    def __init__(self, days: dict, run_start: dt.datetime, hourly):
        self.days, self.start, self.hourly = days, run_start, hourly

    def __call__(self, t_run_s: float) -> float:
        when = self.start + dt.timedelta(seconds=t_run_s)
        rec = self.days.get(when.date())
        if rec is None:
            raise IndexError(f"the station has no day for {when.date()}")
        return float(self.hourly(rec, when.hour + when.minute / 60.0))


def swatplus_inputs(project: Path, station: str, run_start: dt.datetime) -> dict[str, DailyToHourly]:
    """Hourly temperature and humidity for a user from a SWAT+ station's daily files, as ``weather-sta.cli`` names
    them. Where SWAT+ generates a variable (``sim``), there is no series: humidity is then left out, and the user
    needs it as a constant input."""
    project = Path(project)
    lines = (project / "weather-sta.cli").read_text().splitlines()
    head = lines[1].split()
    row = next((ln.split() for ln in lines[2:] if ln.split() and ln.split()[0] == station), None)
    if row is None:
        raise KeyError(f"no weather station {station!r} in weather-sta.cli")
    files = dict(zip(head, row))

    def read(kind):
        name = files.get(kind, "sim")
        if name in ("sim", "null"):
            return None
        out = {}
        for ln in (project / name).read_text().splitlines()[3:]:
            f = ln.split()
            if f:
                out[dt.date(int(f[0]), 1, 1) + dt.timedelta(days=int(f[1]) - 1)] = [float(x) for x in f[2:]]
        return out

    tmp, hmd = read("tmp"), read("hmd")
    if tmp is None:
        raise ValueError(f"station {station!r}: SWAT+ generates its temperature; a user needs a measured series")

    def temperature(r, hour):
        mean, amp = 0.5 * (r[0] + r[1]), 0.5 * (r[0] - r[1])
        return mean + amp * np.cos(2 * np.pi * (hour - 15.0) / 24.0)

    out = {TEMPERATURE: DailyToHourly(tmp, run_start, temperature)}
    if hmd is not None:
        out[HUMIDITY] = DailyToHourly(hmd, run_start, lambda r, hour: min(100.0 * r[0], 100.0))
    return out
