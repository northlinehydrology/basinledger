"""Weather written into SWAT+'s gage records (``[weather]``; swatplus-bmi 0.6's ``set_weather``).

A series replaces a gage's daily records from the first day SWAT+ hasn't read, so the run is the one a weather file
with those values would give: SWAT+ makes the mean temperature, PET, heat units and the rest from it, as it always
does. A change (``[weather.change]``) scales or shifts the series, or the project's own records where no series is
given, for a scenario: wetter, drier, warmer.

    [weather]                                 # optional, SWAT+ only
    precipitation = "weather/pcp.csv"         # each a CSV: a date column, then one column per gage, named as the
    temperature_max = "weather/tmax.csv"      #   gage's file in the project (Imsilpcp.pcp), in the file's units:
    temperature_min = "weather/tmin.csv"      #   mm/day, degC, MJ m-2 day-1, a fraction (humidity), m s-1
    # solar, humidity, wind

    [weather.change]                          # optional: a scenario's change
    precipitation = 1.10                      # a factor (solar, humidity and wind too), or 12 monthly factors
    temperature = 2.0                         # degC added to the maximum and minimum, or 12 monthly values

Where a series comes from is a study's business (station records, ERA5, a climate model). ``from_grid`` maps
gridded daily values to gages, by the nearest cell to each gage or, for each gage, the area-weighted mean of the
cells nearest the HRUs that use it (``hru_points``).

Limits (swatplus-bmi 0.6):
- SWAT+ reads a day's weather before the pause that comes before it, and precipitation a day earlier still. So the
  first day of a run (two for precipitation) comes from the files, and a case started from a saved warm-up takes
  its own weather from the day after its first (two days after, for precipitation). ``apply`` returns the first day
  set for each variable.
- Only inside the period each gage's file covers (SWAT+ generates the days outside it), and only daily gages.
- Gages are per file, not per HRU: a series changes a gage's values, not which gage a station uses. Finer weather
  needs more gages in the project's files, one per cell for gridded data.
"""

from __future__ import annotations

import csv
import datetime as dt
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PCP = "atmosphere_water__precipitation_leq-volume_flux"
TMAX = "atmosphere_bottom_air__max_of_temperature"
TMIN = "atmosphere_bottom_air__min_of_temperature"
SLR = "land_surface_radiation~incoming~shortwave__energy_flux"
HMD = "atmosphere_bottom_air_water~vapor__relative_saturation"
WND = "atmosphere_bottom_air_flow__speed"

# config key: the BMI variable, the factor from the file's unit to the variable's, the gage file's kind, the column
# in that file (temperature files hold the maximum and the minimum)
KINDS = {
    "precipitation": (PCP, 1.0 / 86_400_000.0, "pcp", 0),    # mm/day -> m s-1
    "temperature_max": (TMAX, 1.0, "tmp", 0),                # degC
    "temperature_min": (TMIN, 1.0, "tmp", 1),
    "solar": (SLR, 1e6 / 86_400.0, "slr", 0),                # MJ m-2 day-1 -> W m-2
    "humidity": (HMD, 1.0, "hmd", 0),                        # a fraction
    "wind": (WND, 1.0, "wnd", 0),                            # m s-1
}
FACTORS = {"precipitation", "solar", "humidity", "wind"}     # a change multiplies these
CHANGES = FACTORS | {"temperature"}                          # and shifts temperature


@dataclass
class Series:
    """Daily values, one row per day from ``start`` and one column per gage, in the gage file's units."""
    start: dt.date
    gages: list[str]
    values: np.ndarray

    @property
    def end(self) -> dt.date:
        return self.start + dt.timedelta(days=len(self.values) - 1)

    def dates(self) -> list[dt.date]:
        return [self.start + dt.timedelta(days=i) for i in range(len(self.values))]

    def between(self, first: dt.date, last: dt.date) -> "Series":
        i, j = (first - self.start).days, (last - self.start).days + 1
        if i < 0 or j > len(self.values):
            raise ValueError(f"the series runs {self.start} to {self.end}; {first} to {last} is needed")
        return Series(first, self.gages, self.values[i:j])

    def columns(self, gages: list[str]) -> "Series":
        missing = [g for g in gages if g not in self.gages]
        if missing:
            raise KeyError(f"no series for gages {missing} (the series has {self.gages})")
        return Series(self.start, list(gages), self.values[:, [self.gages.index(g) for g in gages]])


def read_csv(path: str | Path) -> Series:
    """A CSV with a ``date`` column (YYYY-MM-DD, every day in order) and one column per gage."""
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    head, body = rows[0], [r for r in rows[1:] if r]
    if head[0] != "date":
        raise ValueError(f"{path}: the first column must be date")
    dates = [dt.date.fromisoformat(r[0]) for r in body]
    if any((b - a).days != 1 for a, b in zip(dates, dates[1:])):
        raise ValueError(f"{path}: dates must be every day, in order")
    return Series(dates[0], head[1:], np.array([[float(x) for x in r[1:]] for r in body]))


def write_csv(path: str | Path, s: Series):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["date", *s.gages])
        for d, row in zip(s.dates(), s.values):
            w.writerow([d.isoformat(), *(repr(float(v)) for v in row)])


def gage_records(project: str | Path, gage: str, column: int = 0) -> Series:
    """A gage file's daily records (year, day of year, values), as SWAT+ 62.0.1 reads them; -99 is NaN (missing,
    which SWAT+ generates)."""
    lines = (Path(project) / gage).read_text().splitlines()[3:]
    rows = [ln.split() for ln in lines if ln.split()]
    start = dt.date(int(rows[0][0]), 1, 1) + dt.timedelta(days=int(rows[0][1]) - 1)
    v = np.array([float(r[2 + column]) for r in rows])
    v[v <= -99.0] = np.nan
    return Series(start, [gage], v[:, None])


def project_series(project: str | Path, kind: str, gages: list[str]) -> Series:
    """The project's own records for a variable's gages, over the days they share."""
    col = KINDS[kind][3]
    each = [gage_records(project, g, col) for g in gages]
    start, end = max(s.start for s in each), min(s.end for s in each)
    return Series(start, list(gages), np.hstack([s.between(start, end).values for s in each]))


def _monthly(spec, dates: list[dt.date]) -> np.ndarray:
    if isinstance(spec, (int, float)):
        return np.full(len(dates), float(spec))
    if len(spec) != 12:
        raise ValueError(f"a change is one value or 12 monthly values, not {len(spec)}")
    return np.array([float(spec[d.month - 1]) for d in dates])


def changed(kind: str, s: Series, change: dict) -> Series:
    """The series with a scenario's change: a factor for precipitation, solar, humidity (capped at 1) and wind; a
    shift for the temperatures."""
    group = "temperature" if kind.startswith("temperature") else kind
    if group not in change:
        return s
    by_day = _monthly(change[group], s.dates())[:, None]
    v = s.values + by_day if group == "temperature" else s.values * by_day
    if kind == "humidity":
        v = np.minimum(v, 1.0)
    return Series(s.start, s.gages, v)


def apply(model, project: str | Path, spec: dict, last: dt.date) -> dict[str, dict]:
    """Write ``spec``'s series and change into a SWAT+ model's gage records, from the first day each can be set to
    ``last`` or the series' end, whichever comes first; returns, for each variable set, the first and last days set
    and the series written (file units), for the users' weather, the run's check that it has weather for each day
    (``Basin.step``) and its summary."""
    names = model.get_output_var_names()
    change = spec.get("change", {})
    out = {}
    for kind, (var, factor, _, _) in KINDS.items():
        if kind not in spec and not (change and ("temperature" if kind.startswith("temperature") else kind) in change):
            continue
        if var not in names:
            raise ValueError(f"[weather] {kind}: the project reads no daily gage file for it (SWAT+ generates it)")
        gages = model.get_grid_object_names(model.get_var_grid(var))
        base = read_csv(spec[kind]) if kind in spec else project_series(project, kind, gages)
        first = model.first_weather_day(var)
        if first > last:  # a restart at the run's end
            continue
        s = changed(kind, base.columns(gages), change).between(first, min(last, base.end))
        if np.isnan(s.values).any():
            raise ValueError(f"[weather] {kind}: missing values (-99) between {first} and {last}; SWAT+ generates "
                             "those days, so give a series that has them")
        model.set_weather(var, first, s.values * factor)
        out[kind] = {"from": first, "to": s.end, "series": s}
    return out


# ---------- gridded data to gages ----------

def _distance(lon1, lat1, lon2, lat2):
    """Great-circle distance (radians) between points in degrees."""
    p1, p2, dl = np.radians(lat1), np.radians(lat2), np.radians(lon2 - lon1)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * np.arcsin(np.sqrt(a))


def nearest(lon, lat, cell_lon, cell_lat) -> np.ndarray:
    """The nearest cell to each point (longitude and latitude in degrees)."""
    lon, lat = np.atleast_1d(lon), np.atleast_1d(lat)
    return np.array([int(np.argmin(_distance(x, y, cell_lon, cell_lat))) for x, y in zip(lon, lat)])


def from_grid(start: dt.date, values: np.ndarray, cell_lon, cell_lat, gages: list[str], gage_lon=None,
              gage_lat=None, points: dict | None = None) -> Series:
    """Gridded daily values (days, cells) as a series per gage. Each gage takes its nearest cell, or, with
    ``points`` (gage -> (lon, lat, area) of the HRUs that use it; ``hru_points``), the area-weighted mean of the
    cells nearest those HRUs."""
    values = np.asarray(values, float)
    cols = []
    for k, g in enumerate(gages):
        if points is None:
            cols.append(values[:, nearest(gage_lon[k], gage_lat[k], cell_lon, cell_lat)[0]])
        else:
            lon, lat, area = (np.asarray(x, float) for x in points[g])
            cols.append(values[:, nearest(lon, lat, cell_lon, cell_lat)] @ (area / area.sum()))
    return Series(start, list(gages), np.column_stack(cols))


def hru_points(project: str | Path, kind: str) -> dict[str, tuple]:
    """For each gage file of a kind ("pcp", "tmp", ...), the HRUs whose station reads it: (lon, lat, area m2), from
    hru.con and weather-sta.cli."""
    project = Path(project)
    lines = (project / "weather-sta.cli").read_text().splitlines()
    head = lines[1].split()
    station = {r[0]: dict(zip(head, r)) for r in (ln.split() for ln in lines[2:]) if r}
    lines = (project / "hru.con").read_text().splitlines()
    head = [h.lower() for h in lines[1].split()]
    out: dict[str, list] = {}
    for r in (ln.split() for ln in lines[2:]):
        if not r:
            continue
        h = dict(zip(head, r))
        gage = station[h["wst"]].get(kind, "sim")
        if gage in ("sim", "null"):
            continue
        out.setdefault(gage, []).append((float(h["lon"]), float(h["lat"]), float(h["area"]) * 1e4))
    return {g: tuple(np.array(c) for c in zip(*rows)) for g, rows in out.items()}
