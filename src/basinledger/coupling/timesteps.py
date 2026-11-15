"""A MODFLOW 6 model's time discretization rewritten to one-day steps, in a run's copy (docs/design.md, "Extending an
existing MODFLOW 6 model": stress periods).

Basin Ledger steps MODFLOW 6 one day at a time, with SWAT+. A model built for a study usually has longer steps:
monthly or seasonal stress periods of one step, or of a few steps growing by a factor. ``daily_steps`` rewrites the
TDIS file of a run's copy so that every period keeps its place and its length and is cut into one-day steps:

- **The stress periods stay.** Their number and order are unchanged, so every package's period blocks, time series,
  observations and output control apply to the same periods as before; only TDIS changes. Each rewritten row keeps
  the original after it (``was ...``), which MODFLOW 6 reads past.
- **A period's length becomes whole days.** Where periods aren't whole days (30.4375-day months, or a model in hours
  with a period of 36 hours), each period's end moves to the nearest day boundary, measured from the simulation's
  start, so no error builds up over the run: no end moves more than 12 hours. How far they moved is reported. A
  period too short to keep a whole day of its own is refused.
- **Adaptive time stepping is left out** (the ATS6 line commented out), since every step is a day.
- **Steady periods** are cut into days like the others: each day solved without storage, with that day's exchange.
- **Output control** keeps its words: ``ALL``, ``FIRST`` and ``LAST`` mean what they did, but ``STEPS`` and
  ``FREQUENCY`` count the new daily steps, and a model saving every step now writes every day.

A model whose steps are all one day already is left untouched. ``start_date`` reads TDIS's START_DATE_TIME, which the
driver compares with SWAT+'s start.
"""

from __future__ import annotations

import datetime as dt
import math
import re
from pathlib import Path

from .modflow6 import DAYS

KNOWN = ", ".join(repr(k) for k in DAYS)


def tdis_file(ws: str | Path) -> Path:
    """The simulation's TDIS file, from mfsim.nam's timing block."""
    ws = Path(ws)
    block = None
    for line in (ws / "mfsim.nam").read_text().splitlines():
        words = line.split()
        if not words or words[0].startswith("#"):
            continue
        key = words[0].upper()
        if key == "BEGIN":
            block = words[1].upper() if len(words) > 1 else None
        elif key == "END":
            block = None
        elif block == "TIMING" and key == "TDIS6" and len(words) > 1:
            return ws / words[1].strip("'\"")
    raise ValueError(f"{ws / 'mfsim.nam'}: no TDIS6 file in its timing block")


def _read(path: Path):
    """TDIS's options (upper-case keys), and its perioddata rows as (line index, perlen, nstp, tsmult)."""
    lines = path.read_text().splitlines()
    options: dict[str, tuple[int, list[str]]] = {}
    rows = []
    block = None
    for i, line in enumerate(lines):
        words = line.split()
        if not words or words[0].startswith(("#", "!", "//")):
            continue
        key = words[0].upper()
        if key == "BEGIN":
            block = words[1].upper() if len(words) > 1 else None
        elif key == "END":
            block = None
        elif block == "OPTIONS":
            options[key] = (i, words[1:])
        elif block == "PERIODDATA":
            if key in ("OPEN/CLOSE", "INTERNAL"):
                raise ValueError(f"{path.name}: period data read from {key}; Basin Ledger rewrites it in the file")
            rows.append((i, float(words[0].replace("d", "e").replace("D", "e")), int(words[1]),
                         float(words[2].replace("d", "e").replace("D", "e"))))
    return lines, options, rows


def _unit(path: Path, options, time_units: str | None) -> str:
    own = options.get("TIME_UNITS", (0, ["unknown"]))[1][0].lower()
    if own in DAYS:
        return own  # a stated unit that differs is refused when MODFLOW 6 starts (Modflow6)
    if own != "unknown" and own != "undefined":
        raise ValueError(f"MODFLOW 6 model in {own} (TDIS TIME_UNITS): a one-day step isn't a whole number of them; "
                         f"Basin Ledger converts {', '.join(DAYS)}")
    if time_units not in DAYS:
        raise ValueError(f"MODFLOW 6 model's time units are undefined (TDIS TIME_UNITS): state them in the config, "
                         f"[modflow6] time_units = one of {KNOWN}")
    return time_units


def daily_steps(ws: str | Path, time_units: str | None = None) -> list[str]:
    """Rewrite the TDIS of the MODFLOW 6 simulation in ``ws`` (a run's copy: it is changed in place) so every time
    step is one day, each stress period kept. ``time_units``: the model's, where its TDIS leaves them undefined.
    Returns what was changed, as notes for the start-up (empty when every step was a day already)."""
    path = tdis_file(ws)
    lines, options, rows = _read(path)
    unit = _unit(path, options, time_units)
    per_day = 1.0 / DAYS[unit]  # the model's time units in a day
    ats = "ATS6" in options
    t, ends, new = 0.0, [0], []
    for k, (_, perlen, nstp, tsmult) in enumerate(rows, start=1):
        if perlen <= 0.0 or nstp < 1:
            raise ValueError(f"{path.name}: period {k} has length {perlen} and {nstp} steps")
        t += perlen / per_day
        end = math.floor(t + 0.5)  # each end to the nearest day from the start, so no error builds up
        days = end - ends[-1]
        if days < 1:
            raise ValueError(f"{path.name}: period {k} is {perlen / per_day:g} days, too short to keep a whole day "
                             f"once periods end on day boundaries; Basin Ledger steps one day at a time")
        ends.append(end)
        new.append(days)
    daily = all(abs(p / per_day - n) < 1e-9 and n == d and (n == 1 or s == 1.0)
                for (_, p, n, s), d in zip(rows, new))
    if daily and not ats:
        return []

    out = list(lines)
    for (i, perlen, nstp, tsmult), days in zip(rows, [] if daily else new):
        out[i] = f"  {days * per_day!r} {days} 1.0  was {' '.join(lines[i].split())}"
    if ats:
        i = options["ATS6"][0]
        out[i] = f"# {lines[i].strip()}  (adaptive time stepping left out by Basin Ledger: every step is a day)"
    head = "# Basin Ledger: every time step one day (coupling/timesteps.py, daily_steps), each period kept"
    path.write_text("\n".join([head] + out) + "\n")

    steps = sum(n for _, _, n, _ in rows)
    growing = sum(1 for _, _, n, s in rows if n > 1 and s != 1.0)
    was = f"{steps:,} steps" + (f", growing in {growing} periods" if growing else "")
    notes = [] if daily else [(f"MODFLOW 6's time steps rewritten to one day each in the run's copy: {len(rows)} "
                               f"stress periods, {ends[-1]:,} days (they had {was}); each period's input stays as "
                               f"given")]
    moved = [abs(e - s) for e, s in zip(ends[1:], _starts(rows, per_day)) if abs(e - s) > 1e-9]
    if moved:
        notes.append(f"MODFLOW 6's stress periods aren't whole days: {len(moved)} of {len(rows)} period ends moved to "
                     f"the nearest day, at most {max(moved) * 24:.3g} hours")
    if ats:
        notes.append("MODFLOW 6's adaptive time stepping (ATS6) left out: every step is a day")
    return notes


def _starts(rows, per_day):
    """Each period's end in days from the simulation's start, as the model has it."""
    t, out = 0.0, []
    for _, perlen, _, _ in rows:
        t += perlen / per_day
        out.append(t)
    return out


def start_date(ws: str | Path) -> dt.date | None:
    """TDIS's START_DATE_TIME as a date, or None where it is missing or isn't an ISO date."""
    _, options, _ = _read(tdis_file(ws))
    words = options.get("START_DATE_TIME", (0, []))[1]
    m = re.match(r"\d{4}-\d{2}-\d{2}", words[0]) if words else None
    try:
        return dt.date.fromisoformat(m.group(0)) if m else None
    except ValueError:
        return None
