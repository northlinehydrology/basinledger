"""Irrigation pumped from MODFLOW 6's wells (docs/design.md, "Large users and irrigation, MODFLOW first").

SWAT+ decides when and how much a field is irrigated (its decision tables and schedules, from the weather and the
soil), and MODFLOW 6 pumps it. Each day, after SWAT+'s day, Basin Ledger reads the water the day's irrigation put on
each HRU (``land_irrigation_water__volume_flow_rate``, swatplus-bmi 0.11.0: what reached the soil and what ran off)
and pumps the same volume that day from the wells ``[irrigation]`` maps to that HRU, in their shares. So the water
is pumped once, from the model that owns groundwater, and lands on SWAT+'s fields, where its return flow and losses
are SWAT+'s; the demand follows the weather and the land use of a scenario.

An HRU mapped to wells must irrigate from an outside source in SWAT+: a decision table's ``irrigate`` whose source
isn't an aquifer, channel or reservoir, ``irr_demand`` with an unlimited source, or a scheduled ``irrm``. Irrigation
from a SWAT+ channel or reservoir is surface water and already comes out of that water body, so pumping it as well
would take it twice; from SWAT+'s aquifer it is refused in two-way mode anyway (``swatplus.aquifer_withdrawals``);
through water allocation (``irr_demand`` with a water allocation source) its source isn't known here. These are
refused for a mapped HRU. An HRU irrigated from an outside source but not mapped keeps SWAT+'s water from outside,
and the ledger books it so (``irrigation:from_outside``).
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

IRRIGATION = "land_irrigation_water__volume_flow_rate"  # m3 s-1 by HRU, the day just done (swatplus-bmi 0.11.0)

OUTSIDE = "an outside source"
SOURCES = {"aqu": "SWAT+'s aquifer", "cha": "a SWAT+ channel", "sdc": "a SWAT+ channel", "res": "a SWAT+ reservoir"}
ALLOCATION = "water allocation"


def read_map(base: Path, spec) -> dict[str, list[tuple[str, float]]]:
    """``[irrigation] wells``: a list of {hru, well, share} tables, or a CSV file with the columns hru, well, share.
    An HRU's shares are positive and add up to 1; a share left out is 1 (the HRU's one well)."""
    if isinstance(spec, str):
        with open(Path(base) / spec, newline="") as f:
            rows = list(csv.DictReader(f))
        where = spec
    else:
        rows = list(spec)
        where = "[irrigation] wells"
    out: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for r in rows:
        if not r.get("hru") or not r.get("well"):
            raise ValueError(f"{where}: each row needs an hru and a well; got {dict(r)}")
        share = float(r.get("share") or 1.0)
        if share <= 0:
            raise ValueError(f"{where}: HRU {r['hru']}'s share for well {r['well']} must be positive")
        out[str(r["hru"])].append((str(r["well"]), share))
    for hru, ws in out.items():
        if len({w for w, _ in ws}) != len(ws):
            raise ValueError(f"{where}: HRU {hru} names a well twice")
        if abs(sum(s for _, s in ws) - 1.0) > 1e-9:
            raise ValueError(f"{where}: HRU {hru}'s shares add up to {sum(s for _, s in ws):g}, not 1")
    if not out:
        raise ValueError(f"{where}: no HRUs")
    return dict(out)


def _rows(path: Path) -> list[list[str]]:
    """A SWAT+ table's rows after its title and header, split on whitespace."""
    if not path.exists():
        return []
    return [ln.split() for ln in path.read_text(errors="replace").splitlines()[2:] if ln.split()]


def _schedules(project: Path) -> dict[str, tuple[list[str], bool]]:
    """management.sch: each schedule's decision tables, and whether it has a scheduled irrigation (irrm)."""
    lines = [ln.split() for ln in (project / "management.sch").read_text(errors="replace").splitlines()[2:]] \
        if (project / "management.sch").exists() else []
    out, i = {}, 0
    while i < len(lines):
        w = lines[i]
        if len(w) >= 3 and w[1].isdigit() and w[2].isdigit():
            n_ops, n_auto = int(w[1]), int(w[2])
            body = [x for x in lines[i + 1:i + 1 + n_auto + n_ops]]
            out[w[0]] = ([x[0] for x in body[:n_auto] if x], any(x and x[0] == "irrm" for x in body[n_auto:]))
            i += 1 + n_auto + n_ops
        else:
            i += 1
    return out


def _tables(project: Path) -> dict[str, set[str]]:
    """The land use decision tables (lum.dtl): each table's irrigation sources."""
    out: dict[str, set[str]] = {}
    f = project / "lum.dtl"
    if not f.exists():
        return out
    lines = f.read_text(errors="replace").splitlines()
    table = None
    for i, ln in enumerate(lines):
        w = ln.split()
        if len(w) >= 2 and w[0].lower() == "name" and w[1].lower() == "conds":  # a table's header; its name next
            table = next((x.split() for x in lines[i + 1:] if x.split()), [None])[0]
            out.setdefault(table, set())
        elif table and len(w) >= 2 and w[0] == "irrigate":
            out[table].add(SOURCES.get(w[1], OUTSIDE))
        elif table and len(w) >= 8 and w[0] == "irr_demand":
            out[table].add(OUTSIDE if w[7] == "unlim" else ALLOCATION)
    return out


def sources(project: str | Path) -> dict[str, set[str]]:
    """Each HRU's irrigation sources, by its name in hru.con: from its land use's management schedule
    (hru.con -> hru-data.hru's lu_mgt -> landuse.lum's mgt -> management.sch), the decision tables it runs and a
    scheduled irrm. An HRU that never irrigates has none."""
    project = Path(project)
    data = {r[0]: r[5] for r in _rows(project / "hru-data.hru") if len(r) > 5}   # id -> land use
    mgt = {r[0]: r[3] for r in _rows(project / "landuse.lum") if len(r) > 3}      # land use -> schedule
    scheds, tables = _schedules(project), _tables(project)
    out = {}
    for r in _rows(project / "hru.con"):
        if len(r) < 8:
            continue
        name, lu = r[1], data.get(r[7])
        auto, irrm = scheds.get(mgt.get(lu, "null"), ([], False))
        got = set().union(*[tables.get(t, set()) for t in auto]) if auto else set()
        if irrm:
            got.add(OUTSIDE)
        out[name] = got
    return out


def problems(project: str | Path, mapped: dict[str, list[tuple[str, float]]], wells: set[str]) -> tuple[list[str], list[str]]:
    """What the start-up check refuses in an irrigation map, and what it notes."""
    bad, notes = [], []
    src = sources(project)
    for hru, ws in mapped.items():
        if hru not in src:
            bad.append(f"[irrigation] HRU {hru}: no such HRU in the SWAT+ project (hru.con)")
            continue
        for w, _ in ws:
            if w not in wells:
                bad.append(f"[irrigation] HRU {hru}: well {w!r} isn't one of the [[wells]]")
        other = src[hru] - {OUTSIDE}
        if other:
            bad.append(f"[irrigation] HRU {hru} irrigates from {', '.join(sorted(other))}; pumping it from MODFLOW 6's "
                       "wells as well would take that water twice. Only irrigation from an outside source in SWAT+ "
                       "(irrigate from null, irr_demand from unlim, a scheduled irrm) is pumped")
        elif not src[hru]:
            notes.append(f"[irrigation] HRU {hru} never irrigates in SWAT+ (its land use's schedule has no irrigation)")
    outside = sorted(h for h, s in src.items() if OUTSIDE in s and h not in mapped)
    if outside:
        notes.append(f"irrigation on {len(outside)} HRU(s) not in [irrigation] ({', '.join(outside[:5])}"
                     f"{', ...' if len(outside) > 5 else ''}) stays SWAT+'s water from outside the basin")
    return bad, notes
