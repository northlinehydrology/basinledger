"""What a scenario did, from its ledger and a baseline's: where a well's water came from, and what a return does to a
stream's concentrations.

Both work on ledger tables (``Ledger.table()``, or ``ledger.csv`` read back with ``read_ledger``) of runs over the
same days.

**Capture** (``capture``). MODFLOW 6's budget closes each day, and every term of it is a ledger account:
recharge received, boundary flows in (``groundwater:boundary:*``), the release from storage, the exchange with the
coupled streams (``exchange:*:to_stream``, rivers and drains) and the wells (``route:well:*:delivered``). So between
a run with a user and the same run without it, the extra pumping must equal the extra release from storage, the
water the streams lose (less discharge, more leakage) and the extra boundary inflow, plus any change in recharge
(none in two-way mode, where SWAT+'s recharge doesn't depend on the aquifer). ``capture`` returns each part and the
residual, which a test holds to zero.

**Concentration change** (``concentration_change``). A user's return carries dissolved mass (``user:<u>:<c>:kg``).
At the reach receiving it, the day's concentration is the stream's own water at its background concentration plus
the return's mass, over the day's flow there in the run with the user. The percent change is against the
background. It needs the background concentration [A], and the receiving reach's flow, which already has the
depletion and the return in it.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


def read_ledger(path: str | Path) -> dict[str, np.ndarray]:
    """``ledger.csv`` back to a table: account -> daily values (m3 unless the account says otherwise)."""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) if r[k] not in ("", "nan") else np.nan for r in rows])
            for k in rows[0] if k not in ("day", "date")}


def _sum(t: dict, start: str = "", end: str = "") -> np.ndarray:
    keys = [k for k in t if k.startswith(start) and k.endswith(end)]
    n = len(next(iter(t.values())))
    return np.nansum([t[k] for k in keys], axis=0) if keys else np.zeros(n)


def groundwater_terms(t: dict) -> dict[str, np.ndarray]:
    """One run's daily groundwater budget from its ledger, m3, positive into the aquifer except the outflows named."""
    return {
        "recharge": _sum(t, "recharge:received"),
        "boundary": _sum(t, "groundwater:boundary:"),
        "storage_release": _sum(t, "groundwater:storage_release"),
        "to_streams": _sum(t, "exchange:", ":to_stream"),
        "pumping": _sum(t, "route:well:", ":delivered"),
    }


def capture(base: dict, case: dict) -> dict:
    """Where a case's extra pumping came from, against a baseline over the same days. Returns daily arrays (m3) and
    totals: ``pumping`` (extra), ``from_storage``, ``from_streams`` (less discharge and more leakage),
    ``from_boundary``, ``from_recharge`` and ``residual`` (pumping less the rest), and each part's share."""
    b, c = groundwater_terms(base), groundwater_terms(case)
    if len(b["pumping"]) != len(c["pumping"]):
        raise ValueError("the runs cover different days")
    parts = {
        "pumping": c["pumping"] - b["pumping"],
        "from_storage": c["storage_release"] - b["storage_release"],
        "from_streams": b["to_streams"] - c["to_streams"],
        "from_boundary": c["boundary"] - b["boundary"],
        "from_recharge": c["recharge"] - b["recharge"],
    }
    parts["residual"] = parts["pumping"] - sum(v for k, v in parts.items() if k.startswith("from_"))
    total = {k: float(np.sum(v)) for k, v in parts.items()}
    share = {k: total[k] / total["pumping"] if total["pumping"] else float("nan") for k in total if k != "pumping"}
    return {"daily": parts, "total": total, "share": share}


def concentration_change(flow_m3: np.ndarray, return_m3: np.ndarray, return_kg: np.ndarray,
                         background_kg_m3: float) -> dict:
    """At the receiving reach, each day: the concentration with the return (kg/m3) and its change against the
    background, percent. ``flow_m3`` is the reach's flow that day in the run with the return (so the return and
    any depletion are in it); the stream's own water, ``flow - return``, is at the background concentration."""
    q = np.asarray(flow_m3, float)
    own = np.clip(q - np.asarray(return_m3, float), 0.0, None)
    conc = np.where(q > 0, (background_kg_m3 * own + np.asarray(return_kg, float)) / np.where(q > 0, q, 1.0), np.nan)
    return {"concentration_kg_m3": conc, "change_percent": 100.0 * (conc / background_kg_m3 - 1.0)}
