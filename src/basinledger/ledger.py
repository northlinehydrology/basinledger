"""The ledger: every volume and mass that crosses between models, by day, and the checks that it balances.

Each transfer is entered twice, as the sending side reported it and as the receiving side took it. A day's check
fails if any pair differs by more than its tolerance, if an account that must stay zero isn't, or if a user's own
balance doesn't close. A pair's receiving side may be several accounts summed (taken + shortfall; delivered + lost).

Tolerances: transfers between float64 models (MODFLOW 6, the users, the chain) hold to a relative 1e-9. SWAT+ keeps
4-byte reals, so a transfer into or out of SWAT+ holds to its rounding, a relative 1e-6.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

FLOAT32 = 1e-6


class Ledger:
    def __init__(self, rtol: float = 1e-9, atol: float = 1e-6):
        self.rtol, self.atol = rtol, atol
        self.rows: dict[str, list[float]] = defaultdict(list)
        self.pairs: list[tuple[str, tuple[str, ...], float, str | None]] = []
        self.zeros: list[tuple[str, float]] = []
        self.days = 0

    def pair(self, sent: str, received: str | list[str], rtol: float | None = None, bound: str | None = None):
        """Declare that ``sent`` equals ``received`` (or the sum of several accounts) every day. ``bound``: an account
        entered each day with how far the two may differ for a known reason (the precision of the files the models
        read, say), added to the tolerance that day."""
        rec = (received,) if isinstance(received, str) else tuple(received)
        self.pairs.append((sent, rec, self.rtol if rtol is None else rtol, bound))

    def zero(self, account: str, atol: float | None = None):
        """Declare an account that must be zero every day (a flux a coupling has switched off)."""
        self.zeros.append((account, self.atol if atol is None else atol))

    def enter(self, account: str, value: float):
        col = self.rows[account]
        while len(col) < self.days:
            col.append(np.nan)
        if len(col) > self.days:  # a second entry the same day adds to the first
            col[-1] += float(value)
        else:
            col.append(float(value))

    def today(self, account: str) -> float:
        col = self.rows.get(account, [])
        return col[-1] if len(col) > self.days else 0.0

    def close_day(self):
        bad = []
        for a, recs, rtol, bound in self.pairs:
            x = self.today(a)
            y = sum(self.today(b) for b in recs)
            extra = abs(self.today(bound)) if bound else 0.0
            if abs(x - y) > self.atol + extra + rtol * max(abs(x), abs(y)):
                bad.append(f"{a}={x!r} {'+'.join(recs)}={y!r}")
        for a, atol in self.zeros:
            if abs(self.today(a)) > atol:
                bad.append(f"{a}={self.today(a)!r} (must be zero)")
        self.days += 1
        for col in self.rows.values():
            while len(col) < self.days:
                col.append(np.nan)
        if bad:
            raise AssertionError(f"day {self.days - 1}: " + "; ".join(bad))

    def table(self) -> dict[str, np.ndarray]:
        return {k: np.asarray(v) for k, v in self.rows.items()}

    def totals(self) -> dict[str, float]:
        return {k: float(np.nansum(v)) for k, v in self.rows.items()}

    def write_csv(self, path: Path, dates: list):
        cols = sorted(self.rows)
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["day", "date", *cols])
            for d in range(self.days):
                w.writerow([d, dates[d], *(repr(self.rows[c][d]) for c in cols)])

    # restarts
    def state(self) -> dict:
        return {"days": self.days, "rows": {k: list(v) for k, v in self.rows.items()}}

    def load(self, state: dict):
        self.days = state["days"]
        self.rows = defaultdict(list, {k: list(v) for k, v in state["rows"].items()})
