"""Scenarios from a shared warm-up (design.md; testing G5).

One warm-up runs the basin with no water users, and its state is saved at the pause after its last day
(``Basin.save``). Each scenario then carries on from that state (``Basin(..., restart=...)``), with or without the
users and with its own changes to the config, so every case starts from the same SWAT+ and MODFLOW 6 states and
differs only from the end of the warm-up. A case with no change matches the unbroken run with no users, bit for bit
(``tests/test_scenarios.py``).

The config the warm-up uses must already hold every case's MODFLOW 6 wells and SWAT+ objects (at zero), since a
restart carries on the same models: a case changes users, rules and rates, not the models' structure.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Callable

from . import config as C
from .driver import Basin


def warm_up(cfg: C.ExampleConfig, days: int, run: str = "warmup", prepare=None, wrap_surface=None) -> Path:
    """Run ``days`` with no users and save the state; returns the state's folder."""
    b = Basin(cfg, run=run, users=False, days=None, prepare=prepare, wrap_surface=wrap_surface)
    try:
        b.run(days)
        state = b.folder / "state"
        b.save(state)
    finally:
        b.finalize()
    return state


def run_case(cfg: C.ExampleConfig, state: Path, name: str, users: bool = True,
             change: Callable[[C.ExampleConfig], None] | None = None, prepare=None, wrap_surface=None) -> Basin:
    """One scenario from the saved warm-up: ``change(cfg)`` edits a copy of the config first (users, rules, rates).
    Returns the finished Basin, with its outputs written; the caller finalizes it."""
    c = copy.deepcopy(cfg)
    if change:
        change(c)
    b = Basin(c, run=name, users=users, restart=state, prepare=prepare, wrap_surface=wrap_surface)
    b.run()
    b.write_outputs()
    return b
