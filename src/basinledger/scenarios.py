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
import json
import warnings
from pathlib import Path
from typing import Callable

from . import config as C
from .balance import soil_water_trend
from .coupling.swatplus import StartWarning
from .driver import Basin
from .temperature import print_daily

SOIL_SETTLED_SHARE = 0.12  # [A] above what weather alone moves a settled year (10.4%, a study), below a fill (15%)


def warm_up(cfg: C.ExampleConfig, days: int, run: str = "warmup", prepare=None, wrap_surface=None,
            soil_settled_share: float = SOIL_SETTLED_SHARE) -> Path:
    """Run ``days`` with no users and save the state; returns the state's folder.

    With SWAT+, the warm-up also reports whether its soils have settled. ``warm_up.json`` in the state's folder holds
    the basin's soil water change over the warm-up's last year (all of it, if shorter) as a share of that time's
    precipitation (``balance.soil_water_trend``). A share above ``soil_settled_share`` is a warning
    (``StartWarning``), not a stop: a soil still filling or draining carries its trend into every case."""
    swat = cfg.surface == "swatplus"

    def prep(model, d):
        if model == "swatplus":
            print_daily(d, ["basin_wb"])
        if prepare:
            prepare(model, d)

    b = Basin(cfg, run=run, users=False, days=None, prepare=prep if swat else prepare, wrap_surface=wrap_surface)
    try:
        b.run(days)
        state = b.folder / "state"
        b.save(state)
    finally:
        b.finalize()
    if swat:
        trend = soil_water_trend(b.folder / "swatplus", min(days, 365))
        trend |= {"limit": soil_settled_share, "settled": bool(trend["share"] <= soil_settled_share)}
        (state / "warm_up.json").write_text(json.dumps(trend, indent=1))
        if not trend["settled"]:
            warnings.warn(f"the warm-up's soils haven't settled: soil water changed {trend['soil_water_change_mm']:+.0f}"
                          f" mm over its last {trend['days']} days, {trend['share']:.0%} of their precipitation "
                          f"(limit {soil_settled_share:.0%}); every case carries the trend, so warm up longer or "
                          "start the soils wetter", StartWarning, stacklevel=2)
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
