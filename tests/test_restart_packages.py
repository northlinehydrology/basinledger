"""Restarts with the model's own advanced packages (docs/design.md, roadmap item 4): a lake's stage (LAK) and a
multi-aquifer well's head (MAW) are state MODFLOW 6 carries from day to day, as the heads are. A restart saves them
and sets them back, so a restarted run is the unbroken run's to the last digit (tier 2): the synthetic foreign model
with its own pond (LAK, rain and evaporation, left as the model's) and a farm's MAW well, idle, screened through both
layers.
"""

from __future__ import annotations

import warnings

import numpy as np
from conftest import RUNS
from test_existing import config


def run(cfg, name, days, restart=None, until=None):
    from basinledger.driver import Basin
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Basin(cfg, run=name, days=days, restart=restart)
        try:
            if until is not None:
                b.run(until)
                b.save(RUNS / f"{name}-state")
            else:
                b.run()
            stage = b.mf6.get("XNEWPAK", b.mf6.model, "LAK_POND").copy()
            well = b.mf6.get("HEAD", b.mf6.model, "MAW_FARM").copy()
            return b.ledger.table(), b.mf6.head(), stage, well
        finally:
            b.finalize()


def test_a_restart_carries_the_models_lake_and_well(swatplus):
    cfg = lambda: config(lake=True, maw=True)  # noqa: E731
    ta, ha, sa, wa = run(cfg(), "restart-packages-unbroken", 60)
    run(cfg(), "restart-packages-warm", 60, until=30)
    tb, hb, sb, wb = run(cfg(), "restart-packages-restarted", 60, restart=RUNS / "restart-packages-warm-state")
    assert np.array_equal(sa, sb) and np.array_equal(wa, wb)
    assert np.array_equal(ha, hb)
    for key in ta:
        assert np.array_equal(ta[key], tb[key], equal_nan=True), key
