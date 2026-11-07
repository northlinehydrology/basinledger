"""Scenarios from a shared warm-up (scenarios.py, testing G5), on the synthetic two-way example."""

from __future__ import annotations

import numpy as np
import pytest

from basinledger import config as C
from basinledger import scenarios as SC
from basinledger.driver import Basin
from conftest import ROOT, RUNS

pytestmark = pytest.mark.slow
N, WARM = 400, 200


def config():
    cfg = C.load(ROOT / "examples" / "synthetic" / "two_way.toml")
    cfg.run_folder = RUNS
    cfg.days = N
    return cfg


def test_cases_start_from_the_same_warm_up(swatplus):
    a = Basin(config(), run="scen-unbroken", users=False, days=N)
    a.run()
    ta, ha = a.ledger.table(), a.mf6.head()
    a.finalize()

    state = SC.warm_up(config(), WARM, run="scen-warmup")
    b = SC.run_case(config(), state, "scen-base", users=False)
    tb, hb = b.ledger.table(), b.mf6.head()
    b.finalize()
    c = SC.run_case(config(), state, "scen-user", users=True)
    tc = c.ledger.table()
    c.finalize()

    # a case with no change is the unbroken run, bit for bit
    assert np.array_equal(ha, hb)
    for k in ta:
        assert np.array_equal(ta[k], tb[k], equal_nan=True), k
    # a case with the user is the same up to the end of the warm-up, and differs after it
    for k in ("recharge:received", "exchange:cha01:to_stream"):
        assert np.array_equal(ta[k][:WARM], tc[k][:WARM])
    assert not np.allclose(ta["exchange:cha01:to_stream"][WARM:], tc["exchange:cha01:to_stream"][WARM:])
    assert np.nansum(tc["route:well:pw1:delivered"][:WARM]) == 0 and np.nansum(tc["route:well:pw1:delivered"]) > 0


def test_a_wrapped_surface_is_what_the_driver_steps(swatplus):
    """wrap_surface: the driver steps the wrapper, so a study can hold or replace a value each day."""
    class Counting:
        def __init__(self, m):
            self._m, self.days = m, 0

        def update(self):
            self.days += 1
            self._m.update()

        def __getattr__(self, name):
            return getattr(self._m, name)

    b = Basin(config(), run="scen-wrapped", users=False, days=30, wrap_surface=Counting)
    b.run()
    assert b.swat.days == 30
    b.finalize()
