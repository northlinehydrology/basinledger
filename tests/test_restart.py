"""Restarts: a run saved after a warm-up and carried on in a new process matches the unbroken run exactly."""

from __future__ import annotations

import numpy as np
import pytest

from basinledger import config as C
from basinledger.driver import Basin
from conftest import ROOT, RUNS

pytestmark = pytest.mark.slow


def config(name):
    cfg = C.load(ROOT / "examples" / "synthetic" / name)
    cfg.run_folder = RUNS
    return cfg


@pytest.mark.parametrize("example", ["surface.toml", "two_way.toml", "sfr.toml"])
def test_restart_matches_the_unbroken_run(example, swatplus):
    n, k = 500, 200
    tag = example.split(".")[0]
    a = Basin(config(example), run=f"{tag}-unbroken", days=n)
    a.run()
    ta, sa, ha = a.ledger.table(), a.series, a.mf6.head()
    a.finalize()

    w = Basin(config(example), run=f"{tag}-warm", days=n)
    w.run(k)
    w.save(RUNS / f"{tag}-state")
    w.finalize()

    b = Basin(config(example), run=f"{tag}-restarted", days=n, restart=RUNS / f"{tag}-state")
    assert b.day == k
    b.run()
    tb, sb, hb = b.ledger.table(), b.series, b.mf6.head()
    b.finalize()

    assert np.array_equal(ha, hb)
    assert set(ta) == set(tb)
    for key in ta:
        assert np.array_equal(ta[key], tb[key], equal_nan=True), key
    for key in sa:
        assert np.array_equal(sa[key], sb[key]), key
