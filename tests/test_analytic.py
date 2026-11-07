"""Glover and Balmer stream depletion (the closed form and its superposition), and Ogata and Banks transport."""

from __future__ import annotations

import math

import numpy as np
import pytest

from basinledger.analytic import glover_depletion, glover_share, ogata_banks


def test_share_limits():
    s = glover_share(np.array([0.0, 1e-6, 1.0, 1e9]), 500.0, 500.0, 0.1)
    assert s[0] == 0.0 and s[1] < 1e-12 and 0 < s[2] < 1 and s[3] == pytest.approx(1.0, abs=1e-3)
    # the closed form at one point: u = sqrt(500^2 0.1 / (4 500 30)) = 0.6455
    assert glover_share(np.array([30.0]), 500.0, 500.0, 0.1)[0] == pytest.approx(math.erfc(math.sqrt(25000 / 60000)))


def test_constant_pumping_is_the_share():
    q = np.full(400, 0.08)
    d = glover_depletion(q, 500.0, 500.0, 0.1)
    np.testing.assert_allclose(d, 0.08 * glover_share(np.arange(1, 401), 500.0, 500.0, 0.1), rtol=1e-12)


def test_superposition_is_linear_and_bounded():
    rng = np.random.default_rng(1)
    q = rng.uniform(0.03, 0.12, 900)
    d = glover_depletion(q, 500.0, 500.0, 0.1)
    assert (d >= -1e-12).all() and (d <= q.max() + 1e-12).all()
    np.testing.assert_allclose(glover_depletion(2 * q, 500.0, 500.0, 0.1), 2 * d, rtol=1e-12)
    # pumping that stops: depletion continues, then decays
    stop = np.r_[np.full(300, 0.1), np.zeros(600)]
    ds = glover_depletion(stop, 500.0, 500.0, 0.1)
    assert ds[310] > 0 and ds[-1] < 0.05 * ds[299]  # the recovery tail decays slowly (about t^-1/2)


def test_ogata_banks_closed_form_and_limits():
    x, v, D = np.array([0.0, 10.0, 100.0, 150.0, 400.0]), 0.5, 5.0
    c = ogata_banks(x, 200.0, v, D)
    direct = [0.5 * (math.erfc((xi - v * 200) / (2 * math.sqrt(D * 200)))
                     + math.exp(v * xi / D) * math.erfc((xi + v * 200) / (2 * math.sqrt(D * 200)))) for xi in x]
    np.testing.assert_allclose(c, direct, rtol=1e-12, atol=1e-300)
    assert c[0] == pytest.approx(1.0) and c[-1] < 1e-9 and (np.diff(c) < 0).all()
    # far ahead of the front at a high Peclet number, where exp(v x / D) alone would overflow: still finite, in [0, 1]
    far = ogata_banks(np.linspace(0, 5000, 101), 4000.0, 1.0, 0.5)
    assert np.isfinite(far).all() and far.min() >= 0 and far.max() <= 1 + 1e-12
    assert far[np.searchsorted(np.linspace(0, 5000, 101), 4000.0)] == pytest.approx(0.5, abs=0.02)
