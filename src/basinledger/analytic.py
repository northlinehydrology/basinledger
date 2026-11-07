"""Analytical solutions, for comparison with the coupled models: stream depletion and one-dimensional transport.

Glover and Balmer (1954): a well at distance d from a straight, fully penetrating stream in a uniform, infinite
aquifer of transmissivity T and storativity S. After pumping at a constant rate starts, the share of the pumping the
stream supplies at time t is

    q(t) / Q = erfc( sqrt( d^2 S / (4 T t) ) ).

A varying pumping series is handled by superposition: each step's pumping change starts its own response. This is
the method nextgen-michigan uses for wells, and the one Michigan's Water Withdrawal Assessment Tool's depletion is
built on. Its assumptions (one straight stream, no streambed resistance, a uniform aquifer, no other boundaries) are
what a coupled groundwater model relaxes.

Ogata and Banks (1961): a solute carried at a uniform velocity v with longitudinal dispersion D through a
semi-infinite column, initially free of it, with the concentration at x = 0 held at C0 from t = 0:

    C(x, t) / C0 = 1/2 [ erfc( (x - v t) / (2 sqrt(D t)) ) + exp(v x / D) erfc( (x + v t) / (2 sqrt(D t)) ) ].

Glover, R.E., Balmer, G.G., 1954. River depletion resulting from pumping a well near a river. Eos Trans. AGU 35(3),
468-470. https://doi.org/10.1029/TR035i003p00468

Ogata, A., Banks, R.B., 1961. A solution of the differential equation of longitudinal dispersion in porous media.
U.S. Geological Survey Professional Paper 411-A. https://doi.org/10.3133/pp411A
"""

from __future__ import annotations

import math

import numpy as np


def glover_share(t_days: np.ndarray, distance_m: float, transmissivity_m2_d: float, storativity: float) -> np.ndarray:
    """The share of a constant pumping rate the stream supplies t days after it starts (0 at t <= 0)."""
    t = np.asarray(t_days, float)
    out = np.zeros_like(t)
    pos = t > 0
    u = np.sqrt(distance_m ** 2 * storativity / (4.0 * transmissivity_m2_d * t[pos]))
    out[pos] = np.array([math.erfc(x) for x in u])
    return out


def glover_depletion(pumping: np.ndarray, distance_m: float, transmissivity_m2_d: float, storativity: float,
                     dt_days: float = 1.0) -> np.ndarray:
    """Daily stream depletion (same units as ``pumping``) for a pumping series of one value per step, by
    superposition of step changes. Each step's pumping is taken as starting at the step's beginning; the depletion
    reported for a step is at its end."""
    p = np.asarray(pumping, float)
    n = len(p)
    share = glover_share(dt_days * np.arange(1, n + 1), distance_m, transmissivity_m2_d, storativity)
    dp = np.diff(np.r_[0.0, p])
    return np.convolve(dp, share)[:n]


def _erfcx(b: float) -> float:
    """exp(b^2) erfc(b), the scaled complementary error function, for b >= 0: directly below 25, and above by its
    asymptotic series (four terms; relative error below 1e-10 there)."""
    if b < 25.0:
        return math.exp(b * b) * math.erfc(b)
    u = 1.0 / (2.0 * b * b)
    return (1.0 - u + 3.0 * u * u - 15.0 * u ** 3) / (b * math.sqrt(math.pi))


def ogata_banks(x_m: np.ndarray, t_days: float, velocity_m_d: float, dispersion_m2_d: float) -> np.ndarray:
    """C / C0 at distances x >= 0 (m) downstream of the held concentration, t days after it starts. The second term
    is computed as exp(v x / D - b^2) erfcx(b), b = (x + v t) / (2 sqrt(D t)), which can't overflow: the exponent is
    -(x - v t)^2 / (4 D t)."""
    x = np.asarray(x_m, float)
    v, D, t = velocity_m_d, dispersion_m2_d, float(t_days)
    s = 2.0 * math.sqrt(D * t)
    first = np.array([math.erfc(a) for a in (x - v * t) / s])
    second = np.array([math.exp(-(xi - v * t) ** 2 / (4.0 * D * t)) * _erfcx((xi + v * t) / s) for xi in x])
    return 0.5 * (first + second)
