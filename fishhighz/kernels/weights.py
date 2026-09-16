"""Array-only cumulative recipe adapted from lyaforecast weights.py (GPLv3).

McDonald & Eisenstein (2007); explicit finite iterations, no optimality claim.
"""

import numpy as np


def _iterate(r, variance, length, pixel, signal, alias, iterations):
    support = r > 0
    w = np.zeros_like(r)
    w[support] = 1 / (1 + pixel * variance[support] / alias)
    changes = np.empty(iterations)
    for t in range(iterations):
        previous = w
        density = np.cumsum(r * previous) * length / pixel
        w = np.zeros_like(r)
        w[support] = 1 / (1 + variance[support] / density[support] / signal)
        changes[t] = np.max(np.abs(w - previous))
    return w, changes


def _integrals(r, w, variance, length, pixel):
    # Trap individual inexact underflows before a later large multiplier or
    # another positive sample can conceal the lost contribution. Exact zeros
    # do not signal underflow. The host adds field context and raises ValueError.
    with np.errstate(under="raise"):
        rw = r * w
        rww = rw * w
        rwwv = rww * variance
    i1 = np.cumsum(rw)
    i2 = np.cumsum(rww)
    i3 = np.cumsum(rwwv)
    # Sequential divisions avoid squaring the total density.
    a = i2[-1] / i1[-1] / i1[-1] / length
    p = i3[-1] / i1[-1] / i1[-1] * pixel / length
    return i1, i2, i3, a, p
