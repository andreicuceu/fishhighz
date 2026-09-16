"""Intrinsic PD2013 P1D, adapted from lyaforecast analytic_p1d_PD2013.py.

Scientific formula and constants: Palanque-Delabrouille et al. (2013), with
lyaforecast's redshift-dependent stationary low-k floor. The reference routine
is distributed under GPLv3; see the repository LICENSE. No license change.
Mathematical validity outside the fit calibration range implies no empirical
accuracy claim. Response, noise and comoving conversion belong to the consumer.
"""

import numpy as np

from .._arrays import real_array, scalar


def _shape(z):
    z = scalar(z, "redshift")
    if z <= -1:
        raise ValueError("require 1+z > 0")
    evolution = np.log1p(z) - np.log(4.0)
    slope = -2.55 - 0.28 * evolution
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        floor = 0.009 * np.exp((-0.5 * slope - 1) / -0.1)
    if not np.isfinite(floor) or floor <= 0:
        raise ValueError("P1D floor is not positive representable float64")
    return evolution, slope, floor


def p1d_floor(z):
    """Return the stationary low-k floor (s/km) at scalar finite z>-1."""
    return float(_shape(z)[2])


def default_p1d(theta_local, z, k_parallel_velocity):
    """Return intrinsic km/s power (node,) for a zero-local-parameter callable.

    theta_local must have shape (0,). Velocity k must be finite real nonnegative
    nonempty 1D, including zero; slices and read-only arrays retain their order.
    The plateau joins with continuous zero first k derivative; the second
    derivative need not be continuous. No external P3D is called or inferred.
    """
    local = real_array(theta_local, "theta_local")
    if local.shape != (0,):
        raise ValueError("default_p1d requires zero local parameters, shape (0,)")
    k = real_array(k_parallel_velocity, "k_parallel_velocity")
    if k.ndim != 1 or not k.size or np.any(k < 0):
        raise ValueError("require nonempty 1D nonnegative velocity wavenumbers")
    evolution, slope, floor = _shape(z)
    # log(q)-log(k0) avoids overflowing q/k0 at large finite k.
    u = np.log(np.maximum(k, floor)) - np.log(0.009)
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        power = (np.pi * 0.064 / 0.009) * np.exp(
            (2 + slope) * u - 0.1 * u**2 + 3.55 * evolution
        )
    if not np.all(np.isfinite(power)) or np.any(power <= 0):
        raise ValueError("P1D output is not positive representable float64")
    return power
