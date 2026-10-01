"""Array-only cumulative recipe adapted from lyaforecast weights.py (GPLv3).

McDonald & Eisenstein (2007); explicit finite iterations, no optimality claim.
"""

import numpy as np


def _iterate(r, variance, length, pixel, signal, alias, iterations):
    """Apply the finite cumulative forest-weight recurrence.

    Parameters
    ----------
    r : ndarray of shape (n_magnitude,)
        Source density times magnitude quadrature, in deg^-2 (km/s)^-1.
    variance : ndarray of shape (n_magnitude,)
        Dimensionless pixel-noise variance.
    length, pixel : float
        Forest length and pixel width, respectively, in km/s.
    signal : float
        Three-dimensional reference power in angular/velocity coordinates,
        in deg^2 km/s.
    alias : float
        One-dimensional reference forest power, in km/s.
    iterations : int
        Number of simultaneous updates after the initial weights.

    Returns
    -------
    weights : ndarray of shape (n_magnitude,)
        Dimensionless weights, zero outside positive density support.
    changes : ndarray of shape (iterations,)
        Maximum absolute weight change at each update.

    Notes
    -----
    Each magnitude uses its cumulative density, evaluated with the preceding
    weights. The finite update count does not imply convergence or optimality.
    """
    support_mask = r > 0
    weights = np.zeros_like(r)
    weights[support_mask] = 1 / (1 + pixel * variance[support_mask] / alias)
    changes = np.empty(iterations)

    for t in range(iterations):
        previous = weights
        density = np.cumsum(r * previous) * length / pixel
        weights = np.zeros_like(r)
        weights[support_mask] = 1 / (
            1 + variance[support_mask] / density[support_mask] / signal
        )
        changes[t] = np.max(np.abs(weights - previous))

    return weights, changes


def _integrals(r, w, variance, length, pixel):
    """Integrate cumulative weight moments and normalize forest-noise terms.

    Parameters
    ----------
    r : ndarray of shape (n_magnitude,)
        Source density times magnitude quadrature, in deg^-2 (km/s)^-1.
    w : ndarray of shape (n_magnitude,)
        Dimensionless forest weights; the argument name is retained for
        keyword compatibility.
    variance : ndarray of shape (n_magnitude,)
        Dimensionless pixel-noise variance.
    length, pixel : float
        Forest length and pixel width, respectively, in km/s.

    Returns
    -------
    first_moment, second_moment, noise_moment : ndarray of shape (n_magnitude,)
        Cumulative density times weight, weight squared, and weight squared
        times variance; each has units deg^-2 (km/s)^-1.
    aliasing_coefficient : float
        Coefficient multiplying one-dimensional forest power, in deg^2.
    pixel_power : float
        Pixel-noise power, in deg^2 km/s.

    Raises
    ------
    FloatingPointError
        If an intermediate weighted density underflows inexactly. Other
        arithmetic errors follow the caller's NumPy error policy.
    """
    # Trap individual inexact underflows before a later large multiplier or
    # another positive sample can conceal the lost contribution. Exact zeros
    # do not signal underflow. The host adds field context and raises ValueError.
    with np.errstate(under="raise"):
        weighted_density = r * w
        squared_weight_density = weighted_density * w
        noise_weighted_density = squared_weight_density * variance

    first_moment = np.cumsum(weighted_density)
    second_moment = np.cumsum(squared_weight_density)
    noise_moment = np.cumsum(noise_weighted_density)

    # Sequential divisions avoid squaring the total density.
    aliasing_coefficient = (
        second_moment[-1] / first_moment[-1] / first_moment[-1] / length
    )
    pixel_power = (
        noise_moment[-1] / first_moment[-1] / first_moment[-1] * pixel / length
    )
    return first_moment, second_moment, noise_moment, aliasing_coefficient, pixel_power
