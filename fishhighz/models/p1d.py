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
    """Compute the PD2013 redshift evolution, slope, and stationary cutoff.

    Parameters
    ----------
    z : float
        Finite dimensionless redshift satisfying z > -1.

    Returns
    -------
    evolution : float
        ln((1+z)/4), dimensionless.
    slope : float
        Redshift-dependent logarithmic slope coefficient.
    floor : float
        Positive stationary low-wavenumber cutoff in s/km.

    Raises
    ------
    ValueError
        If redshift or the represented cutoff is invalid.
    """
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
    """Return the stationary low-wavenumber cutoff of the PD2013 prescription.

    Parameters
    ----------
    z : float
        Finite dimensionless redshift satisfying z > -1.

    Returns
    -------
    cutoff : float
        Velocity wavenumber in s/km below which the spectrum is constant.

    Raises
    ------
    ValueError
        If redshift or the represented cutoff is invalid.
    """
    return float(_shape(z)[2])


def default_p1d(theta_local, z, k_parallel_velocity):
    """Evaluate intrinsic one-dimensional forest power with the adopted cutoff.

    Parameters
    ----------
    theta_local : array_like of shape (0,)
        Empty parameter vector; this prescription has no free local parameters.
    z : float
        Finite dimensionless forest redshift satisfying z > -1.
    k_parallel_velocity : array_like of shape (n_node,)
        Finite nonnegative velocity wavenumbers in s/km, including zero.

    Returns
    -------
    power : ndarray of shape (n_node,)
        Positive intrinsic power in km/s, preserving query order.

    Raises
    ------
    ValueError
        If inputs violate the callable contract or output power is not positive
        and representable.

    Notes
    -----
    The low-k plateau joins with a continuous zero first derivative; the
    second derivative need not be continuous. Instrumental response, noise,
    and comoving conversion belong to the consumer. No P3D model is inferred.
    """
    local = real_array(theta_local, "theta_local")
    if local.shape != (0,):
        raise ValueError("default_p1d requires zero local parameters, shape (0,)")
    velocity_wavenumber = real_array(k_parallel_velocity, "k_parallel_velocity")
    if (
        velocity_wavenumber.ndim != 1
        or not velocity_wavenumber.size
        or np.any(velocity_wavenumber < 0)
    ):
        raise ValueError("require nonempty 1D nonnegative velocity wavenumbers")
    evolution, slope, floor = _shape(z)
    # log(q)-log(k0) avoids overflowing q/k0 at large finite k.
    log_wavenumber_ratio = np.log(np.maximum(velocity_wavenumber, floor)) - np.log(
        0.009
    )
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        power = (np.pi * 0.064 / 0.009) * np.exp(
            (2 + slope) * log_wavenumber_ratio
            - 0.1 * log_wavenumber_ratio**2
            + 3.55 * evolution
        )
    if not np.all(np.isfinite(power)) or np.any(power <= 0):
        raise ValueError("P1D output is not positive representable float64")
    return power
