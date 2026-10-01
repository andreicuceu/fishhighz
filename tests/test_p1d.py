"""Independent original power-law PD2013 oracle and stationary join limits."""

import math

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fishhighz.models.external import BoundParameters, evaluate_p1d
from fishhighz.models.p1d import default_p1d, p1d_floor
from fishhighz.parameters import Parameter, ParameterRegistry


def scalar_original(z, k):
    """Evaluate the original low-k-regulated P1D formula in scalar arithmetic.

    Parameters
    ----------
    z : float
        Dimensionless forest redshift.
    k : float
        Velocity wavenumber in s/km; the original formula applies its low-k
        floor.

    Returns
    -------
    power : float
        One-dimensional forest power in km/s.
    """
    slope = -2.55 - 0.28 * math.log((1 + z) / 4)
    floor = 0.009 * math.exp((-0.5 * slope - 1) / -0.1)
    scaled_wavenumber = max(k, floor) / 0.009
    return (
        math.pi
        * 0.064
        / 0.009
        * scaled_wavenumber ** (3 + slope - 0.1 * math.log(scaled_wavenumber) - 1)
        * ((1 + z) / 4) ** 3.55
    )


@pytest.mark.parametrize("z", [2.0, 3.0, 4.0])
def test_original_formula_and_floor(z):
    """Check original formula and floor.

    Parameters
    ----------
    z : float
        Dimensionless redshift test input, supplied by pytest parametrization.
    """
    floor = p1d_floor(z)
    velocity_k = np.array([0, floor / 2, floor, floor * 1.5, 0.001, 0.009, 0.03])
    power_1d = default_p1d([], z, velocity_k)
    assert_allclose(power_1d, [scalar_original(z, v) for v in velocity_k], rtol=5e-14)
    assert np.all(power_1d > 0) and power_1d[0] == power_1d[1] == power_1d[2]
    assert_allclose(
        default_p1d([], z, velocity_k[::-2]), power_1d[::-2], rtol=0, atol=0
    )
    assert p1d_floor(2) > p1d_floor(3) > p1d_floor(4)
    assert_allclose(default_p1d([], 3, [0.009]), [np.pi * 0.064 / 0.009], rtol=1e-15)


@pytest.mark.parametrize("z", [2.0, 3.0, 4.0])
def test_join_derivative_limits(z):
    """Check join derivative limits.

    Parameters
    ----------
    z : float
        Dimensionless redshift test input, supplied by pytest parametrization.
    """
    floor = p1d_floor(z)
    peak = default_p1d([], z, [floor])[0]
    deviations, slopes = [], []
    for epsilon in (1e-2, 1e-3, 1e-4):
        below, at, above = default_p1d(
            [], z, floor * np.array([1 - epsilon, 1, 1 + epsilon])
        )
        assert below == at
        deviations.append(abs(above / at - 1))
        # Dimensionless forward slope goes to zero linearly; cancellation
        # dominates below these offsets, so no machine-epsilon derivative claim.
        slopes.append(abs((above - at) / (epsilon * at)))
    assert deviations[2] < deviations[1] / 90 < deviations[0] / 8100
    assert slopes[2] < slopes[1] / 9 < slopes[0] / 81
    assert slopes[-1] < 1.1e-5
    velocity_k = 2 * floor
    step = 1e-5 * velocity_k
    numerical = np.diff(default_p1d([], z, [velocity_k - step, velocity_k + step]))[
        0
    ] / (2 * step)
    slope = -2.55 - 0.28 * np.log((1 + z) / 4)
    analytic = (
        default_p1d([], z, [velocity_k])[0]
        * ((2 + slope) - 0.2 * np.log(velocity_k / 0.009))
        / velocity_k
    )
    assert_allclose(numerical, analytic, rtol=3e-9)
    # Right curvature tends to -0.2 P/k_floor^2, left curvature is zero.
    eps = 1e-3
    power_1d = default_p1d([], z, floor * np.array([1, 1 + eps, 1 + 2 * eps]))
    assert_allclose(
        (power_1d[2] - 2 * power_1d[1] + power_1d[0]) / (eps**2 * peak),
        -0.2,
        rtol=0.004,
    )


def test_empty_binding_and_external_independence():
    """Check empty binding and external independence."""
    registry = ParameterRegistry([Parameter("a", 2.0, "target")])
    binding = BoundParameters(registry, (), {})
    velocity_k = np.array([0.0, 0.001, 0.03, 0.009])
    velocity_k.flags.writeable = False
    saved = velocity_k.copy()
    power_1d = evaluate_p1d(default_p1d, binding, [2], 3, velocity_k)
    assert_allclose(
        power_1d, evaluate_p1d(default_p1d, binding, [8], 3, velocity_k), rtol=0, atol=0
    )
    other = BoundParameters(registry, ("amp",), {"amp": "a"})
    calls = []

    def external(t, z, q):
        """Evaluate an independent synthetic one-dimensional forest spectrum.

        Parameters
        ----------
        t : ndarray of shape (n_parameters,)
            Local model parameters in the provider binding order.
        z : float or ndarray
            Dimensionless redshift.
        q : ndarray of shape (n_nodes,)
            Line-of-sight velocity wavenumbers in s/km.

        Returns
        -------
        power : ndarray of shape (n_nodes,)
            Intrinsic one-dimensional power in km/s.

        Notes
        -----
        Appends to the enclosing test call log so provider dispatch can be checked.
        """
        calls.append(1)
        return t[0] / (1 + q)

    assert_allclose(
        evaluate_p1d(external, other, [2], 3, velocity_k), 2 / (1 + velocity_k)
    )
    assert len(calls) == 1
    assert_allclose(velocity_k, saved, rtol=0, atol=0)


@pytest.mark.parametrize(
    "local,z,k",
    [
        ([1], 3, [0]),
        (0, 3, [0]),
        ([[]], 3, [0]),
        ([], [-1], [0]),
        ([], -1, [0]),
        ([], np.nan, [0]),
        ([], np.inf, [0]),
        ([], 3, []),
        ([], 3, [-1]),
        ([], 3, [np.nan]),
        ([], 3, [1j]),
        ([], 3, [[1]]),
        ([], 3, [1e308]),
        ([], 1e308, [0.009]),
    ],
)
def test_invalid(local, z, k):
    """Check invalid.

    Parameters
    ----------
    local : int or list
        Local parameter vector or mapping, supplied by pytest parametrization.
    z : int or float or list
        Dimensionless redshift test input, supplied by pytest parametrization.
    k : list
        Wavenumber test input in the convention stated by the tested function,
        supplied by pytest parametrization.
    """
    with pytest.raises(ValueError):
        default_p1d(local, z, k)


def test_mathematically_valid_negative_redshift():
    """Check mathematically valid negative redshift."""
    assert_allclose(
        default_p1d([], -0.5, [0, 0.009]),
        [scalar_original(-0.5, k) for k in [0, 0.009]],
        rtol=5e-14,
    )
