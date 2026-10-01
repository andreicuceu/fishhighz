"""Plain provider integration and strict output boundaries."""

import numpy as np
import pytest

from fishhighz.fields import ObservedField, PairSelection
from fishhighz.grids import gauss_legendre_grid
from fishhighz.models.protocols import (
    P1D,
    P3D,
    P3DJacobian,
    validate_p1d,
    validate_p3d,
    validate_p3d_jacobian,
)
from fishhighz.parameters import (
    Parameter,
    ParameterBinding,
    ParameterRegistry,
    gather_local,
    map_jacobian,
)


def test_plain_providers_join_contracts_and_fixed_grid():
    """Check plain providers join contracts and fixed grid."""
    selection = PairSelection(
        [
            ObservedField("A", "galaxy", "shared"),
            ObservedField("B", "galaxy", "shared"),
        ],
        [("B", "A")],
    )
    grid = gauss_legendre_grid([0.1, 0.4], k_order=2, mu_order=3, h_fid=0.7)
    registry = ParameterRegistry(
        [Parameter("global", 2, "target"), Parameter("unused", 0, "nuisance")]
    )
    binding = ParameterBinding(registry, ["x", "y"], {"x": "global", "y": "global"})

    def power(theta, z, k, mu, pairs):
        """Evaluate the synthetic spectrum used by the enclosing regression test.

        Parameters
        ----------
        theta : ndarray of shape (n_parameters,)
            Local model parameters in the provider binding order.
        z : float or ndarray
            Dimensionless redshift.
        k : ndarray of shape (n_nodes,)
            Comoving wavenumbers in h/Mpc.
        mu : ndarray of shape (n_nodes,)
            Dimensionless line-of-sight direction cosines.
        pairs : ndarray of int, shape (n_pairs, 2)
            Observed-field indices defining the requested spectra.

        Returns
        -------
        power : nested list with shape (n_nodes, n_pairs)
            Synthetic intrinsic power in (Mpc/h)^3 before response and noise.
        """
        return [
            [
                (2 * theta[0] + 3 * theta[1]) * (kk + mm + z) * (1 if i == j else -1)
                for i, j in pairs
            ]
            for kk, mm in zip(k, mu)
        ]

    def derivative(theta, z, k, mu, pairs):
        """Evaluate analytic local-parameter derivatives of the synthetic spectrum.

        Parameters
        ----------
        theta : ndarray of shape (n_parameters,)
            Local model parameters in the provider binding order.
        z : float or ndarray
            Dimensionless redshift.
        k : ndarray of shape (n_nodes,)
            Comoving wavenumbers in h/Mpc.
        mu : ndarray of shape (n_nodes,)
            Dimensionless line-of-sight direction cosines.
        pairs : ndarray of int, shape (n_pairs, 2)
            Observed-field indices defining the requested spectra.

        Returns
        -------
        jacobian : nested list with shape (n_nodes, n_pairs, n_parameters)
            Derivatives of synthetic intrinsic power with respect to local
            parameters; units are power divided by parameter units.
        """
        return [
            [
                [(kk + mm + z) * factor * (1 if i == j else -1) for factor in (2, 3)]
                for i, j in pairs
            ]
            for kk, mm in zip(k, mu)
        ]

    def one_dimensional(theta, z, velocity):
        """Evaluate an independent synthetic one-dimensional forest spectrum.

        Parameters
        ----------
        theta : ndarray of shape (n_parameters,)
            Local model parameters in the provider binding order.
        z : float or ndarray
            Dimensionless redshift.
        velocity : array_like of shape (n_nodes,)
            Line-of-sight velocity wavenumbers in s/km.

        Returns
        -------
        power : nested list with shape (n_nodes,)
            Intrinsic one-dimensional power in km/s.
        """
        return [theta[0] + z + v for v in velocity]

    p3d: P3D = power
    jacobian: P3DJacobian = derivative
    p1d: P1D = one_dimensional
    snapshot = {
        name: value.copy() if isinstance(value, np.ndarray) else value
        for name, value in vars(grid).items()
    }
    local = gather_local(registry.fiducials, binding.local_to_global)
    args = (local, 2.0, grid.k_flat, grid.mu_flat, selection.required_pairs)
    result = validate_p3d(p3d(*args), 6, 3)
    jac = validate_p3d_jacobian(jacobian(*args), 6, 3, 2)
    global_jac = map_jacobian(jac, binding.local_to_global, 2)
    np.testing.assert_allclose(
        global_jac[:, 1, 0], -5 * (grid.k_flat + grid.mu_flat + 2)
    )
    assert np.all(global_jac[:, :, 1] == 0)
    assert np.all(result[:, selection.selected_to_required] < 0)
    velocity = [0.001, 0.01, 0.02, 0.1]
    np.testing.assert_allclose(
        validate_p1d(p1d(local, 2.0, velocity), 4), np.array(velocity) + 4
    )
    changed = gather_local([3.0, 0.0], binding.local_to_global)
    assert not np.array_equal(result, p3d(changed, *args[1:]))
    for name, value in snapshot.items():
        np.testing.assert_array_equal(getattr(grid, name), value)
    # Providers may evaluate arbitrary paired points outside tensor quadrature.
    assert validate_p3d(
        p3d(local, 2.0, [0.2, 0.3], [0.8, 0.1], [[0, 1]]), 2, 1
    ).shape == (2, 1)


@pytest.mark.parametrize(
    "bad",
    [
        np.zeros((3, 2)),
        np.zeros((2, 1)),
        1.0,
        [[True] * 3] * 2,
        [["1"] * 3] * 2,
        np.ones((2, 3), dtype=complex),
        np.ones((2, 3), dtype=object),
        np.full((2, 3), np.nan),
        np.full((2, 3), np.inf),
    ],
)
def test_p3d_invalid(bad):
    """Check p3d invalid.

    Parameters
    ----------
    bad : int or float or ndarray or list
        Invalid input exercising the specified rejection path, supplied by
        pytest parametrization.
    """
    with pytest.raises(ValueError):
        validate_p3d(bad, 2, 3)


@pytest.mark.parametrize(
    "validator,shape",
    [(validate_p1d, (5,)), (validate_p3d, (5, 3)), (validate_p3d_jacobian, (5, 3, 2))],
)
def test_all_output_types_shapes_and_normalization(validator, shape):
    """Check all output types shapes and normalization.

    Parameters
    ----------
    validator : callable
        Output-validation callable, supplied by pytest parametrization.
    shape : tuple
        Expected or deliberately invalid array shape, supplied by pytest
        parametrization.
    """
    raw = np.full(shape, -2, dtype=np.float32, order="F")
    out = validator(raw, *shape)
    raw.flat[0] = 9
    assert np.all(out == -2)
    assert out.dtype == np.float64 and out.flags.c_contiguous and out.flags.owndata
    for bad in (
        np.zeros((*shape, 1)),
        np.ones(shape, dtype=bool),
        np.ones(shape, dtype=complex),
        np.full(shape, "1"),
        np.ones(shape, dtype=object),
        np.full(shape, np.nan),
        np.zeros(shape).T,
    ):
        if bad.shape == shape and bad.dtype.kind == "f" and np.all(np.isfinite(bad)):
            continue
        with pytest.raises(ValueError):
            validator(bad, *shape)
    assert validator(out.tolist(), *shape).shape == shape


def test_empty_jacobian_and_invalid_dimensions():
    """Check empty jacobian and invalid dimensions."""
    assert validate_p3d_jacobian(np.empty((2, 3, 0)), 2, 3, 0).shape == (2, 3, 0)
    for dimension in (True, 1.5, 0, -1):
        with pytest.raises(ValueError):
            validate_p1d([1], dimension)
