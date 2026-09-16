"""Analytic integration and fixed-coordinate ownership checks."""

import numpy as np
import pytest

from fishhighz.grids import IntegrationGrid, gauss_legendre_grid


def test_polynomial_moments_mode_volume_and_order():
    edges = [0.1, 0.23, 0.8]
    grid = gauss_legendre_grid(edges, k_order=3, mu_order=4, h_fid=0.67)
    assert (grid.n_k, grid.n_mu) == (6, 4)
    for degree in range(6):
        expected = (0.8 ** (degree + 1) - 0.1 ** (degree + 1)) / (degree + 1)
        np.testing.assert_allclose(grid.w_k @ grid.k**degree, expected, rtol=2e-14)
    for degree in range(8):
        np.testing.assert_allclose(
            grid.w_mu @ grid.mu**degree, 1 / (degree + 1), rtol=2e-14
        )
    assert np.all((grid.k[:3] > 0.1) & (grid.k[:3] < 0.23))
    assert np.all((grid.k[3:] > 0.23) & (grid.k[3:] < 0.8))
    assert np.all((grid.mu > 0) & (grid.mu < 1))
    expected_modes = (0.8**3 - 0.1**3) / (6 * np.pi**2)
    np.testing.assert_allclose(grid.q_mode.sum(), expected_modes, rtol=2e-14)
    np.testing.assert_allclose((1200 * grid.q_mode).sum(), 1200 * expected_modes)
    toy = (10 * grid.k_flat + grid.mu_flat**2).reshape(grid.n_k, grid.n_mu)
    np.testing.assert_array_equal(toy, 10 * grid.k[:, None] + grid.mu[None, :] ** 2)
    np.testing.assert_array_equal(
        grid.weights.reshape(6, 4), np.outer(grid.w_k, grid.w_mu)
    )
    edges[0] = 0.05
    assert grid.k_min == 0.1


def custom(**changes):
    args = dict(
        k=[1.0, 2.0],
        w_k=[0.5, 0.5],
        mu=[0.0, 0.5, 1.0],
        w_mu=[1 / 6, 4 / 6, 1 / 6],
        k_min=1,
        k_max=2,
        h_fid=0.7,
    )
    args.update(changes)
    return IntegrationGrid(**args)


def test_custom_rule_ownership_and_reference_h():
    k = np.array([1.0, 2.0])
    weights = [0.5, 0.5]
    grid = custom(k=k, w_k=weights)
    k[0], weights[0] = 1.2, 0.9
    assert grid.k.tolist() == [1, 2]
    assert grid.w_k.tolist() == [0.5, 0.5]
    np.testing.assert_allclose(grid.k @ grid.w_k, 1.5)
    for degree in range(4):
        np.testing.assert_allclose(grid.mu**degree @ grid.w_mu, 1 / (degree + 1))
    for value in vars(grid).values():
        if isinstance(value, np.ndarray):
            assert (
                value.dtype == np.float64
                and value.flags.owndata
                and value.flags.c_contiguous
            )
            with pytest.raises(ValueError):
                value.flat[0] = 0
    other_h = custom(h_fid=0.5)
    np.testing.assert_array_equal(grid.k, other_h.k)
    np.testing.assert_array_equal(grid.q_mode, other_h.q_mode)


@pytest.mark.parametrize(
    "changes",
    [
        dict(k=[]),
        dict(k=[1, 1]),
        dict(k=[2, 1]),
        dict(k=[[1, 2]]),
        dict(k=[1, np.nan]),
        dict(k=[0.9, 2]),
        dict(k=[1, 2.1]),
        dict(mu=[-0.1, 0.5, 1]),
        dict(mu=[0, 0.5, 1.1]),
        dict(mu=[0, 0, 1]),
        dict(w_k=[1]),
        dict(w_k=[-0.5, 1.5]),
        dict(w_k=[0, 1]),
        dict(w_k=[0.4, 0.5]),
        dict(w_mu=[1, 1, 1]),
        dict(w_mu=[np.inf, 1, 1]),
        dict(w_k=[[0.5, 0.5]]),
        dict(k_min=0),
        dict(k_max=1),
        dict(h_fid=0),
        dict(h_fid=np.nan),
        dict(h_fid=True),
    ],
)
def test_custom_errors(changes):
    with pytest.raises(ValueError):
        custom(**changes)


@pytest.mark.parametrize(
    "edges,ko,mo",
    [
        ([0, 1], 2, 2),
        ([1], 2, 2),
        ([1, 1, 2], 2, 2),
        ([2, 1], 2, 2),
        ([1, np.inf], 2, 2),
        ([1, 2], 0, 2),
        ([1, 2], 2, -1),
        ([1, 2], True, 2),
        ([1, 2], 2, 2.0),
    ],
)
def test_factory_errors(edges, ko, mo):
    with pytest.raises(ValueError):
        gauss_legendre_grid(edges, k_order=ko, mu_order=mo, h_fid=0.7)


def test_unrepresentable_interior_nodes():
    with pytest.raises(ValueError, match="interior"):
        gauss_legendre_grid(
            [1.0, np.nextafter(1.0, 2.0)], k_order=1, mu_order=1, h_fid=0.7
        )
