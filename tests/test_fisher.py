"""Fixed-covariance Fisher assembly against analytic and dense-solve oracles."""

import numpy as np
import pytest

import fishhighz.fisher as assembly
from fishhighz.covariance import gaussian_covariance
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.fisher import factor_covariance, fisher_from_factors, fisher_matrix
from fishhighz.grids import gauss_legendre_grid
from fishhighz.kernels.fisher import _forward_substitute
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.results import FisherResult


def selection(n=1, selected=None):
    return PairSelection(
        [ObservedField(str(i), "galaxy", "model") for i in range(n)], selected
    )


def test_amplitude_normalization_and_volume():
    grid = gauss_legendre_grid([0.1, 0.2, 0.4], k_order=2, mu_order=3, h_fid=0.7)
    n = 300 * grid.q_mode  # includes fractional mode counts
    total = 2 + grid.k_flat
    derivative = 1 - grid.mu_flat
    covariance = gaussian_covariance(total[:, None], n, selection())
    jacobian = derivative[:, None, None]
    fisher = fisher_matrix(jacobian, covariance)
    expected = np.sum(n * derivative**2 / (2 * total**2))
    np.testing.assert_allclose(fisher, [[expected]], rtol=2e-15)
    doubled = fisher_matrix(
        jacobian, gaussian_covariance(total[:, None], 2 * n, selection())
    )
    np.testing.assert_allclose(doubled, 2 * fisher, rtol=2e-15)
    reg = ParameterRegistry([Parameter("A", 1, "target")])
    np.testing.assert_allclose(
        FisherResult(reg, doubled).marginalized_errors(),
        FisherResult(reg, fisher).marginalized_errors() / np.sqrt(2),
    )


def test_dense_oracle_signed_derivatives_and_triangular_solve():
    rng = np.random.default_rng(183)
    raw = rng.normal(size=(5, 4, 4))
    covariance = raw @ raw.transpose(0, 2, 1) + np.eye(4)
    jacobian = rng.normal(size=(5, 4, 3))
    expected = sum(j.T @ np.linalg.solve(c, j) for j, c in zip(jacobian, covariance))
    np.testing.assert_allclose(
        fisher_matrix(jacobian, covariance), expected, rtol=2e-14
    )
    lower = np.linalg.cholesky(covariance[0])
    actual = np.empty((4, 3))
    assert _forward_substitute(lower, jacobian[0], actual) is None
    np.testing.assert_allclose(
        actual, np.linalg.solve(lower, jacobian[0]), rtol=2e-14, atol=1e-15
    )


@pytest.mark.parametrize("chosen", [None, [(3, 1), (0, 0), (4, 2)], [(0, 1)]])
def test_five_fields_selected_order(chosen):
    rng = np.random.default_rng(76)
    fields = selection(5, chosen)
    raw = rng.normal(size=(3, 5, 5))
    total = raw @ raw.transpose(0, 2, 1) + 2 * np.eye(5)
    i, j = fields.required_pairs.T
    covariance = gaussian_covariance(total[:, i, j], [0.3, 2, 4], fields)
    jacobian = rng.normal(size=(3, len(fields.selected_pairs), 4))
    oracle = sum(j.T @ np.linalg.solve(c, j) for j, c in zip(jacobian, covariance))
    np.testing.assert_allclose(
        fisher_matrix(jacobian, covariance), oracle, rtol=2e-14, atol=1e-14
    )


def test_reuse_batching_ownership_and_no_refactorization(monkeypatch):
    covariance = np.tile([[4.0, 1.0], [1.0, 2.0]], (8, 1, 1))[::2]
    jacobian = np.arange(48.0).reshape(8, 2, 3)[::2] - 10
    covariance.flags.writeable = jacobian.flags.writeable = False
    before_c, before_j = covariance.copy(), jacobian.copy()
    factors = factor_covariance(covariance)
    expected = [fisher_matrix(jacobian * scale, covariance) for scale in (1, 2)]
    factors.flags.writeable = False

    def forbidden(*args):
        raise AssertionError("factorization attempted on reuse path")

    monkeypatch.setattr(assembly, "_cholesky", forbidden)
    monkeypatch.setattr(np.linalg, "cholesky", forbidden)
    for scale, oracle in zip((1, 2), expected):
        actual = fisher_from_factors(jacobian * scale, factors)
        np.testing.assert_allclose(actual, oracle)
        assert (
            actual.flags.owndata
            and actual.flags.c_contiguous
            and actual.dtype == np.float64
        )
    partial = sum(
        fisher_from_factors(jacobian[s], factors[s])
        for s in (slice(0, 1), slice(1, None))
    )
    np.testing.assert_allclose(partial, expected[0], rtol=2e-15)
    np.testing.assert_array_equal(covariance, before_c)
    np.testing.assert_array_equal(jacobian, before_j)
    np.testing.assert_allclose(factors @ factors.transpose(0, 2, 1), covariance)
    assert factors.flags.owndata and factors.flags.c_contiguous


def test_selected_covariance_not_parent_rank():
    full = gaussian_covariance([[4.0, -6.0, 9.0]], [2.0], selection(2))
    with pytest.raises(ValueError, match="node 1.*rank"):
        factor_covariance(np.concatenate([np.eye(3)[None], full]))
    cross = gaussian_covariance([[4.0, -6.0, 9.0]], [2.0], selection(2, [(0, 1)]))
    np.testing.assert_allclose(fisher_matrix([[[3.0]]], cross), [[0.25]])


@pytest.mark.parametrize("scale", [np.ones(2), np.array([1e-80, 1e80])])
def test_covariance_rank_tolerance_symmetry_and_units(scale):
    for delta in (1e-15, 1e-12):
        normalized = np.array([[1.0, 1 - delta], [1 - delta, 1.0]])
        covariance = normalized * scale[:, None] * scale[None, :]
        if delta == 1e-15:
            with pytest.raises(ValueError, match="numerically unresolved"):
                factor_covariance(covariance[None])
        else:
            factors = factor_covariance(covariance[None])
            np.testing.assert_allclose(
                (factors @ factors.transpose(0, 2, 1))[0], covariance, rtol=2e-15
            )
    for bad in (
        [[0.0, 0.0], [0.0, 1.0]],
        [[1.0, 2.0], [2.0, 1.0]],
        [[1.0, 0.1], [0.2, 1.0]],
    ):
        with pytest.raises(ValueError, match="node 0"):
            factor_covariance((np.array(bad) * scale[:, None] * scale[None, :])[None])


def test_roundoff_asymmetry_local_copy_and_observable_units():
    covariance = np.array([[[2.0, 0.5 + 1e-15], [0.5, 3.0]]])
    saved = covariance.copy()
    factors = factor_covariance(covariance)
    np.testing.assert_allclose(
        factors @ factors.transpose(0, 2, 1),
        (covariance + covariance.transpose(0, 2, 1)) / 2,
    )
    np.testing.assert_array_equal(covariance, saved)
    jacobian = np.array([[[1.0, -2.0], [3.0, 1.0]]])
    scaling = np.array([1e-70, 1e70])
    scaled = covariance * scaling[None, :, None] * scaling[None, None, :]
    np.testing.assert_allclose(
        fisher_matrix(jacobian * scaling[None, :, None], scaled),
        fisher_matrix(jacobian, covariance),
        rtol=2e-14,
    )


@pytest.mark.parametrize(
    "bad",
    [
        [],
        np.empty((0, 2, 2)),
        np.empty((1, 0, 0)),
        [[1.0]],
        np.ones((1, 2, 3)),
        np.ones((1, 1, 1), dtype=bool),
        np.ones((1, 1, 1), dtype=complex),
        np.ones((1, 1, 1), dtype=object),
        [[["1"]]],
        [[[np.nan]]],
        [[[np.inf]]],
    ],
)
def test_bad_covariance_and_factor_arrays(bad):
    with pytest.raises(ValueError):
        factor_covariance(bad)
    with pytest.raises(ValueError):
        fisher_from_factors([[[1.0]]], bad)


@pytest.mark.parametrize(
    "bad",
    [
        [],
        np.empty((1, 1, 0)),
        np.ones((2, 1, 1)),
        [[1.0]],
        [[[True]]],
        [[[1j]]],
        [[["1"]]],
        np.ones((1, 1, 1), dtype=object),
        [[[np.nan]]],
        [[[np.inf]]],
    ],
)
def test_bad_jacobians(bad):
    with pytest.raises(ValueError):
        fisher_from_factors(bad, [[[1.0]]])


@pytest.mark.parametrize(
    "lower",
    [[[1.0, 1e-300], [0.0, 1.0]], [[0.0, 0.0], [0.0, 1.0]], [[-1.0, 0.0], [0.0, 1.0]]],
)
def test_malformed_factors(lower):
    with pytest.raises(ValueError, match="node 0.*lower triangle"):
        fisher_from_factors(np.ones((1, 2, 1)), [lower])


def test_zero_columns_and_nonfinite_arithmetic():
    np.testing.assert_array_equal(
        fisher_matrix(np.zeros((3, 2, 4)), np.tile(np.eye(2), (3, 1, 1))),
        np.zeros((4, 4)),
    )
    with pytest.raises(ValueError, match="solve at node 0.*column 0"):
        fisher_from_factors([[[1.0]]], [[[1e-320]]])
    with pytest.raises(ValueError, match="accumulation at node 0"):
        fisher_from_factors([[[1e200]]], [[[1.0]]])
    with pytest.raises(ValueError, match="accumulation at node 1"):
        fisher_from_factors([[[1e154]], [[1e154]]], [[[1.0]], [[1.0]]])


def test_integer_inputs_normalize_to_float64():
    factors = factor_covariance([[[4, 0], [0, 9]]])
    assert factors.dtype == np.float64
    np.testing.assert_array_equal(factors, [[[2, 0], [0, 3]]])
    np.testing.assert_allclose(fisher_from_factors([[[2], [3]]], factors), [[2.0]])
