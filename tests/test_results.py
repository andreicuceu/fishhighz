"""Named information, priors once, and analytic uncertainty distinctions."""

import numpy as np
import pytest

from fishhighz.fisher import fisher_matrix
from fishhighz.parameters import (
    Parameter,
    ParameterBinding,
    ParameterRegistry,
    map_jacobian,
)
from fishhighz.results import FisherResult, combine_results, diagonal_prior


def registry(ids=("theta", "eta")):
    return ParameterRegistry(
        [
            Parameter(name, 0.0, "target" if i == 0 else "nuisance", (-10.0, 10.0), 0.1)
            for i, name in enumerate(ids)
        ]
    )


def test_analytic_uncertainties_order_fixing_and_correlation():
    result = FisherResult(registry(), [[4.0, 1.0], [1.0, 1.0]])
    expected = np.array([[1.0, -1.0], [-1.0, 4.0]]) / 3
    np.testing.assert_allclose(result.marginalized_covariance(), expected)
    np.testing.assert_allclose(result.conditional_errors(), [0.5, 1.0])
    np.testing.assert_allclose(
        result.marginalized_errors(), [1 / np.sqrt(3), 2 / np.sqrt(3)]
    )
    np.testing.assert_allclose(
        result.marginalized_covariance(["eta", "theta"]), expected[::-1, ::-1]
    )
    np.testing.assert_allclose(result.conditional_errors(["eta", "theta"]), [1.0, 0.5])
    np.testing.assert_allclose(result.correlations(), [[1.0, -0.5], [-0.5, 1.0]])
    np.testing.assert_allclose(result.marginalized_covariance(["theta"]), [[1 / 3]])
    fixed = result.fix_except(["theta"])
    np.testing.assert_allclose(fixed.marginalized_covariance(), [[0.25]])
    assert fixed.registry.parameters == (result.registry.parameters[0],)
    assert result.diagnostics.rank == 2 and np.isfinite(result.diagnostics.condition)


def test_independent_bins_bindings_and_shared_nuisance():
    reg = registry(("theta", "eta(bin1)", "eta(bin2)"))
    contributions = []
    for name in ("eta(bin1)", "eta(bin2)"):
        binding = ParameterBinding(reg, ["t", "n"], {"t": "theta", "n": name})
        jac = map_jacobian(
            np.array([[[1.0, 1.0], [1.0, -1.0]]]), binding.local_to_global, 3
        )
        result = FisherResult(reg, fisher_matrix(jac, np.eye(2)[None]))
        assert result.diagnostics.rank == 2
        with pytest.raises(ValueError, match="singular Fisher"):
            result.marginalized_errors()
        contributions.append(result)
    combined = combine_results(contributions)
    np.testing.assert_allclose(combined.total_fisher, np.diag([4.0, 2.0, 2.0]))
    np.testing.assert_allclose(
        combined.marginalized_errors(), [0.5, 1 / np.sqrt(2), 1 / np.sqrt(2)]
    )
    reg = registry()
    plus = FisherResult(reg, [[1.0, 1.0], [1.0, 1.0]])
    minus = FisherResult(reg, [[1.0, -1.0], [-1.0, 1.0]])
    for contribution in (plus, minus):
        with pytest.raises(ValueError, match="rank 1/2"):
            contribution.marginalized_errors(["theta"])
    # Each individual Schur complement is zero; marginalizing first loses
    # complementary information about the shared nuisance.
    assert 1.0 - 1.0 * 1.0 / 1.0 == 0
    np.testing.assert_allclose(
        combine_results([plus, minus]).total_fisher, 2 * np.eye(2)
    )


def test_priors_once_correlated_and_resolution():
    reg = registry()
    data = FisherResult(reg, [[1.0, 1.0], [1.0, 1.0]])
    prior = diagonal_prior(reg, {"eta": 2.0})
    np.testing.assert_array_equal(prior, [[0.0, 0.0], [0.0, 0.25]])
    posterior = combine_results([data, data], prior_fisher=prior)
    np.testing.assert_array_equal(posterior.data_fisher, 2 * data.data_fisher)
    np.testing.assert_array_equal(posterior.prior_fisher, prior)
    np.testing.assert_allclose(
        posterior.marginalized_covariance(),
        np.linalg.solve([[2.0, 2.0], [2.0, 2.25]], np.eye(2)),
    )
    one_sum = FisherResult(reg, 2 * data.data_fisher, prior_fisher=prior)
    np.testing.assert_array_equal(one_sum.total_fisher, posterior.total_fisher)
    assert posterior.diagnostics.rank == 2
    with pytest.raises(ValueError, match="data-only"):
        combine_results([posterior, data])
    correlated = np.array([[1.0, -0.5], [-0.5, 1.0]])
    result = FisherResult(reg, data.data_fisher, prior_fisher=correlated)
    np.testing.assert_allclose(
        result.marginalized_covariance(),
        np.linalg.solve([[2.0, 0.5], [0.5, 2.0]], np.eye(2)),
    )
    rank_one_prior = FisherResult(
        reg, np.zeros((2, 2)), prior_fisher=[[1.0, 1.0], [1.0, 1.0]]
    )
    assert rank_one_prior.diagnostics.rank == 1
    np.testing.assert_array_equal(
        FisherResult(reg, np.zeros((2, 2))).prior_fisher, np.zeros((2, 2))
    )  # bounds do not add priors
    fixed = result.fix_except(["eta"])
    assert fixed.data_fisher[0, 0] == 1 and fixed.prior_fisher[0, 0] == 1


@pytest.mark.parametrize(
    "matrix,rank",
    [(np.zeros((2, 2)), 0), (np.diag([4.0, 0.0]), 1), (np.ones((2, 2)), 1)],
)
def test_singular_information_null_space_and_resolution(matrix, rank):
    result = FisherResult(registry(), matrix)
    diagnostic = result.diagnostics
    assert diagnostic.rank == rank and np.isinf(diagnostic.condition)
    assert diagnostic.ids == ("theta", "eta")
    directions = diagnostic.null_directions
    np.testing.assert_allclose(directions @ directions.T, np.eye(2 - rank), atol=1e-15)
    physical = directions / diagnostic.scales[None, :]
    np.testing.assert_allclose(matrix @ physical.T, 0, atol=1e-14)
    with pytest.raises(ValueError, match="null/near-null.*scales"):
        result.marginalized_covariance()
    expected = np.full(2, np.inf)
    nonzero = matrix.diagonal() > 0
    expected[nonzero] = 1 / np.sqrt(matrix.diagonal()[nonzero])
    np.testing.assert_array_equal(result.conditional_errors(), expected)
    regular = FisherResult(registry(), matrix, prior_fisher=np.eye(2))
    assert regular.diagnostics.rank == 2
    assert np.all(np.isfinite(regular.marginalized_errors()))
    if matrix[0, 0] > 0:
        assert result.fix_except(["theta"]).diagnostics.rank == 1


def test_parameter_units_rank_and_null_coordinate_conversion():
    reg = registry()
    base = np.array([[4.0, 1.0], [1.0, 1.0]])
    scales = np.array([1e-70, 1e70])  # theta=S phi
    transformed = base * scales[:, None] * scales[None, :]
    original, result = FisherResult(reg, base), FisherResult(reg, transformed)
    assert result.diagnostics.rank == 2
    np.testing.assert_allclose(
        result.diagnostics.condition, original.diagnostics.condition
    )
    np.testing.assert_allclose(
        result.marginalized_covariance() * scales[:, None] * scales[None, :],
        original.marginalized_covariance(),
        rtol=2e-14,
    )
    jac = np.array([[[1.0, 2.0], [-1.0, 1.0]]])
    data = fisher_matrix(jac, np.eye(2)[None])
    np.testing.assert_allclose(
        fisher_matrix(jac * scales[None, None, :], np.eye(2)[None]),
        data * scales[:, None] * scales[None, :],
        rtol=2e-14,
    )
    # A raw eigenvalue cutoff would incorrectly label this diagonal rank one.
    assert FisherResult(reg, np.diag([1e-200, 1e200])).diagnostics.rank == 2
    singular = FisherResult(reg, np.ones((2, 2)) * scales[:, None] * scales[None, :])
    physical = (
        singular.diagnostics.null_directions / singular.diagnostics.scales[None, :]
    )
    np.testing.assert_allclose(
        (singular.total_fisher @ physical.T) / singular.diagnostics.scales[:, None],
        0,
        atol=2e-15,
    )
    assert "x coordinates" in singular.diagnostics.basis


def test_threshold_and_owned_arrays():
    for delta, rank in ((1e-15, 1), (1e-12, 2)):
        result = FisherResult(registry(), [[1.0, 1 - delta], [1 - delta, 1.0]])
        assert result.diagnostics.rank == rank
        if rank == 1:
            with pytest.raises(ValueError):
                result.marginalized_errors()
        else:
            assert np.all(np.isfinite(result.marginalized_errors()))
    data, prior = np.diag([3.0, 1.0]), np.eye(2)
    result = FisherResult(registry(), data, prior_fisher=prior)
    data[:] = prior[:] = 0
    np.testing.assert_array_equal(result.total_fisher, np.diag([4.0, 2.0]))
    for array in (
        result.data_fisher,
        result.prior_fisher,
        result.total_fisher,
        result.diagnostics.scales,
        result.diagnostics.eigenvalues,
    ):
        assert (
            array.dtype == np.float64
            and array.flags.owndata
            and array.flags.c_contiguous
        )
        with pytest.raises(ValueError):
            array.flat[0] = 0
    null = FisherResult(registry(), np.ones((2, 2))).diagnostics.null_directions
    with pytest.raises(ValueError):
        null[0, 0] = 0
    for array in (
        result.marginalized_covariance(),
        result.marginalized_errors(),
        result.conditional_errors(),
        result.correlations(),
    ):
        assert (
            array.flags.owndata
            and array.flags.c_contiguous
            and array.dtype == np.float64
        )


@pytest.mark.parametrize("ids", [[], ["theta", "theta"], ["missing"], "theta", [True]])
def test_bad_subsets(ids):
    result = FisherResult(registry(), np.eye(2))
    for method in (
        result.fix_except,
        result.conditional_errors,
        result.marginalized_covariance,
        result.marginalized_errors,
        result.correlations,
    ):
        with pytest.raises(ValueError):
            method(ids)


@pytest.mark.parametrize(
    "width", [0.0, -1.0, np.nan, np.inf, True, 1j, "1", 1e-200, 1e200]
)
def test_bad_prior_widths(width):
    with pytest.raises(ValueError):
        diagonal_prior(registry(), {"eta": width})


def test_metadata_mismatches_unknown_prior_and_empty_combination():
    base = registry()
    result = FisherResult(base, np.eye(2))
    for changes in (
        {"fiducial": 1.0},
        {"role": "nuisance"},
        {"bounds": (-20.0, 20.0)},
        {"step": 0.2},
    ):
        args = dict(
            id="theta", fiducial=0.0, role="target", bounds=(-10.0, 10.0), step=0.1
        )
        args.update(changes)
        altered = ParameterRegistry([Parameter(**args), base.parameters[1]])
        with pytest.raises(ValueError, match="metadata"):
            combine_results([result, FisherResult(altered, np.eye(2))])
    with pytest.raises(ValueError, match="metadata"):
        combine_results(
            [result, FisherResult(ParameterRegistry(base.parameters[::-1]), np.eye(2))]
        )
    with pytest.raises(ValueError):
        combine_results([])
    with pytest.raises(ValueError):
        diagonal_prior(base, {"unknown": 1.0})
    with pytest.raises(ValueError):
        result.fix_except(None)


@pytest.mark.parametrize(
    "bad",
    [
        [],
        [[1.0]],
        np.ones((2, 2, 1)),
        [[True, False], [False, True]],
        np.eye(2, dtype=complex),
        np.eye(2, dtype=object),
        [["1", "0"], ["0", "1"]],
        [[np.nan, 0], [0, 1]],
        [[np.inf, 0], [0, 1]],
        [[-1, 0], [0, 1]],
        [[0, 1], [1, 1]],
        [[1, 2], [2, 1]],
        [[1, 0.1], [0.2, 1]],
    ],
)
def test_invalid_information(bad):
    with pytest.raises(ValueError):
        FisherResult(registry(), bad)
    with pytest.raises(ValueError):
        FisherResult(registry(), np.eye(2), prior_fisher=bad)


def test_roundoff_symmetry_and_result_overflow():
    raw = np.array([[2.0, 0.5 + 1e-15], [0.5, 1.0]])
    before = raw.copy()
    result = FisherResult(registry(), raw)
    np.testing.assert_array_equal(raw, before)
    np.testing.assert_array_equal(result.total_fisher, result.total_fisher.T)
    reg = registry(("only",))
    with pytest.raises(ValueError, match="finite"):
        FisherResult(reg, [[1e308]], prior_fisher=[[1e308]])
    big = FisherResult(reg, [[1e308]])
    with pytest.raises(ValueError, match="nonfinite combined"):
        combine_results([big, big])
    tiny = FisherResult(reg, [[1e-320]])
    with pytest.raises(ValueError, match="nonfinite marginalized"):
        tiny.marginalized_covariance()


def test_three_parameter_psd_and_no_eager_solves(monkeypatch):
    import fishhighz.results as results_module

    def forbidden(*args):
        raise AssertionError("uncertainty solve performed eagerly")

    monkeypatch.setattr(results_module, "_forward_substitute", forbidden)
    reg = registry(("a", "b", "c"))
    data = FisherResult(reg, np.eye(3))
    assert combine_results([data, data]).diagnostics.rank == 3
    np.testing.assert_array_equal(data.conditional_errors(), np.ones(3))
    for scales in (np.ones(3), np.array([1e-50, 1.0, 1e50])):
        bad = np.array([[1.0, 0.9, 0.9], [0.9, 1.0, -0.9], [0.9, -0.9, 1.0]])
        with pytest.raises(ValueError, match="indefinite"):
            FisherResult(reg, bad * scales[:, None] * scales[None, :])


def test_symmetric_subnormal_information_is_not_averaged_away():
    smallest = np.nextafter(0.0, 1.0)
    data = np.array([[1e-320, smallest], [smallest, 1e-320]])
    result = FisherResult(registry(), data)
    np.testing.assert_array_equal(result.data_fisher, data)
    np.testing.assert_array_equal(result.total_fisher, data)
    assert result.diagnostics.rank == 2
