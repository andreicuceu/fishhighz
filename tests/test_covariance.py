"""Analytic normalization, physical validity, and independent Wick checks."""

import numpy as np
import pytest

from fishhighz.covariance import combine_observed_power, gaussian_covariance
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.grids import gauss_legendre_grid
from fishhighz.kernels.covariance import _gaussian_covariance_kernel


def selection(n=2, chosen=None):
    return PairSelection(
        [ObservedField(f"field_({i})", "galaxy", "shared") for i in range(n)], chosen
    )


def pack(dense, selected):
    i, j = selected.required_pairs.T
    return dense[:, i, j]


def test_one_field_noise_scaling_and_grid():
    chosen = selection(1)
    signal, noise = np.array([[2.0], [0.0], [3.0]]), np.array([[1.0], [4.0], [0.5]])
    counts = np.array([0.25, 2.0, 7.5])
    total = combine_observed_power(signal, noise)
    result = gaussian_covariance(total, counts, chosen)
    np.testing.assert_allclose(
        result[:, 0, 0], 2 * (signal[:, 0] + noise[:, 0]) ** 2 / counts
    )
    np.testing.assert_allclose(
        gaussian_covariance(total, 3 * counts, chosen), result / 3
    )
    np.testing.assert_allclose(
        gaussian_covariance(5 * total, counts, chosen), result * 25
    )
    grid = gauss_legendre_grid([0.1, 0.3], k_order=2, mu_order=3, h_fid=0.7)
    modes = 1000 * grid.q_mode
    result = gaussian_covariance(np.full((6, 1), 4.0), modes, chosen)
    np.testing.assert_allclose(result[:, 0, 0], 32 / modes)


def test_two_field_analytic_signed_and_uncorrelated():
    total = np.array([[4.0, -3.0, 9.0], [4.0, 0.0, 9.0]])
    n = np.array([0.5, 3.0])
    expected = (
        np.array(
            [
                [[32, -24, 18], [-24, 45, -54], [18, -54, 162]],
                [[32, 0, 0], [0, 36, 0], [0, 0, 162]],
            ]
        )
        / n[:, None, None]
    )
    np.testing.assert_allclose(gaussian_covariance(total, n, selection()), expected)
    cross = gaussian_covariance(total, n, selection(chosen=[(0, 1)]))
    np.testing.assert_allclose(cross[:, 0, 0], [45 / 0.5, 36 / 3])


def test_selection_order_five_fields_and_unused_field():
    rng = np.random.default_rng(14)
    lower = rng.normal(size=(3, 5, 5))
    dense = lower @ lower.transpose(0, 2, 1)
    full = selection(5)
    counts = np.array([1.0, 2.0, 3.0])
    covariance = gaussian_covariance(pack(dense, full), counts, full)
    assert covariance.shape == (3, 15, 15)
    subset = selection(5, [(3, 1), (0, 0), (2, 1)])  # field 4 unused
    indices = [
        list(map(tuple, full.selected_pairs)).index(tuple(p))
        for p in subset.selected_pairs
    ]
    actual = gaussian_covariance(pack(dense, subset), counts, subset)
    np.testing.assert_array_equal(actual, covariance[:, indices][:, :, indices])
    assert np.max(subset.required_pairs) == 3
    for wrong in (np.ones((3, 9)), np.ones((3, 11))):
        with pytest.raises(ValueError, match="required columns"):
            gaussian_covariance(wrong, counts, subset)


def test_supplied_observed_noise_only_added_once():
    # Supplied signal already contains any desired response; no amplitude
    # adjustment is appropriate. Negative noise cross powers are explicit.
    signal = np.array([[4.0, -1.0, 9.0], [2.0, 0.5, 3.0]])
    noise = np.array([[1.0, -0.2, 2.0], [0.5, 0.0, 0.8]])
    expected_total = np.array([[5.0, -1.2, 11.0], [2.5, 0.5, 3.8]])
    actual = combine_observed_power(signal, noise)
    np.testing.assert_array_equal(actual, expected_total)
    np.testing.assert_array_equal(
        gaussian_covariance(actual, [2, 3], selection()),
        gaussian_covariance(expected_total, [2, 3], selection()),
    )
    # Only the total is subject to physical validation.
    np.testing.assert_array_equal(combine_observed_power([[-1.0]], [[2.0]]), [[1.0]])


@pytest.mark.parametrize(
    "total", [[[4.0, -6.0, 9.0]], [[0.0, 0.0, 9.0]], [[0.0, 0.0, 0.0]]]
)
def test_valid_singular_matrices(total):
    result = gaussian_covariance(total, [2.0], selection())
    assert np.all(np.isfinite(result))
    assert np.linalg.matrix_rank(result[0]) < 3
    assert np.linalg.eigvalsh(result[0])[0] >= -1e-12
    if total[0][0]:
        cross = gaussian_covariance(total, [2.0], selection(chosen=[(0, 1)]))
        assert cross[0, 0, 0] == 36
        autos = gaussian_covariance(total, [2.0], selection(chosen=[(0, 0), (1, 1)]))
        assert autos[0, 0, 0] > 0


@pytest.mark.parametrize(
    "packed,message",
    [
        ([-1.0, 0.0, 1.0], "negative auto"),
        ([0.0, 1e-200, 1.0], "zero auto"),
        ([1.0, 1.1, 1.0], "normalized correlation"),
    ],
)
def test_invalid_two_fields_node_and_id(packed, message):
    with pytest.raises(ValueError, match=f"node 1.*field.*{message}") as error:
        gaussian_covariance([[1.0, 0.0, 1.0], packed], [1.0, 1.0], selection())
    assert "field_(" in str(error.value)


@pytest.mark.parametrize("scales", [[1.0, 1.0, 1.0], [1e-60, 1.0, 1e60]])
def test_full_psd_check_with_disparate_amplitudes(scales):
    scales = np.array(scales)
    valid = np.array([[1.0, 0.2, -0.3], [0.2, 1.0, 0.1], [-0.3, 0.1, 1.0]])
    invalid = np.array([[1.0, 0.9, 0.9], [0.9, 1.0, -0.9], [0.9, -0.9, 1.0]])
    scale = scales[:, None] * scales[None, :]
    chosen = selection(3)
    result = gaussian_covariance(pack((valid * scale)[None], chosen), [1.0], chosen)
    assert np.all(np.isfinite(result))
    with pytest.raises(ValueError, match="node 1.*not PSD.*eigenvalue.*tolerance"):
        gaussian_covariance(
            pack(np.array([valid * scale, invalid * scale]), chosen), [1.0, 1.0], chosen
        )
    # An invalid weak two-field sub-block must not hide behind a strong auto.
    invalid[0, 1] = invalid[1, 0] = 1.2
    with pytest.raises(ValueError, match="node 0.*normalized correlation"):
        gaussian_covariance(pack((invalid * scale)[None], chosen), [1.0], chosen)


def test_roundoff_tolerance_does_not_regularize():
    total = [[1.0, 1.0 + 1e-15, 1.0]]
    result = gaussian_covariance(total, [1.0], selection())
    assert result[0, 0, 1] == 2 * total[0][1]  # input cross power preserved
    with pytest.raises(ValueError):
        gaussian_covariance([[1.0, 1.0 + 1e-10, 1.0]], [1.0], selection())


def test_independent_dense_wick_expression():
    rng = np.random.default_rng(771)
    for rank in (2, 4):
        factors = rng.normal(size=(4, 4, rank))
        dense = factors @ factors.transpose(0, 2, 1)
        chosen = selection(4)
        counts = np.array([0.1, 1.0, 2.0, 7.0])
        actual = gaussian_covariance(pack(dense, chosen), counts, chosen)
        # Independent full fourth-moment tensor, then select its matrix entries.
        wick = (
            np.einsum("nim,njl->nijml", dense, dense)
            + np.einsum("nil,njm->nijml", dense, dense)
        ) / counts[:, None, None, None, None]
        i, j = chosen.selected_pairs.T
        expected = wick[:, i[:, None], j[:, None], i[None, :], j[None, :]]
        np.testing.assert_allclose(actual, expected, rtol=2e-15, atol=2e-14)
        np.testing.assert_array_equal(actual, actual.transpose(0, 2, 1))
        eigenvalues = np.linalg.eigvalsh(actual)
        # Dense eigensolver roundoff scales with matrix size and spectral norm.
        tolerance = (
            64 * np.finfo(float).eps * len(i) * np.max(np.abs(eigenvalues), axis=1)
        )
        assert np.all(eigenvalues[:, 0] >= -tolerance)


def test_ownership_noncontiguous_kernel_and_batches():
    chosen = selection()
    original = np.tile([4.0, -3.0, 9.0, 0.0, 0.0, 0.0], (8, 1))
    powers = original[::2, :3]
    counts = np.arange(1.0, 9.0)[::2]
    before = original.copy()
    maps = {
        name: getattr(chosen, name).copy()
        for name in ("im", "jn", "in_", "jm", "required_pairs", "selected_pairs")
    }
    actual = gaussian_covariance(powers, counts, chosen)
    assert (
        actual.dtype == np.float64
        and actual.flags.owndata
        and actual.flags.c_contiguous
    )
    output = np.empty_like(actual)
    assert (
        _gaussian_covariance_kernel(
            powers, counts, chosen.im, chosen.jn, chosen.in_, chosen.jm, output
        )
        is None
    )
    np.testing.assert_array_equal(output, actual)
    batches = [
        gaussian_covariance(powers[:1], counts[:1], chosen),
        gaussian_covariance(powers[1:], counts[1:], chosen),
    ]
    np.testing.assert_array_equal(np.concatenate(batches), actual)
    combined = combine_observed_power(powers, np.zeros_like(powers))
    assert combined.flags.owndata and combined.flags.c_contiguous
    combined[:] = 0
    np.testing.assert_array_equal(original, before)
    for name, value in maps.items():
        np.testing.assert_array_equal(getattr(chosen, name), value)


@pytest.mark.parametrize(
    "bad",
    [
        1.0,
        [1.0],
        [],
        np.empty((0, 1)),
        np.ones((2, 1, 1)),
        [[True]],
        [[1j]],
        [["1"]],
        np.ones((1, 1), dtype=object),
        [[np.nan]],
        [[np.inf]],
    ],
)
def test_invalid_power_arrays(bad):
    with pytest.raises(ValueError):
        gaussian_covariance(bad, [1.0], selection(1))
    with pytest.raises(ValueError):
        combine_observed_power(bad, [[1.0]])
    with pytest.raises(ValueError):
        combine_observed_power([[1.0]], bad)


@pytest.mark.parametrize(
    "bad",
    [
        1.0,
        [[1.0]],
        [],
        [1.0, 2.0],
        [0.0],
        [-1.0],
        [np.nan],
        [np.inf],
        [True],
        [1j],
        ["1"],
        np.array([1.0], dtype=object),
    ],
)
def test_invalid_counts(bad):
    with pytest.raises(ValueError):
        gaussian_covariance([[1.0]], bad, selection(1))


def test_shapes_completeness_overflow_and_no_broadcasting():
    with pytest.raises(ValueError, match="matching shapes"):
        combine_observed_power(np.ones((2, 3)), np.ones((2, 1)))
    with pytest.raises(ValueError, match="selection"):
        gaussian_covariance([[1.0]], [1.0], None)
    with pytest.raises(ValueError, match="nonfinite observed-power sum.*node 0"):
        combine_observed_power([[1e308]], [[1e308]])
    for power, counts in (([[1e200]], [1.0]), ([[1.0]], [1e-320])):
        with pytest.raises(
            ValueError, match="nonfinite covariance arithmetic at node 0"
        ):
            gaussian_covariance(power, counts, selection(1))
    broken = selection()
    object.__setattr__(broken, "required_pairs", np.array([[0, 0], [1, 1]]))
    with pytest.raises(ValueError, match="every canonical pair"):
        gaussian_covariance([[1.0, 1.0]], [1.0], broken)


def test_extreme_auto_scales_cross_only():
    # Normalizing via a product of autos would lose range in related cases;
    # sqrt/axis divisions retain a meaningful dimensionless correlation.
    chosen = selection(chosen=[(0, 1)])
    actual = gaussian_covariance([[1e-300, 0.5, 1e300]], [1.0], chosen)
    np.testing.assert_allclose(actual, [[[1.25]]])
    with pytest.raises(ValueError, match="normalized correlation"):
        gaussian_covariance([[1e-300, 2.0, 1e300]], [1.0], chosen)
