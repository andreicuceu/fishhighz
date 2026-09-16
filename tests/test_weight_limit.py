"""Mathematical controls for the diagnostic cumulative-weight representation."""

from decimal import Decimal, localcontext

import numpy as np
import pytest

from fishhighz.kernels.weights import _integrals, _iterate
from fishhighz.validation.weight_limit import (
    batch_payload,
    coefficients_from_log_weights,
    continuum_linearized_coefficients,
    fixed_mesh_screen,
    initial_log_weights,
    linearized_fixed_grid_limit,
    linearized_spectrum,
    relative_log_change,
    scalar_weight,
    signed_log_ratio,
    snapshot,
    trajectory,
    update_log_weights,
    validate_batch_payload,
)


def decimal_recurrence(weights, masses, variance, length, pixel, signal, count, prec):
    with localcontext() as context:
        context.prec = prec
        w = [Decimal.from_float(float(value)) for value in weights]
        r = [Decimal.from_float(float(value)) for value in masses]
        v = [Decimal.from_float(float(value)) for value in variance]
        length_d = Decimal.from_float(float(length))
        pixel_d = Decimal.from_float(float(pixel))
        signal_d = Decimal.from_float(float(signal))
        for _ in range(count):
            cumulative = Decimal(0)
            updated = []
            for mass, weight, noise in zip(r, w, v):
                cumulative += mass * weight
                if mass == 0:
                    updated.append(Decimal(0))
                elif noise == 0:
                    updated.append(Decimal(1))
                else:
                    density = cumulative * length_d / pixel_d
                    updated.append(density / (density + noise / signal_d))
            w = updated
        i1 = sum(mass * weight for mass, weight in zip(r, w))
        i2 = sum(mass * weight * weight for mass, weight in zip(r, w))
        i3 = sum(mass * weight * weight * noise for mass, weight, noise in zip(r, w, v))
        a = i2 / (length_d * i1 * i1)
        p = pixel_d * i3 / (length_d * i1 * i1)
        return w, a, p


@pytest.mark.parametrize("count", [0, 1, 3, 6, 12, 24])
@pytest.mark.parametrize("scale", [1.0, 1e-80])
def test_finite_recurrence_matches_direct_and_decimal(count, scale):
    masses = np.array([0.2, 0.7, 0.4])
    variance = np.array([0.3, 2.0, 0.8])
    length, pixel, signal, alias = 5.0, 0.7, 1.3, 2.1
    initial = np.exp(initial_log_weights(masses, variance, pixel, alias)) * scale
    log_weights = np.log(initial)
    for _ in range(count):
        log_weights = update_log_weights(
            log_weights, masses, variance, length, pixel, signal
        )
    decimal_weights, decimal_a, decimal_p = decimal_recurrence(
        initial, masses, variance, length, pixel, signal, count, 100
    )
    np.testing.assert_allclose(
        log_weights,
        [float(value.ln()) for value in decimal_weights],
        rtol=5e-12,
        atol=5e-13,
    )
    coefficients = coefficients_from_log_weights(
        masses, variance, log_weights, length, pixel
    )
    assert coefficients.A == pytest.approx(float(decimal_a), rel=5e-12)
    assert coefficients.P_pixel == pytest.approx(float(decimal_p), rel=5e-12)
    if scale == 1:
        direct, _ = _iterate(masses, variance, length, pixel, signal, alias, count)
        np.testing.assert_allclose(np.exp(log_weights), direct, rtol=5e-12, atol=0)


@pytest.mark.parametrize("d", [0.3, 1.0, 1.7])
@pytest.mark.parametrize("count", [0, 1, 3, 12, 24])
def test_scalar_closed_form_and_finite_coefficients(d, count):
    initial = 0.4
    expected = scalar_weight(count, initial, d)
    weight = initial
    for _ in range(count):
        weight = d * weight / (1 + d * weight)
    assert weight == pytest.approx(expected, rel=2e-14)
    masses = np.array([0.7])
    variance = np.array([2.0])
    state = snapshot(
        np.array([np.log(weight)]), masses, variance, 5.0, 0.8, [20.0], count
    )
    assert state.coefficients.A == pytest.approx(1 / (5.0 * 0.7))
    assert state.coefficients.P_pixel == pytest.approx(0.8 * 2 / (5.0 * 0.7))


def test_coefficient_scale_invariance_but_update_is_not_invariant():
    masses = np.array([0.4, 0.6])
    variance = np.array([0.5, 2.0])
    logs = np.log([0.2, 0.7])
    first = coefficients_from_log_weights(masses, variance, logs, 3, 0.8)
    second = coefficients_from_log_weights(masses, variance, logs - 500, 3, 0.8)
    assert first.log_A == pytest.approx(second.log_A, abs=2e-13)
    assert first.log_P_pixel == pytest.approx(second.log_P_pixel, abs=2e-13)
    update = update_log_weights(logs, masses, variance, 3, 0.8, 1.2)
    scaled_update = update_log_weights(logs - 5, masses, variance, 3, 0.8, 1.2)
    assert not np.allclose(update - update.max(), scaled_update - scaled_update.max())


def test_false_normalized_only_recurrence_is_detected():
    masses = np.array([0.3, 0.7])
    variance = np.array([0.2, 3.0])
    logs = initial_log_weights(masses, variance, 0.8, 1.5)
    exact = logs.copy()
    false = logs.copy()
    for _ in range(6):
        exact = update_log_weights(exact, masses, variance, 4, 0.8, 1.1)
        false = update_log_weights(false - false.max(), masses, variance, 4, 0.8, 1.1)
    assert np.max(np.abs((exact - exact.max()) - (false - false.max()))) > 1e-2


def test_log_coefficients_continue_past_both_historical_underflows():
    masses = np.array([1.0, 1.0])
    variance = np.array([1e-150, 1e300])
    logs = np.log([1.0, 1e-200])
    with pytest.raises(FloatingPointError):
        _integrals(masses, np.exp(logs), variance, 1, 1)
    coefficients = coefficients_from_log_weights(masses, variance, logs, 1, 1)
    assert coefficients.P_pixel == pytest.approx(1e-100, rel=5e-12)

    masses = np.array([4.9881e-297])
    logs = np.log([1.6280e-269])
    with pytest.raises(FloatingPointError):
        _integrals(masses, np.exp(logs), [1.0], 1, 1)
    coefficients = coefficients_from_log_weights(masses, [1.0], logs, 1, 1)
    assert coefficients.A == pytest.approx(1 / masses[0], rel=5e-12)
    assert coefficients.P_pixel == pytest.approx(1 / masses[0], rel=5e-12)


def test_unrepresentable_coefficients_keep_finite_logs():
    coefficients = coefficients_from_log_weights([1e-320], [1.0], [np.log(0.5)], 1, 1)
    assert np.isfinite(coefficients.log_A)
    assert not coefficients.A_available and coefficients.A is None
    assert np.isfinite(coefficients.log_P_pixel)
    assert not coefficients.P_pixel_available and coefficients.P_pixel is None


def test_coefficients_cancel_a_very_small_common_amplitude_before_summing():
    masses = np.geomspace(1e-30, 1, 1000)
    variance = np.linspace(0.2, 3, 1000)
    shape = np.linspace(-800, 0, 1000)
    first = coefficients_from_log_weights(masses, variance, shape, 4, 0.8)
    second = coefficients_from_log_weights(masses, variance, shape - 1e6, 4, 0.8)
    # The input vector has already lost a few ulps of relative shape when the
    # common offset was applied; cancellation must not amplify that loss.
    assert first.log_A == pytest.approx(second.log_A, abs=2e-10)
    assert first.log_P_pixel == pytest.approx(second.log_P_pixel, abs=2e-10)


def test_zero_support_and_zero_variance_branches():
    result = trajectory(
        [0, 0.5, 0.5],
        [9, 0, 2],
        4,
        0.8,
        1.1,
        1.5,
        checkpoints=(0, 1, 3),
        coordinates=[19, 20, 21],
    )
    for state in result.snapshots:
        assert np.isneginf(state.log_shape[0])
        assert state.log_shape[1] == 0
        assert state.coefficients.P_pixel > 0
    assert result.has_supported_zero_variance
    zero_noise = coefficients_from_log_weights([1, 2], [0, 0], [0, -1], 3, 1)
    assert zero_noise.P_pixel_exact_zero
    assert zero_noise.P_pixel_available and zero_noise.P_pixel == 0


@pytest.mark.parametrize(
    "masses,variance",
    [
        ([0, 0], [1, 2]),
        ([-1, 2], [1, 2]),
        ([1, 2], [-1, 2]),
        ([1, np.nan], [1, 2]),
        ([1], [1, 2]),
    ],
)
def test_invalid_domain_rejected(masses, variance):
    with pytest.raises(ValueError):
        trajectory(masses, variance, 2, 1, 1, 1, checkpoints=(0, 1))


def test_ordered_coordinate_and_checkpoint_validation():
    with pytest.raises(ValueError, match="checkpoints"):
        trajectory([1], [1], 2, 1, 1, 1, checkpoints=(0, 3, 3))
    with pytest.raises(ValueError, match="coordinates"):
        trajectory([1, 1], [1, 1], 2, 1, 1, 1, coordinates=[1])
    with pytest.raises(ValueError, match="increasing"):
        trajectory([1, 1], [1, 1], 2, 1, 1, 1, coordinates=[2, 1])


@pytest.mark.parametrize("count", [0, 1, 3, 6, 12, 24])
def test_constant_linearized_continuum_control(count):
    nodes, weights = np.polynomial.legendre.leggauss(max(32, count + 2))
    x = (nodes + 1) / 2
    quadrature = weights / 2
    shape = x**count
    total_mass, length, pixel, variance = 2.3, 4.1, 0.7, 1.8
    i1 = total_mass * np.sum(quadrature * shape)
    i2 = total_mass * np.sum(quadrature * shape**2)
    numeric_a = i2 / (length * i1**2)
    analytic_a, analytic_p = continuum_linearized_coefficients(
        count,
        length=length,
        total_mass=total_mass,
        pixel=pixel,
        variance=variance,
    )
    assert numeric_a == pytest.approx(analytic_a, rel=2e-13)
    assert analytic_p == pytest.approx(pixel * variance * analytic_a)


def test_continuum_control_diverges_despite_bounded_shape():
    values = [
        continuum_linearized_coefficients(
            count, length=4, total_mass=2, pixel=1, variance=1
        )[0]
        for count in (8, 16, 32, 64)
    ]
    assert all(upper > lower for lower, upper in zip(values[:-1], values[1:]))
    assert values[-1] / values[-2] > 1.9


def test_finite_grid_repeated_and_near_degenerate_spectra():
    exact = linearized_spectrum([1, 1], [2, 2], 1, 1, 1)
    assert exact[2] == 2 and exact[3] == 0
    near = linearized_spectrum([1, 1], [2, 2 * (1 + 1e-10)], 1, 1, 1)
    assert near[2] == 1
    assert near[3] == pytest.approx(1e-10, rel=2e-6)
    assert linearized_fixed_grid_limit([1, 1], [2, 2], 1, 1, 0.5) is None


def test_unique_fixed_grid_limit_matches_linearized_power_iteration():
    masses = np.array([0.2, 0.3, 0.5])
    variance = np.array([1.0, 0.8, 3.0])
    result = linearized_fixed_grid_limit(
        masses, variance, 1, 1, 0.5, coordinates=[20, 21, 22]
    )
    assert result is not None
    dominant, state = result
    assert dominant == 1
    matrix = np.tril(np.ones((3, 3))) * (0.5 / variance)[:, None] * masses
    vector = np.ones(3)
    for _ in range(1000):
        vector = matrix @ vector
        vector /= vector.max()
    assert np.exp(state.log_shape[0]) == 0
    np.testing.assert_allclose(
        np.exp(state.log_shape[1:]), vector[1:], rtol=2e-12, atol=0
    )


def test_small_nonlinear_shape_approaches_fixed_grid_limit_with_decimal_control():
    masses = np.array([0.2, 0.3, 0.5])
    variance = np.array([1.0, 0.8, 3.0])
    initial = np.exp(initial_log_weights(masses, variance, 1, 1))
    result = trajectory(
        masses,
        variance,
        1,
        1,
        0.5,
        1,
        checkpoints=(128,),
        coordinates=[20, 21, 22],
    )
    low = decimal_recurrence(initial, masses, variance, 1, 1, 0.5, 128, 100)
    high = decimal_recurrence(initial, masses, variance, 1, 1, 0.5, 128, 180)
    for first, second in zip(low[0], high[0]):
        assert abs(first / second - 1) < Decimal("1e-90")
    decimal_logs = np.array([float(value.ln()) for value in high[0]])
    state = result.snapshots[0]
    np.testing.assert_allclose(
        state.log_shape,
        decimal_logs - np.max(decimal_logs),
        rtol=0,
        atol=3e-13,
    )
    fixed = linearized_fixed_grid_limit(
        masses, variance, 1, 1, 0.5, coordinates=[20, 21, 22]
    )[1]
    np.testing.assert_allclose(state.distribution, fixed.distribution, atol=2e-15)


def test_smooth_measure_subdivision_and_iteration_limits_are_separate():
    rows = []
    for cells in (32, 64, 128):
        x = (np.arange(cells) + 0.5) / cells
        mass = np.full(cells, 1 / cells)
        log_weights = np.log(x**8)
        state = snapshot(log_weights, mass, np.ones(cells), 1, 1, x, 8)
        rows.append(state.coefficients.A)
    assert abs(rows[-1] / rows[-2] - 1) < 2e-3
    later = continuum_linearized_coefficients(
        64, length=1, total_mass=1, pixel=1, variance=1
    )[0]
    assert later > 3 * rows[-1]


def test_fixed_mesh_screen_requires_three_doublings_and_rejects_drift():
    result = trajectory(
        np.full(24, 1 / 24),
        np.linspace(1, 1.001, 24),
        1,
        1,
        0.4,
        1,
        checkpoints=(128, 256, 512, 1024),
    )
    screen = fixed_mesh_screen(result)
    assert len(screen["comparisons"]) == 3
    assert not screen["passed"]
    with pytest.raises(ValueError, match="requires checkpoints"):
        fixed_mesh_screen(
            trajectory([1], [1], 1, 1, 0.4, 1, checkpoints=(128, 256, 512))
        )


def test_relative_log_change_handles_zero_and_extreme_values():
    assert relative_log_change(-np.inf, -np.inf) == 0
    assert np.isinf(relative_log_change(-np.inf, 0))
    assert relative_log_change(-1000, -1000 + np.log1p(1e-5)) == pytest.approx(
        1e-5, rel=1e-9
    )


def test_signed_log_ratio_retains_refinement_direction():
    assert signed_log_ratio(np.log(2), np.log(3)) == pytest.approx(0.5)
    assert signed_log_ratio(np.log(3), np.log(2)) == pytest.approx(-1 / 3)
    assert signed_log_ratio(-np.inf, -np.inf) == 0


@pytest.mark.parametrize(
    "mass,variance,length,pixel",
    [
        (1e300, 1e-300, 1e-300, 1e300),
        (1e-300, 1e300, 1e300, 1e-300),
    ],
)
def test_one_cell_extreme_scalar_ratios_remain_representable(
    mass, variance, length, pixel
):
    result = trajectory(
        [mass],
        [variance],
        length,
        pixel,
        1,
        1,
        checkpoints=(0, 1),
        coordinates=[20],
    )
    assert np.exp(result.snapshots[0].log_amplitude) == pytest.approx(1 / 2)
    assert np.exp(result.snapshots[1].log_amplitude) == pytest.approx(1 / 3)
    assert result.spectral_radius == pytest.approx(1)
    for state in result.snapshots:
        assert state.coefficients.A == pytest.approx(1)
        assert state.coefficients.P_pixel == pytest.approx(1)


def test_extreme_tail_bound_remains_in_log_range():
    result = trajectory(
        [1e300],
        [1e-300],
        1e-300,
        1e300,
        0.5,
        1,
        checkpoints=(0, 1),
        coordinates=[20],
    )
    _, report = batch_payload(
        result,
        magnitudes=[20],
        masses=[1e300],
        variance=[1e-300],
        length=1e-300,
        pixel=1e300,
        signal=0.5,
        alias=1,
        bin_index=0,
        field="lya(qso)",
        order=1,
        source={"arrays": "fixture.npz", "sha256": "a" * 64},
    )
    tail = report["fixed_grid_tail"]
    assert tail["A_cauchy_lower_bound"] == pytest.approx(1)
    assert tail["P_pixel_lower_bound"] == pytest.approx(1)
    assert tail["A_cauchy_lower_bound_available"]
    assert tail["P_pixel_lower_bound_available"]


def test_diagnostic_matches_decimal_after_float64_weight_underflow():
    masses = np.array([1e-120, 2e-120, 4e-120])
    variance = np.array([1.0, 1.3, 2.0])
    length, pixel, signal, alias = 3.0, 1.0, 0.7, 1e-180
    result = trajectory(
        masses,
        variance,
        length,
        pixel,
        signal,
        alias,
        checkpoints=(24,),
        coordinates=[20, 21, 22],
    )
    initial = np.exp(initial_log_weights(masses, variance, pixel, alias))
    assert np.all(initial > 0)
    low = decimal_recurrence(initial, masses, variance, length, pixel, signal, 24, 160)
    high = decimal_recurrence(initial, masses, variance, length, pixel, signal, 24, 260)
    state = result.snapshots[0]
    for first, second in zip(low[0], high[0]):
        assert abs(first / second - 1) < Decimal("1e-145")
    decimal_logs = np.array([float(value.ln()) for value in high[0]])
    assert np.all(np.exp(decimal_logs) == 0)
    np.testing.assert_allclose(
        state.log_shape,
        decimal_logs - np.max(decimal_logs),
        rtol=0,
        atol=2e-11,
    )
    assert state.log_amplitude == pytest.approx(np.max(decimal_logs), abs=2e-11)
    assert state.coefficients.log_A == pytest.approx(float(high[1].ln()), abs=2e-11)
    assert state.coefficients.log_P_pixel == pytest.approx(
        float(high[2].ln()), abs=2e-11
    )


def test_refinement_screens_do_not_prove_limit_outcomes():
    rising_but_bounded = [n / (1 + n / 10000) for n in (16, 32, 64)]
    assert all(
        signed_log_ratio(np.log(lower), np.log(upper)) > 0.9
        for lower, upper in zip(rising_but_bounded[:-1], rising_but_bounded[1:])
    )
    assert all(value < 10000 for value in rising_but_bounded)

    false_plateau = {128: 1.0, 256: 1.0, 512: 1.0, 1024: 1.0, 2048: 2048.0}
    assert all(
        relative_log_change(np.log(false_plateau[lower]), np.log(false_plateau[upper]))
        == 0
        for lower, upper in ((128, 256), (256, 512), (512, 1024))
    )
    assert false_plateau[2048] > false_plateau[1024]


def test_constant_volterra_pixel_relation_requires_constant_variance():
    state = snapshot(
        np.log([0.25, 1.0]),
        [0.5, 0.5],
        [1.0, 4.0],
        1,
        1,
        [0.25, 0.75],
        1,
    )
    assert state.coefficients.P_pixel != pytest.approx(state.coefficients.A)


def test_decimal_precision_refinement_is_stable_after_float_underflow():
    masses = [1e-120, 2e-120, 4e-120]
    variance = [1.0, 1.3, 2.0]
    initial = [1e-180, 2e-180, 3e-180]
    low = decimal_recurrence(initial, masses, variance, 3, 1, 0.7, 24, 80)
    high = decimal_recurrence(initial, masses, variance, 3, 1, 0.7, 24, 160)
    for first, second in zip(low[0], high[0]):
        assert abs(first / second - 1) < Decimal("1e-70")
    assert abs(low[1] / high[1] - 1) < Decimal("1e-70")
    assert abs(low[2] / high[2] - 1) < Decimal("1e-70")


def tiny_batch():
    result = trajectory(
        [0.2, 0.3, 0.5],
        [0.5, 1.0, 2.0],
        4,
        0.8,
        0.4,
        1.2,
        checkpoints=(128, 256, 512, 1024),
        coordinates=[20, 21, 22],
    )
    return batch_payload(
        result,
        magnitudes=[20, 21, 22],
        masses=[0.2, 0.3, 0.5],
        variance=[0.5, 1.0, 2.0],
        length=4,
        pixel=0.8,
        signal=0.4,
        alias=1.2,
        bin_index=2,
        field="lya(qso)",
        order=16,
        source={"arrays": "fixture.npz", "sha256": "a" * 64},
    )


def tiny_expectation(arrays, report):
    return dict(
        bin=2,
        field="lya(qso)",
        order=16,
        inputs={
            name: np.array(value, copy=True)
            for name, value in arrays.items()
            if name
            in (
                "magnitudes",
                "masses",
                "variance",
                "length",
                "pixel",
                "signal",
                "alias",
            )
        },
        source=dict(report["source"]),
        attempts=[dict(row) for row in report["attempts"]],
        capped_at=report["capped_at"],
    )


def test_batch_payload_recomputes_every_checkpoint():
    arrays, report = tiny_batch()
    assert validate_batch_payload(
        arrays, report, expected=tiny_expectation(arrays, report)
    )


def test_batch_payload_rejects_self_consistent_replacement_source():
    arrays, report = tiny_batch()
    expected = tiny_expectation(arrays, report)
    replacement = trajectory(
        [0.202, 0.303, 0.505],
        [0.5, 1.0, 2.0],
        4,
        0.8,
        0.4,
        1.2,
        checkpoints=(128, 256, 512, 1024),
        coordinates=[20, 21, 22],
    )
    replaced_arrays, replaced_report = batch_payload(
        replacement,
        magnitudes=[20, 21, 22],
        masses=[0.202, 0.303, 0.505],
        variance=[0.5, 1.0, 2.0],
        length=4,
        pixel=0.8,
        signal=0.4,
        alias=1.2,
        bin_index=2,
        field="lya(qso)",
        order=16,
        source=report["source"],
    )
    with pytest.raises(ValueError, match="source masses"):
        validate_batch_payload(replaced_arrays, replaced_report, expected=expected)


@pytest.mark.parametrize("outcome", ["failed", "capped"])
def test_batch_payload_retains_truthful_incomplete_attempt(outcome):
    arrays, report = tiny_batch()
    report["attempts"].append(
        dict(iterations=2048, outcome=outcome, error="bounded fixture outcome")
    )
    report["capped_at"] = 2048
    expected = tiny_expectation(arrays, report)
    assert validate_batch_payload(arrays, report, expected=expected)


def test_identical_source_arrays_remain_valid_for_distinct_declared_population():
    arrays, report = tiny_batch()
    report["field"] = "lya(lbg)"
    expected = tiny_expectation(arrays, report)
    expected["field"] = "lya(lbg)"
    assert validate_batch_payload(arrays, report, expected=expected)


@pytest.mark.parametrize(
    "quantity,reason",
    [
        ("source", "source provenance"),
        ("tail", "tail diagnostic"),
        ("remainder", "remainder diagnostic"),
    ],
)
def test_batch_payload_rejects_mutated_derived_or_provenance_content(quantity, reason):
    arrays, report = tiny_batch()
    expected = tiny_expectation(arrays, report)
    if quantity == "source":
        report["source"]["sha256"] = "b" * 64
    elif quantity == "tail":
        report["fixed_grid_tail"]["probability"] *= 0.99
    else:
        report["nonlinear_remainder"]["maximum_log_linearized_update"] += 0.1
    with pytest.raises(ValueError, match=reason):
        validate_batch_payload(arrays, report, expected=expected)


@pytest.mark.parametrize(
    "mutation,reason",
    [
        ("field", "field identity"),
        ("order", "order identity"),
        ("coefficient", "log_A"),
        ("shape", "log_shape"),
        ("missing_attempt", "attempts"),
        ("false_verdict", "verdict"),
        ("operand", "input hash"),
    ],
)
def test_batch_payload_rejects_detached_or_mislabeled_evidence(mutation, reason):
    arrays, report = tiny_batch()
    if mutation == "field":
        report["field"] = "lya(lbg)"
    elif mutation == "order":
        report["order"] = 32
    elif mutation == "coefficient":
        arrays["log_A"][0] += 0.1
    elif mutation == "shape":
        arrays["log_shape"][0, 0] += 0.1
    elif mutation == "missing_attempt":
        report["attempts"].pop()
    elif mutation == "false_verdict":
        report["fixed_mesh_screen"]["passed"] = not report["fixed_mesh_screen"][
            "passed"
        ]
    else:
        arrays["masses"][0] *= 1.1
    with pytest.raises(ValueError, match=reason):
        validate_batch_payload(
            arrays, report, expected={"bin": 2, "field": "lya(qso)", "order": 16}
        )
