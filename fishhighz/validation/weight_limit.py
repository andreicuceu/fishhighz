"""Range-safe diagnosis of the existing cumulative forest-weight recurrence.

This module is validation-only.  It neither changes nor supplies a fallback for
the production weighting implementation in :mod:`fishhighz.weights`.
"""

import hashlib
from dataclasses import dataclass

import numpy as np

_LOG_MAX = float(np.log(np.finfo(np.float64).max))
_LOG_MIN_SUBNORMAL = float(np.log(np.nextafter(np.float64(0), np.float64(1))))


def _vector(value, name, *, nonnegative=False):
    """Validate a one-dimensional diagnostic input.

    Parameters
    ----------
    value : array_like, shape (n,)
        Input numerical vector in the units required by its diagnostic.
    name : str
        Quantity or record label used in diagnostics.
    nonnegative : bool
        Whether to reject negative entries. Default is ``False``.

    Returns
    -------
    array : ndarray, shape (n,)
        Finite nonempty float64 view or converted array in the input units.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    array = np.asarray(value, dtype=np.float64)
    if array.ndim != 1 or not array.size or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a nonempty finite 1D array")
    if nonnegative and np.any(array < 0):
        raise ValueError(f"{name} must be nonnegative")
    return array


def _positive(value, name):
    """Validate a positive scalar physical input.

    Parameters
    ----------
    value : float
        Scalar to check, in the units required by the diagnostic.
    name : str
        Quantity or record label used in diagnostics.

    Returns
    -------
    value : float
        Finite strictly positive value in the input units.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _log_ratio(numerator, denominator):
    """Return log(numerator / denominator) without forming the ratio.

    Parameters
    ----------
    numerator : float
        Positive numerator.
    denominator : float
        Positive denominator in the same units.

    Returns
    -------
    log_ratio : float
        Natural logarithm of the dimensionless ratio.
    """
    return float(np.log(numerator) - np.log(denominator))


def _logsumexp(log_values):
    """Sum positive terms directly in logarithmic arithmetic.

    Parameters
    ----------
    log_values : array_like, shape (n,)
        Natural logarithms of terms in a positive sum.

    Returns
    -------
    log_sum : float
        Natural logarithm of the sum; minus infinity for empty input.
    """
    if not len(log_values):
        return float("-inf")
    return float(np.logaddexp.reduce(log_values))


def _ordinary_value(log_value, *, exact_zero=False):
    """Convert a logarithmic quantity only when float64 can represent it.

    Parameters
    ----------
    log_value : float
        Natural logarithm of the positive quantity to represent.
    exact_zero : bool
        Whether the quantity is known to be exactly zero, independently of its
        logarithm. Default is ``False``.

    Returns
    -------
    value : float or None
        Exponentiated quantity, or None if it is not representable.
    available : bool
        Whether a finite positive value, or a declared exact zero, is
        represented.
    """
    if exact_zero:
        return 0.0, True
    if not np.isfinite(log_value) or not _LOG_MIN_SUBNORMAL <= log_value <= _LOG_MAX:
        return None, False
    value = float(np.exp(log_value))
    return (value, value > 0 and np.isfinite(value))


@dataclass(frozen=True)
class NoiseCoefficients:
    """Scale-cancelled coefficients and their range-safe logarithms."""

    log_A: float
    A: float | None
    A_available: bool
    log_P_pixel: float
    P_pixel: float | None
    P_pixel_available: bool
    P_pixel_exact_zero: bool


@dataclass(frozen=True)
class WeightSnapshot:
    """Amplitude, shape and coefficient diagnostics at one update count."""

    iterations: int
    log_amplitude: float
    log_shape: np.ndarray
    distribution: np.ndarray
    coefficients: NoiseCoefficients
    mean_coordinate: float
    sigma_coordinate: float
    log_concentration: float
    log_effective_measure: float
    log_dropped_distribution_bound: float


@dataclass(frozen=True)
class WeightTrajectory:
    """Finite trajectory and fixed-grid linearized spectrum diagnostics."""

    snapshots: tuple[WeightSnapshot, ...]
    spectral_radius: float
    log_spectral_radius: float
    dominant_multiplicity: int
    relative_spectral_gap: float
    has_supported_zero_variance: bool


def initial_log_weights(masses, variance, pixel, alias):
    """Return the exact log representation of the production initial weights.

    Parameters
    ----------
    masses : array_like, shape (n_magnitude,)
        Nonnegative source-density quadrature masses, in deg^-2 (km/s)^-1.
    variance : array_like, shape (n_magnitude,)
        Nonnegative dimensionless pixel-noise variance at each magnitude node.
    pixel : float
        Positive pixel width in km/s.
    alias : float
        Positive one-dimensional auxiliary forest power in km/s.

    Returns
    -------
    log_weights : ndarray, shape (n_magnitude,)
        Logarithmic initial weights, with minus infinity outside positive source
        support.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    masses = _vector(masses, "masses", nonnegative=True)
    variance = _vector(variance, "variance", nonnegative=True)
    if variance.shape != masses.shape or not np.any(masses > 0):
        raise ValueError("masses and variance must match with positive support")
    pixel = _positive(pixel, "pixel")
    alias = _positive(alias, "alias")

    support = masses > 0
    log_weights = np.full(masses.shape, -np.inf)
    zero_variance = support & (variance == 0)
    positive_variance = support & (variance > 0)
    log_weights[zero_variance] = 0.0
    log_ratio = np.log(pixel) + np.log(variance[positive_variance]) - np.log(alias)
    log_weights[positive_variance] = -np.logaddexp(0.0, log_ratio)
    return log_weights


def update_log_weights(log_weights, masses, variance, length, pixel, signal):
    """Apply one algebraically equivalent nonlinear update in log arithmetic.

    Parameters
    ----------
    log_weights : array_like, shape (n_magnitude,)
        Natural logarithms of dimensionless weights; unsupported entries are
        minus infinity.
    masses : array_like, shape (n_magnitude,)
        Nonnegative source-density quadrature masses, in deg^-2 (km/s)^-1.
    variance : array_like, shape (n_magnitude,)
        Nonnegative dimensionless pixel-noise variance at each magnitude node.
    length : float
        Positive forest length in km/s.
    pixel : float
        Positive pixel width in km/s.
    signal : float
        Intrinsic auxiliary three-dimensional power in deg^2 km/s.

    Returns
    -------
    log_weights : ndarray, shape (n_magnitude,)
        Natural logarithms of the updated dimensionless weights.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    masses = _vector(masses, "masses", nonnegative=True)
    variance = _vector(variance, "variance", nonnegative=True)
    log_weights = np.asarray(log_weights, dtype=np.float64)
    if (
        variance.shape != masses.shape
        or log_weights.shape != masses.shape
        or np.any(np.isnan(log_weights))
        or np.any(log_weights > 0)
        or not np.any(masses > 0)
    ):
        raise ValueError("invalid log weights, masses or variance")
    length = _positive(length, "length")
    pixel = _positive(pixel, "pixel")
    signal = _positive(signal, "signal")

    support = masses > 0
    if np.any(~np.isfinite(log_weights[support])):
        raise ValueError("positive supported weights must have finite logarithms")
    log_amplitude = float(np.max(log_weights[support]))
    log_shape = np.full(masses.shape, -np.inf)
    log_shape[support] = log_weights[support] - log_amplitude
    log_amplitude, log_shape = _update_amplitude_shape(
        log_amplitude, log_shape, masses, variance, length, pixel, signal
    )
    return log_amplitude + log_shape


def _update_amplitude_shape(
    log_amplitude, log_shape, masses, variance, length, pixel, signal
):
    """Update retained amplitude and relative shape without recombining them.

    Parameters
    ----------
    log_amplitude : float
        Natural logarithm of the maximum dimensionless weight.
    log_shape : ndarray, shape (n_magnitude,)
        Natural logarithms of relative weights, with maximum zero on positive
        source support.
    masses : array_like, shape (n_magnitude,)
        Nonnegative source-density quadrature masses, in deg^-2 (km/s)^-1.
    variance : array_like, shape (n_magnitude,)
        Nonnegative dimensionless pixel-noise variance at each magnitude node.
    length : float
        Positive forest length in km/s.
    pixel : float
        Positive pixel width in km/s.
    signal : float
        Intrinsic auxiliary three-dimensional power in deg^2 km/s.

    Returns
    -------
    log_amplitude : float
        Updated logarithmic maximum weight.
    log_shape : ndarray, shape (n_magnitude,)
        Updated logarithmic relative weights, normalized to a maximum of zero.
    """

    support = masses > 0
    log_mass = np.full(masses.shape, -np.inf)
    log_mass[support] = np.log(masses[support])
    log_h = _log_ratio(length, pixel) + np.logaddexp.accumulate(log_mass + log_shape)
    zero_variance = support & (variance == 0)
    positive_variance = support & (variance > 0)
    if np.any(zero_variance):
        result = np.full(masses.shape, -np.inf)
        result[zero_variance] = 0.0
        log_a = np.log(variance[positive_variance]) - np.log(signal)
        log_sh = log_amplitude + log_h[positive_variance]
        result[positive_variance] = log_sh - np.logaddexp(log_a, log_sh)
        new_amplitude = float(np.max(result[support]))
        result[support] -= new_amplitude
        return new_amplitude, result
    log_a = np.log(variance[positive_variance]) - np.log(signal)
    log_b = log_h[positive_variance] - np.logaddexp(
        log_a, log_amplitude + log_h[positive_variance]
    )
    log_b_max = float(np.max(log_b))
    result = np.full(masses.shape, -np.inf)
    result[positive_variance] = log_b - log_b_max
    return log_amplitude + log_b_max, result


def coefficients_from_log_weights(masses, variance, log_weights, length, pixel):
    """Evaluate A and P_pixel after exact cancellation of common weight scale.

    Parameters
    ----------
    masses : array_like, shape (n_magnitude,)
        Nonnegative source-density quadrature masses, in deg^-2 (km/s)^-1.
    variance : array_like, shape (n_magnitude,)
        Nonnegative dimensionless pixel-noise variance at each magnitude node.
    log_weights : array_like, shape (n_magnitude,)
        Natural logarithms of dimensionless weights; unsupported entries are
        minus infinity.
    length : float
        Positive forest length in km/s.
    pixel : float
        Positive pixel width in km/s.

    Returns
    -------
    coefficients : NoiseCoefficients
        Aliasing coefficient A in deg^2 and pixel-noise power in deg^2 km/s,
        their logarithms and explicit representability flags.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    masses = _vector(masses, "masses", nonnegative=True)
    variance = _vector(variance, "variance", nonnegative=True)
    log_weights = np.asarray(log_weights, dtype=np.float64)
    if variance.shape != masses.shape or log_weights.shape != masses.shape:
        raise ValueError("masses, variance and log_weights must match")
    length = _positive(length, "length")
    pixel = _positive(pixel, "pixel")

    support = masses > 0
    supported_logs = log_weights[support]
    if (
        not np.any(support)
        or np.any(np.isnan(supported_logs))
        or np.any(np.isposinf(supported_logs))
        or not np.any(np.isfinite(supported_logs))
    ):
        raise ValueError("positive support requires represented nonzero weights")

    log_mass = np.log(masses[support])
    # Remove the common amplitude before forming the moments.  Subtracting it
    # only after log(I2)-2*log(I1) would lose coefficient precision once the
    # finite recurrence has driven the amplitude through many decades.
    log_weight = log_weights[support] - np.max(log_weights[support])
    log_i1 = _logsumexp(log_mass + log_weight)
    log_i2 = _logsumexp(log_mass + 2 * log_weight)
    positive_noise = variance[support] > 0
    if np.any(positive_noise):
        log_i3 = _logsumexp(
            log_mass[positive_noise]
            + 2 * log_weight[positive_noise]
            + np.log(variance[support][positive_noise])
        )
        log_p = np.log(pixel) + log_i3 - 2 * log_i1 - np.log(length)
        pixel_power, p_available = _ordinary_value(log_p)
        exact_zero = False
    else:
        log_p = float("-inf")
        pixel_power, p_available = _ordinary_value(log_p, exact_zero=True)
        exact_zero = True

    log_a = log_i2 - 2 * log_i1 - np.log(length)
    aliasing_coefficient, a_available = _ordinary_value(log_a)
    return NoiseCoefficients(
        log_A=log_a,
        A=aliasing_coefficient,
        A_available=a_available,
        log_P_pixel=log_p,
        P_pixel=pixel_power,
        P_pixel_available=p_available,
        P_pixel_exact_zero=exact_zero,
    )


def _snapshot_from_shape(
    log_amplitude, log_shape, masses, variance, length, pixel, coordinates, iterations
):
    """Summarize one retained amplitude/shape state.

    Parameters
    ----------
    log_amplitude : float
        Natural logarithm of the maximum dimensionless weight.
    log_shape : ndarray, shape (n_magnitude,)
        Natural logarithms of relative weights, with maximum zero on positive
        source support.
    masses : array_like, shape (n_magnitude,)
        Nonnegative source-density quadrature masses, in deg^-2 (km/s)^-1.
    variance : array_like, shape (n_magnitude,)
        Nonnegative dimensionless pixel-noise variance at each magnitude node.
    length : float
        Positive forest length in km/s.
    pixel : float
        Positive pixel width in km/s.
    coordinates : array_like, shape (n_magnitude,), or None
        Ordered source coordinates, usually apparent magnitudes in mag; None
        uses consecutive dimensionless indices.
    iterations : int
        Number of completed nonlinear weight updates.

    Returns
    -------
    state : WeightSnapshot
        Normalized source distribution, coordinate moments and scale-cancelled
        noise coefficients.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.

    Notes
    -----
    The returned log-shape and distribution arrays are marked read-only.
    """
    masses = _vector(masses, "masses", nonnegative=True)
    variance = _vector(variance, "variance", nonnegative=True)
    coordinates = _vector(coordinates, "coordinates")
    log_shape = np.asarray(log_shape, dtype=np.float64)
    if not (variance.shape == coordinates.shape == log_shape.shape == masses.shape):
        raise ValueError("snapshot arrays must have identical shapes")

    support = masses > 0
    supported_logs = log_shape[support]
    if (
        np.any(np.isnan(supported_logs))
        or np.any(np.isposinf(supported_logs))
        or not np.any(np.isfinite(supported_logs))
    ):
        raise ValueError("positive support requires represented nonzero weights")
    log_amplitude = float(log_amplitude)
    if not np.isfinite(log_amplitude) or not np.isclose(
        np.max(log_shape[support]), 0, rtol=0, atol=5e-15
    ):
        raise ValueError("shape must be normalized with a finite amplitude")
    log_mass = np.log(masses[support])
    log_normalization = _logsumexp(log_mass + log_shape[support])
    log_distribution = log_mass + log_shape[support] - log_normalization
    distribution = np.zeros(masses.shape)
    distribution[support] = np.exp(log_distribution)
    dropped = log_distribution[distribution[support] == 0]
    dropped_bound = _logsumexp(dropped)
    if dropped_bound > np.log(1e-12):
        raise ValueError("relative-shape output loses more than 1e-12 mass")

    # Coordinate moments describe the normalized source measure, not weight amplitude.
    mean = float(np.sum(distribution * coordinates))
    sigma = float(np.sqrt(np.sum(distribution * (coordinates - mean) ** 2)))
    log_concentration = _logsumexp(2 * log_distribution - log_mass)
    coefficients = coefficients_from_log_weights(
        masses, variance, log_shape, length, pixel
    )
    if not np.isclose(
        log_concentration,
        np.log(_positive(length, "length")) + coefficients.log_A,
        rtol=0,
        atol=5e-12,
    ):
        raise ValueError("concentration and A identity disagree")
    for array in (log_shape, distribution):
        array.setflags(write=False)
    return WeightSnapshot(
        iterations=int(iterations),
        log_amplitude=log_amplitude,
        log_shape=log_shape,
        distribution=distribution,
        coefficients=coefficients,
        mean_coordinate=mean,
        sigma_coordinate=sigma,
        log_concentration=log_concentration,
        log_effective_measure=-log_concentration,
        log_dropped_distribution_bound=dropped_bound,
    )


def snapshot(log_weights, masses, variance, length, pixel, coordinates, iterations):
    """Summarize log weights, cancelling their represented common amplitude.

    Parameters
    ----------
    log_weights : array_like, shape (n_magnitude,)
        Natural logarithms of dimensionless weights; unsupported entries are
        minus infinity.
    masses : array_like, shape (n_magnitude,)
        Nonnegative source-density quadrature masses, in deg^-2 (km/s)^-1.
    variance : array_like, shape (n_magnitude,)
        Nonnegative dimensionless pixel-noise variance at each magnitude node.
    length : float
        Positive forest length in km/s.
    pixel : float
        Positive pixel width in km/s.
    coordinates : array_like, shape (n_magnitude,), or None
        Ordered source coordinates, usually apparent magnitudes in mag; None
        uses consecutive dimensionless indices.
    iterations : int
        Number of completed nonlinear weight updates.

    Returns
    -------
    state : WeightSnapshot
        Retained logarithmic amplitude, normalized source distribution and noise
        coefficients.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    log_weights = np.asarray(log_weights, dtype=np.float64)
    masses_array = _vector(masses, "masses", nonnegative=True)
    if log_weights.shape != masses_array.shape:
        raise ValueError("snapshot arrays must have identical shapes")
    support = masses_array > 0
    if not np.any(support) or np.any(~np.isfinite(log_weights[support])):
        raise ValueError("positive support requires finite log weights")
    log_amplitude = float(np.max(log_weights[support]))
    log_shape = np.full(masses_array.shape, -np.inf)
    log_shape[support] = log_weights[support] - log_amplitude
    return _snapshot_from_shape(
        log_amplitude,
        log_shape,
        masses_array,
        variance,
        length,
        pixel,
        coordinates,
        iterations,
    )


def linearized_spectrum(masses, variance, length, pixel, signal):
    """Return diagonal-spectrum diagnostics of the zero-weight Jacobian.

    Parameters
    ----------
    masses : array_like, shape (n_magnitude,)
        Nonnegative source-density quadrature masses, in deg^-2 (km/s)^-1.
    variance : array_like, shape (n_magnitude,)
        Nonnegative dimensionless pixel-noise variance at each magnitude node.
    length : float
        Positive forest length in km/s.
    pixel : float
        Positive pixel width in km/s.
    signal : float
        Intrinsic auxiliary three-dimensional power in deg^2 km/s.

    Returns
    -------
    spectral_radius : float
        Largest eigenvalue of the zero-weight Jacobian, or infinity if
        unavailable.
    log_spectral_radius : float
        Natural logarithm of the spectral radius.
    multiplicity : int
        Number of numerically coincident dominant eigenvalues.
    relative_gap : float
        Dimensionless gap to the next smaller eigenvalue.
    has_zero_variance : bool
        Whether positive source support contains zero variance.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    masses = _vector(masses, "masses", nonnegative=True)
    variance = _vector(variance, "variance", nonnegative=True)
    if masses.shape != variance.shape or not np.any(masses > 0):
        raise ValueError("masses and variance must match with positive support")
    length = _positive(length, "length")
    pixel = _positive(pixel, "pixel")
    signal = _positive(signal, "signal")

    support = masses > 0
    if np.any((variance == 0) & support):
        return float("inf"), float("inf"), 0, 0.0, True
    log_diagonal = (
        _log_ratio(length, pixel)
        + np.log(masses[support])
        + np.log(signal)
        - np.log(variance[support])
    )
    order = np.sort(log_diagonal)
    largest = float(order[-1])
    close = np.abs(log_diagonal - largest) <= 64 * np.finfo(float).eps * max(
        1.0, abs(largest)
    )
    multiplicity = int(np.count_nonzero(close))
    if len(order) == 1 or multiplicity == len(order):
        gap = 0.0
    else:
        next_log = float(order[-multiplicity - 1])
        gap = float(-np.expm1(next_log - largest))
    radius, available = _ordinary_value(largest)
    return (
        radius if available else float("inf"),
        largest,
        multiplicity,
        gap,
        False,
    )


def linearized_fixed_grid_limit(
    masses, variance, length, pixel, signal, *, coordinates=None
):
    """Return the unique-dominant eigenvector limit of the discrete linear map.

    Parameters
    ----------
    masses : array_like, shape (n_magnitude,)
        Nonnegative source-density quadrature masses, in deg^-2 (km/s)^-1.
    variance : array_like, shape (n_magnitude,)
        Nonnegative dimensionless pixel-noise variance at each magnitude node.
    length : float
        Positive forest length in km/s.
    pixel : float
        Positive pixel width in km/s.
    signal : float
        Intrinsic auxiliary three-dimensional power in deg^2 km/s.
    coordinates : array_like, shape (n_magnitude,), or None
        Ordered source coordinates, usually apparent magnitudes in mag; None
        uses consecutive dimensionless indices. Default is ``None``.

    Returns
    -------
    limit : tuple of (int, WeightSnapshot) or None
        Dominant source index and normalized discrete eigenvector snapshot; None
        if the stated conditions do not establish this limit.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.

    Notes
    -----
    The result is conditional on a strictly positive variance, a unique dominant
    diagonal entry and spectral radius below one.  It is a fixed-grid result,
    not evidence that a continuum limit exists.
    """
    masses = _vector(masses, "masses", nonnegative=True)
    variance = _vector(variance, "variance", nonnegative=True)
    if masses.shape != variance.shape:
        raise ValueError("masses and variance must match")
    if coordinates is None:
        coordinates = np.arange(len(masses), dtype=np.float64)
    else:
        coordinates = _vector(coordinates, "coordinates")
        if coordinates.shape != masses.shape or np.any(np.diff(coordinates) <= 0):
            raise ValueError("coordinates must match masses in increasing order")
    radius, log_radius, multiplicity, _, has_zero = linearized_spectrum(
        masses, variance, length, pixel, signal
    )
    if has_zero or multiplicity != 1 or not radius < 1:
        return None
    support_indices = np.flatnonzero(masses > 0)
    log_diagonal = (
        _log_ratio(length, pixel)
        + np.log(masses[support_indices])
        + np.log(signal)
        - np.log(variance[support_indices])
    )
    dominant_support_index = int(np.argmax(log_diagonal))
    dominant_index = int(support_indices[dominant_support_index])
    log_vector = np.full(masses.shape, -np.inf)
    log_vector[dominant_index] = 0.0
    log_prefix = np.log(masses[dominant_index])
    for index in support_indices[dominant_support_index + 1 :]:
        log_d = (
            _log_ratio(length, pixel)
            + np.log(masses[index])
            + np.log(signal)
            - np.log(variance[index])
        )
        ratio = np.exp(log_d - log_radius)
        log_difference = log_radius + np.log1p(-ratio)
        log_c = _log_ratio(length, pixel) + np.log(signal) - np.log(variance[index])
        log_vector[index] = log_c + log_prefix - log_difference
        log_prefix = np.logaddexp(log_prefix, np.log(masses[index]) + log_vector[index])
    finite = np.isfinite(log_vector)
    log_vector[finite] -= np.max(log_vector[finite])
    state = _snapshot_from_shape(
        0.0,
        log_vector,
        masses,
        variance,
        length,
        pixel,
        coordinates,
        -1,
    )
    return dominant_index, state


def trajectory(
    masses,
    variance,
    length,
    pixel,
    signal,
    alias,
    *,
    checkpoints=(0, 3, 6, 12, 24, 48, 96, 192, 384, 768, 1024),
    coordinates=None,
):
    """Follow the exact finite recurrence at explicit bounded checkpoints.

    Parameters
    ----------
    masses : array_like, shape (n_magnitude,)
        Nonnegative source-density quadrature masses, in deg^-2 (km/s)^-1.
    variance : array_like, shape (n_magnitude,)
        Nonnegative dimensionless pixel-noise variance at each magnitude node.
    length : float
        Positive forest length in km/s.
    pixel : float
        Positive pixel width in km/s.
    signal : float
        Intrinsic auxiliary three-dimensional power in deg^2 km/s.
    alias : float
        Positive one-dimensional auxiliary forest power in km/s.
    checkpoints : sequence of int
        Strictly increasing nonnegative iteration counts to retain. Default is
        ``(0, 3, 6, 12, 24, 48, 96, 192, 384, 768, 1024)``.
    coordinates : array_like, shape (n_magnitude,), or None
        Ordered source coordinates, usually apparent magnitudes in mag; None
        uses consecutive dimensionless indices. Default is ``None``.

    Returns
    -------
    trajectory : WeightTrajectory
        Saved finite-count snapshots and fixed-grid linearized spectral
        diagnostics.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    masses = _vector(masses, "masses", nonnegative=True)
    variance = _vector(variance, "variance", nonnegative=True)
    if masses.shape != variance.shape:
        raise ValueError("masses and variance must match")
    if coordinates is None:
        coordinates = np.arange(len(masses), dtype=np.float64)
    else:
        coordinates = _vector(coordinates, "coordinates")
        if coordinates.shape != masses.shape or np.any(np.diff(coordinates) <= 0):
            raise ValueError("coordinates must match masses in increasing order")
    requested = tuple(int(value) for value in checkpoints)
    if (
        not requested
        or requested[0] < 0
        or any(value < 0 for value in requested)
        or any(b <= a for a, b in zip(requested[:-1], requested[1:]))
    ):
        raise ValueError("checkpoints must be strictly increasing nonnegative integers")
    log_weights = initial_log_weights(masses, variance, pixel, alias)

    support = masses > 0
    log_amplitude = float(np.max(log_weights[support]))
    log_shape = np.full(masses.shape, -np.inf)
    log_shape[support] = log_weights[support] - log_amplitude
    states = []
    next_index = 0
    for count in range(requested[-1] + 1):
        if count == requested[next_index]:
            states.append(
                _snapshot_from_shape(
                    log_amplitude,
                    log_shape,
                    masses,
                    variance,
                    length,
                    pixel,
                    coordinates,
                    count,
                )
            )
            next_index += 1
            if next_index == len(requested):
                break
        log_amplitude, log_shape = _update_amplitude_shape(
            log_amplitude,
            log_shape,
            masses,
            variance,
            length,
            pixel,
            signal,
        )
    radius, log_radius, multiplicity, gap, has_zero = linearized_spectrum(
        masses, variance, length, pixel, signal
    )
    return WeightTrajectory(
        snapshots=tuple(states),
        spectral_radius=radius,
        log_spectral_radius=log_radius,
        dominant_multiplicity=multiplicity,
        relative_spectral_gap=gap,
        has_supported_zero_variance=has_zero,
    )


def relative_log_change(log_lower, log_upper):
    """Return |upper/lower - 1| without first exponentiating either value.

    Parameters
    ----------
    log_lower : float
        Natural logarithm of the reference quantity.
    log_upper : float
        Natural logarithm of the comparison quantity, in the same units.

    Returns
    -------
    relative_change : float
        Absolute dimensionless upper/lower minus one, preserving overflow as
        infinity.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    values = np.asarray([log_lower, log_upper], dtype=np.float64)
    if np.any(np.isnan(values)):
        raise ValueError("log coefficients must not be NaN")
    if np.all(np.isneginf(values)):
        return 0.0
    if np.any(~np.isfinite(values)):
        return float("inf")
    delta = float(log_upper - log_lower)
    if delta > _LOG_MAX:
        return float("inf")
    return float(abs(np.expm1(delta)))


def signed_log_ratio(log_lower, log_upper):
    """Return upper/lower - 1 while retaining its direction and log range.

    Parameters
    ----------
    log_lower : float
        Natural logarithm of the reference quantity.
    log_upper : float
        Natural logarithm of the comparison quantity, in the same units.

    Returns
    -------
    relative_change : float
        Signed dimensionless upper/lower minus one. Infinite ratios retain
        infinity; comparing two positive infinities returns NaN.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    values = np.asarray([log_lower, log_upper], dtype=np.float64)
    if np.any(np.isnan(values)):
        raise ValueError("log coefficients must not be NaN")
    if np.all(np.isneginf(values)):
        return 0.0
    if np.isneginf(log_lower):
        return float("inf")
    if np.isneginf(log_upper):
        return -1.0
    if np.any(np.isposinf(values)):
        return float("nan") if log_lower == log_upper else float("inf")
    delta = float(log_upper - log_lower)
    if delta > _LOG_MAX:
        return float("inf")
    if delta < _LOG_MIN_SUBNORMAL:
        return -1.0
    return float(np.expm1(delta))


def fixed_mesh_screen(result, *, threshold=1e-4):
    """Screen three successive late doublings without asserting a true limit.

    Parameters
    ----------
    result : WeightTrajectory
        Trajectory containing the required late-doubling snapshots.
    threshold : float
        Dimensionless relative-change threshold. Default is ``0.0001``.

    Returns
    -------
    screen : dict
        Late-doubling coefficient and distribution comparisons, threshold and
        finite-grid screening status.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    threshold = _positive(threshold, "threshold")
    by_count = {state.iterations: state for state in result.snapshots}
    required = (128, 256, 512, 1024)
    if any(count not in by_count for count in required):
        raise ValueError("fixed-mesh screen requires checkpoints 128/256/512/1024")
    comparisons = []
    for lower, upper in zip(required[:-1], required[1:]):
        aliasing_change = relative_log_change(
            by_count[lower].coefficients.log_A,
            by_count[upper].coefficients.log_A,
        )
        pixel_power_change = relative_log_change(
            by_count[lower].coefficients.log_P_pixel,
            by_count[upper].coefficients.log_P_pixel,
        )
        concentration = relative_log_change(
            by_count[lower].log_concentration,
            by_count[upper].log_concentration,
        )
        comparisons.append(
            dict(
                iterations=[lower, upper],
                A_relative=aliasing_change,
                P_pixel_relative=pixel_power_change,
                concentration_relative=concentration,
                passed=aliasing_change <= threshold and pixel_power_change <= threshold,
            )
        )
    return dict(
        threshold=threshold,
        comparisons=comparisons,
        passed=all(row["passed"] for row in comparisons),
        interpretation="bounded fixed-mesh screen, not a continuum-limit verdict",
    )


def continuum_linearized_coefficients(
    iterations, *, length, total_mass, pixel, variance
):
    """Analytic constant-coefficient continuum control from the Step 13 plan.

    Parameters
    ----------
    iterations : int
        Number of completed nonlinear weight updates.
    length : float
        Positive forest length in km/s.
    total_mass : float
        Total positive source measure in deg^-2 (km/s)^-1.
    pixel : float
        Positive pixel width in km/s.
    variance : float
        Constant positive dimensionless pixel variance.

    Returns
    -------
    aliasing : float
        Continuum aliasing coefficient in deg^2.
    pixel_power : float
        Continuum pixel-noise power in deg^2 km/s.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    OverflowError :
        If the requested coefficients are outside the representable float64
        range.
    """
    iterations = int(iterations)
    if iterations < 0:
        raise ValueError("iterations must be nonnegative")
    length = _positive(length, "length")
    total_mass = _positive(total_mass, "total_mass")
    pixel = _positive(pixel, "pixel")
    variance = _positive(variance, "variance")
    log_a = (
        2 * np.log(iterations + 1)
        - np.log(2 * iterations + 1)
        - np.log(length)
        - np.log(total_mass)
    )
    aliasing_coefficient, a_available = _ordinary_value(log_a)
    log_p = np.log(pixel) + np.log(variance) + log_a
    pixel_power, p_available = _ordinary_value(log_p)
    if not a_available or not p_available:
        raise OverflowError("analytic coefficients are outside float64 range")
    return aliasing_coefficient, pixel_power


def scalar_weight(iterations, initial_weight, d):
    """Closed-form one-cell recurrence, including the marginal d=1 case.

    Parameters
    ----------
    iterations : int
        Number of completed nonlinear weight updates.
    initial_weight : float
        Initial dimensionless weight in (0, 1].
    d : float
        Positive dimensionless linearized recurrence multiplier.

    Returns
    -------
    weight : float
        Dimensionless scalar weight after the requested transitions.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    iterations = int(iterations)
    initial_weight = _positive(initial_weight, "initial_weight")
    d = _positive(d, "d")
    if iterations < 0 or initial_weight > 1:
        raise ValueError("require nonnegative iterations and initial_weight <= 1")
    log_weight = float(np.log(initial_weight))
    log_d = float(np.log(d))
    for _ in range(iterations):
        log_product = log_d + log_weight
        log_weight = log_product - np.logaddexp(0.0, log_product)
    return float(np.exp(log_weight))


def _array_identity(arrays, names):
    """Hash exact diagnostic operands and their array identities.

    Parameters
    ----------
    arrays : dict of str to ndarray
        Numerical evidence arrays; Fourier-cell axes and pair order follow the
        declared task. Powers use (Mpc/h)^3 and volumes use (Mpc/h)^3 unless
        separately labeled.
    names : sequence of str
        Ordered names of the quantities to process.

    Returns
    -------
    digest : str
        SHA-256 of the ordered names, dtypes, shapes and contiguous array bytes.
    """
    digest = hashlib.sha256()
    for name in names:
        array = np.ascontiguousarray(arrays[name])
        digest.update(name.encode())
        digest.update(array.dtype.str.encode())
        digest.update(str(array.shape).encode())
        digest.update(array.tobytes())
    return digest.hexdigest()


def _snapshot_record(state):
    """Serialize one weight snapshot without inventing unavailable values.

    Parameters
    ----------
    state : WeightSnapshot
        Retained amplitude, normalized shape and noise-coefficient diagnostics
        at one update count.

    Returns
    -------
    record : dict
        JSON-compatible snapshot scalars and availability flags.
    """
    coefficients = state.coefficients

    def finite(value):
        """Represent a nonfinite diagnostic scalar as unavailable.

        Parameters
        ----------
        value : float
            Diagnostic scalar, which may be nonfinite.

        Returns
        -------
        value : float or None
            Finite scalar, or None for a nonfinite input.
        """
        return float(value) if np.isfinite(value) else None

    return dict(
        iterations=state.iterations,
        log_amplitude=state.log_amplitude,
        log_A=coefficients.log_A,
        A=coefficients.A,
        A_available=coefficients.A_available,
        log_P_pixel=finite(coefficients.log_P_pixel),
        P_pixel=coefficients.P_pixel,
        P_pixel_available=coefficients.P_pixel_available,
        P_pixel_exact_zero=coefficients.P_pixel_exact_zero,
        mean_coordinate=state.mean_coordinate,
        sigma_coordinate=state.sigma_coordinate,
        log_concentration=state.log_concentration,
        log_effective_measure=state.log_effective_measure,
        log_dropped_distribution_bound=finite(state.log_dropped_distribution_bound),
    )


def _fixed_grid_tail_record(state, masses, variance, coordinates, length, pixel):
    """Describe concentration in the faintest tenth of the sampled support.

    Parameters
    ----------
    state : WeightSnapshot
        Retained amplitude, normalized shape and noise-coefficient diagnostics
        at one update count.
    masses : array_like, shape (n_magnitude,)
        Nonnegative source-density quadrature masses, in deg^-2 (km/s)^-1.
    variance : array_like, shape (n_magnitude,)
        Nonnegative dimensionless pixel-noise variance at each magnitude node.
    coordinates : array_like, shape (n_magnitude,), or None
        Ordered source coordinates, usually apparent magnitudes in mag; None
        uses consecutive dimensionless indices.
    length : float
        Positive forest length in km/s.
    pixel : float
        Positive pixel width in km/s.

    Returns
    -------
    record : dict
        Tail probability, source measure and Cauchy lower bounds on A and pixel-
        noise power.
    """
    masses = np.asarray(masses, dtype=np.float64)
    variance = np.asarray(variance, dtype=np.float64)
    coordinates = np.asarray(coordinates, dtype=np.float64)

    support = masses > 0
    lower = float(
        np.max(coordinates[support])
        - 0.1 * (np.max(coordinates[support]) - np.min(coordinates[support]))
    )
    tail = support & (coordinates >= lower)
    log_mass = np.log(masses[tail])
    log_weighted_mass = log_mass + state.log_shape[tail]
    log_total_weighted_mass = _logsumexp(
        np.log(masses[support]) + state.log_shape[support]
    )
    log_probability = _logsumexp(log_weighted_mass) - log_total_weighted_mass
    probability, probability_available = _ordinary_value(log_probability)
    log_measure = _logsumexp(log_mass)
    measure, measure_available = _ordinary_value(log_measure)
    log_a_lower = 2 * log_probability - np.log(length) - log_measure
    a_lower, a_lower_available = _ordinary_value(log_a_lower)
    positive_variance = variance[tail & (variance > 0)]
    minimum_variance = (
        float(np.min(positive_variance)) if positive_variance.size else 0.0
    )
    if minimum_variance > 0:
        log_p_lower = np.log(pixel) + np.log(minimum_variance) + log_a_lower
        p_lower, p_lower_available = _ordinary_value(log_p_lower)
    else:
        log_p_lower = float("-inf")
        p_lower, p_lower_available = _ordinary_value(log_p_lower, exact_zero=True)
    return dict(
        definition="faintest 10 percent of the sampled magnitude interval",
        lower_coordinate=lower,
        upper_coordinate=float(np.max(coordinates[support])),
        probability=probability,
        probability_available=probability_available,
        log_probability=log_probability,
        source_measure=measure,
        source_measure_available=measure_available,
        log_source_measure=log_measure,
        minimum_positive_variance=minimum_variance,
        A_cauchy_lower_bound=a_lower,
        A_cauchy_lower_bound_available=a_lower_available,
        log_A_cauchy_lower_bound=log_a_lower,
        P_pixel_lower_bound=p_lower,
        P_pixel_lower_bound_available=p_lower_available,
        log_P_pixel_lower_bound=log_p_lower,
        interpretation=(
            "fixed-grid distribution diagnostic; finite refinements do not "
            "establish a continuum tail limit"
        ),
    )


def _one_step_remainder_record(state, masses, variance, length, pixel, signal):
    """Quantify the current nonlinear correction to the linearized map.

    Parameters
    ----------
    state : WeightSnapshot
        Retained amplitude, normalized shape and noise-coefficient diagnostics
        at one update count.
    masses : array_like, shape (n_magnitude,)
        Nonnegative source-density quadrature masses, in deg^-2 (km/s)^-1.
    variance : array_like, shape (n_magnitude,)
        Nonnegative dimensionless pixel-noise variance at each magnitude node.
    length : float
        Positive forest length in km/s.
    pixel : float
        Positive pixel width in km/s.
    signal : float
        Intrinsic auxiliary three-dimensional power in deg^2 km/s.

    Returns
    -------
    record : dict
        Maximum linearized update and the exact one-step nonlinear remainder
        identity.
    """
    masses = np.asarray(masses, dtype=np.float64)
    variance = np.asarray(variance, dtype=np.float64)
    support = (masses > 0) & (variance > 0)
    log_mass = np.log(masses[support])
    log_prefix = np.logaddexp.accumulate(log_mass + state.log_shape[support])
    log_x = (
        _log_ratio(length, pixel)
        + np.log(signal)
        - np.log(variance[support])
        + state.log_amplitude
        + log_prefix
    )
    maximum_log_x = float(np.max(log_x))
    maximum_x, available = _ordinary_value(maximum_log_x)
    return dict(
        iterations=state.iterations,
        maximum_log_linearized_update=maximum_log_x,
        maximum_linearized_update=maximum_x,
        maximum_linearized_update_available=available,
        componentwise_identity="M w - F(w) = (M w)^2 / (1 + M w)",
        accumulated_future_remainder_bound=None,
        interpretation=(
            "one-step finite-trajectory diagnostic only; the saved checkpoint "
            "is not used as an error-bounded approximation to the asymptote"
        ),
    )


def batch_payload(
    result,
    *,
    magnitudes,
    masses,
    variance,
    length,
    pixel,
    signal,
    alias,
    bin_index,
    field,
    order,
    source,
):
    """Create a hash-bound array/report pair for one diagnostic batch.

    Parameters
    ----------
    result : WeightTrajectory
        Finite trajectory with a stable unique discrete linearized limit.
    magnitudes : array_like, shape (n_magnitude,)
        Apparent-magnitude quadrature nodes in mag.
    masses : array_like, shape (n_magnitude,)
        Nonnegative source-density quadrature masses, in deg^-2 (km/s)^-1.
    variance : array_like, shape (n_magnitude,)
        Nonnegative dimensionless pixel-noise variance at each magnitude node.
    length : float
        Positive forest length in km/s.
    pixel : float
        Positive pixel width in km/s.
    signal : float
        Intrinsic auxiliary three-dimensional power in deg^2 km/s.
    alias : float
        Positive one-dimensional auxiliary forest power in km/s.
    bin_index : int
        Zero-based redshift-bin index.
    field : str
        Observed forest field identifier.
    order : int
        Gauss-Legendre quadrature order within each magnitude interval.
    source : dict
        Provenance of the physical source inputs.

    Returns
    -------
    arrays : dict of str to ndarray
        Diagnostic inputs and checkpoint states.
    report : dict
        Hash-bound provenance, spectral, tail and finite-trajectory diagnostics.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    if not isinstance(result, WeightTrajectory):
        raise ValueError("result must be WeightTrajectory")
    fixed_limit = linearized_fixed_grid_limit(
        masses,
        variance,
        length,
        pixel,
        signal,
        coordinates=magnitudes,
    )
    if fixed_limit is None:
        raise ValueError(
            "saved batch lacks a stable unique fixed-grid linearized limit"
        )
    dominant_index, limit_state = fixed_limit
    arrays = dict(
        magnitudes=_vector(magnitudes, "magnitudes"),
        masses=_vector(masses, "masses", nonnegative=True),
        variance=_vector(variance, "variance", nonnegative=True),
        length=np.asarray([_positive(length, "length")]),
        pixel=np.asarray([_positive(pixel, "pixel")]),
        signal=np.asarray([_positive(signal, "signal")]),
        alias=np.asarray([_positive(alias, "alias")]),
        iterations=np.asarray([state.iterations for state in result.snapshots]),
        log_amplitude=np.asarray([state.log_amplitude for state in result.snapshots]),
        log_shape=np.stack([state.log_shape for state in result.snapshots]),
        distribution=np.stack([state.distribution for state in result.snapshots]),
        log_A=np.asarray([state.coefficients.log_A for state in result.snapshots]),
        log_P_pixel=np.asarray(
            [state.coefficients.log_P_pixel for state in result.snapshots]
        ),
        mean_coordinate=np.asarray(
            [state.mean_coordinate for state in result.snapshots]
        ),
        sigma_coordinate=np.asarray(
            [state.sigma_coordinate for state in result.snapshots]
        ),
        log_concentration=np.asarray(
            [state.log_concentration for state in result.snapshots]
        ),
        linearized_log_shape=limit_state.log_shape,
        linearized_distribution=limit_state.distribution,
    )
    report = dict(
        schema=2,
        kind="cumulative_weight_limit_diagnostic",
        bin=int(bin_index),
        field=str(field),
        order=int(order),
        source=dict(source),
        input_sha256=_array_identity(
            arrays,
            ("magnitudes", "masses", "variance", "length", "pixel", "signal", "alias"),
        ),
        array_inventory={name: list(value.shape) for name, value in arrays.items()},
        checkpoints=[_snapshot_record(state) for state in result.snapshots],
        attempts=[
            dict(iterations=state.iterations, outcome="completed")
            for state in result.snapshots
        ],
        capped_at=int(arrays["iterations"][-1]),
        linearized=dict(
            spectral_radius=result.spectral_radius,
            log_spectral_radius=result.log_spectral_radius,
            dominant_multiplicity=result.dominant_multiplicity,
            relative_spectral_gap=result.relative_spectral_gap,
            has_supported_zero_variance=result.has_supported_zero_variance,
            dominant_index=dominant_index,
            dominant_coordinate=float(np.asarray(magnitudes)[dominant_index]),
            fixed_grid_limit=_snapshot_record(limit_state),
        ),
        fixed_grid_tail=_fixed_grid_tail_record(
            limit_state, masses, variance, magnitudes, length, pixel
        ),
        nonlinear_remainder=_one_step_remainder_record(
            result.snapshots[-1], masses, variance, length, pixel, signal
        ),
    )
    if all(count in arrays["iterations"] for count in (128, 256, 512, 1024)):
        report["fixed_mesh_screen"] = fixed_mesh_screen(result)
    return arrays, report


def _close(actual, expected, name):
    """Check diagnostic operands at the saved numerical tolerance.

    Parameters
    ----------
    actual : array_like
        Saved quantity.
    expected : array_like
        Recomputed quantity with the same shape and units.
    name : str
        Quantity or record label used in diagnostics.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    if not np.allclose(actual, expected, rtol=5e-12, atol=5e-14, equal_nan=False):
        raise ValueError(f"{name}: inconsistent numerical content")


def _source_close(actual, expected, name):
    """Compare one physical source quantity without an unrelated absolute floor.

    Parameters
    ----------
    actual : array_like
        Saved quantity.
    expected : array_like
        Recomputed quantity with the same shape and units.
    name : str
        Quantity or record label used in diagnostics.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    if actual.shape != expected.shape or not np.allclose(
        actual, expected, rtol=5e-12, atol=0, equal_nan=False
    ):
        raise ValueError(f"source {name}: inconsistent numerical content")


def validate_batch_payload(arrays, report, *, expected=None):
    """Recompute a saved batch and bind every checkpoint to its own operands.

    Parameters
    ----------
    arrays : dict of str to ndarray
        Numerical evidence arrays; Fourier-cell axes and pair order follow the
        declared task. Powers use (Mpc/h)^3 and volumes use (Mpc/h)^3 unless
        separately labeled.
    report : dict
        Scientific settings, provenance, array inventory and validation outcomes
        associated with the numerical evidence.
    expected : dict or None
        Optional externally declared bin, field, order, source operands and
        attempted iteration counts. Default is ``None``.

    Returns
    -------
    valid : bool
        True after all saved identities and recomputed diagnostics agree.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    required = (
        "magnitudes",
        "masses",
        "variance",
        "length",
        "pixel",
        "signal",
        "alias",
        "iterations",
        "log_amplitude",
        "log_shape",
        "distribution",
        "log_A",
        "log_P_pixel",
        "mean_coordinate",
        "sigma_coordinate",
        "log_concentration",
        "linearized_log_shape",
        "linearized_distribution",
    )
    if set(arrays) != set(required):
        raise ValueError("batch array inventory is incomplete or contains extras")
    inventory = {name: list(np.asarray(value).shape) for name, value in arrays.items()}
    if report.get("array_inventory") != inventory:
        raise ValueError("batch array inventory disagrees with report")
    if (
        report.get("schema") != 2
        or report.get("kind") != "cumulative_weight_limit_diagnostic"
    ):
        raise ValueError("unsupported batch report")
    if expected is not None:
        for name in ("bin", "field", "order"):
            if report.get(name) != expected.get(name):
                raise ValueError(f"batch {name} identity mismatch")
        for name, value in expected.get("inputs", {}).items():
            _source_close(arrays[name], value, name)
        if "source" in expected and report.get("source") != expected["source"]:
            raise ValueError("batch source provenance mismatch")
    identity = _array_identity(
        arrays,
        ("magnitudes", "masses", "variance", "length", "pixel", "signal", "alias"),
    )
    if report.get("input_sha256") != identity:
        raise ValueError("batch input hash mismatch")
    counts = tuple(int(value) for value in np.asarray(arrays["iterations"]))
    attempts = report.get("attempts")
    if not isinstance(attempts, list) or not attempts:
        raise ValueError("batch attempts must be a nonempty list")
    attempted_counts = []
    completed_counts = []
    for attempt in attempts:
        if not isinstance(attempt, dict) or attempt.get("outcome") not in {
            "completed",
            "failed",
            "capped",
        }:
            raise ValueError("batch attempts contain an invalid outcome")
        count = attempt.get("iterations")
        if not isinstance(count, int) or count < 0:
            raise ValueError("batch attempts contain an invalid iteration")
        attempted_counts.append(count)
        if attempt["outcome"] == "completed":
            completed_counts.append(count)
        elif not isinstance(attempt.get("error"), str) or not attempt["error"]:
            raise ValueError("failed or capped attempts require an error reason")
    if any(
        upper <= lower for lower, upper in zip(attempted_counts, attempted_counts[1:])
    ):
        raise ValueError("batch attempts must be strictly ordered")
    if tuple(completed_counts) != counts:
        raise ValueError("completed attempts do not bind every checkpoint")
    target_attempts = [dict(iterations=value, outcome="completed") for value in counts]
    if expected is not None and "attempts" in expected:
        target_attempts = expected["attempts"]
    if attempts != target_attempts:
        raise ValueError("batch attempts do not bind every checkpoint")
    target_cap = counts[-1]
    if expected is not None and "capped_at" in expected:
        target_cap = expected["capped_at"]
    if report.get("capped_at") != target_cap:
        raise ValueError("batch cap does not match the final attempt")
    if target_cap != attempted_counts[-1]:
        raise ValueError("batch cap does not match attempted scope")
    recomputed = trajectory(
        arrays["masses"],
        arrays["variance"],
        float(arrays["length"][0]),
        float(arrays["pixel"][0]),
        float(arrays["signal"][0]),
        float(arrays["alias"][0]),
        checkpoints=counts,
        coordinates=arrays["magnitudes"],
    )
    records = [_snapshot_record(state) for state in recomputed.snapshots]
    if len(report.get("checkpoints", ())) != len(records):
        raise ValueError("batch checkpoint report inventory mismatch")
    for index, (actual, target) in enumerate(zip(report["checkpoints"], records)):
        if actual.keys() != target.keys():
            raise ValueError(f"checkpoint {index}: report schema mismatch")
        for name in actual:
            if isinstance(target[name], bool) or target[name] is None:
                if actual[name] != target[name]:
                    raise ValueError(f"checkpoint {index} {name}: inconsistent report")
            else:
                _close(actual[name], target[name], f"checkpoint {index} {name}")

    # Reconstruct each checkpoint from its own source operands before comparison.
    comparisons = dict(
        log_amplitude=[state.log_amplitude for state in recomputed.snapshots],
        log_shape=np.stack([state.log_shape for state in recomputed.snapshots]),
        distribution=np.stack([state.distribution for state in recomputed.snapshots]),
        log_A=[state.coefficients.log_A for state in recomputed.snapshots],
        log_P_pixel=[state.coefficients.log_P_pixel for state in recomputed.snapshots],
        mean_coordinate=[state.mean_coordinate for state in recomputed.snapshots],
        sigma_coordinate=[state.sigma_coordinate for state in recomputed.snapshots],
        log_concentration=[state.log_concentration for state in recomputed.snapshots],
    )
    for name, target in comparisons.items():
        _close(arrays[name], target, name)
    fixed_limit = linearized_fixed_grid_limit(
        arrays["masses"],
        arrays["variance"],
        float(arrays["length"][0]),
        float(arrays["pixel"][0]),
        float(arrays["signal"][0]),
        coordinates=arrays["magnitudes"],
    )
    if fixed_limit is None:
        raise ValueError("batch no longer has a unique fixed-grid limit")
    dominant_index, limit_state = fixed_limit
    _close(arrays["linearized_log_shape"], limit_state.log_shape, "linearized shape")
    _close(
        arrays["linearized_distribution"],
        limit_state.distribution,
        "linearized distribution",
    )
    linearized = report.get("linearized", {})
    target_linearized = dict(
        spectral_radius=recomputed.spectral_radius,
        log_spectral_radius=recomputed.log_spectral_radius,
        dominant_multiplicity=recomputed.dominant_multiplicity,
        relative_spectral_gap=recomputed.relative_spectral_gap,
        has_supported_zero_variance=recomputed.has_supported_zero_variance,
        dominant_index=dominant_index,
        dominant_coordinate=float(arrays["magnitudes"][dominant_index]),
        fixed_grid_limit=_snapshot_record(limit_state),
    )
    if linearized.keys() != target_linearized.keys():
        raise ValueError("linearized diagnostic schema mismatch")
    for name, target in target_linearized.items():
        if isinstance(target, dict):
            if linearized[name].keys() != target.keys():
                raise ValueError("linearized fixed-grid limit schema mismatch")
            for quantity, value in target.items():
                if value is None or isinstance(value, bool):
                    if linearized[name][quantity] != value:
                        raise ValueError(f"linearized fixed-grid {quantity} mismatch")
                else:
                    _close(
                        linearized[name][quantity],
                        value,
                        f"linearized fixed-grid {quantity}",
                    )
        elif isinstance(target, bool) or isinstance(target, int):
            if linearized[name] != target:
                raise ValueError(f"linearized {name} mismatch")
        else:
            _close(linearized[name], target, f"linearized {name}")
    if "fixed_mesh_screen" in report:
        if report["fixed_mesh_screen"] != fixed_mesh_screen(recomputed):
            raise ValueError("fixed-mesh convergence verdict is inconsistent")
    target_tail = _fixed_grid_tail_record(
        limit_state,
        arrays["masses"],
        arrays["variance"],
        arrays["magnitudes"],
        float(arrays["length"][0]),
        float(arrays["pixel"][0]),
    )
    if report.get("fixed_grid_tail") != target_tail:
        raise ValueError("fixed-grid tail diagnostic is inconsistent")
    target_remainder = _one_step_remainder_record(
        recomputed.snapshots[-1],
        arrays["masses"],
        arrays["variance"],
        float(arrays["length"][0]),
        float(arrays["pixel"][0]),
        float(arrays["signal"][0]),
    )
    if report.get("nonlinear_remainder") != target_remainder:
        raise ValueError("nonlinear remainder diagnostic is inconsistent")
    return True
