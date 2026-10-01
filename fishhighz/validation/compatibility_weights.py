"""Signed, validation-only compatibility weighting comparisons.

Literal arithmetic adapted from lyaforecast weights.py/covariance.py (GPLv3).
No production preparation, interpolation, response or estimator policy is changed.
"""

from dataclasses import dataclass

import numpy as np

from ..kernels import full_sum_weights as full_sum
from ..models.p1d import default_p1d
from ..response import velocity_response
from .numerics import wick

VARIANTS = (
    "prefix_intrinsic",
    "sum_intrinsic",
    "prefix_aliasing",
    "sum_aliasing",
    "sum_historical",
)


@dataclass(frozen=True)
class WeightInputs:
    """One saved pair context in angular/velocity units, including signed inputs."""

    magnitudes: np.ndarray
    density: np.ndarray
    quadrature: np.ndarray
    variance: np.ndarray
    length: float
    pixel: float
    signal: float
    p1d: float

    @classmethod
    def from_pair(cls, row):
        """Preserve each pair's own uniform compatibility magnitude measure.

        Parameters
        ----------
        row : dict
            Captured pair context containing magnitude, density, variance and
            auxiliary-signal arrays.

        Returns
        -------
        inputs : WeightInputs
            Captured angular/velocity source inputs with the literal uniform
            magnitude measure.
        """
        magnitude_grid = np.asarray(row["magnitudes"], dtype=float)
        return cls(
            magnitude_grid,
            np.asarray(row["density"], dtype=float),
            np.full(magnitude_grid.shape, magnitude_grid[1] - magnitude_grid[0]),
            np.asarray(row["variance"], dtype=float),
            float(row["forest_length"]),
            float(row["_pix_kms"]),
            float(row["auxiliary_signal"]),
            float(row["auxiliary_p1d"]),
        )


def seed(inputs):
    """Common legacy seed, preserving the saved division order.

    Parameters
    ----------
    inputs : WeightInputs
        Magnitude-dependent source density, integration measure, pixel variance
        and auxiliary signal in angular/velocity units.

    Returns
    -------
    weights : ndarray, shape (n_magnitude,)
        Dimensionless initial inverse-variance weights.
    """
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        power = inputs.p1d / inputs.pixel
        return power / (power + inputs.variance)


def moments(inputs, weights):
    """Cumulative I1/I2/I3 with literal legacy multiplication/reduction order.

    Parameters
    ----------
    inputs : WeightInputs
        Magnitude-dependent source density, integration measure, pixel variance
        and auxiliary signal in angular/velocity units.
    weights : array_like, shape (n_magnitude,)
        Dimensionless source weights, retaining their common amplitude.

    Returns
    -------
    first_moment : ndarray, shape (n_magnitude,)
        Cumulative density-weight integral I1 in deg^-2 (km/s)^-1.
    second_moment : ndarray, shape (n_magnitude,)
        Cumulative squared-weight integral I2 in the same units.
    noise_moment : ndarray, shape (n_magnitude,)
        Cumulative variance-weighted integral I3 in the same units.
    """
    density, magnitude_measure, pixel_variance = (
        inputs.density,
        inputs.quadrature,
        inputs.variance,
    )
    return (
        np.cumsum(density * weights * magnitude_measure),
        np.cumsum(density * weights**2 * magnitude_measure),
        np.cumsum(density * weights**2 * pixel_variance * magnitude_measure),
    )


def update(inputs, weights, variant):
    """One simultaneous update, retaining amplitude and signed arithmetic.

    Parameters
    ----------
    inputs : WeightInputs
        Magnitude-dependent source density, integration measure, pixel variance
        and auxiliary signal in angular/velocity units.
    weights : array_like, shape (n_magnitude,)
        Dimensionless source weights, retaining their common amplitude.
    variant : str
        Named compatibility recurrence, selecting prefix or full-sample moments
        and the auxiliary signal prescription.

    Returns
    -------
    weights : ndarray, shape (n_magnitude,)
        Updated dimensionless weights, including their common amplitude.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    FloatingPointError :
        If the weighting arithmetic produces an invalid or nonfinite state.

    Notes
    -----
    Floating-point failures propagate; callers must report the last finite state.
    Full-sample sums use the same cumulative reduction as the prefix controls.
    """
    # Full-sample production kernels retain their own literal reduction order.
    if variant in full_sum.VARIANTS:
        return full_sum.update(inputs, weights, variant)
    if variant not in VARIANTS:
        raise ValueError(f"unknown compatibility weighting variant: {variant}")
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        first_moment, second_moment, _ = moments(inputs, weights)
        if variant.startswith("sum_"):
            first_moment, second_moment = first_moment[-1], second_moment[-1]

        signal = inputs.signal
        if variant.endswith("aliasing"):
            signal = signal + inputs.p1d * second_moment / (
                first_moment**2 * inputs.length
            )
        elif variant == "sum_historical":
            signal = signal + inputs.p1d / (first_moment * inputs.length)

        noise = inputs.variance / (first_moment * (inputs.length / inputs.pixel))
        result = signal / (signal + noise)
    if not np.all(np.isfinite(result)):
        raise FloatingPointError("nonfinite updated weights")
    return result


def coefficients(inputs, weights):
    """Full-sample I1, I2, I3, A and pixel power; no positivity substitution.

    Parameters
    ----------
    inputs : WeightInputs
        Magnitude-dependent source density, integration measure, pixel variance
        and auxiliary signal in angular/velocity units.
    weights : array_like, shape (n_magnitude,)
        Dimensionless source weights, retaining their common amplitude.

    Returns
    -------
    coefficients : ndarray, shape (5,)
        I1, I2, I3 in deg^-2 (km/s)^-1, aliasing coefficient A in deg^2, and
        pixel power in deg^2 km/s.

    Raises
    ------
    FloatingPointError :
        If the weighting arithmetic produces an invalid or nonfinite state.
    """
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        first_moment, second_moment, noise_moment = (
            a[-1] for a in moments(inputs, weights)
        )
        denominator = first_moment**2 * inputs.length
        result = np.array(
            [
                first_moment,
                second_moment,
                noise_moment,
                second_moment / denominator,
                noise_moment * inputs.pixel / denominator,
            ]
        )
    if not np.all(np.isfinite(result)):
        raise FloatingPointError("nonfinite final moments or coefficients")
    return result


def fixed_weights(inputs, variant="prefix_intrinsic", updates=3):
    """Return weights after exactly updates transitions following the seed.

    Parameters
    ----------
    inputs : WeightInputs
        Magnitude-dependent source density, integration measure, pixel variance
        and auxiliary signal in angular/velocity units.
    variant : str
        Named compatibility recurrence, selecting prefix or full-sample moments
        and the auxiliary signal prescription. Default is
        ``'prefix_intrinsic'``.
    updates : int
        Number of weight transitions after the initial seed. Default is ``3``.

    Returns
    -------
    weights : ndarray, shape (n_magnitude,)
        Dimensionless weights at the requested finite update count.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    if isinstance(updates, bool) or not isinstance(updates, int) or updates < 0:
        raise ValueError("updates must be a nonnegative integer")
    if variant not in VARIANTS:
        raise ValueError("unknown compatibility weighting variant")
    weights = seed(inputs)
    for _ in range(updates):
        weights = update(inputs, weights, variant)
    return weights


def prepare_contexts(pair_inputs):
    """Retain all forest-related pair contexts, including mixed-pair signals.

    Parameters
    ----------
    pair_inputs : dict
        Captured source and auxiliary-power inputs keyed by field-pair name.

    Returns
    -------
    contexts : dict of str to WeightInputs
        Forest-related inputs keyed by captured pair name.
    """
    return {
        name: WeightInputs.from_pair(row)
        for name, row in pair_inputs.items()
        if "density" in row
    }


def forest_noise(row, k, mu, aliasing, pixel_power):
    """Observed auto noise in (Mpc/h)^3; P1D response applied exactly once.

    Parameters
    ----------
    row : dict
        Captured forest auto context with redshift and angular/velocity
        conversion factors.
    k : array_like
        Comoving Fourier wavenumbers in h/Mpc; array shape follows the model or
        paired grid.
    mu : array_like
        Dimensionless line-of-sight direction cosines aligned with the Fourier
        grid.
    aliasing : float
        Aliasing coefficient A in deg^2.
    pixel_power : float
        Pixel-noise power in deg^2 km/s.

    Returns
    -------
    noise : ndarray
        Observed forest auto noise in (Mpc/h)^3, with the broadcast shape of k
        and mu.
    """
    velocity = row["_distance_to_velocity"]
    velocity_wavenumber = np.asarray(k) * np.asarray(mu) / velocity
    response = velocity_response(
        velocity_wavenumber,
        pixel_width_velocity=row["_pix_kms"],
        gaussian_sigma_velocity=row["_res_kms"],
    )
    return (
        (
            aliasing
            * default_p1d([], row["_z_mean"], velocity_wavenumber)
            * response**2
            + pixel_power
        )
        * row["_angle_to_distance"] ** 2
        / velocity
    )


def reassemble(arrays, task, pair_inputs, variant="prefix_intrinsic", updates=3):
    """Replace forest auto noise; retain saved intrinsic/cross/galaxy powers.

    Parameters
    ----------
    arrays : dict of str to ndarray
        Numerical evidence arrays; Fourier-cell axes and pair order follow the
        declared task. Powers use (Mpc/h)^3 and volumes use (Mpc/h)^3 unless
        separately labeled.
    task : dict
        Declared case, bin, selected field pairs, parameter order and validation
        thresholds.
    pair_inputs : dict
        Captured source and auxiliary-power inputs keyed by field-pair name.
    variant : str
        Named compatibility recurrence, selecting prefix or full-sample moments
        and the auxiliary signal prescription. Default is
        ``'prefix_intrinsic'``.
    updates : int
        Number of weight transitions after the initial seed. Default is ``3``.

    Returns
    -------
    total : ndarray, shape (n_cell, n_required_pair)
        Total powers with replaced forest-auto noise in (Mpc/h)^3.
    covariance : ndarray, shape (n_cell, n_selected_pair, n_selected_pair)
        Reconstructed Wick covariance in (Mpc/h)^6.
    states : dict
        Weights and coefficients for each captured forest context.
    noise_checks : dict
        Intrinsic powers and baseline/replacement noise for each forest auto.

    Notes
    -----
    Baseline noise is independently reconstructed from captured coefficients.
    The caller must check its intrinsic residual against saved mean powers before
    interpreting variants. Pair-specific auxiliary contexts are retained, while
    only auto-context coefficients enter the captured Wick covariance recipe.
    """
    contexts = prepare_contexts(pair_inputs)
    states = {}
    for name, inputs in contexts.items():
        weights = fixed_weights(inputs, variant, updates)
        states[name] = dict(weights=weights, coefficients=coefficients(inputs, weights))

    # Only forest autos change; all cross powers retain their captured values.
    total = arrays["total"].copy()
    noise_checks = {}
    for col, (i, j) in enumerate(task["required_pairs"]):
        field = task["fields"][i]
        if i != j or field["kind"] != "forest":
            continue
        name = field["id"] + "_" + field["id"]
        row = pair_inputs[name]
        baseline_noise = forest_noise(
            row,
            arrays["k"],
            arrays["mu"],
            row["_aliasing_weights"][-1],
            row["_effective_noise_power"][-1],
        )
        aliasing_coefficient, pixel_noise_power = states[name]["coefficients"][-2:]
        changed_noise = forest_noise(
            row, arrays["k"], arrays["mu"], aliasing_coefficient, pixel_noise_power
        )
        intrinsic = arrays["total"][:, col] - baseline_noise
        total[:, col] = intrinsic + changed_noise
        noise_checks[name] = dict(
            intrinsic=intrinsic, baseline_noise=baseline_noise, noise=changed_noise
        )

    covariance = wick(
        total,
        arrays["modes"],
        task["required_pairs"],
        task["selected_pairs"],
        len(task["fields"]),
    )
    return total, covariance, states, noise_checks
