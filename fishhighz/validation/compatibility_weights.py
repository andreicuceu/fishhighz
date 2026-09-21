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
        """Preserve each pair's own uniform compatibility magnitude measure."""
        m = np.asarray(row["magnitudes"], dtype=float)
        return cls(
            m,
            np.asarray(row["density"], dtype=float),
            np.full(m.shape, m[1] - m[0]),
            np.asarray(row["variance"], dtype=float),
            float(row["forest_length"]),
            float(row["_pix_kms"]),
            float(row["auxiliary_signal"]),
            float(row["auxiliary_p1d"]),
        )


def seed(inputs):
    """Common legacy seed, preserving the saved division order."""
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        power = inputs.p1d / inputs.pixel
        return power / (power + inputs.variance)


def moments(inputs, weights):
    """Cumulative I1/I2/I3 with literal legacy multiplication/reduction order."""
    d, q, v = inputs.density, inputs.quadrature, inputs.variance
    return (
        np.cumsum(d * weights * q),
        np.cumsum(d * weights**2 * q),
        np.cumsum(d * weights**2 * v * q),
    )


def update(inputs, weights, variant):
    """One simultaneous update, retaining amplitude and signed arithmetic.

    Floating-point failures propagate; callers must report the last finite state.
    Full-sample sums use the same cumulative reduction as the prefix controls.
    """
    if variant in full_sum.VARIANTS:
        return full_sum.update(inputs, weights, variant)
    if variant not in VARIANTS:
        raise ValueError(f"unknown compatibility weighting variant: {variant}")
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        j1, j2, _ = moments(inputs, weights)
        if variant.startswith("sum_"):
            j1, j2 = j1[-1], j2[-1]
        signal = inputs.signal
        if variant.endswith("aliasing"):
            signal = signal + inputs.p1d * j2 / (j1**2 * inputs.length)
        elif variant == "sum_historical":
            signal = signal + inputs.p1d / (j1 * inputs.length)
        noise = inputs.variance / (j1 * (inputs.length / inputs.pixel))
        result = signal / (signal + noise)
    if not np.all(np.isfinite(result)):
        raise FloatingPointError("nonfinite updated weights")
    return result


def coefficients(inputs, weights):
    """Full-sample I1, I2, I3, A and pixel power; no positivity substitution."""
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        i1, i2, i3 = (a[-1] for a in moments(inputs, weights))
        denominator = i1**2 * inputs.length
        result = np.array(
            [i1, i2, i3, i2 / denominator, i3 * inputs.pixel / denominator]
        )
    if not np.all(np.isfinite(result)):
        raise FloatingPointError("nonfinite final moments or coefficients")
    return result


def fixed_weights(inputs, variant="prefix_intrinsic", updates=3):
    """Return weights after exactly updates transitions following the seed."""
    if isinstance(updates, bool) or not isinstance(updates, int) or updates < 0:
        raise ValueError("updates must be a nonnegative integer")
    if variant not in VARIANTS:
        raise ValueError("unknown compatibility weighting variant")
    weights = seed(inputs)
    for _ in range(updates):
        weights = update(inputs, weights, variant)
    return weights


def prepare_contexts(pair_inputs):
    """Retain all forest-related pair contexts, including mixed-pair signals."""
    return {
        name: WeightInputs.from_pair(row)
        for name, row in pair_inputs.items()
        if "density" in row
    }


def forest_noise(row, k, mu, aliasing, pixel_power):
    """Observed auto noise in (Mpc/h)^3; P1D response applied exactly once."""
    velocity = row["_distance_to_velocity"]
    q = np.asarray(k) * np.asarray(mu) / velocity
    response = velocity_response(
        q, pixel_width_velocity=row["_pix_kms"], gaussian_sigma_velocity=row["_res_kms"]
    )
    return (
        (aliasing * default_p1d([], row["_z_mean"], q) * response**2 + pixel_power)
        * row["_angle_to_distance"] ** 2
        / velocity
    )


def reassemble(arrays, task, pair_inputs, variant="prefix_intrinsic", updates=3):
    """Replace forest auto noise; retain saved intrinsic/cross/galaxy powers.

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
        a, p = states[name]["coefficients"][-2:]
        changed_noise = forest_noise(row, arrays["k"], arrays["mu"], a, p)
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
