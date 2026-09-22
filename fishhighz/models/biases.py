"""Small NumPy-only bias and redshift-space evolution utilities.

The constants and powers are the established lyaforecast analytic prescription.
The optional tabulated route is a linear interpolator with linear extrapolation,
matching ``scipy.interpolate.interp1d(..., fill_value='extrapolate')`` without
making SciPy a FishHighz runtime dependency.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

OPTIONS = ("lya", "qso", "elgqso", "lbg", "lae")


def _z_array(redshift):
    array = np.asarray(redshift)
    if array.dtype.kind not in "iuf" or not np.all(np.isfinite(array)):
        raise ValueError("redshift must contain finite real values")
    return np.asarray(array, dtype=np.float64)


def _scalar_or_array(value, input_value):
    array = np.asarray(value, dtype=np.float64)
    return float(array) if np.asarray(input_value).ndim == 0 else array


@dataclass(frozen=True)
class LinearTabulatedBias:
    """Linear interpolation and endpoint-slope extrapolation of ``b(z)``."""

    redshifts: np.ndarray
    values: np.ndarray

    def __post_init__(self):
        z = np.asarray(self.redshifts, dtype=np.float64)
        values = np.asarray(self.values, dtype=np.float64)
        if z.ndim != 1 or values.shape != z.shape or len(z) < 2:
            raise ValueError(
                "tabulated bias requires matching 1D arrays with >=2 points"
            )
        if not np.all(np.isfinite(z)) or not np.all(np.isfinite(values)):
            raise ValueError("tabulated bias values must be finite")
        if np.any(np.diff(z) <= 0):
            raise ValueError("tabulated bias redshifts must be strictly increasing")
        # Keep source arrays owned and read-only; no sorting or duplicate merging.
        object.__setattr__(
            self, "redshifts", np.frombuffer(z.tobytes(), dtype=np.float64)
        )
        object.__setattr__(
            self, "values", np.frombuffer(values.tobytes(), dtype=np.float64)
        )

    def __call__(self, redshift):
        query = _z_array(redshift)
        flat = query.reshape(-1)
        result = np.interp(flat, self.redshifts, self.values)
        left_slope = (self.values[1] - self.values[0]) / (
            self.redshifts[1] - self.redshifts[0]
        )
        right_slope = (self.values[-1] - self.values[-2]) / (
            self.redshifts[-1] - self.redshifts[-2]
        )
        left = flat < self.redshifts[0]
        right = flat > self.redshifts[-1]
        result[left] = self.values[0] + left_slope * (flat[left] - self.redshifts[0])
        result[right] = self.values[-1] + right_slope * (
            flat[right] - self.redshifts[-1]
        )
        result = result.reshape(query.shape)
        return _scalar_or_array(result, redshift)


def linear_tabulated_bias(redshifts, values):
    """Return a callable reproducing lyaforecast's linear bias interpolation."""
    return LinearTabulatedBias(redshifts, values)


def _tracer_parameters(tracer):
    if not isinstance(tracer, str):
        raise ValueError(f"invalid biasing: {tracer!r}, select from: {OPTIONS}")
    if tracer == "lya":
        return 2.9, -0.1352, 2.33
    if tracer in ("qso", "elgqso"):
        return 1.44, 3.54, 2.33
    if tracer == "lbg":
        return 1.44, 3.48, 2.9
    if tracer == "lae":
        return 1.44, 2.2, 2.9
    raise ValueError(f"invalid biasing: {tracer}, select from: {OPTIONS}")


def analytic_density_bias(redshift, tracer):
    """Return the established analytic linear density bias ``b(z)``."""
    z = _z_array(redshift)
    alpha, bias_zref, zref = _tracer_parameters(tracer)
    result = bias_zref * ((1 + z) / (1 + zref)) ** alpha
    return _scalar_or_array(result, redshift)


def analytic_beta_rsd(redshift, tracer, growth_rate=None):
    """Return the established analytic Kaiser ``beta(z)`` prescription.

    Ly-alpha uses the fixed ``1.45`` normalization.  Other tracers use the
    supplied growth rate divided by their analytic density bias, as in
    lyaforecast's ``AnalyticBias`` implementation.
    """
    z = _z_array(redshift)
    if tracer == "lya":
        result = 1.45 * ((1 + z) / (1 + 2.33)) ** 0.0
    else:
        if growth_rate is None:
            raise ValueError("growth_rate is required for non-lya beta evolution")
        if callable(growth_rate):
            growth = growth_rate(redshift)
        elif hasattr(growth_rate, "growth_rate_at"):
            # Exact-redshift CAMB backgrounds intentionally reject interpolation.
            growth = np.asarray(
                [growth_rate.growth_rate_at(value) for value in z.reshape(-1)]
            )
            growth = growth.reshape(z.shape)
        else:
            growth = growth_rate
        growth = np.asarray(growth, dtype=np.float64)
        if growth.shape not in ((), z.shape):
            raise ValueError(
                "growth_rate must be scalar, callable, or match redshift shape"
            )
        result = growth / np.asarray(analytic_density_bias(redshift, tracer))
    if not np.all(np.isfinite(result)):
        raise ValueError("bias evolution produced nonfinite beta")
    return _scalar_or_array(result, redshift)


# Concise aliases for callers that use the quantity names directly.
density_bias = analytic_density_bias
beta_rsd = analytic_beta_rsd
analytic_beta = analytic_beta_rsd
get_density_bias = analytic_density_bias
get_beta_rsd = analytic_beta_rsd
tabulated_bias = linear_tabulated_bias


@dataclass
class AnalyticBias:
    """Stateful compatibility wrapper with optional tabulated density biases."""

    growth_rate: object = None
    _density_bias_functions: dict = field(default_factory=dict, init=False, repr=False)

    def set_density_bias_func(self, tracer, bias_func):
        if tracer not in OPTIONS:
            raise ValueError(f"tracer name must be in {OPTIONS}")
        if not callable(bias_func):
            raise ValueError("bias_func must be callable")
        self._density_bias_functions[tracer] = bias_func

    def density_bias(self, redshift, tracer):
        function = self._density_bias_functions.get(tracer)
        return (
            function(redshift)
            if function is not None
            else analytic_density_bias(redshift, tracer)
        )

    def beta_rsd(self, redshift, tracer):
        if tracer == "lya":
            return analytic_beta_rsd(redshift, tracer)
        function = self._density_bias_functions.get(tracer)
        if function is not None:
            bias = function(redshift)
        else:
            bias = analytic_density_bias(redshift, tracer)
        growth = self.growth_rate
        if growth is None:
            raise ValueError("growth_rate is required for non-lya beta evolution")
        if callable(growth):
            growth = growth(redshift)
        elif hasattr(growth, "growth_rate_at"):
            values = _z_array(redshift)
            growth = np.asarray(
                [growth.growth_rate_at(value) for value in values.reshape(-1)]
            ).reshape(values.shape)
        return _scalar_or_array(np.asarray(growth) / np.asarray(bias), redshift)

    # Names retained for direct migration from the reference implementation.
    get_density_bias = density_bias
    get_beta_rsd = beta_rsd
    _get_density_bias = density_bias
    _get_beta_rsd = beta_rsd

    def compute_bias(self, redshift, k_hmpc, mu, corr, linear=True):
        del k_hmpc, linear
        tracers = corr.split("_")
        if len(tracers) != 2:
            raise ValueError("corr must be of the form tracer1_tracer2")
        result = 1.0
        for tracer in tracers:
            tracer = "lya" if tracer.startswith("lya") else tracer
            result = (
                result
                * self.density_bias(redshift, tracer)
                * (1 + self.beta_rsd(redshift, tracer) * np.asarray(mu) ** 2)
            )
        return result


__all__ = [
    "OPTIONS",
    "AnalyticBias",
    "LinearTabulatedBias",
    "analytic_beta_rsd",
    "analytic_beta",
    "analytic_density_bias",
    "beta_rsd",
    "density_bias",
    "get_beta_rsd",
    "get_density_bias",
    "linear_tabulated_bias",
    "tabulated_bias",
]
