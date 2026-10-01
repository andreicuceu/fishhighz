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
    """Validate finite real redshifts without changing their shape.

    Parameters
    ----------
    redshift : float or array_like
        Dimensionless evaluation redshift; scalar or arbitrary array shape.

    Returns
    -------
    redshifts : ndarray
        Float64 redshifts with the input shape.

    Raises
    ------
    ValueError
        If any redshift is nonfinite or not real numeric data.
    """
    array = np.asarray(redshift)
    if array.dtype.kind not in "iuf" or not np.all(np.isfinite(array)):
        raise ValueError("redshift must contain finite real values")
    return np.asarray(array, dtype=np.float64)


def _scalar_or_array(value, input_value):
    """Preserve the scalar-versus-array convention of an evaluation input.

    Parameters
    ----------
    value : array_like
        Numeric result, retaining its physical units.
    input_value : array_like
        Original query used to decide whether the result is scalar.

    Returns
    -------
    result : float or ndarray
        Python float for scalar input, otherwise a float64 array.
    """
    array = np.asarray(value, dtype=np.float64)
    return float(array) if np.asarray(input_value).ndim == 0 else array


@dataclass(frozen=True)
class LinearTabulatedBias:
    """Linear interpolation and endpoint-slope extrapolation of ``b(z)``."""

    redshifts: np.ndarray
    values: np.ndarray

    def __post_init__(self):
        """Validate and freeze tabulated redshifts and density biases.

        Returns
        -------
        None
            Store immutable float64 copies of both one-dimensional arrays.

        Raises
        ------
        ValueError
            If arrays differ in shape, contain nonfinite entries, have fewer than
            two nodes, or redshifts are not strictly increasing.
        """
        redshift_grid = np.asarray(self.redshifts, dtype=np.float64)
        values = np.asarray(self.values, dtype=np.float64)
        if (
            redshift_grid.ndim != 1
            or values.shape != redshift_grid.shape
            or len(redshift_grid) < 2
        ):
            raise ValueError(
                "tabulated bias requires matching 1D arrays with >=2 points"
            )
        if not np.all(np.isfinite(redshift_grid)) or not np.all(np.isfinite(values)):
            raise ValueError("tabulated bias values must be finite")
        if np.any(np.diff(redshift_grid) <= 0):
            raise ValueError("tabulated bias redshifts must be strictly increasing")
        # Keep source arrays owned and read-only; no sorting or duplicate merging.
        object.__setattr__(
            self, "redshifts", np.frombuffer(redshift_grid.tobytes(), dtype=np.float64)
        )
        object.__setattr__(
            self, "values", np.frombuffer(values.tobytes(), dtype=np.float64)
        )

    def __call__(self, redshift):
        """Interpolate density bias and extrapolate using the endpoint slopes.

        Parameters
        ----------
        redshift : float or array_like
            Dimensionless evaluation redshift; scalar or arbitrary array shape.

        Returns
        -------
        bias : float or ndarray
            Dimensionless bias, preserving scalar input or the redshift array shape.

        Raises
        ------
        ValueError
            If query redshifts are not finite real values.
        """
        query = _z_array(redshift)
        flat = query.reshape(-1)
        result = np.interp(flat, self.redshifts, self.values)
        left_slope = (self.values[1] - self.values[0]) / (
            self.redshifts[1] - self.redshifts[0]
        )
        right_slope = (self.values[-1] - self.values[-2]) / (
            self.redshifts[-1] - self.redshifts[-2]
        )
        left_mask = flat < self.redshifts[0]
        right_mask = flat > self.redshifts[-1]
        result[left_mask] = self.values[0] + left_slope * (
            flat[left_mask] - self.redshifts[0]
        )
        result[right_mask] = self.values[-1] + right_slope * (
            flat[right_mask] - self.redshifts[-1]
        )
        result = result.reshape(query.shape)
        return _scalar_or_array(result, redshift)


def linear_tabulated_bias(redshifts, values):
    """Prepare linear density-bias interpolation with endpoint extrapolation.

    Parameters
    ----------
    redshifts : array_like of shape (n_redshift,)
        Increasing finite redshift nodes, with at least two entries.
    values : array_like of shape (n_redshift,)
        Dimensionless density biases at the corresponding redshifts.

    Returns
    -------
    bias : LinearTabulatedBias
        Callable preserving scalar or array query shape.

    Raises
    ------
    ValueError
        If the table does not contain matching finite arrays on increasing
        nodes.
    """
    return LinearTabulatedBias(redshifts, values)


def _tracer_parameters(tracer):
    """Select the reference bias and power-law redshift evolution.

    Parameters
    ----------
    tracer : {'lya', 'qso', 'elgqso', 'lbg', 'lae'}
        Physical tracer whose bias prescription is evaluated.

    Returns
    -------
    exponent : float
        Dimensionless exponent of (1+z)/(1+z_ref).
    bias_reference : float
        Dimensionless density bias at the reference redshift.
    redshift_reference : float
        Dimensionless reference redshift.

    Raises
    ------
    ValueError
        If the tracer is unknown or not a string.
    """
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
    """Evaluate the adopted power-law density-bias evolution.

    Parameters
    ----------
    redshift : float or array_like
        Dimensionless evaluation redshift; scalar or arbitrary array shape.
    tracer : {'lya', 'qso', 'elgqso', 'lbg', 'lae'}
        Physical tracer whose bias prescription is evaluated.

    Returns
    -------
    bias : float or ndarray
        Dimensionless bias, preserving scalar input or the redshift array shape.

    Raises
    ------
    ValueError
        If the tracer or redshift input is invalid.

    Notes
    -----
    Constants and exponents follow the established lyaforecast prescription.
    """
    redshift_grid = _z_array(redshift)
    evolution_exponent, reference_bias, reference_redshift = _tracer_parameters(tracer)
    result = (
        reference_bias
        * ((1 + redshift_grid) / (1 + reference_redshift)) ** evolution_exponent
    )
    return _scalar_or_array(result, redshift)


def analytic_beta_rsd(redshift, tracer, growth_rate=None):
    """Evaluate the adopted Kaiser redshift-space distortion parameter.

    Parameters
    ----------
    redshift : float or array_like
        Dimensionless evaluation redshift; scalar or arbitrary array shape.
    tracer : {'lya', 'qso', 'elgqso', 'lbg', 'lae'}
        Physical tracer whose bias prescription is evaluated.
    growth_rate : float, array_like, callable, or CAMBBackground, optional
        Dimensionless f(z), supplied directly, through a callable, or through
        exact prepared redshift lookup. Default None; required for non-Ly-alpha
        tracers.

    Returns
    -------
    beta : float or ndarray
        Dimensionless beta with the redshift input shape, or a float for scalar
        input.

    Raises
    ------
    ValueError
        If growth is absent for a non-Ly-alpha tracer, has incompatible shape,
        or gives nonfinite beta.

    Notes
    -----
    Ly-alpha uses the fixed normalization 1.45. Other tracers use f(z)/b(z).
    Prepared CAMB growth values require an exact requested redshift.
    """
    redshift_grid = _z_array(redshift)
    if tracer == "lya":
        result = 1.45 * ((1 + redshift_grid) / (1 + 2.33)) ** 0.0
    else:
        if growth_rate is None:
            raise ValueError("growth_rate is required for non-lya beta evolution")
        if callable(growth_rate):
            growth = growth_rate(redshift)
        elif hasattr(growth_rate, "growth_rate_at"):
            # Exact-redshift CAMB backgrounds intentionally reject interpolation.
            growth = np.asarray(
                [
                    growth_rate.growth_rate_at(value)
                    for value in redshift_grid.reshape(-1)
                ]
            )
            growth = growth.reshape(redshift_grid.shape)
        else:
            growth = growth_rate
        growth = np.asarray(growth, dtype=np.float64)
        if growth.shape not in ((), redshift_grid.shape):
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
        """Register a density-bias callable for one physical tracer.

        Parameters
        ----------
        tracer : {'lya', 'qso', 'elgqso', 'lbg', 'lae'}
            Physical tracer whose bias prescription is evaluated.
        bias_func : callable
            Function of scalar or array redshift returning dimensionless density
            bias.

        Returns
        -------
        None
            Replace the stored callable for this tracer.

        Raises
        ------
        ValueError
            If the tracer is unknown or the supplied object is not callable.
        """
        if tracer not in OPTIONS:
            raise ValueError(f"tracer name must be in {OPTIONS}")
        if not callable(bias_func):
            raise ValueError("bias_func must be callable")
        self._density_bias_functions[tracer] = bias_func

    def density_bias(self, redshift, tracer):
        """Evaluate a registered density-bias function or the analytic prescription.

        Parameters
        ----------
        redshift : float or array_like
            Dimensionless evaluation redshift; scalar or arbitrary array shape.
        tracer : {'lya', 'qso', 'elgqso', 'lbg', 'lae'}
            Physical tracer whose bias prescription is evaluated.

        Returns
        -------
        bias : float or ndarray
            Dimensionless bias, preserving scalar input or the redshift array shape.
        """
        function = self._density_bias_functions.get(tracer)
        return (
            function(redshift)
            if function is not None
            else analytic_density_bias(redshift, tracer)
        )

    def beta_rsd(self, redshift, tracer):
        """Evaluate beta using the stored growth rate and selected density bias.

        Parameters
        ----------
        redshift : float or array_like
            Dimensionless evaluation redshift; scalar or arbitrary array shape.
        tracer : {'lya', 'qso', 'elgqso', 'lbg', 'lae'}
            Physical tracer whose bias prescription is evaluated.

        Returns
        -------
        beta : float or ndarray
            Dimensionless beta with the redshift input shape, or a float for scalar
            input.

        Raises
        ------
        ValueError
            If the tracer/redshift is invalid or required growth is absent.

        Notes
        -----
        Ly-alpha retains beta=1.45; other tracers use the registered density bias
        when present and the stored f(z).
        """
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
        """Multiply the two tracer Kaiser factors for a named correlation.

        Parameters
        ----------
        redshift : float or array_like
            Dimensionless evaluation redshift; scalar or arbitrary array shape.
        k_hmpc : array_like
            Wavenumbers in h/Mpc; accepted for compatibility and unused.
        mu : float or array_like
            Dimensionless direction cosine, broadcast against the redshift-dependent
            biases.
        corr : str
            Two tracer names joined by an underscore; names beginning with lya use
            the Ly-alpha prescription.
        linear : bool, default=True
            Compatibility argument; this method always uses linear Kaiser factors.

        Returns
        -------
        factor : float or ndarray
            Dimensionless pair factor b_i*(1+beta_i*mu**2)*b_j*(1+beta_j*mu**2),
            with the broadcast input shape.

        Raises
        ------
        ValueError
            If corr does not identify two tracers or required bias inputs are
            invalid.
        """
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
