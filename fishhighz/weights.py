"""Fixed per-field forest preparation from normalized arrays; no survey readers.

Cumulative formulas adapted from lyaforecast/weights.py, GPLv3, implementing
McDonald & Eisenstein (2007). No raw assets or normalization policies are copied.
"""

import os
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np

from ._arrays import integer, real_array, scalar
from .fields import ObservedField, PairSelection
from .geometry import SPEED_LIGHT_KMS, BinGeometry, _immutable, _positive
from .kernels.full_sum_weights import METHODS, fixed_weights, solve
from .kernels.integrated_weights import (
    fixed_weights_integrated,
    integrated_moments,
    solve_integrated,
)
from .kernels.weights import _integrals, _iterate
from .models.external import P3DProvider, PreparedP3D, evaluate_p1d, evaluate_p3d
from .response import InstrumentResponse, velocity_response


def _nonnegative(value, name):
    """Validate finite nonnegative real array values.

    Parameters
    ----------
    value : array_like
        Numeric input of arbitrary shape, in the named quantity's units.
    name : str
        Quantity name used in diagnostics.

    Returns
    -------
    array : ndarray
        Owned float64 copy with unchanged shape and units.

    Raises
    ------
    ValueError
        If values are nonnumeric, nonfinite, complex, Boolean, or negative.
    """
    array = real_array(value, name)
    if np.any(array < 0):
        raise ValueError(f"{name} must be nonnegative")
    return array


def _nonnegative_transient(value, name):
    """Validate finite nonnegative real array values for use within one call.

    Parameters
    ----------
    value : array_like
        Numeric input of arbitrary shape, in the named quantity's units.
    name : str
        Quantity name used in diagnostics.

    Returns
    -------
    array : ndarray
        The input itself when it is a C-contiguous float64 array (the
        (n_pixel, n_magnitude) arrays of an integrated source hold 1e7
        elements and the result is only read, never stored), otherwise an
        owned float64 copy; unchanged shape and units.

    Raises
    ------
    ValueError
        If values are nonnumeric, nonfinite, complex, Boolean, or negative.
    """
    if (
        isinstance(value, np.ndarray)
        and value.dtype == np.float64
        and value.flags.c_contiguous
        and value.size
    ):
        if not np.all(np.isfinite(value)):
            raise ValueError(f"{name} must be finite")
        if value.min() < 0:
            raise ValueError(f"{name} must be nonnegative")
        return value
    return _nonnegative(value, name)


def _redshift(value, name):
    """Validate a finite nonnegative scalar redshift.

    Parameters
    ----------
    value : float
        Dimensionless redshift.
    name : str
        Redshift label used in diagnostics.

    Returns
    -------
    redshift : float
        Validated dimensionless redshift.

    Raises
    ------
    ValueError
        If the value is not a finite nonnegative real scalar.
    """
    redshift = scalar(value, name)
    if redshift < 0:
        raise ValueError(f"{name} must be nonnegative")
    return redshift


def density_per_velocity(dndzdm, *, z_source):
    """Convert normalized dN/(dz dm deg²) to deg⁻² (km/s)⁻¹ mag⁻¹.

    Parameters
    ----------
    dndzdm : array_like of shape (n_magnitude,)
        Nonnegative source density dN/(dz dm deg^2).
    z_source : float
        Finite nonnegative dimensionless source redshift.

    Returns
    -------
    rho : ndarray of shape (n_magnitude,)
        Immutable source density in deg^-2 (km/s)^-1 mag^-1, preserving
        magnitude order.

    Raises
    ------
    ValueError
        If density or source redshift is invalid, or conversion is not
        representable.

    Notes
    -----
    No target normalization, integration, sorting, or area factor is applied.
    """
    row = _nonnegative(dndzdm, "dndzdm")
    if row.ndim != 1 or not row.size:
        raise ValueError("dndzdm must be nonempty 1D")
    source_redshift = _redshift(z_source, "z_source")
    with np.errstate(over="raise", invalid="raise", under="ignore"):
        try:
            result = row * ((1 + source_redshift) / SPEED_LIGHT_KMS)
        except FloatingPointError as error:
            raise ValueError("density conversion is not representable") from error
    if np.any((row > 0) & (result == 0)):
        raise ValueError("density conversion is not representable")
    return _immutable(result)


def _context(field, geometry, response):
    """Identify the field, geometry, and response defining forest-weight reuse.

    Parameters
    ----------
    field : ObservedField
        Forest sample identity, including its physical tracer and background
        population.
    geometry : BinGeometry
        Fixed forest evaluation geometry and fiducial coordinate conversion.
    response : InstrumentResponse
        Fixed instrumental widths in km/s, with positive forest pixel width.

    Returns
    -------
    context : tuple
        Field identity, z_eval, h_fid, a_v, d_deg, speed of light, and response
        in that order; retains the original objects and conventions.

    Raises
    ------
    ValueError
        If types are invalid or forest pixel width is not positive.
    """
    if not isinstance(field, ObservedField) or field.kind != "forest":
        raise ValueError("require an observed forest field")
    if not isinstance(geometry, BinGeometry):
        raise ValueError("require BinGeometry")
    if not isinstance(response, InstrumentResponse):
        raise ValueError("require InstrumentResponse")
    _positive(response.pixel_width_velocity, "forest pixel width")
    return (
        field,
        geometry.z_eval,
        geometry.h_fid,
        geometry.a_v,
        geometry.d_deg,
        geometry.speed_light_kms,
        response,
    )


@dataclass(frozen=True, init=False, eq=False)
class AuxiliarySamples:
    """Immutable auto-only fiducial S (deg² km/s), B (km/s), states and mode."""

    context: tuple
    k_t_deg: float
    k_p_velocity: float
    k: float
    mu: float
    signal: float
    alias: float
    theta_p3d: np.ndarray
    theta_p1d: np.ndarray


def sample_auxiliary(
    field,
    geometry,
    response,
    prepared,
    theta,
    *,
    p1d_model,
    p1d_parameters,
    theta_p1d,
    k_t_deg,
    k_p_velocity,
):
    """Query one matching auto and independent P1D once, outside all weight loops.

    Parameters
    ----------
    field : ObservedField
        Forest sample identity, including its physical tracer and background
        population.
    geometry : BinGeometry
        Fixed forest evaluation geometry and fiducial coordinate conversion.
    response : InstrumentResponse
        Fixed instrumental widths in km/s, with positive forest pixel width.
    prepared : PreparedP3D
        Intrinsic model routes containing the forest auto-spectrum.
    theta : array_like of shape (n_global_p3d,)
        P3D fiducial parameter values in their registry units.
    p1d_model : callable
        Independent intrinsic one-dimensional forest-power model.
    p1d_parameters : BoundParameters
        Explicit local bindings for the P1D model.
    theta_p1d : array_like of shape (n_global_p1d,)
        P1D fiducial parameter values in their own registry units.
    k_t_deg : float
        Nonnegative transverse angular wavenumber in deg^-1.
    k_p_velocity : float
        Nonnegative line-of-sight velocity wavenumber in s/km; both auxiliary
        components cannot vanish.

    Returns
    -------
    samples : AuxiliarySamples
        Immutable fiducial state with response-smoothed signal S in deg^2 km/s,
        alias B in km/s, and corresponding observed k and mu.

    Raises
    ------
    ValueError
        If the preparation context, auxiliary mode, model domain, or positive
        signal/alias requirement fails.

    Notes
    -----
    Preserve original field indices/bindings, routing only this auto. S is the
    caller's full intrinsic prediction, not a promise of legacy 'smooth' power.
    Explicit legacy examples are (2.4, .00035) for BAO and (7, .001) for P1D.
    Auxiliary nodes add no forecast modes and must be within provider domains.
    """
    context = _context(field, geometry, response)
    transverse_angular_wavenumber = _redshift(k_t_deg, "k_t_deg")
    parallel_velocity_wavenumber = _redshift(k_p_velocity, "k_p_velocity")
    if transverse_angular_wavenumber == parallel_velocity_wavenumber == 0:
        raise ValueError(f"{field.id}: auxiliary mode cannot be zero")
    if not isinstance(prepared, PreparedP3D) or field not in prepared.selection.fields:
        raise ValueError(f"{field.id}: require matching PreparedP3D field")
    i = prepared.selection.fields.index(field)
    owners = [
        r.provider for r in prepared.routes if (r.pairs == (i, i)).all(axis=1).any()
    ]
    if len(owners) != 1:
        raise ValueError(f"{field.id}: missing auxiliary auto provider")
    owner = owners[0]
    auto = PreparedP3D(
        prepared.registry,
        PairSelection(prepared.selection.fields, [(i, i)]),
        [P3DProvider(owner.label, owner.model, owner.parameters, [(i, i)])],
    )
    with np.errstate(over="ignore", invalid="ignore"):
        parallel = geometry.a_v * parallel_velocity_wavenumber
        wavenumber = np.hypot(parallel, transverse_angular_wavenumber / geometry.d_deg)
        direction_cosine = parallel / wavenumber
    try:
        amplitude_response = velocity_response(
            [parallel_velocity_wavenumber],
            pixel_width_velocity=response.pixel_width_velocity,
            gaussian_sigma_velocity=response.gaussian_sigma_velocity,
        )[0]
        signal_power = evaluate_p3d(
            auto, theta, geometry.z_eval, [wavenumber], [direction_cosine]
        )[0, 0]
        alias_power = evaluate_p1d(
            p1d_model,
            p1d_parameters,
            theta_p1d,
            geometry.z_eval,
            [parallel_velocity_wavenumber],
        )[0]
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            signal_power = _positive(
                signal_power * amplitude_response**2 * geometry.a_v / geometry.d_deg**2,
                "auxiliary S",
            )
            alias_power = _positive(alias_power * amplitude_response**2, "auxiliary B")
    except (ValueError, OverflowError) as error:
        raise ValueError(
            f"{field.id}, auxiliary ({transverse_angular_wavenumber}, {parallel_velocity_wavenumber}), provider {owner.label}: {error}"
        ) from error
    result = object.__new__(AuxiliarySamples)
    for name, value in dict(
        context=context,
        k_t_deg=transverse_angular_wavenumber,
        k_p_velocity=parallel_velocity_wavenumber,
        k=float(wavenumber),
        mu=float(direction_cosine),
        signal=signal_power,
        alias=alias_power,
        theta_p3d=_immutable(theta),
        theta_p1d=_immutable(theta_p1d),
    ).items():
        object.__setattr__(result, name, value)
    return result


class _ForestContextMixin:
    """Shared reuse check for prepared forest weights carrying ``context``."""

    def validate_context(self, field, geometry, response):
        """Reject reuse with a different field, evaluation geometry, or response.

        Parameters
        ----------
        field : ObservedField
            Forest sample identity, including its physical tracer and background
            population.
        geometry : BinGeometry
            Fixed forest evaluation geometry and fiducial coordinate conversion.
        response : InstrumentResponse
            Fixed instrumental widths in km/s, with positive forest pixel width.

        Returns
        -------
        None
            Validate exact equality with the stored preparation context.

        Raises
        ------
        ValueError
            If the field, fiducial units, local geometry, or response differs.
        """
        if self.context != _context(field, geometry, response):
            raise ValueError("forest preparation field/geometry/response mismatch")


@dataclass(frozen=True, init=False, eq=False)
class ForestWeights(_ForestContextMixin):
    """Fixed immutable input arrays, prefix integrals and noise coefficients.

    A is deg²; P_pixel is deg² km/s. I1/I2 are density per source velocity;
    I3 adds dimensionless pixel variance. weight_changes is max absolute change
    per update, not a convergence assertion. Construct with prepare_forest_weights.
    """

    context: tuple
    z_source: float
    magnitudes: np.ndarray
    quadrature: np.ndarray
    rho: np.ndarray
    variance: np.ndarray
    length_velocity: float
    method: str
    iterations: int | None
    auxiliary: AuxiliarySamples | None
    signal: float | None
    alias: float | None
    weights: np.ndarray
    weight_changes: np.ndarray
    I1: np.ndarray
    I2: np.ndarray
    I3: np.ndarray
    A: float
    P_pixel: float
    convergence: object


def prepare_forest_weights(
    field,
    geometry,
    response,
    *,
    z_source,
    magnitudes,
    quadrature,
    rho,
    variance,
    length_velocity,
    method,
    weights=None,
    iterations=None,
    signal=None,
    alias=None,
    auxiliary=None,
    rtol=1e-4,
    min_updates=3,
    stable_steps=3,
    max_updates=96,
):
    """Prepare from same-shaped 1D m, positive dm weights, rho and pixel variance.

    Parameters
    ----------
    field : ObservedField
        Forest sample identity, including its physical tracer and background
        population.
    geometry : BinGeometry
        Fixed forest evaluation geometry and fiducial coordinate conversion.
    response : InstrumentResponse
        Fixed instrumental widths in km/s, with positive forest pixel width.
    z_source : float
        Dimensionless source redshift, strictly above the forest evaluation
        redshift.
    magnitudes : array_like of shape (n_magnitude,)
        Finite strictly increasing source magnitudes.
    quadrature : array_like of shape (n_magnitude,)
        Positive integration weights in magnitudes.
    rho : array_like of shape (n_magnitude,)
        Nonnegative source density in deg^-2 (km/s)^-1 mag^-1, with positive
        support.
    variance : array_like of shape (n_magnitude,)
        Nonnegative dimensionless pixel-noise variance.
    length_velocity : float
        Positive forest length in km/s.
    method : str
        Explicit supplied, legacy, inverse_variance, early_lyaforecast, or
        mcdonald weighting prescription.
    weights : array_like of shape (n_magnitude,), optional
        Dimensionless supplied weights; default None. Required only for
        method='supplied'.
    iterations : int, optional
        Nonnegative fixed update count. Default None requests adaptive stopping
        for the full-sample methods; legacy requires an explicit count.
    signal : float, optional
        Positive response-smoothed reference S in deg^2 km/s. Default None; used
        only by iterative prescriptions.
    alias : float, optional
        Positive response-smoothed reference B or B_star in km/s. Default None.
    auxiliary : AuxiliarySamples, optional
        Matching fiducial samples supplying S/B instead of explicit signal and
        alias. Default None.
    rtol : float, default=1e-4
        Relative tolerance for adaptive weight and noise-coefficient changes.
    min_updates : int, default=3
        Minimum adaptive update count before nomination of convergence.
    stable_steps : int, default=3
        Consecutive stable updates required before confirmation.
    max_updates : int, default=96
        Maximum adaptive update count.

    Returns
    -------
    prepared : ForestWeights
        Immutable input arrays, dimensionless weights, cumulative I1/I2/I3 in
        deg^-2 (km/s)^-1, A in deg^2, pixel power in deg^2 km/s, and convergence
        metadata.

    Raises
    ------
    ValueError
        If inputs or method settings conflict, arithmetic is unrepresentable,
        weighted support is absent, or adaptive convergence is not confirmed.

    Notes
    -----
    'early_lyaforecast' and 'mcdonald' use full-sample moments and require
    positive signal/alias or matching fiducial auxiliary samples. With no
    iterations they require confirmed adaptive convergence (rtol=1e-4, at
    least three updates and three stable transitions, doubled-count
    confirmation, cap 96). Failure raises; no last iterate becomes a forecast.
    Explicit iterations selects a labeled fixed-count diagnostic.

    Require explicit method='supplied' with weights, or 'legacy' with iterations
    and positive S/B (signal/alias), optionally from sample_auxiliary.
    'inverse_variance' requires only alias=B_star, the positive finite P1D
    already response-smoothed at a fixed weighting mode, and computes
    B_star/(B_star+pixel_width*variance) once on positive density support.
    Its weights, iterations, signal and auxiliary arguments must be None.
    z_source exceeds z_eval. No model is queried here. Final noise uses its
    own mode-dependent intrinsic P1D and response; B_star only sets the weights.
    """
    context = _context(field, geometry, response)
    z_source = _redshift(z_source, "z_source")
    if z_source <= geometry.z_eval:
        raise ValueError("z_source must exceed z_eval")

    magnitude_grid = real_array(magnitudes, "magnitudes")
    magnitude_weights = _nonnegative(quadrature, "quadrature")
    source_density = _nonnegative(rho, "rho")
    pixel_variance = _nonnegative(variance, "variance")
    if (
        magnitude_grid.ndim != 1
        or not magnitude_grid.size
        or any(
            sample_array.shape != magnitude_grid.shape
            for sample_array in (
                magnitude_weights,
                source_density,
                pixel_variance,
            )
        )
        or np.any(magnitude_grid[1:] <= magnitude_grid[:-1])
        or np.any(magnitude_weights <= 0)
        or not np.any(source_density > 0)
    ):
        raise ValueError(
            "require ordered 1D magnitudes, positive quadrature and density support"
        )
    length = _positive(length_velocity, "length_velocity")
    pixel = response.pixel_width_velocity

    # Resolve the chosen prescription before any nonlinear weight update.
    if method == "supplied":
        if weights is None or any(
            x is not None for x in (iterations, signal, alias, auxiliary)
        ):
            raise ValueError("supplied weights conflict with legacy settings")
        forest_weights = _nonnegative(weights, "weights")
        if forest_weights.shape != magnitude_grid.shape:
            raise ValueError("weights shape mismatch")
        changes = np.empty(0)
    elif method == "legacy" or method in METHODS:
        if weights is not None:
            raise ValueError("legacy method conflicts with supplied weights")
        if method == "legacy" or iterations is not None:
            iterations = integer(iterations, "iterations")
        if auxiliary is not None:
            if (
                not isinstance(auxiliary, AuxiliarySamples)
                or auxiliary.context != context
                or signal is not None
                or alias is not None
            ):
                raise ValueError("auxiliary context/settings mismatch")
            signal, alias = auxiliary.signal, auxiliary.alias
        signal, alias = _positive(signal, "S"), _positive(alias, "B")
    elif method == "inverse_variance":
        if any(x is not None for x in (weights, iterations, signal, auxiliary)):
            raise ValueError(
                "inverse_variance requires weights, iterations, signal and auxiliary "
                "to be None"
            )
        alias = _positive(alias, "B_star (alias)")
        changes = np.empty(0)
    else:
        raise ValueError(
            "method must be supplied, legacy, inverse_variance, early_lyaforecast or mcdonald"
        )

    convergence = None
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            masses = source_density * magnitude_weights
            if np.any((source_density > 0) & (masses == 0)):
                raise ValueError("quadrature masses are not representable")

            if method in METHODS:
                inputs = SimpleNamespace(
                    density=source_density,
                    quadrature=magnitude_weights,
                    variance=pixel_variance,
                    length=length,
                    pixel=pixel,
                    signal=signal,
                    p1d=alias,
                )
                if iterations is None:
                    convergence = solve(
                        inputs,
                        METHODS[method],
                        rtol=rtol,
                        min_updates=min_updates,
                        stable_steps=stable_steps,
                        max_updates=max_updates,
                    )
                    if convergence["status"] != "converged":
                        raise ValueError(
                            f"{field.id}: {method} {convergence['status']}: "
                            f"{convergence['reason']}"
                        )
                    forest_weights = convergence["weights"]
                else:
                    forest_weights = fixed_weights(inputs, METHODS[method], iterations)
                    convergence = dict(status="fixed_count", updates=iterations)
                changes = np.empty(0)
            elif method == "legacy":
                forest_weights, changes = _iterate(
                    masses, pixel_variance, length, pixel, signal, alias, iterations
                )
                if np.any((source_density > 0) & (forest_weights == 0)):
                    raise ValueError("legacy weights are not representable")
            elif method == "inverse_variance":
                support_mask = source_density > 0
                forest_weights = np.zeros_like(source_density)
                with np.errstate(under="raise"):
                    instrumental_power = pixel * pixel_variance[support_mask]
                forest_weights[support_mask] = alias / (alias + instrumental_power)
                if np.any(
                    (forest_weights[support_mask] <= 0)
                    | ~np.isfinite(forest_weights[support_mask])
                ):
                    raise ValueError("inverse_variance weights are not representable")

            # Normalize moments only after the weights have been determined.
            (
                first_moment,
                second_moment,
                noise_moment,
                aliasing_coefficient,
                pixel_power,
            ) = _integrals(masses, forest_weights, pixel_variance, length, pixel)
            if (
                not np.isfinite(aliasing_coefficient)
                or aliasing_coefficient <= 0
                or not np.isfinite(pixel_power)
                or pixel_power < 0
                or (
                    np.any((masses > 0) & (forest_weights > 0) & (pixel_variance > 0))
                    and pixel_power == 0
                )
                or first_moment[-1] <= 0
                or second_moment[-1] <= 0
            ):
                raise ValueError(
                    "integrals/coefficients are not representable or lack weighted support"
                )
    except FloatingPointError as error:
        raise ValueError(
            f"{field.id}: weighting integrals/coefficients are not representable "
            f"or lack support: {error}"
        ) from error

    result = object.__new__(ForestWeights)
    values = dict(
        context=context,
        z_source=z_source,
        magnitudes=magnitude_grid,
        quadrature=magnitude_weights,
        rho=source_density,
        variance=pixel_variance,
        length_velocity=length,
        method=method,
        iterations=iterations,
        auxiliary=auxiliary,
        signal=signal,
        alias=alias,
        weights=forest_weights,
        weight_changes=changes,
        I1=first_moment,
        I2=second_moment,
        I3=noise_moment,
        A=float(aliasing_coefficient),
        P_pixel=float(pixel_power),
        convergence=convergence,
    )
    for name, value in values.items():
        object.__setattr__(
            result, name, _immutable(value) if isinstance(value, np.ndarray) else value
        )
    return result


# Attributes an integrated forest source must expose (duck-typed so that this
# module does not import fishhighz.forest_integration).
_INTEGRATED_SOURCE_ATTRIBUTES = (
    "nodes",
    "density",
    "magnitudes",
    "quadrature",
    "variance",
    "pixel_width_velocity",
    "info",
    "measure",
)
INTEGRATED_WEIGHT_OPTIONS = frozenset(
    {
        "method",
        "iterations",
        "signal",
        "alias",
        "auxiliary",
        "rtol",
        "min_updates",
        "stable_steps",
        "max_updates",
    }
)
_INTEGRATED_NODE_ATTRIBUTES = (
    "geom",
    "y_index",
    "z_q",
    "lam_obs",
    "z_pix",
    "zq_nodes",
    "info",
)


def _require_integrated_source(source):
    """Check that an object exposes the integrated forest-source attributes.

    Parameters
    ----------
    source : object
        Candidate integrated forest source, such as
        ``fishhighz.forest_integration.IntegratedForestSource``.

    Returns
    -------
    source : object
        The unchanged input.

    Raises
    ------
    ValueError
        If a required attribute of the source or of ``source.nodes`` is absent.
    """
    missing = [a for a in _INTEGRATED_SOURCE_ATTRIBUTES if not hasattr(source, a)]
    if not missing:
        missing = [
            f"nodes.{a}"
            for a in _INTEGRATED_NODE_ATTRIBUTES
            if not hasattr(source.nodes, a)
        ]
    if missing:
        raise ValueError(f"integrated forest source lacks attributes {missing}")
    return source


@dataclass(frozen=True, init=False, eq=False)
class IntegratedForestWeights(_ForestContextMixin):
    """Fixed weights and noise coefficients of an integrated forest source.

    The source redshift z_q and forest pixel are integrated inside the bin with
    measure mu (n_pixel, n_magnitude) in deg^-2. N1 = sum(mu w) and
    N2 = sum(mu w^2) are in deg^-2, and N3 = sum(mu w^2 v) is in deg^-2 times
    the dimensionless pixel variance v (the weights w are dimensionless);
    A = N2/N1^2 is in deg^2 and P_pixel = pixel_width N3/N1^2 in deg^2 km/s, so
    ``forest_noise`` consumes them exactly as for ForestWeights. z_eff = sum(z_pix mu w)/N1 is the
    weight-averaged pixel redshift. ``weights`` is dimensionless with shape
    (n_pixel, n_magnitude); ``convergence`` has the same keys as the central
    solver record, with its 'weights' entry being the same (n_pixel,
    n_magnitude) array. Construct with prepare_integrated_forest_weights.
    """

    context: tuple
    method: str
    iterations: int | None
    auxiliary: AuxiliarySamples | None
    signal: float
    alias: float
    convergence: object
    A: float
    P_pixel: float
    N1: float
    N2: float
    N3: float
    z_eff: float
    weights: np.ndarray


def prepare_integrated_forest_weights(
    field,
    geometry,
    response,
    source,
    *,
    method,
    rtol=1e-4,
    min_updates=3,
    stable_steps=3,
    max_updates=96,
    auxiliary=None,
    signal=None,
    alias=None,
    iterations=None,
    backend=None,
):
    """Prepare early-lyaforecast weights for a source integrated in z_q and pixel.

    Parameters
    ----------
    field : ObservedField
        Forest sample identity, including its physical tracer and background
        population.
    geometry : BinGeometry
        Fixed forest evaluation geometry and fiducial coordinate conversion.
    response : InstrumentResponse
        Fixed instrumental widths in km/s, with positive forest pixel width.
    source : IntegratedForestSource
        Duck-typed integrated source exposing ``measure`` (n_pixel,
        n_magnitude) in deg^-2, ``variance`` (n_pixel, n_magnitude) of the
        dimensionless pixel noise, ``magnitudes`` and ``quadrature``
        (n_magnitude,), ``nodes.z_pix`` (n_pixel,) and ``pixel_width_velocity``
        in km/s, which must equal the response pixel width.
    method : str
        Must be 'early_lyaforecast'.
    rtol : float, default=1e-4
        Relative tolerance for adaptive weight and noise-coefficient changes.
    min_updates : int, default=3
        Minimum adaptive update count before nomination of convergence.
    stable_steps : int, default=3
        Consecutive stable updates required before confirmation.
    max_updates : int, default=96
        Maximum adaptive update count.
    auxiliary : AuxiliarySamples, optional
        Matching fiducial samples supplying S/B instead of explicit signal and
        alias. Default None.
    signal : float, optional
        Positive response-smoothed reference S in deg^2 km/s. Default None.
    alias : float, optional
        Positive response-smoothed reference B in km/s. Default None.
    iterations : int, optional
        Nonnegative fixed update count. Default None requests adaptive stopping
        with confirmed convergence.
    backend : {'reference', 'numpy', 'numba'}, optional
        Implementation of the recurrence. Default None reads
        ``FISHHIGHZ_INTEGRATED_BACKEND`` and falls back to 'numba' (which
        itself falls back to 'numpy' when Numba is unavailable). 'reference'
        runs the central kernels ``full_sum_weights.solve``/``fixed_weights``
        unchanged on the flattened arrays; 'numpy' and 'numba' run the
        dedicated streaming recurrence of ``kernels.integrated_weights``, which
        reproduces it to floating-point summation order (relative 1e-12) with
        identical update counts and stopping decisions, except when ``rtol``
        coincides with a stopping metric to within its rounding (~1e-16
        absolute), where the backends may stop at a different state that is
        still within ``rtol``.

    Returns
    -------
    prepared : IntegratedForestWeights
        Immutable weights (n_pixel, n_magnitude), N1, N2, N3, A in deg^2,
        pixel power in deg^2 km/s, z_eff and convergence metadata.

    Raises
    ------
    ValueError
        If the method is not 'early_lyaforecast', the source or the S/B inputs
        are invalid, arithmetic is unrepresentable, weighted support is absent,
        or adaptive convergence is not confirmed.

    Notes
    -----
    The 'sum_historical' recurrence of ``kernels.full_sum_weights`` is applied
    to the flattened arrays with density := measure, quadrature := 1,
    length := 1 and variance := variance. The central moments L*sum(rho q w)
    then become N1 = sum(mu w), so that S = P + B/N1, the per-pixel noise is
    pixel_width v/N1, A = N2/N1^2 and P_pixel = pixel_width N3/N1^2. The
    recurrence uses only full sums, never magnitude prefixes, so the (pixel,
    magnitude) ordering of the flattened arrays affects only floating-point
    summation. Arrays are flattened explicitly (C order, as views) and weights
    are reshaped afterwards.

    The default backends run the dedicated streaming recurrence of
    ``kernels.integrated_weights``, which reproduces the central kernel
    (same update count, candidate, stopping decisions and record keys) with
    full sums that are more accurate than the central kernel's cumulative
    prefix sums (about 1e-16 against n*eps relative, i.e. up to about 1e-12
    at the 1e7 nodes of a production source); ``backend='reference'`` runs the
    central kernels themselves.
    """
    if method != "early_lyaforecast":
        raise ValueError(
            "integrated forest weights support only method='early_lyaforecast'"
        )
    context = _context(field, geometry, response)
    _require_integrated_source(source)

    measure = _nonnegative_transient(source.measure, "integrated measure")
    pixel_variance = _nonnegative_transient(source.variance, "integrated variance")
    magnitude_grid = real_array(source.magnitudes, "magnitudes")
    magnitude_weights = _nonnegative(source.quadrature, "quadrature")
    pixel_redshift = real_array(source.nodes.z_pix, "z_pix")
    if (
        measure.ndim != 2
        or not measure.size
        or pixel_variance.shape != measure.shape
        or magnitude_grid.shape != (measure.shape[1],)
        or magnitude_weights.shape != magnitude_grid.shape
        or pixel_redshift.shape != (measure.shape[0],)
        or np.any(magnitude_grid[1:] <= magnitude_grid[:-1])
        or np.any(magnitude_weights <= 0)
        or not np.any(measure > 0)
    ):
        raise ValueError(
            "require 2D measure/variance of shape (n_pixel, n_magnitude), "
            "ordered magnitudes, positive quadrature and measure support"
        )
    pixel = _positive(source.pixel_width_velocity, "integrated pixel width")
    if not np.isclose(pixel, response.pixel_width_velocity, rtol=1e-12, atol=0):
        raise ValueError("source pixel width differs from the response pixel width")

    if iterations is not None:
        iterations = integer(iterations, "iterations")
    if auxiliary is not None:
        if (
            not isinstance(auxiliary, AuxiliarySamples)
            or auxiliary.context != context
            or signal is not None
            or alias is not None
        ):
            raise ValueError("auxiliary context/settings mismatch")
        signal, alias = auxiliary.signal, auxiliary.alias
    signal, alias = _positive(signal, "S"), _positive(alias, "B")

    if backend is None:
        backend = os.environ.get("FISHHIGHZ_INTEGRATED_BACKEND", "numba")
    if backend not in ("reference", "numpy", "numba"):
        raise ValueError(
            "FISHHIGHZ_INTEGRATED_BACKEND must be reference, numpy or numba"
        )

    # Explicit C-order flattening (views; the validated arrays are contiguous).
    flat_measure = measure.reshape(-1)
    flat_variance = pixel_variance.reshape(-1)

    try:
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            if backend == "reference":
                # Central kernels, unchanged, on the flattened arrays.
                inputs = SimpleNamespace(
                    density=flat_measure,
                    quadrature=np.ones_like(flat_measure),
                    variance=flat_variance,
                    length=1.0,  # N1 already carries the forest length
                    pixel=pixel,
                    signal=signal,
                    p1d=alias,
                )
                if iterations is None:
                    solution = solve(
                        inputs,
                        METHODS[method],
                        rtol=rtol,
                        min_updates=min_updates,
                        stable_steps=stable_steps,
                        max_updates=max_updates,
                    )
                else:
                    flat_weights = fixed_weights(inputs, METHODS[method], iterations)
            elif iterations is None:
                solution = solve_integrated(
                    flat_measure,
                    flat_variance,
                    pixel=pixel,
                    signal=signal,
                    alias=alias,
                    rtol=rtol,
                    min_updates=min_updates,
                    stable_steps=stable_steps,
                    max_updates=max_updates,
                    backend=backend,
                )
            else:
                flat_weights = fixed_weights_integrated(
                    flat_measure,
                    flat_variance,
                    pixel=pixel,
                    signal=signal,
                    alias=alias,
                    updates=iterations,
                    backend=backend,
                )

            if iterations is None:
                if solution["status"] != "converged":
                    raise ValueError(
                        f"{field.id}: {method} {solution['status']}: "
                        f"{solution['reason']}"
                    )
                flat_weights = solution["weights"]
                convergence = dict(solution)
            else:
                convergence = dict(status="fixed_count", updates=iterations)

            # Normalize moments only after the weights have been determined.
            if backend == "reference":
                (
                    first_moment,
                    second_moment,
                    noise_moment,
                    aliasing_coefficient,
                    pixel_power,
                ) = _integrals(flat_measure, flat_weights, flat_variance, 1.0, pixel)
                total_first, total_second, total_noise = (
                    first_moment[-1],
                    second_moment[-1],
                    noise_moment[-1],
                )
            else:
                (
                    total_first,
                    total_second,
                    total_noise,
                    aliasing_coefficient,
                    pixel_power,
                ) = integrated_moments(flat_measure, flat_weights, flat_variance, pixel)
            if (
                not np.isfinite(aliasing_coefficient)
                or aliasing_coefficient <= 0
                or not np.isfinite(pixel_power)
                or pixel_power < 0
                or (
                    pixel_power == 0
                    and np.any(
                        (flat_measure > 0) & (flat_weights > 0) & (flat_variance > 0)
                    )
                )
                or total_first <= 0
                or total_second <= 0
            ):
                raise ValueError(
                    "integrals/coefficients are not representable or lack weighted support"
                )

            weights = flat_weights.reshape(measure.shape)

            # sum_p z_p sum_j mu_pj w_pj without a full-size temporary.
            effective_redshift = (
                np.einsum("pj,pj->p", measure, weights) @ pixel_redshift / total_first
            )
    except FloatingPointError as error:
        raise ValueError(
            f"{field.id}: weighting integrals/coefficients are not representable "
            f"or lack support: {error}"
        ) from error

    weights = _immutable(weights)
    if iterations is None:
        # Share the read-only (n_pixel, n_magnitude) array instead of keeping a
        # second flattened copy of the converged weights.
        convergence["weights"] = weights

    result = object.__new__(IntegratedForestWeights)
    values = dict(
        context=context,
        method=method,
        iterations=iterations,
        auxiliary=auxiliary,
        signal=signal,
        alias=alias,
        convergence=convergence,
        A=float(aliasing_coefficient),
        P_pixel=float(pixel_power),
        N1=float(total_first),
        N2=float(total_second),
        N3=float(total_noise),
        z_eff=float(effective_redshift),
        weights=weights,
    )
    for name, value in values.items():
        object.__setattr__(result, name, value)
    return result
