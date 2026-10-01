"""Fixed per-field forest preparation from normalized arrays; no survey readers.

Cumulative formulas adapted from lyaforecast/weights.py, GPLv3, implementing
McDonald & Eisenstein (2007). No raw assets or normalization policies are copied.
"""

from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np

from ._arrays import integer, real_array, scalar
from .fields import ObservedField, PairSelection
from .geometry import SPEED_LIGHT_KMS, BinGeometry, _immutable, _positive
from .kernels.full_sum_weights import METHODS, fixed_weights, solve
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


@dataclass(frozen=True, init=False, eq=False)
class ForestWeights:
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
