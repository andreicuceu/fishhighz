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
    a = real_array(value, name)
    if np.any(a < 0):
        raise ValueError(f"{name} must be nonnegative")
    return a


def _redshift(value, name):
    z = scalar(value, name)
    if z < 0:
        raise ValueError(f"{name} must be nonnegative")
    return z


def density_per_velocity(dndzdm, *, z_source):
    """Convert normalized dN/(dz dm deg²) to deg⁻² (km/s)⁻¹ mag⁻¹.

    No target normalization, integration, sorting, or area factor is applied.
    """
    row = _nonnegative(dndzdm, "dndzdm")
    if row.ndim != 1 or not row.size:
        raise ValueError("dndzdm must be nonempty 1D")
    z = _redshift(z_source, "z_source")
    with np.errstate(over="raise", invalid="raise", under="ignore"):
        try:
            result = row * ((1 + z) / SPEED_LIGHT_KMS)
        except FloatingPointError as error:
            raise ValueError("density conversion is not representable") from error
    if np.any((row > 0) & (result == 0)):
        raise ValueError("density conversion is not representable")
    return _immutable(result)


def _context(field, geometry, response):
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

    Preserve original field indices/bindings, routing only this auto. S is the
    caller's full intrinsic prediction, not a promise of legacy 'smooth' power.
    Explicit legacy examples are (2.4, .00035) for BAO and (7, .001) for P1D.
    Auxiliary nodes add no forecast modes and must be within provider domains.
    """
    context = _context(field, geometry, response)
    kt = _redshift(k_t_deg, "k_t_deg")
    kp = _redshift(k_p_velocity, "k_p_velocity")
    if kt == kp == 0:
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
        parallel = geometry.a_v * kp
        k = np.hypot(parallel, kt / geometry.d_deg)
        mu = parallel / k
    try:
        w = velocity_response(
            [kp],
            pixel_width_velocity=response.pixel_width_velocity,
            gaussian_sigma_velocity=response.gaussian_sigma_velocity,
        )[0]
        s = evaluate_p3d(auto, theta, geometry.z_eval, [k], [mu])[0, 0]
        b = evaluate_p1d(p1d_model, p1d_parameters, theta_p1d, geometry.z_eval, [kp])[0]
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            s = _positive(s * w**2 * geometry.a_v / geometry.d_deg**2, "auxiliary S")
            b = _positive(b * w**2, "auxiliary B")
    except (ValueError, OverflowError) as error:
        raise ValueError(
            f"{field.id}, auxiliary ({kt}, {kp}), provider {owner.label}: {error}"
        ) from error
    result = object.__new__(AuxiliarySamples)
    for name, value in dict(
        context=context,
        k_t_deg=kt,
        k_p_velocity=kp,
        k=float(k),
        mu=float(mu),
        signal=s,
        alias=b,
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
        """Reject reuse across field, evaluation geometry, h or response settings."""
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
    m = real_array(magnitudes, "magnitudes")
    q = _nonnegative(quadrature, "quadrature")
    r = _nonnegative(rho, "rho")
    v = _nonnegative(variance, "variance")
    if (
        m.ndim != 1
        or not m.size
        or any(a.shape != m.shape for a in (q, r, v))
        or np.any(m[1:] <= m[:-1])
        or np.any(q <= 0)
        or not np.any(r > 0)
    ):
        raise ValueError(
            "require ordered 1D magnitudes, positive quadrature and density support"
        )
    length = _positive(length_velocity, "length_velocity")
    pixel = response.pixel_width_velocity
    if method == "supplied":
        if weights is None or any(
            x is not None for x in (iterations, signal, alias, auxiliary)
        ):
            raise ValueError("supplied weights conflict with legacy settings")
        w = _nonnegative(weights, "weights")
        if w.shape != m.shape:
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
            masses = r * q
            if np.any((r > 0) & (masses == 0)):
                raise ValueError("quadrature masses are not representable")
            if method in METHODS:
                inputs = SimpleNamespace(
                    density=r,
                    quadrature=q,
                    variance=v,
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
                    w = convergence["weights"]
                else:
                    w = fixed_weights(inputs, METHODS[method], iterations)
                    convergence = dict(status="fixed_count", updates=iterations)
                changes = np.empty(0)
            elif method == "legacy":
                w, changes = _iterate(
                    masses, v, length, pixel, signal, alias, iterations
                )
                if np.any((r > 0) & (w == 0)):
                    raise ValueError("legacy weights are not representable")
            elif method == "inverse_variance":
                support = r > 0
                w = np.zeros_like(r)
                with np.errstate(under="raise"):
                    instrumental_power = pixel * v[support]
                w[support] = alias / (alias + instrumental_power)
                if np.any((w[support] <= 0) | ~np.isfinite(w[support])):
                    raise ValueError("inverse_variance weights are not representable")
            i1, i2, i3, a, p = _integrals(masses, w, v, length, pixel)
            if (
                not np.isfinite(a)
                or a <= 0
                or not np.isfinite(p)
                or p < 0
                or (np.any((masses > 0) & (w > 0) & (v > 0)) and p == 0)
                or i1[-1] <= 0
                or i2[-1] <= 0
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
        magnitudes=m,
        quadrature=q,
        rho=r,
        variance=v,
        length_velocity=length,
        method=method,
        iterations=iterations,
        auxiliary=auxiliary,
        signal=signal,
        alias=alias,
        weights=w,
        weight_changes=changes,
        I1=i1,
        I2=i2,
        I3=i3,
        A=float(a),
        P_pixel=float(p),
        convergence=convergence,
    )
    for name, value in values.items():
        object.__setattr__(
            result, name, _immutable(value) if isinstance(value, np.ndarray) else value
        )
    return result
