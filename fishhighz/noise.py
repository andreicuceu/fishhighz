"""Fixed sampling noise in (Mpc/h_fid)^3; response applies only to aliasing.

Forest coefficients follow lyaforecast/covariance.py (GPLv3). Full supplied
noise replaces generated noise and is validated independently of the signal.
"""

from dataclasses import dataclass

import numpy as np

from ._arrays import integer, real_array
from .covariance import _validate_field_power
from .fields import PairSelection
from .geometry import BinGeometry, _immutable, _positive
from .response import velocity_response
from .weights import ForestWeights, _nonnegative, density_per_velocity


def local_galaxy_density(dndzdm, quadrature, geometry):
    """Convert a local surface-density distribution to comoving number density.

    Parameters
    ----------
    dndzdm : array_like of shape (n_magnitude,)
        Source density dN/(dz dm deg^2) at geometry.z_eval.
    quadrature : array_like of shape (n_magnitude,)
        Positive magnitude integration weights, in magnitudes.
    geometry : BinGeometry
        Fixed background conversions at the bin evaluation redshift.

    Returns
    -------
    n_bar : float
        Positive comoving density in (h_fid/Mpc)^3.

    Raises
    ------
    ValueError
        If geometry, density, quadrature, or the converted density is invalid.

    Notes
    -----
    This local approximation uses neither bin-integrated counts nor survey
    area. Supply n_bar directly when another density convention is required.
    """
    if not isinstance(geometry, BinGeometry):
        raise ValueError("require BinGeometry")
    rho = density_per_velocity(dndzdm, z_source=geometry.z_eval)
    magnitude_weights = _nonnegative(quadrature, "quadrature")
    if magnitude_weights.shape != rho.shape or np.any(magnitude_weights <= 0):
        raise ValueError("quadrature must match density with positive weights")
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        number_density = (
            np.sum(magnitude_weights * rho) * geometry.a_v / geometry.d_deg**2
        )
    return _positive(number_density, "local galaxy density")


def galaxy_noise(n_bar):
    """Return unsmoothed Poisson auto-noise for a discrete tracer.

    Parameters
    ----------
    n_bar : float
        Positive comoving number density in (h_fid/Mpc)^3.

    Returns
    -------
    power : float
        Positive shot-noise power 1/n_bar in (Mpc/h_fid)^3, without area or
        volume rescaling.

    Raises
    ------
    ValueError
        If the density or its reciprocal is not positive and representable.
    """
    number_density = _positive(n_bar, "n_bar")
    with np.errstate(over="ignore", under="ignore"):
        return _positive(np.float64(1) / number_density, "galaxy noise")


@dataclass(frozen=True, init=False, eq=False)
class ForestNoise:
    """Immutable aliasing, unsmoothed pixel, and total arrays (node,) in P3D units."""

    aliasing: np.ndarray
    pixel: np.ndarray
    total: np.ndarray


def forest_noise(prepared, field, geometry, response, k, mu, p1d):
    """Evaluate forest aliasing and pixel-noise power at observed Fourier nodes.

    Parameters
    ----------
    prepared : ForestWeights
        Fixed magnitude weights and integrated noise coefficients.
    field : ObservedField
        Forest identity matching the weight preparation.
    geometry : BinGeometry
        Geometry matching the preparation context.
    response : InstrumentResponse
        Pixel and Gaussian widths matching the preparation context.
    k : array_like of shape (n_node,)
        Nonnegative observed wavenumbers in h_fid/Mpc.
    mu : array_like of shape (n_node,)
        Paired direction cosines on [0, 1].
    p1d : array_like of shape (n_node,)
        Intrinsic one-dimensional power in km/s at k*mu/a_v.

    Returns
    -------
    noise : ForestNoise
        Immutable aliasing, pixel, and total arrays, each of shape (n_node,) in
        (Mpc/h_fid)^3.

    Raises
    ------
    ValueError
        If contexts or arrays disagree, or coordinate/noise conversion is not
        representable.

    Notes
    -----
    Both coefficients receive d_deg**2/a_v once. Only aliasing receives the
    squared instrumental response; supplied pixel noise is not smoothed.
    Zero P1D and sinc-null aliasing are valid.
    """
    if not isinstance(prepared, ForestWeights):
        raise ValueError("require ForestWeights")
    prepared.validate_context(field, geometry, response)
    k, mu = _nonnegative(k, "k"), _nonnegative(mu, "mu")
    p1d = _nonnegative(p1d, "P1D")
    if (
        k.ndim != 1
        or not k.size
        or mu.shape != k.shape
        or p1d.shape != k.shape
        or np.any(mu > 1)
    ):
        raise ValueError("require matching nonempty 1D k, mu in [0,1], P1D")
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        velocity_wavenumber = k * mu / geometry.a_v
    if np.any((k > 0) & (mu > 0) & (velocity_wavenumber == 0)):
        raise ValueError("noise coordinates are not representable")
    amplitude_response = velocity_response(
        velocity_wavenumber,
        pixel_width_velocity=response.pixel_width_velocity,
        gaussian_sigma_velocity=response.gaussian_sigma_velocity,
    )
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        conversion = geometry.d_deg**2 / geometry.a_v
        aliasing = prepared.A * p1d * amplitude_response**2 * conversion
        pixel = np.full(k.shape, prepared.P_pixel * conversion)
        total = aliasing + pixel
    if (
        not all(np.all(np.isfinite(x)) for x in (aliasing, pixel, total))
        or np.any((p1d > 0) & (amplitude_response != 0) & (aliasing == 0))
        or (prepared.P_pixel > 0 and np.any(pixel == 0))
    ):
        raise ValueError(f"{field.id}: noise conversion is not representable")
    result = object.__new__(ForestNoise)
    for name, value in dict(aliasing=aliasing, pixel=pixel, total=total).items():
        object.__setattr__(result, name, _immutable(value))
    return result


def prepare_noise(
    selection, n_node, *, diagonal=None, independent_sampling=None, full=None
):
    """Prepare the complete known-noise contribution in required-pair order.

    Parameters
    ----------
    selection : PairSelection
        Selected observables and their covariance dependencies.
    n_node : int
        Positive number of Fourier nodes.
    diagonal : mapping, optional
        Auto-noise arrays of shape (n_node,) in (Mpc/h_fid)^3 for exactly the
        active field IDs. Default None.
    independent_sampling : bool, optional
        Must explicitly be True for generated diagonal noise; default None.
    full : array_like of shape (n_node, n_required_pair), optional
        Complete replacement noise in (Mpc/h_fid)^3. Default None selects
        generated noise.

    Returns
    -------
    noise : ndarray of shape (n_node, n_required_pair)
        Immutable float64 noise, including covariance-required spectra.

    Raises
    ------
    ValueError
        If noise inputs conflict, have invalid shape, or violate the positive-
        semidefinite field-noise contract.

    Notes
    -----
    Full replacement noise excludes diagonal and independence arguments.
    Signed cross noise, singular positive-semidefinite matrices, and zero noise
    are accepted without jitter under the normalized 64*eps64 convention.
    """
    if not isinstance(selection, PairSelection):
        raise ValueError("require PairSelection")
    n_node = integer(n_node, "n_node", minimum=1)
    shape = (n_node, len(selection.required_pairs))
    if full is not None:
        if diagonal is not None or independent_sampling is not None:
            raise ValueError(
                "full noise replaces generated inputs; conflicting settings"
            )
        out = real_array(full, "full noise")
        if out.shape != shape:
            raise ValueError(f"full noise must have exact shape {shape}")
    else:
        if independent_sampling is not True:
            raise ValueError(
                "generated noise requires explicit independent_sampling=True"
            )
        active = np.unique(selection.selected_pairs)
        ids = {selection.fields[i].id for i in active}
        if not hasattr(diagonal, "keys") or set(diagonal) != ids:
            raise ValueError("diagonal noise must cover exactly active field IDs")
        out = np.zeros(shape)
        for col, (i, j) in enumerate(selection.required_pairs):
            if i == j:
                row = _nonnegative(diagonal[selection.fields[i].id], "diagonal noise")
                if row.shape != (n_node,):
                    raise ValueError("diagonal noise must have shape (n_node,)")
                out[:, col] = row
    try:
        _validate_field_power(out, selection)
    except ValueError as error:
        raise ValueError(f"known noise: {error}") from error
    return _immutable(out)
