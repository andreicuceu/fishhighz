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
    """Local n_bar in (h_fid/Mpc)^3 from normalized dN/(dz dm deg²) at z_eval.

    This is a local-density approximation, not a bin-integrated count/volume.
    No area or forest length is used. Supply n_bar directly for other conventions.
    """
    if not isinstance(geometry, BinGeometry):
        raise ValueError("require BinGeometry")
    rho = density_per_velocity(dndzdm, z_source=geometry.z_eval)
    q = _nonnegative(quadrature, "quadrature")
    if q.shape != rho.shape or np.any(q <= 0):
        raise ValueError("quadrature must match density with positive weights")
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        n = np.sum(q * rho) * geometry.a_v / geometry.d_deg**2
    return _positive(n, "local galaxy density")


def galaxy_noise(n_bar):
    """Poisson auto-noise 1/n_bar, no smoothing, volume or area rescaling."""
    n = _positive(n_bar, "n_bar")
    with np.errstate(over="ignore", under="ignore"):
        return _positive(np.float64(1) / n, "galaxy noise")


@dataclass(frozen=True, init=False, eq=False)
class ForestNoise:
    """Immutable aliasing, unsmoothed pixel, and total arrays (node,) in P3D units."""

    aliasing: np.ndarray
    pixel: np.ndarray
    total: np.ndarray


def forest_noise(prepared, field, geometry, response, k, mu, p1d):
    """Evaluate noise at paired observed k (h_fid/Mpc), mu; intrinsic P1D is km/s.

    Supply P1D at q=k*mu/a_v, e.g. via evaluate_p1d with its independent binding.
    Both coefficients get d_deg²/a_v once; only aliasing gets W(q)². Zero P1D
    and sinc-null aliasing are valid. Preparation context must match exactly.
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
        q = k * mu / geometry.a_v
    if np.any((k > 0) & (mu > 0) & (q == 0)):
        raise ValueError("noise coordinates are not representable")
    w = velocity_response(
        q,
        pixel_width_velocity=response.pixel_width_velocity,
        gaussian_sigma_velocity=response.gaussian_sigma_velocity,
    )
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        conversion = geometry.d_deg**2 / geometry.a_v
        aliasing = prepared.A * p1d * w**2 * conversion
        pixel = np.full(k.shape, prepared.P_pixel * conversion)
        total = aliasing + pixel
    if (
        not all(np.all(np.isfinite(x)) for x in (aliasing, pixel, total))
        or np.any((p1d > 0) & (w != 0) & (aliasing == 0))
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
    """Return immutable known noise (node, required_pair) in required-pair order.

    Generated path: diagonal maps exactly active field IDs to nonnegative (node,)
    arrays, and independent_sampling=True must be explicit. Full path: full is
    the entire packed replacement, with neither diagonal nor independence set.
    Signed cross terms, singular and zero PSD matrices are allowed. Validation
    uses the covariance normalized 64*eps64 convention without jitter.
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
