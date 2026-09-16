"""Observed-coordinate field transfers; intrinsic models and supplied noise untouched."""

from dataclasses import dataclass

import numpy as np

from ._arrays import real_array, scalar
from .fields import PairSelection
from .geometry import SPEED_LIGHT_KMS, _immutable, _positive
from .kernels.response import _transfer

_FWHM_PER_SIGMA = 2 * np.sqrt(2 * np.log(2))


def _width(value, name):
    value = scalar(value, name)
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def _width_conversion(width, factor, name):
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        result = np.float64(width) * factor
    if not np.isfinite(result) or (width > 0 and result <= 0):
        raise ValueError(f"{name} conversion is not representable")
    return float(result)


def pixel_width_angstrom_to_velocity(pixel_width_angstrom, *, lambda_obs_angstrom):
    """Full pixel width to km/s; local narrow-width wavelength approximation."""
    width = _width(pixel_width_angstrom, "pixel_width_angstrom")
    wavelength = _positive(lambda_obs_angstrom, "lambda_obs_angstrom")
    return _width_conversion(width, SPEED_LIGHT_KMS / wavelength, "pixel width")


def gaussian_sigma_angstrom_to_velocity(sigma_angstrom, *, lambda_obs_angstrom):
    """Gaussian one-sigma wavelength width to km/s, using c*sigma/lambda_obs."""
    width = _width(sigma_angstrom, "sigma_angstrom")
    wavelength = _positive(lambda_obs_angstrom, "lambda_obs_angstrom")
    return _width_conversion(width, SPEED_LIGHT_KMS / wavelength, "Gaussian sigma")


def gaussian_fwhm_velocity_to_sigma(fwhm_velocity):
    """Gaussian velocity FWHM (km/s) to one-sigma (km/s), exactly once."""
    return _width_conversion(
        _width(fwhm_velocity, "fwhm_velocity"), 1 / _FWHM_PER_SIGMA, "Gaussian FWHM"
    )


def resolving_power_fwhm_to_sigma(resolving_power_fwhm):
    """R=lambda/FWHM_lambda to Gaussian sigma in km/s."""
    resolving = _positive(resolving_power_fwhm, "resolving_power_fwhm")
    return _positive((SPEED_LIGHT_KMS / resolving) / _FWHM_PER_SIGMA, "sigma_velocity")


def legacy_resolving_power_to_sigma(resolving_power_legacy):
    """Legacy compatibility only: sigma=c/R, using the new c=299792.458 km/s.

    lyaforecast consumes res_kms as sigma and uses c=299800 km/s. This helper
    preserves that interpretation, with an explicitly different light constant.
    """
    resolving = _positive(resolving_power_legacy, "resolving_power_legacy")
    return _positive(SPEED_LIGHT_KMS / resolving, "legacy sigma_velocity")


@dataclass(frozen=True)
class InstrumentResponse:
    """Explicit full pixel width and Gaussian one-sigma in km/s per field/bin.

    InstrumentResponse(0, 0) deliberately selects identity, including galaxies.
    """

    pixel_width_velocity: float
    gaussian_sigma_velocity: float

    def __post_init__(self):
        for name in ("pixel_width_velocity", "gaussian_sigma_velocity"):
            object.__setattr__(self, name, _width(getattr(self, name), name))


def velocity_response(
    k_parallel_velocity, *, pixel_width_velocity, gaussian_sigma_velocity
):
    """Return signed W(node,) for nonnegative paired velocity wavenumbers s/km.

    W=sinc(q*Delta_v/(2*pi))*exp(-(q*sigma_v)^2/2), with exact W(0)=1.
    Strong attenuation may underflow to zero. Inputs are never mutated.
    """
    q = real_array(k_parallel_velocity, "k_parallel_velocity")
    if q.ndim != 1 or not q.size or np.any(q < 0):
        raise ValueError("require nonempty 1D nonnegative velocity wavenumbers")
    pixel = _width(pixel_width_velocity, "pixel_width_velocity")
    sigma = _width(gaussian_sigma_velocity, "gaussian_sigma_velocity")
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        result, x, gaussian = _transfer(q, np.array([pixel]), np.array([sigma]))
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(gaussian)):
        raise ValueError("response arguments are not representable")
    if not np.all(np.isfinite(result)):
        raise ValueError("response result is not finite")
    return result[:, 0]


def prepare_response(fields, k, mu, *, a_v, settings):
    """Return immutable W(node,field) in explicit ObservedField order.

    k (h_fid/Mpc), mu in [0,1] are paired nonempty 1D observed nodes, including
    zero k. settings maps every field ID to InstrumentResponse, without defaults.
    Physical labels never merge settings; no AP mapping is applied here.
    """
    fields = PairSelection(fields).fields
    if not hasattr(settings, "keys") or set(settings) != {f.id for f in fields}:
        raise ValueError("settings must cover exactly all field IDs")
    if any(not isinstance(settings[f.id], InstrumentResponse) for f in fields):
        raise ValueError("each field requires an explicit InstrumentResponse")
    k, mu = real_array(k, "k"), real_array(mu, "mu")
    if (
        k.ndim != 1
        or not k.size
        or k.shape != mu.shape
        or np.any(k < 0)
        or np.any((mu < 0) | (mu > 1))
    ):
        raise ValueError("require paired 1D k>=0 and mu in [0,1]")
    a_v = _positive(a_v, "a_v")
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        q = (k * mu) / a_v
    if not np.all(np.isfinite(q)) or np.any((k > 0) & (mu > 0) & (q == 0)):
        raise ValueError("response velocity coordinates are not representable")
    pixel = np.array([settings[f.id].pixel_width_velocity for f in fields])
    sigma = np.array([settings[f.id].gaussian_sigma_velocity for f in fields])
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        result, x, gaussian = _transfer(q, pixel, sigma)
    if (
        not np.all(np.isfinite(x))
        or not np.all(np.isfinite(gaussian))
        or not np.all(np.isfinite(result))
    ):
        raise ValueError("response arguments/results are not representable")
    return _immutable(result)


def pair_response(response, selection):
    """Gather W_i*W_j in original required-pair order, shape (node,required_pair).

    Multiply intrinsic power and Jacobians by this result, then gather selected
    means. Supplied noise gets no automatic W; P1D consumers explicitly apply W^2.
    """
    if not isinstance(selection, PairSelection):
        raise ValueError("require PairSelection")
    response = real_array(response, "response")
    if (
        response.ndim != 2
        or response.shape[0] == 0
        or response.shape[1] != len(selection.fields)
    ):
        raise ValueError("response must have shape (node,field)")
    i, j = selection.required_pairs.T
    with np.errstate(over="ignore", invalid="ignore"):
        result = response[:, i] * response[:, j]
    if not np.all(np.isfinite(result)):
        raise ValueError("pair response is not finite")
    return result
