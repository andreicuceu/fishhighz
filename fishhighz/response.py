"""Observed-coordinate field transfers; intrinsic models and supplied noise untouched."""

from dataclasses import dataclass

import numpy as np

from ._arrays import real_array, scalar
from .fields import PairSelection
from .geometry import SPEED_LIGHT_KMS, _immutable, _positive
from .kernels.response import _transfer

_FWHM_PER_SIGMA = 2 * np.sqrt(2 * np.log(2))


def _width(value, name):
    """Validate a finite nonnegative instrumental width.

    Parameters
    ----------
    value : float
        Width in the units of the named quantity.
    name : str
        Quantity name used in validation errors.

    Returns
    -------
    width : float
        Validated width in unchanged units.

    Raises
    ------
    ValueError
        If a width is negative or nonfinite, or the conversion is not
        representable.
    """
    value = scalar(value, name)
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def _width_conversion(width, factor, name):
    """Apply a scalar width conversion while checking representability.

    Parameters
    ----------
    width : float
        Input width, already validated as finite and nonnegative.
    factor : float
        Multiplicative conversion factor in output units per input unit.
    name : str
        Width name used in errors.

    Returns
    -------
    width : float
        Converted width in output units.

    Raises
    ------
    ValueError
        If a width is negative or nonfinite, or the conversion is not
        representable.
    """
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        result = np.float64(width) * factor
    if not np.isfinite(result) or (width > 0 and result <= 0):
        raise ValueError(f"{name} conversion is not representable")
    return float(result)


def pixel_width_angstrom_to_velocity(pixel_width_angstrom, *, lambda_obs_angstrom):
    """Convert a full pixel width from wavelength to velocity units.

    Parameters
    ----------
    pixel_width_angstrom : float
        Nonnegative full pixel width in Angstrom.
    lambda_obs_angstrom : float
        Positive observed wavelength in Angstrom.

    Returns
    -------
    width : float
        Converted full pixel width in km/s.

    Raises
    ------
    ValueError
        If a width is negative or nonfinite, or the conversion is not
        representable.

    Notes
    -----
    Use the local narrow-width approximation c*width/lambda_obs.
    """
    width = _width(pixel_width_angstrom, "pixel_width_angstrom")
    wavelength = _positive(lambda_obs_angstrom, "lambda_obs_angstrom")
    return _width_conversion(width, SPEED_LIGHT_KMS / wavelength, "pixel width")


def gaussian_sigma_angstrom_to_velocity(sigma_angstrom, *, lambda_obs_angstrom):
    """Convert a Gaussian standard deviation from wavelength to velocity units.

    Parameters
    ----------
    sigma_angstrom : float
        Nonnegative Gaussian standard deviation in Angstrom.
    lambda_obs_angstrom : float
        Positive observed wavelength in Angstrom.

    Returns
    -------
    width : float
        Converted Gaussian standard deviation in km/s.

    Raises
    ------
    ValueError
        If a width is negative or nonfinite, or the conversion is not
        representable.

    Notes
    -----
    Use the local narrow-width approximation c*width/lambda_obs.
    """
    width = _width(sigma_angstrom, "sigma_angstrom")
    wavelength = _positive(lambda_obs_angstrom, "lambda_obs_angstrom")
    return _width_conversion(width, SPEED_LIGHT_KMS / wavelength, "Gaussian sigma")


def gaussian_fwhm_velocity_to_sigma(fwhm_velocity):
    """Convert a Gaussian velocity FWHM to its standard deviation.

    Parameters
    ----------
    fwhm_velocity : float
        Nonnegative full width at half maximum in km/s.

    Returns
    -------
    sigma : float
        Gaussian standard deviation in km/s.

    Raises
    ------
    ValueError
        If a width is negative or nonfinite, or the conversion is not
        representable.
    """
    return _width_conversion(
        _width(fwhm_velocity, "fwhm_velocity"), 1 / _FWHM_PER_SIGMA, "Gaussian FWHM"
    )


def resolving_power_fwhm_to_sigma(resolving_power_fwhm):
    """Convert FWHM-based resolving power to a Gaussian velocity dispersion.

    Parameters
    ----------
    resolving_power_fwhm : float
        Positive dimensionless R=lambda/FWHM_lambda.

    Returns
    -------
    sigma : float
        Gaussian standard deviation c/(R*sqrt(8*ln(2))) in km/s.

    Raises
    ------
    ValueError
        If resolving power or the converted dispersion is invalid.
    """
    resolving = _positive(resolving_power_fwhm, "resolving_power_fwhm")
    return _positive((SPEED_LIGHT_KMS / resolving) / _FWHM_PER_SIGMA, "sigma_velocity")


def legacy_resolving_power_to_sigma(resolving_power_legacy):
    """Convert legacy resolving power using the explicit sigma=c/R convention.

    Parameters
    ----------
    resolving_power_legacy : float
        Positive dimensionless legacy resolving power.

    Returns
    -------
    sigma : float
        Gaussian standard deviation in km/s.

    Raises
    ------
    ValueError
        If resolving power or the converted dispersion is invalid.

    Notes
    -----
    The legacy reference treats res_kms as sigma. This function retains that
    interpretation but uses c=299792.458 km/s instead of 299800 km/s.
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
        """Validate and normalize the pixel width and Gaussian standard deviation.

        Returns
        -------
        None
            Store finite nonnegative widths in km/s as Python floats.

        Raises
        ------
        ValueError
            If a width is negative or nonfinite, or the conversion is not
            representable.
        """
        for name in ("pixel_width_velocity", "gaussian_sigma_velocity"):
            object.__setattr__(self, name, _width(getattr(self, name), name))


def velocity_response(
    k_parallel_velocity, *, pixel_width_velocity, gaussian_sigma_velocity
):
    """Evaluate the signed instrumental amplitude response in velocity units.

    Parameters
    ----------
    k_parallel_velocity : array_like of shape (n_node,)
        Nonnegative line-of-sight wavenumbers in s/km, including zero.
    pixel_width_velocity : float
        Nonnegative full pixel width in km/s.
    gaussian_sigma_velocity : float
        Nonnegative Gaussian standard deviation in km/s.

    Returns
    -------
    response : ndarray of shape (n_node,)
        Dimensionless signed sinc-times-Gaussian response, with W(0)=1.

    Raises
    ------
    ValueError
        If inputs or intermediate response arguments are invalid.

    Notes
    -----
    The Gaussian multiplies the amplitude. Auto power uses W**2. Strong
    attenuation may underflow to zero; input arrays are not mutated.
    """
    velocity_wavenumber = real_array(k_parallel_velocity, "k_parallel_velocity")
    if (
        velocity_wavenumber.ndim != 1
        or not velocity_wavenumber.size
        or np.any(velocity_wavenumber < 0)
    ):
        raise ValueError("require nonempty 1D nonnegative velocity wavenumbers")
    pixel = _width(pixel_width_velocity, "pixel_width_velocity")
    sigma = _width(gaussian_sigma_velocity, "gaussian_sigma_velocity")
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        result, pixel_phase, gaussian_phase = _transfer(
            velocity_wavenumber, np.array([pixel]), np.array([sigma])
        )
    if not np.all(np.isfinite(pixel_phase)) or not np.all(np.isfinite(gaussian_phase)):
        raise ValueError("response arguments are not representable")
    if not np.all(np.isfinite(result)):
        raise ValueError("response result is not finite")
    return result[:, 0]


def prepare_response(fields, k, mu, *, a_v, settings):
    """Prepare immutable instrumental amplitudes in observed field order.

    Parameters
    ----------
    fields : iterable of ObservedField
        Ordered field identities; physical tracer labels do not merge settings.
    k : array_like of shape (n_node,)
        Nonnegative observed wavenumbers in h_fid/Mpc.
    mu : array_like of shape (n_node,)
        Paired direction cosines on [0, 1].
    a_v : float
        Positive H(z)/((1+z)*h_fid) in (km/s)/(Mpc/h_fid).
    settings : mapping of str to InstrumentResponse
        Explicit pixel and Gaussian widths for every field ID.

    Returns
    -------
    response : ndarray of shape (n_node, n_field)
        Immutable dimensionless amplitude responses.

    Raises
    ------
    ValueError
        If settings, coordinates, or converted response values are invalid.

    Notes
    -----
    No Alcock–Paczynski mapping is applied to the instrumental response.
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
        velocity_wavenumber = (k * mu) / a_v
    if not np.all(np.isfinite(velocity_wavenumber)) or np.any(
        (k > 0) & (mu > 0) & (velocity_wavenumber == 0)
    ):
        raise ValueError("response velocity coordinates are not representable")
    pixel = np.array([settings[f.id].pixel_width_velocity for f in fields])
    sigma = np.array([settings[f.id].gaussian_sigma_velocity for f in fields])
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        result, pixel_phase, gaussian_phase = _transfer(
            velocity_wavenumber, pixel, sigma
        )
    if (
        not np.all(np.isfinite(pixel_phase))
        or not np.all(np.isfinite(gaussian_phase))
        or not np.all(np.isfinite(result))
    ):
        raise ValueError("response arguments/results are not representable")
    return _immutable(result)


def pair_response(response, selection):
    """Multiply field amplitudes for each covariance-required spectrum.

    Parameters
    ----------
    response : array_like of shape (n_node, n_field)
        Finite dimensionless instrumental amplitudes in field order.
    selection : PairSelection
        Pair definitions including spectra needed to form covariance.

    Returns
    -------
    products : ndarray of shape (n_node, n_required_pair)
        Dimensionless W_i*W_j in required-pair order.

    Raises
    ------
    ValueError
        If the selection, array shape, or resulting products are invalid.

    Notes
    -----
    Apply this factor to intrinsic power and its derivatives. Supplied noise
    receives no automatic response; a P1D consumer explicitly applies W**2.
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
