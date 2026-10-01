"""Fixed single-bin geometry in physical background and fiducial h units."""

from dataclasses import dataclass

import numpy as np

from ._arrays import integer, real_array, scalar
from .grids import IntegrationGrid

SPEED_LIGHT_KMS = 299792.458
LYA_REST_ANGSTROM = 1215.67


def _immutable(value):
    # A bytes-backed view cannot have WRITEABLE re-enabled by a consumer.
    """Copy float64 values into an array backed by immutable bytes.

    Parameters
    ----------
    value : array_like
        Numeric values of arbitrary shape; physical units are retained.

    Returns
    -------
    array : ndarray
        Float64 copy with unchanged shape and units; writeability cannot be re-
        enabled.
    """
    array = np.asarray(value, dtype=np.float64)
    return np.frombuffer(array.tobytes(), dtype=np.float64).reshape(array.shape)


def _positive(value, name):
    """Validate a finite, strictly positive real scalar.

    Parameters
    ----------
    value : float
        Scalar to validate, in the units of the named quantity.
    name : str
        Quantity name used in errors.

    Returns
    -------
    value : float
        Positive scalar in the input units.

    Raises
    ------
    ValueError
        If the value is not a finite positive real scalar.
    """
    value = scalar(value, name)
    if value <= 0:
        raise ValueError(f"{name} must be positive and representable")
    return value


def _convert(value, a_v, *, inverse, name):
    """Apply a fixed conversion between comoving and velocity coordinates.

    Parameters
    ----------
    value : array_like
        Finite real quantity, with arbitrary shape.
    a_v : float
        Positive fixed velocity conversion H(z)/((1+z)*h_fid), in
        (km/s)/(Mpc/h_fid).
    inverse : bool
        Multiply by a_v if True; divide by it otherwise.
    name : str
        Quantity name used in errors.

    Returns
    -------
    result : numpy.float64 or ndarray
        Converted float64 values; scalar input returns a NumPy scalar, and array input retains its shape.

    Raises
    ------
    ValueError
        If inputs are invalid or a nonzero result overflows or underflows to
        zero.
    """
    value = real_array(value, name)
    a_v = _positive(a_v, "a_v")
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        result = value * a_v if inverse else value / a_v
    if not np.all(np.isfinite(result)) or np.any((value != 0) & (result == 0)):
        raise ValueError(f"{name} conversion is not representable")
    return result


def wavenumber_comoving_to_velocity(k_parallel_comoving, *, a_v):
    """Convert line-of-sight wavenumber from h_fid/Mpc to s/km.

    Parameters
    ----------
    k_parallel_comoving : array_like
        Line-of-sight wavenumber in h_fid/Mpc, of arbitrary shape.
    a_v : float
        Positive fixed velocity conversion H(z)/((1+z)*h_fid), in
        (km/s)/(Mpc/h_fid).

    Returns
    -------
    result : numpy.float64 or ndarray
        Converted values in s/km, with a NumPy scalar returned for scalar input and the input shape retained otherwise.

    Raises
    ------
    ValueError
        If inputs are invalid or the converted values cannot be represented.

    Notes
    -----
    The conversion changes units only; any smoothing or width convention
    is supplied by the caller.
    """
    return _convert(k_parallel_comoving, a_v, inverse=False, name="wavenumber")


def wavenumber_velocity_to_comoving(k_parallel_velocity, *, a_v):
    """Convert line-of-sight wavenumber from s/km to h_fid/Mpc.

    Parameters
    ----------
    k_parallel_velocity : array_like
        Line-of-sight wavenumber in s/km, of arbitrary shape.
    a_v : float
        Positive fixed velocity conversion H(z)/((1+z)*h_fid), in
        (km/s)/(Mpc/h_fid).

    Returns
    -------
    result : numpy.float64 or ndarray
        Converted values in h_fid/Mpc, with a NumPy scalar returned for scalar input and the input shape retained otherwise.

    Raises
    ------
    ValueError
        If inputs are invalid or the converted values cannot be represented.

    Notes
    -----
    The conversion changes units only; any smoothing or width convention
    is supplied by the caller.
    """
    return _convert(k_parallel_velocity, a_v, inverse=True, name="wavenumber")


def p1d_velocity_to_comoving(power_velocity, *, a_v):
    """Convert one-dimensional power from km/s to Mpc/h_fid.

    Parameters
    ----------
    power_velocity : array_like
        One-dimensional power in km/s, of arbitrary shape.
    a_v : float
        Positive fixed velocity conversion H(z)/((1+z)*h_fid), in
        (km/s)/(Mpc/h_fid).

    Returns
    -------
    result : numpy.float64 or ndarray
        Converted values in Mpc/h_fid, with a NumPy scalar returned for scalar input and the input shape retained otherwise.

    Raises
    ------
    ValueError
        If inputs are invalid or the converted values cannot be represented.

    Notes
    -----
    The conversion changes units only; any smoothing or width convention
    is supplied by the caller.
    """
    return _convert(power_velocity, a_v, inverse=False, name="P1D")


def p1d_comoving_to_velocity(power_comoving, *, a_v):
    """Convert one-dimensional power from Mpc/h_fid to km/s.

    Parameters
    ----------
    power_comoving : array_like
        One-dimensional power in Mpc/h_fid, of arbitrary shape.
    a_v : float
        Positive fixed velocity conversion H(z)/((1+z)*h_fid), in
        (km/s)/(Mpc/h_fid).

    Returns
    -------
    result : numpy.float64 or ndarray
        Converted values in km/s, with a NumPy scalar returned for scalar input and the input shape retained otherwise.

    Raises
    ------
    ValueError
        If inputs are invalid or the converted values cannot be represented.

    Notes
    -----
    The conversion changes units only; any smoothing or width convention
    is supplied by the caller.
    """
    return _convert(power_comoving, a_v, inverse=True, name="P1D")


def width_velocity_to_comoving(width_velocity, *, a_v):
    """Convert nonnegative length or smoothing width from km/s to Mpc/h_fid.

    Parameters
    ----------
    width_velocity : array_like
        Nonnegative length or smoothing width in km/s, of arbitrary shape.
    a_v : float
        Positive fixed velocity conversion H(z)/((1+z)*h_fid), in
        (km/s)/(Mpc/h_fid).

    Returns
    -------
    result : numpy.float64 or ndarray
        Converted values in Mpc/h_fid, with a NumPy scalar returned for scalar input and the input shape retained otherwise.

    Raises
    ------
    ValueError
        If inputs are invalid or the converted values cannot be represented.

    Notes
    -----
    The conversion changes units only; any smoothing or width convention
    is supplied by the caller.
    """
    width = real_array(width_velocity, "width")
    if np.any(width < 0):
        raise ValueError("width must be nonnegative")
    return _convert(width, a_v, inverse=False, name="width")


def width_comoving_to_velocity(width_comoving, *, a_v):
    """Convert nonnegative length or smoothing width from Mpc/h_fid to km/s.

    Parameters
    ----------
    width_comoving : array_like
        Nonnegative length or smoothing width in Mpc/h_fid, of arbitrary shape.
    a_v : float
        Positive fixed velocity conversion H(z)/((1+z)*h_fid), in
        (km/s)/(Mpc/h_fid).

    Returns
    -------
    result : numpy.float64 or ndarray
        Converted values in km/s, with a NumPy scalar returned for scalar input and the input shape retained otherwise.

    Raises
    ------
    ValueError
        If inputs are invalid or the converted values cannot be represented.

    Notes
    -----
    The conversion changes units only; any smoothing or width convention
    is supplied by the caller.
    """
    width = real_array(width_comoving, "width")
    if np.any(width < 0):
        raise ValueError("width must be nonnegative")
    return _convert(width, a_v, inverse=True, name="width")


@dataclass(frozen=True, init=False, eq=False)
class BinGeometry:
    """Immutable preparation result; construct with prepare_geometry.

    z_nodes/w_z are interior redshifts/dz weights. hubble_nodes and
    transverse_distance_nodes are km/s/Mpc and Mpc. volume is (Mpc/h_fid)^3;
    a_v is (km/s)/(Mpc/h_fid), d_deg is (Mpc/h_fid)/degree. No background
    callable or cosmology object is retained. z_order controls volume accuracy
    independently of the fixed Fourier integration grid.
    """

    z_min: float
    z_max: float
    z_eval: float
    area_deg2: float
    solid_angle: float
    h_fid: float
    z_order: int
    z_nodes: np.ndarray
    w_z: np.ndarray
    hubble_nodes: np.ndarray
    transverse_distance_nodes: np.ndarray
    hubble_eval: float
    transverse_distance_eval: float
    volume: float
    a_v: float
    d_deg: float
    speed_light_kms: float


def _background(function, nodes, name):
    """Evaluate a background callable at immutable redshift nodes.

    Parameters
    ----------
    function : callable
        Batched background evaluator returning H in km/s/Mpc or transverse
        comoving distance in Mpc.
    nodes : ndarray of shape (n_redshift,)
        Dimensionless redshift queries.
    name : str
        Background quantity name used in errors.

    Returns
    -------
    values : ndarray of shape (n_redshift,)
        Positive finite background values, with the callable's physical units.

    Raises
    ------
    ValueError
        If evaluation fails or returns nonpositive, nonfinite, or incorrectly
        shaped data.
    """
    if not callable(function):
        raise ValueError(f"{name} must be callable")
    try:
        values = real_array(function(_immutable(nodes)), name)
        if values.shape != nodes.shape or np.any(values <= 0):
            raise ValueError("require same-shaped positive 1D output")
        return values
    except Exception as error:
        raise ValueError(f"{name} background at z={nodes.tolist()}: {error}") from error


def prepare_geometry(
    z_min, z_max, *, z_eval, area_deg2, h_fid, hubble, transverse_distance, z_order
):
    """Integrate the comoving bin volume with explicit redshift quadrature.

    Parameters
    ----------
    z_min, z_max : float
        Dimensionless redshift-bin bounds, with 0 <= z_min < z_max.
    z_eval : float
        Explicit positive evaluation redshift within the bin.
    area_deg2 : float
        Common survey area in square degrees, at most the full sky.
    h_fid : float
        Positive dimensionless reference Hubble parameter, independent of the
        physical cosmology.
    z_order : int
        Positive Gauss–Legendre order for the volume integral.
    hubble : callable
        Batched H(z) in km/s/Mpc, accepting a one-dimensional redshift array.
    transverse_distance : callable
        Batched transverse comoving D_M(z) in Mpc with the same array contract.

    Returns
    -------
    geometry : BinGeometry
        Immutable quadrature and background state, volume in (Mpc/h_fid)^3, a_v
        in (km/s)/(Mpc/h_fid), and d_deg in (Mpc/h_fid)/degree.

    Raises
    ------
    ValueError
        If bin bounds, area, background values, or represented quadrature and
        conversion factors are invalid.

    Notes
    -----
    The volume is Omega*h_fid**3*integral(c*D_M**2/H dz). Both callables
    are evaluated at quadrature nodes and z_eval during preparation only.
    No flatness, background consistency, interpolation, or convergence is inferred.
    """
    z_min, z_max, z_eval = (
        scalar(v, n)
        for v, n in ((z_min, "z_min"), (z_max, "z_max"), (z_eval, "z_eval"))
    )
    if not 0 <= z_min < z_max or not z_min <= z_eval <= z_max or z_eval <= 0:
        raise ValueError("require 0 <= z_min < z_max and positive z_eval in bin")
    area_deg2 = _positive(area_deg2, "area_deg2")
    h_fid = _positive(h_fid, "h_fid")
    full_sky_deg2 = 4 * np.pi / (np.pi / 180) ** 2
    if area_deg2 > full_sky_deg2:
        raise ValueError("area exceeds full sky")
    # Recognize the exact public full-sky endpoint before degree/radian
    # roundoff can move it one ulp outside the closed solid-angle interval.
    solid_angle = _positive(
        4 * np.pi if area_deg2 == full_sky_deg2 else area_deg2 * (np.pi / 180) ** 2,
        "solid_angle",
    )
    z_order = integer(z_order, "z_order", 1)
    legendre_nodes, legendre_weights = np.polynomial.legendre.leggauss(z_order)
    redshift_half_width = (z_max - z_min) / 2
    redshift_nodes, redshift_weights = (
        z_min + redshift_half_width * (legendre_nodes + 1),
        redshift_half_width * legendre_weights,
    )
    if (
        np.any(redshift_nodes <= z_min)
        or np.any(redshift_nodes >= z_max)
        or np.any(np.diff(redshift_nodes) <= 0)
        or not np.all(np.isfinite(redshift_weights))
        or np.any(redshift_weights <= 0)
    ):
        raise ValueError("redshift quadrature nodes/weights are not representable")
    hubble_nodes = _background(hubble, redshift_nodes, "H")
    distance_nodes = _background(transverse_distance, redshift_nodes, "D_M")
    hubble_eval = float(_background(hubble, np.array([z_eval]), "H(z_eval)")[0])
    distance_eval = float(
        _background(transverse_distance, np.array([z_eval]), "D_M(z_eval)")[0]
    )
    with np.errstate(over="ignore", under="ignore", invalid="ignore", divide="ignore"):
        integrand = SPEED_LIGHT_KMS * distance_nodes**2 / hubble_nodes
        volume = (
            solid_angle * np.float64(h_fid) ** 3 * np.sum(redshift_weights * integrand)
        )
        a_v = np.float64(hubble_eval) / (1 + z_eval) / h_fid
        d_deg = np.float64(h_fid) * distance_eval * (np.pi / 180)
    if not np.all(np.isfinite(integrand)) or np.any(integrand <= 0):
        raise ValueError("volume integrand is not positive representable float64")
    result = object.__new__(BinGeometry)
    values = dict(
        z_min=z_min,
        z_max=z_max,
        z_eval=z_eval,
        area_deg2=area_deg2,
        solid_angle=solid_angle,
        h_fid=h_fid,
        z_order=z_order,
        z_nodes=_immutable(redshift_nodes),
        w_z=_immutable(redshift_weights),
        hubble_nodes=_immutable(hubble_nodes),
        transverse_distance_nodes=_immutable(distance_nodes),
        hubble_eval=hubble_eval,
        transverse_distance_eval=distance_eval,
        volume=_positive(volume, "volume"),
        a_v=_positive(a_v, "a_v"),
        d_deg=_positive(d_deg, "d_deg"),
        speed_light_kms=SPEED_LIGHT_KMS,
    )
    for name, value in values.items():
        object.__setattr__(result, name, value)
    return result


def mode_counts(geometry, grid):
    """Multiply the fiducial bin volume by the Fourier mode-density weights.

    Parameters
    ----------
    geometry : BinGeometry
        Prepared bin geometry with volume in (Mpc/h_fid)^3.
    grid : IntegrationGrid
        Fixed Fourier quadrature with q_mode in (h_fid/Mpc)^3.

    Returns
    -------
    modes : ndarray of shape (n_node,)
        Positive dimensionless mode counts in C order with mu fastest.

    Raises
    ------
    ValueError
        If input types or h_fid values disagree, or mode counts are not positive
        finite float64 values.
    """
    if not isinstance(geometry, BinGeometry) or not isinstance(grid, IntegrationGrid):
        raise ValueError("require BinGeometry and IntegrationGrid")
    if geometry.h_fid != grid.h_fid:
        raise ValueError("geometry/grid h_fid mismatch")
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        modes = geometry.volume * grid.q_mode
    if not np.all(np.isfinite(modes)) or np.any(modes <= 0):
        raise ValueError("mode counts must be positive representable float64")
    return modes


def prepare_astropy_geometry(
    cosmology, z_min, z_max, *, z_eval, area_deg2, h_fid, z_order
):
    """Prepare bin geometry from a caller-supplied Astropy FLRW cosmology.

    Parameters
    ----------
    cosmology : astropy.cosmology.FLRW
        Physical background cosmology; cosmology.h does not replace h_fid.
    z_min, z_max : float
        Dimensionless redshift-bin bounds, with 0 <= z_min < z_max.
    z_eval : float
        Explicit positive evaluation redshift within the bin.
    area_deg2 : float
        Common survey area in square degrees, at most the full sky.
    h_fid : float
        Positive dimensionless reference Hubble parameter, independent of the
        physical cosmology.
    z_order : int
        Positive Gauss–Legendre order for the volume integral.

    Returns
    -------
    geometry : BinGeometry
        Prepared volume and background conversions; see prepare_geometry for
        units.

    Raises
    ------
    ImportError
        If the optional Astropy/SciPy cosmology dependencies are unavailable.
    ValueError
        If the cosmology type or geometry inputs are invalid.

    Notes
    -----
    Curved backgrounds use transverse comoving distance. Astropy quantities
    remain confined to preparation; the forecast receives numeric arrays.
    """
    try:
        import scipy  # noqa: F401
        from astropy import units as u
        from astropy.cosmology import FLRW
    except ImportError as error:
        raise ImportError(
            "Astropy geometry requires pip install 'fishhighz[cosmology]'"
        ) from error
    if not isinstance(cosmology, FLRW):
        raise ValueError("cosmology must be a caller-created Astropy FLRW")
    return prepare_geometry(
        z_min,
        z_max,
        z_eval=z_eval,
        area_deg2=area_deg2,
        h_fid=h_fid,
        z_order=z_order,
        hubble=lambda z: cosmology.H(z).to_value(u.km / u.s / u.Mpc),
        transverse_distance=lambda z: cosmology.comoving_transverse_distance(
            z
        ).to_value(u.Mpc),
    )
