"""Fixed single-bin geometry in physical background and fiducial h units."""

from dataclasses import dataclass

import numpy as np

from ._arrays import integer, real_array, scalar
from .grids import IntegrationGrid

SPEED_LIGHT_KMS = 299792.458
LYA_REST_ANGSTROM = 1215.67


def _immutable(value):
    # A bytes-backed view cannot have WRITEABLE re-enabled by a consumer.
    array = np.asarray(value, dtype=np.float64)
    return np.frombuffer(array.tobytes(), dtype=np.float64).reshape(array.shape)


def _positive(value, name):
    value = scalar(value, name)
    if value <= 0:
        raise ValueError(f"{name} must be positive and representable")
    return value


def _convert(value, a_v, *, inverse, name):
    value = real_array(value, name)
    a_v = _positive(a_v, "a_v")
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        result = value * a_v if inverse else value / a_v
    if not np.all(np.isfinite(result)) or np.any((value != 0) & (result == 0)):
        raise ValueError(f"{name} conversion is not representable")
    return result


def wavenumber_comoving_to_velocity(k_parallel_comoving, *, a_v):
    """Convert h_fid/Mpc to s/km using fixed a_v in (km/s)/(Mpc/h_fid)."""
    return _convert(k_parallel_comoving, a_v, inverse=False, name="wavenumber")


def wavenumber_velocity_to_comoving(k_parallel_velocity, *, a_v):
    """Convert s/km to h_fid/Mpc."""
    return _convert(k_parallel_velocity, a_v, inverse=True, name="wavenumber")


def p1d_velocity_to_comoving(power_velocity, *, a_v):
    """Convert intrinsic or explicitly smoothed km/s power to Mpc/h_fid."""
    return _convert(power_velocity, a_v, inverse=False, name="P1D")


def p1d_comoving_to_velocity(power_comoving, *, a_v):
    """Convert Mpc/h_fid power to km/s."""
    return _convert(power_comoving, a_v, inverse=True, name="P1D")


def width_velocity_to_comoving(width_velocity, *, a_v):
    """Convert km/s widths to Mpc/h_fid, without changing sigma/full-width meaning."""
    width = real_array(width_velocity, "width")
    if np.any(width < 0):
        raise ValueError("width must be nonnegative")
    return _convert(width, a_v, inverse=False, name="width")


def width_comoving_to_velocity(width_comoving, *, a_v):
    """Convert nonnegative Mpc/h_fid widths to km/s."""
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
    """Integrate Omega*h_fid^3*c*D_M^2/H with an explicit Gauss–Legendre rule.

    Supply independent ordinary batched callables H(z) in km/s/Mpc and
    transverse comoving D_M(z) in Mpc. Inputs are immutable 1D snapshots;
    outputs must have exactly matching shapes. Both callables run at quadrature
    nodes and the explicit evaluation redshift during preparation only.
    area_deg2 is one common area, at most the full sky. No flatness, background
    consistency, interpolation, fiducial cosmology, or convergence is inferred.
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
    x, w = np.polynomial.legendre.leggauss(z_order)
    half = (z_max - z_min) / 2
    nodes, weights = z_min + half * (x + 1), half * w
    if (
        np.any(nodes <= z_min)
        or np.any(nodes >= z_max)
        or np.any(np.diff(nodes) <= 0)
        or not np.all(np.isfinite(weights))
        or np.any(weights <= 0)
    ):
        raise ValueError("redshift quadrature nodes/weights are not representable")
    h_nodes = _background(hubble, nodes, "H")
    dm_nodes = _background(transverse_distance, nodes, "D_M")
    h_eval = float(_background(hubble, np.array([z_eval]), "H(z_eval)")[0])
    dm_eval = float(
        _background(transverse_distance, np.array([z_eval]), "D_M(z_eval)")[0]
    )
    with np.errstate(over="ignore", under="ignore", invalid="ignore", divide="ignore"):
        integrand = SPEED_LIGHT_KMS * dm_nodes**2 / h_nodes
        volume = solid_angle * np.float64(h_fid) ** 3 * np.sum(weights * integrand)
        a_v = np.float64(h_eval) / (1 + z_eval) / h_fid
        d_deg = np.float64(h_fid) * dm_eval * (np.pi / 180)
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
        z_nodes=_immutable(nodes),
        w_z=_immutable(weights),
        hubble_nodes=_immutable(h_nodes),
        transverse_distance_nodes=_immutable(dm_nodes),
        hubble_eval=h_eval,
        transverse_distance_eval=dm_eval,
        volume=_positive(volume, "volume"),
        a_v=_positive(a_v, "a_v"),
        d_deg=_positive(d_deg, "d_deg"),
        speed_light_kms=SPEED_LIGHT_KMS,
    )
    for name, value in values.items():
        object.__setattr__(result, name, value)
    return result


def mode_counts(geometry, grid):
    """Return V*q_mode (node,), rejecting a geometry/grid h_fid mismatch."""
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
    """Prepare from a caller-created Astropy FLRW; never use cosmology.h as h_fid.

    Requires fishhighz[cosmology]. Quantities and cosmology stay in preparation;
    curved backgrounds use comoving_transverse_distance, not radial distance.
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
