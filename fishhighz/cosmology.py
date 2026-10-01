"""Optional CAMB preparation for explicit FishHighz background quantities.

The package's numerical kernels do not need CAMB.  This module is a small
preparation boundary for callers that choose the bundled Planck18 background:
CAMB is imported only by :func:`prepare_camb`.  The complete requested set is
solved in one strictly decreasing-redshift CAMB run; its returned redshift
metadata is validated before growth arrays are mapped back to caller order.

The transfer-function solve is only needed for sigma8(z) and f(z)sigma8(z).
Those arrays are cached on disk, keyed by the CAMB ini bytes, the requested
redshift set and the CAMB version; a cache hit recomputes the background (H,
D_M), which CAMB reproduces bitwise without transfer functions.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType

import numpy as np

from .resources import bundled_path


def _immutable(values):
    """Copy values into an immutable float64 array.

    Parameters
    ----------
    values : array_like
        Numeric values of arbitrary shape.

    Returns
    -------
    array : ndarray
        Read-only float64 copy with the input shape and units.
    """
    array = np.asarray(values, dtype=np.float64)
    return np.frombuffer(array.tobytes(), dtype=np.float64).reshape(array.shape)


def _redshift(value, name):
    """Validate a scalar redshift.

    Parameters
    ----------
    value : float
        Candidate dimensionless redshift.
    name : str
        Quantity label for validation errors.

    Returns
    -------
    redshift : float
        Finite nonnegative redshift.

    Raises
    ------
    ValueError
        If the value is not a finite nonnegative real scalar.
    """
    array = np.asarray(value)
    if array.ndim != 0 or array.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a finite scalar redshift")
    result = float(array)
    if not np.isfinite(result) or result < 0:
        raise ValueError(f"{name} must be a finite nonnegative redshift")
    return result


def _redshift_sequence(values, name):
    """Validate an explicitly ordered redshift sequence.

    Parameters
    ----------
    values : array_like
        Dimensionless redshifts, shape (n_redshift,).
    name : str
        Quantity label for validation errors.

    Returns
    -------
    redshifts : tuple of float
        Unique nonnegative redshifts in input order.

    Raises
    ------
    ValueError
        If the sequence is empty, malformed, nonfinite or duplicated.
    """
    raw = np.asarray(values)
    if raw.ndim != 1 or not len(raw) or raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a nonempty one-dimensional redshift set")
    result = tuple(_redshift(value, name) for value in raw)
    if len(set(result)) != len(result):
        raise ValueError(f"{name} must contain unique redshifts")
    return result


def _single_result(value, name, redshift):
    """Extract one finite scalar from a CAMB result.

    Parameters
    ----------
    value : array_like
        CAMB output containing one value.
    name : str
        Physical quantity name for diagnostics.
    redshift : float
        Dimensionless evaluation redshift.

    Returns
    -------
    result : float
        Scalar in the units of the supplied CAMB quantity.

    Raises
    ------
    ValueError
        If the output does not contain exactly one finite value.
    """
    array = np.asarray(value, dtype=np.float64)
    if array.size != 1 or not np.all(np.isfinite(array)):
        raise ValueError(f"CAMB {name} at z={redshift:g} must be one finite value")
    return float(array.reshape(-1)[0])


def _background_value(results, method, redshift, name):
    """Evaluate one CAMB background quantity.

    Parameters
    ----------
    results : object
        Prepared CAMB results.
    method : str
        Name of the background evaluation method.
    redshift : float
        Dimensionless evaluation redshift.
    name : str
        Physical quantity name for diagnostics.

    Returns
    -------
    value : float
        Finite background value in the CAMB method units.

    Raises
    ------
    ValueError
        If the method is absent or its output is not one finite value.
    """
    function = getattr(results, method, None)
    if function is None:
        raise ValueError(f"CAMB results do not provide {method} for {name}")
    return _single_result(function(redshift), name, redshift)


def _growth_values(results, method, redshifts):
    """Read growth values in a validated CAMB redshift order.

    Parameters
    ----------
    results : object
        Prepared CAMB results.
    method : str
        Growth accessor, such as get_sigma8 or get_fsigma8.
    redshifts : sequence of float
        Previously validated dimensionless redshift order.

    Returns
    -------
    values : ndarray
        Dimensionless growth values, shape (n_redshift,).

    Raises
    ------
    ValueError
        If the accessor is absent or the output shape or values are invalid.
    """
    function = getattr(results, method, None)
    if function is None:
        raise ValueError(f"CAMB results do not provide {method}")
    value = np.asarray(function(), dtype=np.float64)
    if value.shape != (len(redshifts),) or not np.all(np.isfinite(value)):
        raise ValueError(
            f"CAMB {method} must return one finite value per validated redshift; "
            "FishHighz does not infer redshift ordering"
        )
    return value


def _set_redshifts(parameters, redshifts):
    """Set CAMB's transfer request to an explicit ordered redshift set.

    Parameters
    ----------
    parameters : object
        CAMB parameters modified in place.
    redshifts : sequence of float
        Dimensionless redshifts in the requested transfer order.

    Returns
    -------
    None
        No value is returned.

    Raises
    ------
    ValueError
        If the parameters expose no Transfer settings.
    """
    transfer = getattr(parameters, "Transfer", None)
    if transfer is None:
        raise ValueError("CAMB parameters do not expose Transfer settings")
    transfer.PK_redshifts = list(redshifts)
    transfer.PK_num_redshifts = len(redshifts)


def _returned_redshifts(results, parameters, requested):
    """Read and validate CAMB's returned transfer-redshift metadata.

    Parameters
    ----------
    results : object
        CAMB results whose parameter metadata are checked first.
    parameters : object
        Requested CAMB parameters used as the final metadata source.
    requested : sequence of float
        Exact decreasing dimensionless redshift order.

    Returns
    -------
    redshifts : tuple of float
        Validated order, identical to the request.

    Raises
    ------
    ValueError
        If metadata are absent, malformed or inconsistent with the request.
    """
    surfaces = [
        getattr(results, "Params", None),
        getattr(results, "params", None),
        parameters,
    ]
    for surface in surfaces:
        transfer = getattr(surface, "Transfer", None)
        values = None if transfer is None else getattr(transfer, "PK_redshifts", None)
        if values is None:
            continue
        values = np.asarray(values, dtype=np.float64)
        if values.shape != (len(requested),) or not np.all(np.isfinite(values)):
            raise ValueError("CAMB returned malformed transfer-redshift metadata")
        if not np.array_equal(values, np.asarray(requested, dtype=np.float64)):
            raise ValueError(
                "CAMB returned transfer redshifts that differ from the requested "
                f"decreasing order: returned={values.tolist()}, requested={list(requested)}"
            )
        return tuple(float(value) for value in values)
    raise ValueError(
        "CAMB results expose no transfer-redshift metadata; refusing to infer "
        "sigma8/fsigma8 ordering"
    )


_CAMB_CACHE_SCHEMA = 1
_CACHE_DISABLED = ("0", "false", "no", "off")


def _camb_cache_dir(cache_dir, camb_module):
    """Resolve the growth cache directory, or None when caching is disabled.

    Parameters
    ----------
    cache_dir : path-like or None
        Explicit cache directory, or None for environment defaults.
    camb_module : module-like or None
        Injected CAMB module; disables implicit caching when supplied.

    Returns
    -------
    directory : pathlib.Path or None
        Resolved directory, or None when caching is disabled.

    Notes
    -----
    ``FISHHIGHZ_CAMB_CACHE=0`` disables caching. An explicit ``cache_dir`` wins;
    otherwise ``$FISHHIGHZ_CACHE_DIR/camb``, ``$XDG_CACHE_HOME/fishhighz/camb``
    or ``~/.cache/fishhighz/camb``. Injected CAMB test surfaces are only
    cached when a directory is given explicitly.
    """
    setting = os.environ.get("FISHHIGHZ_CAMB_CACHE", "1").strip().lower()
    if setting in _CACHE_DISABLED:
        return None
    if cache_dir is not None:
        return Path(cache_dir).expanduser()
    if camb_module is not None:
        return None
    root = os.environ.get("FISHHIGHZ_CACHE_DIR")
    if root:
        return Path(root).expanduser() / "camb"
    xdg_cache_home = os.environ.get("XDG_CACHE_HOME")
    base = (
        Path(xdg_cache_home).expanduser() if xdg_cache_home else Path.home() / ".cache"
    )
    return base / "fishhighz" / "camb"


def _camb_cache_identity(ini_bytes, camb_order, camb):
    """Hash the CAMB inputs that determine the cached growth arrays.

    Parameters
    ----------
    ini_bytes : bytes
        Exact CAMB configuration contents.
    camb_order : sequence of float
        Dimensionless redshifts in CAMB order.
    camb : module-like
        CAMB module supplying its version string.

    Returns
    -------
    key : str
        SHA256 digest used as the cache filename.
    identity : str
        JSON representation of the cache identity.
    """
    identity = dict(
        schema=_CAMB_CACHE_SCHEMA,
        ini_sha256=hashlib.sha256(ini_bytes).hexdigest(),
        redshifts=[float(z).hex() for z in camb_order],
        camb_version=str(getattr(camb, "__version__", "unknown")),
    )
    text = json.dumps(identity, sort_keys=True)
    return hashlib.sha256(text.encode()).hexdigest(), text


def _read_growth_cache(path, identity, camb_order):
    """Read cached growth arrays only when their identity and order match.

    Parameters
    ----------
    path : path-like
        Existing or candidate NPZ cache file.
    identity : str
        Expected serialized cache identity.
    camb_order : sequence of float
        Expected dimensionless redshifts in CAMB order.

    Returns
    -------
    growth : tuple of ndarray or None
        Dimensionless sigma8 and f*sigma8 arrays, each shape (n_redshift,), or
        None for a missing or unusable cache.
    """
    try:
        with np.load(path, allow_pickle=False) as stored:
            if str(stored["identity"]) != identity:
                return None
            redshifts = np.asarray(stored["redshifts"], dtype=np.float64)
            sigma8 = np.array(stored["sigma8"], dtype=np.float64)
            fsigma8 = np.array(stored["fsigma8"], dtype=np.float64)
    except (OSError, KeyError, ValueError):
        return None
    shape = (len(camb_order),)
    if (
        not np.array_equal(redshifts, np.asarray(camb_order, dtype=np.float64))
        or sigma8.shape != shape
        or fsigma8.shape != shape
        or not np.all(np.isfinite(sigma8))
        or not np.all(np.isfinite(fsigma8))
    ):
        return None
    return sigma8, fsigma8


def _write_growth_cache(path, identity, camb_order, sigma8, fsigma8):
    """Atomically store growth arrays when the cache is writable.

    Parameters
    ----------
    path : pathlib.Path
        Destination NPZ cache file; parent directories are created.
    identity : str
        Serialized cache identity.
    camb_order : sequence of float
        Dimensionless redshifts in CAMB order.
    sigma8 : array_like
        Dimensionless sigma8 values, shape (n_redshift,).
    fsigma8 : array_like
        Dimensionless f*sigma8 values, shape (n_redshift,).

    Returns
    -------
    stored : bool
        True after a successful atomic replacement; False for an I/O failure.
    """
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            dir=path.parent, prefix=".tmp-", suffix=".npz", delete=False
        ) as handle:
            np.savez(
                handle,
                identity=np.asarray(identity),
                redshifts=np.asarray(camb_order, dtype=np.float64),
                sigma8=np.asarray(sigma8, dtype=np.float64),
                fsigma8=np.asarray(fsigma8, dtype=np.float64),
            )
        os.replace(handle.name, path)
    except OSError:
        try:
            os.unlink(handle.name)
        except (OSError, NameError):
            pass
        return False
    return True


def _load_camb(camb_module):
    """Load the optional CAMB dependency only when needed.

    Parameters
    ----------
    camb_module : module-like or None
        Explicit module to retain, or None to import CAMB.

    Returns
    -------
    camb : module-like
        Injected or imported CAMB module.

    Raises
    ------
    ImportError
        If CAMB is unavailable and no module was supplied.
    """
    if camb_module is not None:
        return camb_module
    try:
        import camb
    except ImportError as error:
        raise ImportError(
            "CAMB preparation requires the optional CAMB dependency; install "
            "'fishhighz[camb]'"
        ) from error
    return camb


@dataclass(frozen=True, eq=False)
class CAMBBackground:
    """Immutable CAMB values at the requested exact redshifts.

    ``redshifts`` and all value arrays have identical order.  ``z_to_index``
    records that mapping explicitly, and the ``*_at`` methods reject redshifts
    that were not prepared.  The background callables used by volume
    quadrature are retained separately because quadrature nodes are not part of
    the exact growth-evaluation set.
    """

    redshifts: np.ndarray
    hubble_values: np.ndarray
    transverse_distance_values: np.ndarray
    sigma8_values: np.ndarray
    growth_rate_values: np.ndarray
    template_growth_redshift: float
    damping_reference_redshift: float
    sigma8_template: float
    sigma8_damping_reference: float
    H0: float
    _geometry_results: object = field(repr=False, compare=False)
    _z_to_index: MappingProxyType = field(repr=False, compare=False)
    cache: MappingProxyType | None = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        """Validate common redshift shapes and immutable background arrays.

        Returns
        -------
        None
            No value is returned.

        Raises
        ------
        ValueError
            If arrays differ in shape or allow writes.
        """
        n_redshifts = len(self.redshifts)
        arrays = (
            self.redshifts,
            self.hubble_values,
            self.transverse_distance_values,
            self.sigma8_values,
            self.growth_rate_values,
        )
        if any(np.asarray(array).shape != (n_redshifts,) for array in arrays):
            raise ValueError("CAMB background arrays must share one redshift shape")
        for value in arrays:
            if np.asarray(value).flags.writeable:
                raise ValueError("CAMB background arrays must be immutable")

    @property
    def z_bins(self):
        """Return the prepared redshift order.

        Returns
        -------
        value : ndarray
            Dimensionless redshifts, shape (n_redshift,).
        """
        return self.redshifts

    @property
    def H(self):
        """Return the prepared Hubble parameter.

        Returns
        -------
        value : ndarray
            H(z) in km/s/Mpc, shape (n_redshift,).
        """
        return self.hubble_values

    @property
    def D_M(self):
        """Return the prepared transverse comoving distance.

        Returns
        -------
        value : ndarray
            D_M(z) in Mpc, shape (n_redshift,).
        """
        return self.transverse_distance_values

    @property
    def sigma8_zbins(self):
        """Return the prepared density fluctuation amplitude.

        Returns
        -------
        value : ndarray
            Dimensionless sigma8(z), shape (n_redshift,).
        """
        return self.sigma8_values

    @property
    def growth_rate_zbins(self):
        """Return the prepared logarithmic growth rate.

        Returns
        -------
        value : ndarray
            Dimensionless f(z), shape (n_redshift,).
        """
        return self.growth_rate_values

    @property
    def f(self):
        """Return the prepared logarithmic growth rate.

        Returns
        -------
        value : ndarray
            Dimensionless f(z), shape (n_redshift,).
        """
        return self.growth_rate_values

    @property
    def H_values(self):
        """Return the prepared Hubble parameter.

        Returns
        -------
        value : ndarray
            H(z) in km/s/Mpc, shape (n_redshift,).
        """
        return self.hubble_values

    @property
    def D_M_values(self):
        """Return the prepared transverse comoving distance.

        Returns
        -------
        value : ndarray
            D_M(z) in Mpc, shape (n_redshift,).
        """
        return self.transverse_distance_values

    @property
    def f_values(self):
        """Return the prepared logarithmic growth rate.

        Returns
        -------
        value : ndarray
            Dimensionless f(z), shape (n_redshift,).
        """
        return self.growth_rate_values

    @property
    def h_fid(self):
        """Return the fiducial reduced Hubble constant.

        Returns
        -------
        value : float
            Dimensionless H0/(100 km/s/Mpc).
        """
        return self.H0 / 100.0

    @property
    def sigma8(self):
        """Return sigma8 at the damping reference redshift.

        Returns
        -------
        value : float
            Dimensionless damping normalization, independent of template
            normalization.
        """
        return self.sigma8_damping_reference

    @property
    def results(self):
        """Return the CAMB results retained for geometry quadrature.

        Returns
        -------
        value : object
            Prepared CAMB background result object.
        """
        return self._geometry_results

    @property
    def z_to_index(self):
        """Return the exact prepared-redshift index mapping.

        Returns
        -------
        value : mappingproxy
            Mapping from dimensionless redshift to array index.
        """
        return self._z_to_index

    def index(self, redshift):
        """Find an exact prepared redshift.

        Parameters
        ----------
        redshift : float
            Dimensionless redshift present in the prepared set.

        Returns
        -------
        value : int
            Index in all prepared background arrays.

        Raises
        ------
        ValueError
            If redshift is invalid or was not prepared exactly.
        """
        redshift = _redshift(redshift, "redshift")
        try:
            return self._z_to_index[redshift]
        except KeyError as error:
            raise ValueError(
                f"CAMB was not prepared at exact redshift {redshift:g}; "
                f"prepared={list(self.redshifts)}"
            ) from error

    def sigma8_at(self, redshift):
        """Read sigma8 at an exact prepared redshift.

        Parameters
        ----------
        redshift : float
            Dimensionless redshift present in the prepared set.

        Returns
        -------
        value : float
            Dimensionless sigma8.

        Raises
        ------
        ValueError
            If redshift is invalid or was not prepared exactly.
        """
        return float(self.sigma8_values[self.index(redshift)])

    def growth_rate_at(self, redshift):
        """Read the growth rate at an exact prepared redshift.

        Parameters
        ----------
        redshift : float
            Dimensionless redshift present in the prepared set.

        Returns
        -------
        value : float
            Dimensionless logarithmic growth rate f.

        Raises
        ------
        ValueError
            If redshift is invalid or was not prepared exactly.
        """
        return float(self.growth_rate_values[self.index(redshift)])

    def growth_rate(self, redshift):
        """Read the growth rate at an exact prepared redshift.

        Parameters
        ----------
        redshift : float
            Dimensionless redshift present in the prepared set.

        Returns
        -------
        value : float
            Dimensionless logarithmic growth rate f.

        Raises
        ------
        ValueError
            If redshift is invalid or was not prepared exactly.
        """
        return self.growth_rate_at(redshift)

    def hubble_at(self, redshift):
        """Read H(z) at an exact prepared redshift.

        Parameters
        ----------
        redshift : float
            Dimensionless redshift present in the prepared set.

        Returns
        -------
        value : float
            Hubble parameter in km/s/Mpc.

        Raises
        ------
        ValueError
            If redshift is invalid or was not prepared exactly.
        """
        return float(self.hubble_values[self.index(redshift)])

    def transverse_distance_at(self, redshift):
        """Read D_M(z) at an exact prepared redshift.

        Parameters
        ----------
        redshift : float
            Dimensionless redshift present in the prepared set.

        Returns
        -------
        value : float
            Transverse comoving distance in Mpc.

        Raises
        ------
        ValueError
            If redshift is invalid or was not prepared exactly.
        """
        return float(self.transverse_distance_values[self.index(redshift)])

    def _call_background(self, method, redshift, name):
        """Evaluate CAMB geometry at scalar or vector redshifts.

        Parameters
        ----------
        method : str
            Name of the CAMB background method.
        redshift : float or array_like
            Dimensionless scalar or shape (n_redshift,) evaluation coordinates.
        name : str
            Quantity label for diagnostics.

        Returns
        -------
        values : float or ndarray
            Result in the CAMB method units, preserving scalar or vector shape.

        Raises
        ------
        ValueError
            If redshifts, method availability, output shape or finiteness are
            invalid.

        Notes
        -----
        A scalar loop supports result objects that do not accept vector queries.
        """
        values = np.asarray(redshift)
        if values.ndim == 0:
            values = np.asarray([_redshift(values, "redshift")])
            scalar = True
        elif values.ndim == 1 and values.dtype.kind in "iuf":
            values = np.asarray([_redshift(value, "redshift") for value in values])
            scalar = False
        else:
            raise ValueError("redshift must be a scalar or one-dimensional array")
        function = getattr(self._geometry_results, method, None)
        if function is None:
            raise ValueError(f"CAMB results do not provide {method} for {name}")
        # CAMB accepts arrays in current releases, while a scalar loop also
        # supports small fake result surfaces used by deterministic tests.
        try:
            output = np.asarray(function(values), dtype=np.float64)
            if output.shape != values.shape:
                raise ValueError
        except (TypeError, ValueError):
            output = np.asarray(
                [function(float(value)) for value in values], dtype=np.float64
            )
        if output.shape != values.shape or not np.all(np.isfinite(output)):
            raise ValueError(f"CAMB {name} returned a nonfinite or wrong-shaped value")
        return float(output[0]) if scalar else output

    def hubble_parameter(self, redshift):
        """Evaluate the background for geometry quadrature.

        Parameters
        ----------
        redshift : float or array_like
            Dimensionless scalar or shape (n_redshift,) coordinates.

        Returns
        -------
        values : float or ndarray
            H(z) in km/s/Mpc, with the input shape.

        Raises
        ------
        ValueError
            If redshifts or CAMB background results are invalid.
        """
        return self._call_background("hubble_parameter", redshift, "H")

    def comoving_radial_distance(self, redshift):
        """Evaluate the background for geometry quadrature.

        Parameters
        ----------
        redshift : float or array_like
            Dimensionless scalar or shape (n_redshift,) coordinates.

        Returns
        -------
        values : float or ndarray
            Transverse comoving distance D_M(z) in Mpc, with the input shape.

        Raises
        ------
        ValueError
            If redshifts or CAMB background results are invalid.
        """
        values = np.asarray(redshift)
        if hasattr(self._geometry_results, "angular_diameter_distance"):
            angular = self._call_background(
                "angular_diameter_distance", redshift, "angular diameter distance"
            )
            return angular * (1 + values)
        if hasattr(self._geometry_results, "comoving_radial_distance"):
            return self._call_background(
                "comoving_radial_distance", redshift, "transverse distance"
            )
        raise ValueError("CAMB results provide neither transverse-distance method")

    def transverse_comoving_distance(self, redshift):
        """Evaluate the background for geometry quadrature.

        Parameters
        ----------
        redshift : float or array_like
            Dimensionless scalar or shape (n_redshift,) coordinates.

        Returns
        -------
        values : float or ndarray
            Transverse comoving distance D_M(z) in Mpc, with the input shape.

        Raises
        ------
        ValueError
            If redshifts or CAMB background results are invalid.
        """
        return self.comoving_radial_distance(redshift)


def prepare_camb(
    ini: str | Path | None = None,
    redshifts=None,
    *,
    template_growth_redshift=None,
    template_redshift=None,
    damping_reference_redshift=2.3,
    camb_module=None,
    cache_dir=None,
):
    """Prepare the native CAMB background at an exact, explicit redshift set.

    Parameters
    ----------
    ini : path-like, optional
        CAMB parameter file.  ``None`` selects the byte-preserved bundled
        ``camb_configs/Planck18.ini`` resource.
    redshifts : sequence of float
        Ordered evaluation redshifts.  The two named normalization redshifts
        are appended when absent, preserving this input order.
    template_growth_redshift, template_redshift : float
        Redshift used for the template power-growth ratio.  The alias is kept
        for callers that use the shorter historical name; supplying both names
        is rejected rather than silently choosing one.
    damping_reference_redshift : float
        Independent redshift used in the damping-width sigma8 normalization.
    camb_module : module-like, optional
        Test surface or already imported CAMB module.  Normal callers should
        leave this unset; importing CAMB remains lazy.
    cache_dir : path-like, optional
        Directory for the sigma8/fsigma8 cache.  By default the user cache
        directory is used for the real CAMB module (see ``_camb_cache_dir``);
        ``FISHHIGHZ_CAMB_CACHE=0`` disables caching.  A hit skips only the
        transfer-function solve: the ini is read, the redshifts are set and the
        background is recomputed with ``camb.get_background``.

    Returns
    -------
    background : CAMBBackground
        Immutable background and growth arrays in caller redshift order, with
        missing normalization redshifts appended. H is in km/s/Mpc, distances
        are in Mpc, and sigma8 and growth rate are dimensionless.

    Raises
    ------
    ValueError
        If redshifts, normalization choices or CAMB results are invalid.
    ImportError
        If CAMB is unavailable and no module was injected.

    Notes
    -----
    The default damping reference redshift is 2.3. Both template-redshift
    arguments default to None; exactly one must be specified. Growth caching
    may create or replace files in the resolved cache directory.
    """
    if template_growth_redshift is not None and template_redshift is not None:
        raise ValueError(
            "provide either template_growth_redshift or template_redshift, not both"
        )
    if template_growth_redshift is None:
        template_growth_redshift = template_redshift
    if template_growth_redshift is None:
        raise ValueError("template_growth_redshift must be explicit")
    template_growth_redshift = _redshift(
        template_growth_redshift, "template_growth_redshift"
    )
    damping_reference_redshift = _redshift(
        damping_reference_redshift, "damping_reference_redshift"
    )
    if redshifts is None:
        raise ValueError("redshifts must be an explicit nonempty sequence")
    caller_order = list(_redshift_sequence(redshifts, "redshifts"))
    ordered = list(caller_order)
    for redshift in (template_growth_redshift, damping_reference_redshift):
        if redshift not in ordered:
            ordered.append(redshift)
    ordered = tuple(ordered)
    camb = _load_camb(camb_module)

    # CAMB 2.x returns growth arrays in increasing-time (decreasing-z) order.
    # Request that order explicitly, validate the returned metadata, and only
    # then map values back to the caller's original order.
    camb_order = tuple(sorted(ordered, reverse=True))

    directory = _camb_cache_dir(cache_dir, camb_module)

    def run_bulk(path):
        """Prepare CAMB geometry and retrieve growth arrays for one INI.

        Parameters
        ----------
        path : path-like
            CAMB parameter file to read.

        Returns
        -------
        parameters : object
            Parsed CAMB parameters with explicit transfer redshifts.
        results : object
            CAMB results retaining background geometry.
        growth : tuple of ndarray
            Dimensionless sigma8 and f*sigma8 arrays in CAMB order.
        record : dict or None
            Cache identity and hit/write status, or None when disabled.

        Notes
        -----
        The enclosing preparation supplies the redshift order and cache directory; a cache miss may write growth arrays.
        """
        parameters = camb.read_ini(str(path))
        _set_redshifts(parameters, camb_order)
        record = None
        if directory is not None:
            key, identity = _camb_cache_identity(
                Path(path).read_bytes(), camb_order, camb
            )
            cache_path = directory / f"{key}.npz"
            record = dict(key=key, path=str(cache_path))
            cached = _read_growth_cache(cache_path, identity, camb_order)
            if cached is not None and hasattr(camb, "get_background"):
                record["status"] = "hit"
                return parameters, camb.get_background(parameters), cached, record
        results = camb.get_results(parameters)
        returned = _returned_redshifts(results, parameters, camb_order)
        growth = (
            _growth_values(results, "get_sigma8", returned),
            _growth_values(results, "get_fsigma8", returned),
        )
        if record is not None:
            written = _write_growth_cache(cache_path, identity, returned, *growth)
            record["status"] = "miss, stored" if written else "miss, unwritable"
        return parameters, results, growth, record

    if ini is None:
        with bundled_path("camb_configs/Planck18.ini") as path:
            parameters, results, growth, record = run_bulk(path)
    else:
        path = Path(ini).expanduser().resolve(strict=True)
        parameters, results, growth, record = run_bulk(path)

    # A cache entry is only accepted for exactly this decreasing request order.
    returned = camb_order
    sigma_returned, fsigma8_returned = growth
    returned_index = {redshift: index for index, redshift in enumerate(returned)}
    sigma_by_z = {
        redshift: sigma_returned[index] for redshift, index in returned_index.items()
    }
    fsigma8_by_z = {
        redshift: fsigma8_returned[index] for redshift, index in returned_index.items()
    }
    values = []
    for redshift in ordered:
        sigma8_value = float(sigma_by_z[redshift])
        fsigma8 = float(fsigma8_by_z[redshift])
        if sigma8_value <= 0:
            raise ValueError(f"CAMB sigma8 at z={redshift:g} must be positive")
        hubble = _background_value(results, "hubble_parameter", redshift, "H")
        if hasattr(results, "angular_diameter_distance"):
            distance = _background_value(
                results,
                "angular_diameter_distance",
                redshift,
                "angular diameter distance",
            )
            distance *= 1 + redshift
        elif hasattr(results, "comoving_radial_distance"):
            distance = _background_value(
                results, "comoving_radial_distance", redshift, "transverse distance"
            )
        else:
            raise ValueError(
                "CAMB results provide neither angular nor comoving distance"
            )
        if hubble <= 0:
            raise ValueError(f"CAMB H at z={redshift:g} must be positive")
        # Normalization references may be at the observer, where D_M(0)=0.
        # Survey geometry retains its own positive-distance requirements.
        if (
            not np.isfinite(distance)
            or distance < 0
            or (redshift > 0 and distance == 0)
        ):
            raise ValueError(
                f"CAMB transverse distance at z={redshift:g} must be finite and "
                "nonnegative, and positive at positive redshift"
            )
        values.append((hubble, distance, sigma8_value, fsigma8 / sigma8_value))

    # Retain one result surface only for the arbitrary quadrature nodes used by
    # geometry preparation.  Exact growth/background arrays above remain the
    # authoritative redshift-indexed values.
    geometry_results = results
    hubble_constant = float(getattr(parameters, "H0"))
    if not np.isfinite(hubble_constant) or hubble_constant <= 0:
        raise ValueError("CAMB H0 must be positive and finite")
    arrays = np.asarray(values, dtype=np.float64).T
    redshift_array = _immutable(ordered)
    background = CAMBBackground(
        redshifts=redshift_array,
        hubble_values=_immutable(arrays[0]),
        transverse_distance_values=_immutable(arrays[1]),
        sigma8_values=_immutable(arrays[2]),
        growth_rate_values=_immutable(arrays[3]),
        template_growth_redshift=template_growth_redshift,
        damping_reference_redshift=damping_reference_redshift,
        sigma8_template=float(arrays[2][ordered.index(template_growth_redshift)]),
        sigma8_damping_reference=float(
            arrays[2][ordered.index(damping_reference_redshift)]
        ),
        H0=hubble_constant,
        _geometry_results=geometry_results,
        _z_to_index=MappingProxyType({z: i for i, z in enumerate(ordered)}),
        cache=None if record is None else MappingProxyType(record),
    )
    return background


# Descriptive alias used by callers that prefer the preparation name.
prepare_camb_background = prepare_camb


__all__ = ["CAMBBackground", "prepare_camb", "prepare_camb_background"]
