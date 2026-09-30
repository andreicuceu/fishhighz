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
    """Return an owned read-only float64 array backed by immutable bytes."""
    array = np.asarray(values, dtype=np.float64)
    return np.frombuffer(array.tobytes(), dtype=np.float64).reshape(array.shape)


def _redshift(value, name):
    array = np.asarray(value)
    if array.ndim != 0 or array.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a finite scalar redshift")
    result = float(array)
    if not np.isfinite(result) or result < 0:
        raise ValueError(f"{name} must be a finite nonnegative redshift")
    return result


def _redshift_sequence(values, name):
    raw = np.asarray(values)
    if raw.ndim != 1 or not len(raw) or raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a nonempty one-dimensional redshift set")
    result = tuple(_redshift(value, name) for value in raw)
    if len(set(result)) != len(result):
        raise ValueError(f"{name} must contain unique redshifts")
    return result


def _single_result(value, name, redshift):
    array = np.asarray(value, dtype=np.float64)
    if array.size != 1 or not np.all(np.isfinite(array)):
        raise ValueError(f"CAMB {name} at z={redshift:g} must be one finite value")
    return float(array.reshape(-1)[0])


def _background_value(results, method, redshift, name):
    function = getattr(results, method, None)
    if function is None:
        raise ValueError(f"CAMB results do not provide {method} for {name}")
    return _single_result(function(redshift), name, redshift)


def _growth_values(results, method, redshifts):
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
    """Set CAMB's transfer request to an explicit ordered redshift set."""
    transfer = getattr(parameters, "Transfer", None)
    if transfer is None:
        raise ValueError("CAMB parameters do not expose Transfer settings")
    transfer.PK_redshifts = list(redshifts)
    transfer.PK_num_redshifts = len(redshifts)


def _returned_redshifts(results, parameters, requested):
    """Read and validate CAMB's returned transfer-redshift metadata."""
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
    xdg = os.environ.get("XDG_CACHE_HOME")
    base = Path(xdg).expanduser() if xdg else Path.home() / ".cache"
    return base / "fishhighz" / "camb"


def _camb_cache_identity(ini_bytes, camb_order, camb):
    identity = dict(
        schema=_CAMB_CACHE_SCHEMA,
        ini_sha256=hashlib.sha256(ini_bytes).hexdigest(),
        redshifts=[float(z).hex() for z in camb_order],
        camb_version=str(getattr(camb, "__version__", "unknown")),
    )
    text = json.dumps(identity, sort_keys=True)
    return hashlib.sha256(text.encode()).hexdigest(), text


def _read_growth_cache(path, identity, camb_order):
    """Return cached (sigma8, fsigma8) in CAMB order, or None if unusable."""
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
    """Atomically store growth arrays; return False when the cache is unwritable."""
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
        n = len(self.redshifts)
        arrays = (
            self.redshifts,
            self.hubble_values,
            self.transverse_distance_values,
            self.sigma8_values,
            self.growth_rate_values,
        )
        if any(np.asarray(array).shape != (n,) for array in arrays):
            raise ValueError("CAMB background arrays must share one redshift shape")
        for value in arrays:
            if np.asarray(value).flags.writeable:
                raise ValueError("CAMB background arrays must be immutable")

    @property
    def z_bins(self):
        """Compatibility alias for the explicitly ordered redshift array."""
        return self.redshifts

    @property
    def H(self):
        """H(z) in km/s/Mpc, in ``redshifts`` order."""
        return self.hubble_values

    @property
    def D_M(self):
        """Transverse comoving distance in Mpc, in redshift order."""
        return self.transverse_distance_values

    @property
    def sigma8_zbins(self):
        return self.sigma8_values

    @property
    def growth_rate_zbins(self):
        return self.growth_rate_values

    @property
    def f(self):
        return self.growth_rate_values

    @property
    def H_values(self):
        return self.hubble_values

    @property
    def D_M_values(self):
        return self.transverse_distance_values

    @property
    def f_values(self):
        return self.growth_rate_values

    @property
    def h_fid(self):
        return self.H0 / 100.0

    @property
    def sigma8(self):
        """Damping-reference sigma8 (not the template-growth normalization)."""
        return self.sigma8_damping_reference

    @property
    def results(self):
        """CAMB result surface retained for geometry quadrature callables."""
        return self._geometry_results

    @property
    def z_to_index(self):
        return self._z_to_index

    def index(self, redshift):
        redshift = _redshift(redshift, "redshift")
        try:
            return self._z_to_index[redshift]
        except KeyError as error:
            raise ValueError(
                f"CAMB was not prepared at exact redshift {redshift:g}; "
                f"prepared={list(self.redshifts)}"
            ) from error

    def sigma8_at(self, redshift):
        return float(self.sigma8_values[self.index(redshift)])

    def growth_rate_at(self, redshift):
        return float(self.growth_rate_values[self.index(redshift)])

    def growth_rate(self, redshift):
        """Evaluate the exact prepared growth rate ``f(z)``."""
        return self.growth_rate_at(redshift)

    def hubble_at(self, redshift):
        return float(self.hubble_values[self.index(redshift)])

    def transverse_distance_at(self, redshift):
        return float(self.transverse_distance_values[self.index(redshift)])

    def _call_background(self, method, redshift, name):
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
        """Evaluate H(z) in km/s/Mpc for geometry quadrature."""
        return self._call_background("hubble_parameter", redshift, "H")

    def comoving_radial_distance(self, redshift):
        """Evaluate transverse comoving distance for geometry quadrature."""
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
        sigma = float(sigma_by_z[redshift])
        fsigma8 = float(fsigma8_by_z[redshift])
        if sigma <= 0:
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
        values.append((hubble, distance, sigma, fsigma8 / sigma))

    # Retain one result surface only for the arbitrary quadrature nodes used by
    # geometry preparation.  Exact growth/background arrays above remain the
    # authoritative redshift-indexed values.
    geometry_results = results
    h0 = float(getattr(parameters, "H0"))
    if not np.isfinite(h0) or h0 <= 0:
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
        H0=h0,
        _geometry_results=geometry_results,
        _z_to_index=MappingProxyType({z: i for i, z in enumerate(ordered)}),
        cache=None if record is None else MappingProxyType(record),
    )
    return background


# Descriptive alias used by callers that prefer the preparation name.
prepare_camb_background = prepare_camb


__all__ = ["CAMBBackground", "prepare_camb", "prepare_camb_background"]
