"""Strict legacy density/SNR readers, adapted from lyaforecast (GPLv3).

Provenance: lyaforecast/tracer.py (_setup_dndzdm_lya/_tracer) and
spectrograph.py (_setup_desi_spectro, get_pixel_rms_noise). Only interpolation
and explicit normalization conventions are retained; boundary fallbacks are
not. No raw assets are distributed. SciPy is imported only on reader setup.
"""

import hashlib
import io
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .._arrays import label as validate_label
from .._arrays import real_array, scalar
from ..geometry import LYA_REST_ANGSTROM, _immutable, _positive
from ..noise import local_galaxy_density
from ..response import pixel_width_angstrom_to_velocity
from ..survey import freeze
from ..weights import density_per_velocity


def _scipy():
    try:
        from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator
        from scipy.ndimage import gaussian_filter1d
    except ImportError as error:
        raise ImportError(
            "raw survey readers require SciPy; install fishhighz[survey]"
        ) from error
    return RectBivariateSpline, RegularGridInterpolator, gaussian_filter1d


def _read(path):
    path = Path(path).resolve(strict=True)
    content = path.read_bytes()
    try:
        data = real_array(np.loadtxt(io.BytesIO(content), ndmin=2), str(path))
    except Exception as error:
        raise ValueError(f"{path}: malformed numeric table: {error}") from error
    return str(path), hashlib.sha256(content).hexdigest(), content.decode("utf-8"), data


def _axis(axis, name, minimum=1):
    axis = real_array(axis, name)
    if axis.ndim != 1 or len(axis) < minimum or np.any(np.diff(axis) <= 0):
        raise ValueError(f"{name}: require ordered unique axis with >= {minimum} nodes")
    return _immutable(axis)


def _query(axis, values, name, context):
    values = real_array(values, name)
    if np.any((values < axis[0]) | (values > axis[-1])):
        raise ValueError(
            f"{context}: {name}={values.tolist()} outside closed domain [{axis[0]}, {axis[-1]}]"
        )
    return values


@dataclass(frozen=True, init=False, eq=False)
class DensityReader:
    """Quadratic density adapter for explicit cell counts per square degree.

    Parameters
    ----------
    path : path-like
        Three-column (z, magnitude, cell_count_per_deg2) rectangular table.
    semantics : str
        Must be 'cell_count_per_deg2'; normalized arrays bypass this reader.
    target_density : float or None
        Explicit positive target or None for no renormalization.
    z_norm_min : float or None
        Strict z>threshold selection, or None for the whole masked raw grid.
    magnitude_bounds : pair, optional
        Inclusive raw mask; either endpoint may be None.
    redshift_widths : array_like, optional
        Positive physical cell widths in sorted redshift-axis order, shape (n_z,).
        These are explicit measures, never inferred from centre spacing or bins.
    width_policy : {'uniform', 'legacy_first_spacing'}, optional
        Omitted with widths selects 'explicit'. Omitted without widths requires
        uniform z. Explicit legacy_first_spacing uses z[1]-z[0] for all rows,
        even on irregular z; this is a compatibility convention, not physical
        width inference. Cannot be combined with redshift_widths.
    label : str
        Caller population label, used only for provenance and errors.

    Notes
    -----
    Uniform-axis tolerance is 64*eps64*max(1,max(abs(axis))), an absolute
    float64 coordinate-roundoff allowance. Magnitudes must remain uniform;
    irregular redshifts require explicit widths or legacy_first_spacing.
    """

    z: np.ndarray
    magnitudes: np.ndarray
    raw_counts: np.ndarray
    redshift_widths: np.ndarray
    density: np.ndarray
    provenance: object
    _spline: object
    _context: str

    def __init__(
        self,
        path,
        *,
        semantics,
        target_density,
        z_norm_min,
        magnitude_bounds=None,
        redshift_widths=None,
        width_policy=None,
        label="density",
    ):
        spline, _, _ = _scipy()
        if semantics != "cell_count_per_deg2":
            raise ValueError("require explicit cell_count_per_deg2 semantics")
        caller = validate_label(label, "density label")
        path, digest, _, data = _read(path)
        context = f"{caller} ({path})"
        if data.shape[1] != 3 or np.any(data[:, 2] < 0):
            raise ValueError(
                f"{context}: require three columns and nonnegative cell counts"
            )
        z, iz = np.unique(data[:, 0], return_inverse=True)
        m, im = np.unique(data[:, 1], return_inverse=True)
        z, m = _axis(z, "redshift", 3), _axis(m, "magnitude", 3)
        if z[0] < 0:
            raise ValueError(f"{context}: negative redshift")
        flat = iz * len(m) + im
        if len(np.unique(flat)) != len(flat) or len(flat) != len(z) * len(m):
            raise ValueError(f"{context}: duplicate or missing grid cells")
        dz, dm = np.diff(z), np.diff(m)
        tolerance = 64 * np.finfo(float).eps * max(1, np.max(np.abs(m)))
        if np.any(np.abs(dm - dm[0]) > tolerance):
            raise ValueError(f"{context}: nonuniform magnitude grid spacing")
        dm = float(dm[0])
        if redshift_widths is not None:
            if width_policy is not None:
                raise ValueError(
                    f"{context}: explicit widths conflict with width_policy"
                )
            widths = real_array(redshift_widths, "redshift_widths")
            if widths.shape != z.shape or np.any(widths <= 0):
                raise ValueError(
                    f"{context}: redshift_widths must be positive (n_z,) in sorted z order"
                )
            policy = "explicit"
        else:
            policy = "uniform" if width_policy is None else width_policy
            if policy not in ("uniform", "legacy_first_spacing"):
                raise ValueError(f"{context}: unknown width_policy {policy!r}")
            tolerance = 64 * np.finfo(float).eps * max(1, np.max(np.abs(z)))
            if policy == "uniform" and np.any(np.abs(dz - dz[0]) > tolerance):
                raise ValueError(
                    f"{context}: nonuniform redshift nodes require redshift_widths or explicit width_policy='legacy_first_spacing'"
                )
            widths = np.full(len(z), dz[0])
        widths = _immutable(widths)
        raw = np.empty((len(z), len(m)))
        raw[iz, im] = data[:, 2]
        mask = np.ones(len(m), dtype=bool)
        bounds = None
        if magnitude_bounds is not None:
            if len(magnitude_bounds) != 2:
                raise ValueError("magnitude_bounds must have two endpoints")
            lo, hi = (
                None if x is None else scalar(x, "magnitude bound")
                for x in magnitude_bounds
            )
            if lo is not None and hi is not None and lo > hi:
                raise ValueError("reversed magnitude bounds")
            mask = (m >= (-np.inf if lo is None else lo)) & (
                m <= (np.inf if hi is None else hi)
            )
            bounds = (lo, hi)
        threshold = None if z_norm_min is None else scalar(z_norm_min, "z_norm_min")
        selected = np.ones(len(z), dtype=bool) if threshold is None else z > threshold
        target = (
            None
            if target_density is None
            else _positive(target_density, "target_density")
        )
        masked = raw * mask
        try:
            with np.errstate(
                over="raise", under="raise", invalid="raise", divide="raise"
            ):
                original = float(raw.sum())
                measure = float(masked[selected].sum())
                if target is not None and measure <= 0:
                    raise ValueError(f"{context}: zero target-normalization support")
                scale = 1.0 if target is None else np.float64(target) / measure
                density = masked * scale / (widths[:, None] * dm)
        except FloatingPointError as error:
            raise ValueError(
                f"{context}: density normalization not representable"
            ) from error
        if not np.all(np.isfinite(density)):
            raise ValueError(f"{context}: nonfinite prepared density")
        provenance = freeze(
            dict(
                path=path,
                sha256=digest,
                label=caller,
                semantics=semantics,
                original_measure=original,
                selected_measure=measure,
                magnitude_bounds=bounds,
                z_norm_min=threshold,
                target_density=target,
                scale=float(scale),
                dz=None if policy == "explicit" else float(widths[0]),
                dm=dm,
                width_policy=policy,
                redshift_widths=widths,
                redshift_axis=z,
                magnitude_axis=m,
                width_order="reconstructed sorted redshift axis",
                units="deg^-2 redshift^-1 mag^-1",
                interpolation="RectBivariateSpline kx=2 ky=2 s=0",
                z_domain=(z[0], z[-1]),
                magnitude_domain=(m[0], m[-1]),
            )
        )
        for name, value in dict(
            z=z,
            magnitudes=m,
            raw_counts=_immutable(raw),
            redshift_widths=widths,
            density=_immutable(density),
            provenance=provenance,
            _context=context,
            _spline=spline(
                z, m, density, kx=2, ky=2, s=0, bbox=[z[0], z[-1], m[0], m[-1]]
            ),
        ).items():
            object.__setattr__(self, name, value)

    def query(self, z, magnitudes):
        """Return owned dndzdm (magnitude,) at one redshift, preserving order."""
        z = scalar(z, "redshift")
        _query(self.z, z, "redshift", self._context)
        m = _query(self.magnitudes, magnitudes, "magnitude", self._context)
        if m.ndim != 1 or not len(m):
            raise ValueError("magnitudes must be nonempty 1D")
        result = self._spline.ev(np.full(m.shape, z), m)
        if np.any(result < 0) or not np.all(np.isfinite(result)):
            raise ValueError(
                f"{self._context}: negative/nonfinite interpolated density at z={z}, magnitudes={m.tolist()}"
            )
        return _immutable(result)

    def local_galaxy_density(self, geometry, magnitudes, quadrature):
        """Explicit local z_eval approximation, no area or volume averaging."""
        return local_galaxy_density(
            self.query(geometry.z_eval, magnitudes), quadrature, geometry
        )


def _header(text, path):
    metadata, redshifts = {}, None
    for line in text.splitlines():
        if not line.lstrip().startswith("#"):
            continue
        line = line.lstrip()[1:].strip()
        for key, value in re.findall(r"\b(BAND|MAG|EXPTIME|NEXP)\s*=\s*([^\s]+)", line):
            if key in metadata:
                raise ValueError(f"{path}: duplicate header {key}")
            metadata[key] = value
        if line.startswith("Wave"):
            tokens = line.split()
            if tokens[0] != "Wave" or len(tokens) < 2 or redshifts is not None:
                raise ValueError(f"{path}: malformed Wave header")
            values = []
            for token in tokens[1:]:
                match = re.fullmatch(r"SN\(z=([^()]+)\)", token)
                if match is None:
                    raise ValueError(f"{path}: malformed Wave header")
                values.append(float(match[1]))
            redshifts = _axis(values, "SNR redshift")
    if set(metadata) != {"BAND", "MAG", "EXPTIME", "NEXP"} or redshifts is None:
        raise ValueError(f"{path}: missing BAND/MAG/EXPTIME/NEXP/Wave header")
    if redshifts[0] < 0:
        raise ValueError(f"{path}: negative source redshift")
    return (
        metadata["BAND"],
        scalar(float(metadata["MAG"]), "MAG"),
        _positive(float(metadata["EXPTIME"]), "EXPTIME"),
        _positive(float(metadata["NEXP"]), "NEXP"),
        redshifts,
    )


@dataclass(frozen=True, init=False, eq=False)
class SNRReader:
    """Header-driven SNR tensor (magnitude,source_redshift,wavelength).

    smoothing must explicitly be 'legacy' (sigma=10 sample indices, reflect,
    truncate=4) or 'none'. Interpolation is linear in all three closed domains.
    """

    magnitudes: np.ndarray
    z: np.ndarray
    wavelength: np.ndarray
    raw_snr: np.ndarray
    smoothed_snr: np.ndarray
    provenance: object
    _interpolator: object
    _context: str

    def __init__(self, paths, *, smoothing, label="SNR"):
        _, interpolator, gaussian = _scipy()
        if smoothing not in ("legacy", "none"):
            raise ValueError("smoothing must be explicitly 'legacy' or 'none'")
        if isinstance(paths, (str, Path)):
            raise ValueError("paths must be a nonempty explicit sequence")
        paths = tuple(paths)
        if not paths:
            raise ValueError("paths must be a nonempty explicit sequence")
        rows = []
        caller = validate_label(label, "SNR label")
        for path in paths:
            path, digest, text, data = _read(path)
            band, mag, time, nexp, z = _header(text, path)
            wave = _axis(data[:, 0], "wavelength")
            if wave[0] <= 0 or data.shape[1] != len(z) + 1 or np.any(data[:, 1:] < 0):
                raise ValueError(
                    f"{path}: SNR column count, wavelength or negative raw SNR"
                )
            if rows:
                first = rows[0]
                if (
                    (band, time, nexp) != first[3:6]
                    or not np.array_equal(z, first[6])
                    or not np.array_equal(wave, first[7])
                ):
                    raise ValueError(f"{path}: inconsistent SNR metadata/grids")
            rows.append((mag, path, digest, band, time, nexp, z, wave, data[:, 1:].T))
        rows.sort(key=lambda r: r[0])
        mags = _axis([r[0] for r in rows], "SNR magnitude")
        band, time, nexp, z, wave = rows[0][3:8]
        raw = np.stack([r[8] for r in rows])
        smooth = (
            raw.copy()
            if smoothing == "none"
            else gaussian(raw, 10, axis=2, mode="reflect", truncate=4)
        )
        if not np.all(np.isfinite(smooth)) or np.any(smooth < 0):
            raise ValueError("smoothed SNR must be finite nonnegative")
        provenance = freeze(
            dict(
                label=caller,
                paths=[r[1] for r in rows],
                sha256=[r[2] for r in rows],
                band=band,
                exposure_time=time,
                exposure_count=nexp,
                smoothing=smoothing,
                sigma_samples=0 if smoothing == "none" else 10,
                mode="reflect",
                truncate=4,
                units="SNR per Angstrom",
                interpolation="linear RegularGridInterpolator",
                magnitudes=mags,
                source_redshifts=z,
                wavelengths=wave,
            )
        )
        for name, value in dict(
            magnitudes=mags,
            z=z,
            wavelength=wave,
            raw_snr=_immutable(raw),
            smoothed_snr=_immutable(smooth),
            provenance=provenance,
            _context=f"{caller} ({[r[1] for r in rows]})",
            _interpolator=interpolator(
                (mags, z, wave), smooth, method="linear", bounds_error=True
            ),
        ).items():
            object.__setattr__(self, name, value)

    def query(self, *, z_source, magnitudes, wavelength):
        """SNR per Angstrom (magnitude,) at explicit source z and observed lambda."""
        z = scalar(z_source, "z_source")
        wave = scalar(wavelength, "wavelength")
        _query(self.z, z, "source redshift", self._context)
        _query(self.wavelength, wave, "wavelength", self._context)
        m = _query(self.magnitudes, magnitudes, "magnitude", self._context)
        if m.ndim != 1 or not len(m):
            raise ValueError("magnitudes must be nonempty 1D")
        values = self._interpolator(
            np.column_stack((m, np.full(len(m), z), np.full(len(m), wave)))
        )
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            raise ValueError(f"{self._context}: nonfinite/negative interpolated SNR")
        return _immutable(values)

    def variance(
        self,
        *,
        z_source,
        magnitudes,
        wavelength,
        pixel_width_angstrom,
        exposure_count,
        exposure_time=None,
    ):
        """Dimensionless delta-flux variance, inverse SNR²/pixel/exposure ratio."""
        pixel = _positive(pixel_width_angstrom, "pixel_width_angstrom")
        nexp = _positive(exposure_count, "exposure_count")
        if (
            exposure_time is not None
            and _positive(exposure_time, "exposure_time")
            != self.provenance["exposure_time"]
        ):
            raise ValueError("incompatible per-exposure EXPTIME; no implicit rescaling")
        snr = self.query(
            z_source=z_source, magnitudes=magnitudes, wavelength=wavelength
        )
        if np.any(snr <= 0):
            raise ValueError(
                f"{self._context}: strictly positive SNR required for variance"
            )
        try:
            with np.errstate(
                over="raise", under="raise", invalid="raise", divide="raise"
            ):
                result = 1 / (
                    snr**2
                    * pixel
                    * (np.float64(nexp) / self.provenance["exposure_count"])
                )
        except FloatingPointError as error:
            raise ValueError(f"{self._context}: variance not representable") from error
        return _immutable(result)


def sample_forest_readers(
    density,
    snr,
    geometry,
    response,
    *,
    z_source,
    magnitudes,
    pixel_width_angstrom,
    exposure_count,
    exposure_time=None,
):
    """Query each adapter once; return immutable normalized rho/variance metadata.

    Pass these arrays to ForestInput with explicit magnitude quadrature, L_v and
    weighting method. The clustering wavelength uses z_eval; both source readers
    use z_source. No reader object survives this boundary.
    """
    if scalar(z_source, "z_source") <= geometry.z_eval:
        raise ValueError("z_source must exceed z_eval")
    wavelength = LYA_REST_ANGSTROM * (1 + geometry.z_eval)
    pixel = pixel_width_angstrom_to_velocity(
        pixel_width_angstrom, lambda_obs_angstrom=wavelength
    )
    if not np.isclose(
        pixel, response.pixel_width_velocity, rtol=8 * np.finfo(float).eps, atol=0
    ):
        raise ValueError(
            "simultaneous wavelength/velocity pixel widths are inconsistent"
        )
    rho = density_per_velocity(density.query(z_source, magnitudes), z_source=z_source)
    variance = snr.variance(
        z_source=z_source,
        magnitudes=magnitudes,
        wavelength=wavelength,
        pixel_width_angstrom=pixel_width_angstrom,
        exposure_count=exposure_count,
        exposure_time=exposure_time,
    )
    return freeze(
        dict(
            rho=rho,
            variance=variance,
            provenance=dict(
                density=density.provenance,
                snr=snr.provenance,
                z_source=z_source,
                z_eval=geometry.z_eval,
                magnitudes=real_array(magnitudes, "magnitudes"),
                wavelength=wavelength,
                pixel_width_angstrom=pixel_width_angstrom,
                pixel_width_velocity=pixel,
                exposure_count=exposure_count,
            ),
        )
    )
