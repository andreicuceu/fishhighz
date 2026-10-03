"""Opt-in input fallbacks adapted from lyaforecast tracer/spectrograph (GPLv3).

Strict readers remain unchanged. Negative-density flooring is a named extension
beyond the reference, which preserves negative spline overshoot. Piecewise-constant
densities follow lyaforecast's cell histogram, with 1e-20 outside the cell edges. Sample records
are owned snapshots; no mutable counters or reader objects enter forecasts.
"""

from dataclasses import dataclass

import numpy as np

from .._arrays import real_array, scalar
from ..geometry import LYA_REST_ANGSTROM, _positive
from ..response import pixel_width_angstrom_to_velocity
from ..survey import freeze
from ..weights import density_per_velocity
from .legacy_inputs import (
    CellHistogram2D,
    DensityReader,
    SNRReader,
    _grid_vector,
    _paired_points,
    _row_chunks,
    _scipy,
    snr_grid_values,
)

DENSITY_FLOOR = 1e-20


def plain(value):
    """Convert immutable provenance into plain serialization values.

    Parameters
    ----------
    value : object
        Nested immutable mappings, arrays, sequences or scalars.

    Returns
    -------
    converted : object
        Plain containers and scalar values preserving metadata.
    """
    from collections.abc import Mapping

    if isinstance(value, Mapping):
        return {k: plain(v) for k, v in value.items()}
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (tuple, list)):
        return [plain(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _mags(value):
    """Validate a nonempty magnitude vector.

    Parameters
    ----------
    value : array_like
        Magnitude nodes, shape (n_magnitude,).

    Returns
    -------
    magnitudes : ndarray
        Validated finite magnitude vector.

    Raises
    ------
    ValueError
        If magnitudes are empty, nonfinite or not one dimensional.
    """
    magnitudes = real_array(value, "magnitudes")
    if magnitudes.ndim != 1 or not len(magnitudes):
        raise ValueError("magnitudes must be nonempty 1D")
    return magnitudes


def _record(values, raw, original, effective, policies, masks, source):
    """Freeze sampled values and explicit fallback diagnostics.

    Parameters
    ----------
    values : ndarray
        Final density or variance values, shape (n_magnitude,).
    raw : ndarray
        Interpolated density or SNR before fallback/scaling, shape
        (n_magnitude,).
    original : mapping
        Requested coordinates.
    effective : mapping
        Coordinates and selections used by interpolation.
    policies : mapping
        Named fallback and scaling conventions.
    masks : mapping of ndarray
        Boolean fallback selections, each shape (n_magnitude,).
    source : mapping
        Original reader provenance.

    Returns
    -------
    record : mappingproxy
        Immutable sampled values and diagnostic metadata, including fallback
        counts.
    """
    return freeze(
        dict(
            values=values,
            raw=raw,
            provenance=dict(
                compatibility=True,
                original=original,
                effective=effective,
                policies=policies,
                masks=masks,
                counts={k: int(np.count_nonzero(v)) for k, v in masks.items()},
                raw_values=raw,
                source=source,
            ),
        )
    )


@dataclass(frozen=True)
class LegacyDensity:
    """Sample a validated DensityReader with explicit negative policy.

    negative_policy must be 'reject' or 'floor_negative'. The latter replaces
    negative in-domain-magnitude interpolants only; zero and positive values
    below 1e-20 are retained. The interpolation follows the reader. For the
    spline, redshift extension is the reference behavior. For piecewise-constant
    cells, coordinates outside the redshift or magnitude cell edges receive 1e-20,
    as in lyaforecast, and in-domain values cannot be negative.
    """

    reader: DensityReader
    negative_policy: str

    def __post_init__(self):
        """Prepare the explicit legacy density normalization and spline.

        Returns
        -------
        None
            No value is returned.

        Raises
        ------
        ValueError
            If the density reader or negative-density policy is invalid.

        Notes
        -----
        Retains the reference reduction order, including the flat zero-masked sum for forest normalization. The source reader is unchanged. The interpolant (spline or cell histogram) matches the reader's interpolation.
        """
        if not isinstance(self.reader, DensityReader):
            raise ValueError("require validated DensityReader")
        if self.negative_policy not in ("reject", "floor_negative"):
            raise ValueError("explicit negative_policy: reject or floor_negative")
        # Reference forest normalization sums the entire masked flat array,
        # including zeroed redshift rows. Preserve its arithmetic ordering here;
        # the strict reader's accepted selected-row reduction stays unchanged.
        reader = self.reader
        values = reader.raw_counts.copy()
        bounds = reader.provenance["magnitude_bounds"]
        if bounds is not None:
            magnitude_min, magnitude_max = bounds
            if magnitude_min is not None:
                values *= reader.magnitudes[None, :] >= magnitude_min
            if magnitude_max is not None:
                values *= reader.magnitudes[None, :] <= magnitude_max
        threshold = reader.provenance["z_norm_min"]
        redshift_mask = (
            np.ones(len(reader.z), dtype=bool)
            if threshold is None
            else reader.z > threshold
        )
        normalization_count = (
            np.sum(values * redshift_mask[:, None])
            if bounds is not None
            else np.sum(values[redshift_mask])
        )
        target = reader.provenance["target_density"]
        with np.errstate(over="raise", under="raise", invalid="raise", divide="raise"):
            if target is not None:
                values *= target / normalization_count
            values /= reader.redshift_widths[:, None] * reader.provenance["dm"]
        object.__setattr__(self, "_normalization_measure", float(normalization_count))
        if reader.interpolation == "piecewise_constant":
            object.__setattr__(self, "_spline", None)
            object.__setattr__(
                self,
                "_histogram",
                CellHistogram2D(
                    reader.z_edges, reader.magnitude_edges, values, DENSITY_FLOOR
                ),
            )
            return
        spline, _, _ = _scipy()
        object.__setattr__(self, "_histogram", None)
        object.__setattr__(
            self,
            "_spline",
            spline(reader.z, reader.magnitudes, values, kx=2, ky=2, s=0),
        )

    def sample(self, z, magnitudes):
        """Evaluate density and retain the selected legacy fallbacks.

        Parameters
        ----------
        z : float
            Nonnegative dimensionless source redshift; spline extension is
            permitted, and piecewise-constant cells return 1e-20 outside their
            redshift edges.
        magnitudes : array_like
            Magnitude nodes, shape (n_magnitude,), in requested order.

        Returns
        -------
        record : mappingproxy
            Density values in deg^-2 redshift^-1 mag^-1, shape (n_magnitude,), with
            raw interpolants and fallback provenance.

        Raises
        ------
        ValueError
            If redshift is invalid or negative density occurs under reject policy.

        Notes
        -----
        Out-of-domain magnitudes receive 1e-20. Negative interpolants receive the same value only under the explicitly selected floor_negative policy.
        """
        z, magnitude_grid = scalar(z, "redshift"), _mags(magnitudes)
        if z < 0:
            raise ValueError("redshift must be nonnegative")
        reader = self.reader
        if self._histogram is not None:
            return self._sample_cells(z, magnitude_grid)

        # Magnitude boundaries and negative spline overshoot have separate policies.
        outside_mask = (magnitude_grid < reader.magnitudes[0]) | (
            magnitude_grid > reader.magnitudes[-1]
        )
        raw = real_array(
            self._spline.ev(np.full(magnitude_grid.shape, z), magnitude_grid),
            "raw density",
        )
        negative_mask = (raw < 0) & ~outside_mask
        if np.any(negative_mask) and self.negative_policy == "reject":
            raise ValueError("negative interpolated density; floor_negative is opt-in")
        values = raw.copy()
        values[outside_mask | negative_mask] = 1e-20
        redshift_extension_mask = np.full(
            magnitude_grid.shape, z < reader.z[0] or z > reader.z[-1]
        )
        return _record(
            values,
            raw,
            dict(z=z, magnitudes=magnitude_grid),
            dict(
                z=float(np.clip(z, reader.z[0], reader.z[-1])),
                magnitudes=np.clip(
                    magnitude_grid, reader.magnitudes[0], reader.magnitudes[-1]
                ),
            ),
            dict(
                magnitude_domain="legacy_floor",
                density_floor=1e-20,
                redshift="legacy_spline_extension",
                negative=self.negative_policy,
                negative_is_reference_extension=self.negative_policy
                == "floor_negative",
                compatibility_normalization_measure=self._normalization_measure,
                normalization_reduction="legacy flat zero-masked forest sum; selected-row galaxy sum",
            ),
            dict(
                density_floor=outside_mask,
                negative_density=negative_mask,
                redshift_extension=redshift_extension_mask,
            ),
            reader.provenance,
        )

    def sample_grid(self, z, magnitudes):
        """Evaluate the cell density on a redshift-magnitude grid with floors.

        Parameters
        ----------
        z : array_like
            Nonnegative dimensionless source redshifts, shape (n_z,). Redshifts
            outside the cell edges receive 1e-20, as in sample.
        magnitudes : array_like
            Magnitude nodes, shape (n_magnitude,), in requested order.

        Returns
        -------
        record : mappingproxy
            values and raw densities in deg^-2 redshift^-1 mag^-1, shape
            (n_z, n_magnitude); element [i, j] equals sample(z[i],
            magnitudes)["values"][j]. Provenance counts are in grid elements.

        Raises
        ------
        ValueError
            If the reader uses spline interpolation (integrated forest sources
            support piecewise-constant cells only), or a redshift is negative.
        """
        if self._histogram is None:
            raise ValueError(
                "sample_grid supports piecewise_constant interpolation only"
            )
        z = _grid_vector(z, "redshift")
        magnitude_grid = _mags(magnitudes)
        if np.any(z < 0):
            raise ValueError("redshift must be nonnegative")
        reader, histogram = self.reader, self._histogram
        z_edges, magnitude_edges = histogram.z_edges, histogram.magnitude_edges

        # Same exterior-floor policy as the scalar cell sampling.
        magnitude_outside = (magnitude_grid < magnitude_edges[0]) | (
            magnitude_grid > magnitude_edges[-1]
        )
        redshift_outside = (z < z_edges[0]) | (z > z_edges[-1])
        outside_mask = redshift_outside[:, None] | magnitude_outside[None, :]
        raw = real_array(histogram(z[:, None], magnitude_grid[None, :]), "raw density")
        values = raw.copy()
        values[outside_mask] = DENSITY_FLOOR
        return freeze(
            dict(
                values=values,
                raw=raw,
                provenance=dict(
                    compatibility=True,
                    policies=dict(
                        magnitude_domain="legacy_floor_outside_cell_edges",
                        density_floor=DENSITY_FLOOR,
                        redshift="floor_outside_cells",
                        interpolation="piecewise_constant_cells",
                        negative=self.negative_policy,
                        negative_is_reference_extension=False,
                        compatibility_normalization_measure=self._normalization_measure,
                    ),
                    counts=dict(
                        density_floor=int(np.count_nonzero(outside_mask)),
                        negative_density=0,
                        redshift_extension=0,
                        redshift_outside_cells=int(
                            np.count_nonzero(
                                np.broadcast_to(redshift_outside[:, None], values.shape)
                            )
                        ),
                    ),
                    source=reader.provenance,
                ),
            )
        )

    def _sample_cells(self, z, magnitude_grid):
        """Evaluate the piecewise-constant density with exterior floors.

        Parameters
        ----------
        z : float
            Nonnegative dimensionless source redshift.
        magnitude_grid : ndarray
            Validated magnitude nodes, shape (n_magnitude,).

        Returns
        -------
        record : mappingproxy
            Density values in deg^-2 redshift^-1 mag^-1, shape (n_magnitude,), with
            cell values and exterior-floor provenance.
        """
        reader, histogram = self.reader, self._histogram
        z_edges, magnitude_edges = histogram.z_edges, histogram.magnitude_edges
        magnitude_outside = (magnitude_grid < magnitude_edges[0]) | (
            magnitude_grid > magnitude_edges[-1]
        )
        redshift_outside = np.full(
            magnitude_grid.shape, z < z_edges[0] or z > z_edges[-1]
        )
        raw = real_array(histogram(z, magnitude_grid), "raw density")
        values = raw.copy()
        values[magnitude_outside | redshift_outside] = DENSITY_FLOOR
        return _record(
            values,
            raw,
            dict(z=z, magnitudes=magnitude_grid),
            dict(
                z=float(np.clip(z, z_edges[0], z_edges[-1])),
                magnitudes=np.clip(
                    magnitude_grid, magnitude_edges[0], magnitude_edges[-1]
                ),
            ),
            dict(
                magnitude_domain="legacy_floor_outside_cell_edges",
                density_floor=DENSITY_FLOOR,
                redshift="floor_outside_cells",
                interpolation="piecewise_constant_cells",
                negative=self.negative_policy,
                negative_is_reference_extension=False,
                compatibility_normalization_measure=self._normalization_measure,
                normalization_reduction="legacy flat zero-masked forest sum; selected-row galaxy sum",
            ),
            dict(
                density_floor=magnitude_outside | redshift_outside,
                negative_density=np.zeros(magnitude_grid.shape, dtype=bool),
                redshift_extension=np.zeros(magnitude_grid.shape, dtype=bool),
                redshift_outside_cells=redshift_outside,
            ),
            reader.provenance,
        )


@dataclass(frozen=True)
class LegacySNR:
    """Legacy bright clamp/sentinel and post-scaling SNR floor, per population."""

    reader: SNRReader

    def __post_init__(self):
        """Prepare linear interpolation of the validated SNR table.

        Returns
        -------
        None
            No value is returned.

        Raises
        ------
        ValueError
            If the supplied object is not a validated SNRReader.
        """
        if not isinstance(self.reader, SNRReader):
            raise ValueError("require validated SNRReader")
        _, interpolator, _ = _scipy()
        reader = self.reader
        object.__setattr__(
            self,
            "_interpolator",
            interpolator(
                (reader.magnitudes, reader.z, reader.wavelength),
                reader.smoothed_snr.copy(),
                method="linear",
                bounds_error=True,
            ),
        )

    def sample(
        self,
        *,
        z_source,
        magnitudes,
        wavelength,
        pixel_width_angstrom,
        exposure_count,
        exposure_time=None,
    ):
        """Evaluate pixel variance with explicit legacy SNR fallbacks.

        Parameters
        ----------
        z_source : float
            Dimensionless redshift of the background source.
        magnitudes : array_like
            Magnitude nodes, shape (n_magnitude,), in requested order.
        wavelength : float
            Observed wavelength in angstrom.
        pixel_width_angstrom : float
            Positive observed pixel width in angstrom.
        exposure_count : float
            Positive number of exposures.
        exposure_time : float or None, optional
            Per-exposure duration in seconds; None uses the table duration, and an
            explicit value must match it.

        Returns
        -------
        record : mappingproxy
            Dimensionless variance and raw SNR, shape (n_magnitude,), with fallback
            and exposure provenance.

        Raises
        ------
        ValueError
            If coordinates or exposure settings are invalid, or scaled variance is
            not representable.

        Notes
        -----
        Bright magnitudes clamp to the brightest node. Other out-of-range coordinates use variance 1e20 independently of exposure. In-domain SNR is scaled before the 1e-10 floor.
        """
        reader = self.reader
        source_redshift, observed_wavelength = (
            scalar(z_source, "z_source"),
            _positive(wavelength, "wavelength"),
        )
        magnitude_grid = _mags(magnitudes)
        if source_redshift < 0:
            raise ValueError("source redshift must be nonnegative")
        pixel_width = _positive(pixel_width_angstrom, "pixel_width_angstrom")
        n_exposures = _positive(exposure_count, "exposure_count")
        if (
            exposure_time is not None
            and _positive(exposure_time, "exposure_time")
            != reader.provenance["exposure_time"]
        ):
            raise ValueError("incompatible per-exposure EXPTIME")

        # Out-of-domain samples retain the exposure-independent sentinel variance.
        outside_mask = (
            (magnitude_grid > reader.magnitudes[-1])
            | (source_redshift < reader.z[0])
            | (source_redshift > reader.z[-1])
            | (observed_wavelength < reader.wavelength[0])
            | (observed_wavelength > reader.wavelength[-1])
        )
        bright_mask = (magnitude_grid < reader.magnitudes[0]) & ~outside_mask
        effective_magnitudes = np.maximum(magnitude_grid, reader.magnitudes[0])
        raw_snr, scaled_snr = (
            np.zeros(magnitude_grid.shape),
            np.zeros(magnitude_grid.shape),
        )
        values = np.full(magnitude_grid.shape, 1e20)

        inside_mask = ~outside_mask
        if np.any(inside_mask):
            raw_snr[inside_mask] = real_array(
                self._interpolator(
                    np.column_stack(
                        (
                            effective_magnitudes[inside_mask],
                            np.full(np.count_nonzero(inside_mask), source_redshift),
                            np.full(np.count_nonzero(inside_mask), observed_wavelength),
                        )
                    )
                ),
                "interpolated SNR",
            )
            if np.any(raw_snr < 0):
                raise ValueError("negative interpolated SNR")
            try:
                with np.errstate(
                    over="raise", under="raise", invalid="raise", divide="raise"
                ):
                    # Apply pixel/exposure scaling before the reference SNR floor.
                    scaled_snr[inside_mask] = (
                        raw_snr[inside_mask]
                        * np.sqrt(pixel_width)
                        * np.sqrt(
                            np.float64(n_exposures)
                            / reader.provenance["exposure_count"]
                        )
                    )
                    values[inside_mask] = (
                        1 / np.maximum(scaled_snr[inside_mask], 1e-10) ** 2
                    )
            except FloatingPointError as error:
                raise ValueError("scaled SNR/variance not representable") from error
        result = _record(
            values,
            raw_snr,
            dict(
                z_source=source_redshift,
                wavelength=observed_wavelength,
                magnitudes=magnitude_grid,
            ),
            dict(
                z_source=source_redshift,
                wavelength=observed_wavelength,
                magnitudes=effective_magnitudes,
                evaluated=inside_mask,
            ),
            dict(
                snr="legacy_floor_clamp",
                snr_floor=1e-10,
                sentinel_variance=1e20,
                pixel_width_angstrom=pixel_width,
                exposure_count=n_exposures,
                scaled_snr=scaled_snr,
            ),
            dict(
                bright_clamp=bright_mask,
                out_of_range=outside_mask,
                snr_floor=inside_mask & (scaled_snr < 1e-10),
            ),
            reader.provenance,
        )
        return result

    def variance_grid(
        self,
        *,
        z_source,
        wavelength,
        magnitudes,
        pixel_width_angstrom,
        exposure_count,
        exposure_time=None,
    ):
        """Evaluate pixel variance on paired points times magnitudes (legacy policy).

        Parameters
        ----------
        z_source : array_like
            Nonnegative dimensionless source redshifts, shape (n_point,).
        wavelength : array_like
            Observed wavelengths in angstrom, shape (n_point,), paired with
            z_source.
        magnitudes : array_like
            Magnitude nodes, shape (n_magnitude,), in requested order.
        pixel_width_angstrom : array_like
            Positive observed pixel widths in angstrom, shape (n_point,).
        exposure_count : float
            Positive number of exposures.
        exposure_time : float or None, optional
            Per-exposure duration in seconds; None uses the table duration, and an
            explicit value must match it.

        Returns
        -------
        record : mappingproxy
            Dimensionless variance values, shape (n_point, n_magnitude), and
            provenance. Element [p, j] equals sample(...)["values"][j] at point p.
            Fallback counts (bright_clamp, out_of_range, snr_floor) count grid
            elements.

        Raises
        ------
        ValueError
            If coordinates or exposure settings are invalid, a raw interpolated
            SNR is negative or nonfinite, or a scaled variance is not
            representable.

        Notes
        -----
        Reproduces the scalar policy per element: out-of-range magnitude (faint),
        redshift or wavelength gives the 1e20 sentinel independently of exposure;
        bright magnitudes clamp to the brightest node; in-domain SNR is scaled by
        sqrt(pixel width) and sqrt(n_exposure/table exposures) before the 1e-10
        floor. The interpolator is called once per chunk of at most 2**20 points.
        """
        reader = self.reader
        z_source, observed_wavelength, pixel_width = _paired_points(
            z_source, wavelength, pixel_width_angstrom
        )
        magnitude_grid = _mags(magnitudes)
        if np.any(z_source < 0):
            raise ValueError("source redshift must be nonnegative")
        n_exposures = _positive(exposure_count, "exposure_count")
        if (
            exposure_time is not None
            and _positive(exposure_time, "exposure_time")
            != reader.provenance["exposure_time"]
        ):
            raise ValueError("incompatible per-exposure EXPTIME")

        # Sentinel selections factorize into a per-point and a per-magnitude part.
        point_inside = (
            (z_source >= reader.z[0])
            & (z_source <= reader.z[-1])
            & (observed_wavelength >= reader.wavelength[0])
            & (observed_wavelength <= reader.wavelength[-1])
        )
        magnitude_inside = magnitude_grid <= reader.magnitudes[-1]
        inside_mask = point_inside[:, None] & magnitude_inside[None, :]
        bright_mask = (magnitude_grid < reader.magnitudes[0]) & magnitude_inside
        effective_magnitudes = np.maximum(magnitude_grid, reader.magnitudes[0])

        # Sentinel variance outside the table domain; the inside block is
        # overwritten below, so only the outside rows and columns are filled.
        values = np.empty(inside_mask.shape)
        values[~point_inside, :] = 1e20
        values[:, ~magnitude_inside] = 1e20
        floor_count = 0
        point_rows = np.flatnonzero(point_inside)
        magnitude_columns = np.flatnonzero(magnitude_inside)
        exposure_ratio = np.float64(n_exposures) / reader.provenance["exposure_count"]
        if len(point_rows) and len(magnitude_columns):
            for rows in _row_chunks(len(point_rows), len(magnitude_columns)):
                chunk = point_rows[rows]
                # Fresh array owned here, so the scaling below is done in place
                # (same operations and order as the elementwise expressions).
                raw_snr = snr_grid_values(
                    self._interpolator,
                    z_source[chunk],
                    observed_wavelength[chunk],
                    effective_magnitudes[magnitude_columns],
                )
                if (
                    raw_snr.dtype != np.float64
                    or not raw_snr.flags.owndata
                    or not raw_snr.flags.writeable
                ):
                    # Not a fresh float64 array of ours: copy before scaling.
                    raw_snr = real_array(raw_snr, "interpolated SNR")
                if not np.all(np.isfinite(raw_snr)):
                    raise ValueError("interpolated SNR must be finite")
                if raw_snr.min() < 0:
                    raise ValueError("negative interpolated SNR")
                try:
                    with np.errstate(
                        over="raise", under="raise", invalid="raise", divide="raise"
                    ):
                        # Apply pixel/exposure scaling before the reference SNR floor.
                        scaled_snr = np.multiply(
                            raw_snr, np.sqrt(pixel_width[chunk, None]), out=raw_snr
                        )
                        np.multiply(scaled_snr, np.sqrt(exposure_ratio), out=scaled_snr)
                        floor_count += int(np.count_nonzero(scaled_snr < 1e-10))
                        np.maximum(scaled_snr, 1e-10, out=scaled_snr)
                        np.square(scaled_snr, out=scaled_snr)
                        variance_block = np.divide(1, scaled_snr, out=scaled_snr)
                except FloatingPointError as error:
                    raise ValueError("scaled SNR/variance not representable") from error
                if len(magnitude_columns) == values.shape[1]:
                    values[chunk] = variance_block
                else:
                    values[np.ix_(chunk, magnitude_columns)] = variance_block
        return freeze(
            dict(
                values=values,
                provenance=dict(
                    compatibility=True,
                    policies=dict(
                        snr="legacy_floor_clamp",
                        snr_floor=1e-10,
                        sentinel_variance=1e20,
                        exposure_count=n_exposures,
                    ),
                    counts=dict(
                        bright_clamp=int(
                            np.count_nonzero(point_inside)
                            * np.count_nonzero(bright_mask)
                        ),
                        out_of_range=int(np.count_nonzero(~inside_mask)),
                        snr_floor=floor_count,
                    ),
                    source=reader.provenance,
                ),
            )
        )


def sample_legacy_forest(
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
    """Sample source density and pixel variance with legacy provenance.

    Parameters
    ----------
    density : object
        Prepared source-density adapter.
    snr : object
        Prepared signal-to-noise adapter.
    geometry : BinGeometry
        Foreground geometry and evaluation redshift.
    response : InstrumentResponse
        Response with the independently specified velocity pixel width.
    z_source : float
        Dimensionless redshift of the background source.
    magnitudes : array_like
        Magnitude nodes, shape (n_magnitude,), in requested order.
    pixel_width_angstrom : float
        Positive observed pixel width in angstrom.
    exposure_count : float
        Positive number of exposures.
    exposure_time : float or None, optional
        Per-exposure duration in seconds; None uses the table duration, and an
        explicit value must match it.

    Returns
    -------
    sample : dict
        rho in deg^-2 (km/s)^-1 mag^-1 and dimensionless variance arrays, shape
        (n_magnitude,), plus plain reader/fallback provenance.

    Raises
    ------
    ValueError
        If source/foreground order or wavelength/velocity pixel widths disagree.

    Notes
    -----
    Density and SNR use source redshift; observed wavelength uses the foreground evaluation redshift.
    """
    if scalar(z_source, "z_source") <= geometry.z_eval:
        raise ValueError("z_source must exceed z_eval")
    observed_wavelength = LYA_REST_ANGSTROM * (1 + geometry.z_eval)
    pixel_width_velocity = pixel_width_angstrom_to_velocity(
        pixel_width_angstrom, lambda_obs_angstrom=observed_wavelength
    )
    if not np.isclose(
        pixel_width_velocity,
        response.pixel_width_velocity,
        rtol=8 * np.finfo(float).eps,
        atol=0,
    ):
        raise ValueError("inconsistent pixel widths")
    density_sample = density.sample(z_source, magnitudes)
    snr_sample = snr.sample(
        z_source=z_source,
        magnitudes=magnitudes,
        wavelength=observed_wavelength,
        pixel_width_angstrom=pixel_width_angstrom,
        exposure_count=exposure_count,
        exposure_time=exposure_time,
    )
    return dict(
        rho=density_per_velocity(density_sample["values"], z_source=z_source),
        variance=snr_sample["values"],
        provenance=plain(
            dict(
                compatibility=True,
                density=density_sample["provenance"],
                snr=snr_sample["provenance"],
            )
        ),
    )
