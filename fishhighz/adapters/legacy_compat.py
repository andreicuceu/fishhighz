"""Opt-in input fallbacks adapted from lyaforecast tracer/spectrograph (GPLv3).

Strict readers remain unchanged. Negative-density flooring is a named extension
beyond the reference, which preserves negative spline overshoot. Sample records
are owned snapshots; no mutable counters or reader objects enter forecasts.
"""

from dataclasses import dataclass

import numpy as np

from .._arrays import real_array, scalar
from ..geometry import LYA_REST_ANGSTROM, _positive
from ..response import pixel_width_angstrom_to_velocity
from ..survey import freeze
from ..weights import density_per_velocity
from .legacy_inputs import DensityReader, SNRReader, _scipy


def plain(value):
    """Convert immutable provenance to plain JSON-compatible data."""
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
    m = real_array(value, "magnitudes")
    if m.ndim != 1 or not len(m):
        raise ValueError("magnitudes must be nonempty 1D")
    return m


def _record(values, raw, original, effective, policies, masks, source):
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
    below 1e-20 are retained. Redshift extension is the reference spline behavior.
    """

    reader: DensityReader
    negative_policy: str

    def __post_init__(self):
        if not isinstance(self.reader, DensityReader):
            raise ValueError("require validated DensityReader")
        if self.negative_policy not in ("reject", "floor_negative"):
            raise ValueError("explicit negative_policy: reject or floor_negative")
        # Reference forest normalization sums the entire masked flat array,
        # including zeroed redshift rows. Preserve its arithmetic ordering here;
        # the strict reader's accepted selected-row reduction stays unchanged.
        r = self.reader
        values = r.raw_counts.copy()
        bounds = r.provenance["magnitude_bounds"]
        if bounds is not None:
            lo, hi = bounds
            if lo is not None:
                values *= r.magnitudes[None, :] >= lo
            if hi is not None:
                values *= r.magnitudes[None, :] <= hi
        threshold = r.provenance["z_norm_min"]
        selected = (
            np.ones(len(r.z), dtype=bool) if threshold is None else r.z > threshold
        )
        measure = (
            np.sum(values * selected[:, None])
            if bounds is not None
            else np.sum(values[selected])
        )
        target = r.provenance["target_density"]
        with np.errstate(over="raise", under="raise", invalid="raise", divide="raise"):
            if target is not None:
                values *= target / measure
            values /= r.redshift_widths[:, None] * r.provenance["dm"]
        spline, _, _ = _scipy()
        object.__setattr__(self, "_normalization_measure", float(measure))
        object.__setattr__(
            self, "_spline", spline(r.z, r.magnitudes, values, kx=2, ky=2, s=0)
        )

    def sample(self, z, magnitudes):
        """Return values, raw interpolants and immutable per-sample diagnostics."""
        z, m = scalar(z, "redshift"), _mags(magnitudes)
        if z < 0:
            raise ValueError("redshift must be nonnegative")
        r = self.reader
        outside = (m < r.magnitudes[0]) | (m > r.magnitudes[-1])
        raw = real_array(self._spline.ev(np.full(m.shape, z), m), "raw density")
        negative = (raw < 0) & ~outside
        if np.any(negative) and self.negative_policy == "reject":
            raise ValueError("negative interpolated density; floor_negative is opt-in")
        values = raw.copy()
        values[outside | negative] = 1e-20
        extended = np.full(m.shape, z < r.z[0] or z > r.z[-1])
        return _record(
            values,
            raw,
            dict(z=z, magnitudes=m),
            dict(
                z=float(np.clip(z, r.z[0], r.z[-1])),
                magnitudes=np.clip(m, r.magnitudes[0], r.magnitudes[-1]),
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
                density_floor=outside,
                negative_density=negative,
                redshift_extension=extended,
            ),
            r.provenance,
        )


@dataclass(frozen=True)
class LegacySNR:
    """Legacy bright clamp/sentinel and post-scaling SNR floor, per population."""

    reader: SNRReader

    def __post_init__(self):
        if not isinstance(self.reader, SNRReader):
            raise ValueError("require validated SNRReader")
        _, interpolator, _ = _scipy()
        r = self.reader
        object.__setattr__(
            self,
            "_interpolator",
            interpolator(
                (r.magnitudes, r.z, r.wavelength),
                r.smoothed_snr.copy(),
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
        """Return variance; sentinel branches equal 1e20 independent of exposure."""
        r = self.reader
        z, wave = scalar(z_source, "z_source"), _positive(wavelength, "wavelength")
        m = _mags(magnitudes)
        if z < 0:
            raise ValueError("source redshift must be nonnegative")
        pixel = _positive(pixel_width_angstrom, "pixel_width_angstrom")
        count = _positive(exposure_count, "exposure_count")
        if (
            exposure_time is not None
            and _positive(exposure_time, "exposure_time")
            != r.provenance["exposure_time"]
        ):
            raise ValueError("incompatible per-exposure EXPTIME")
        outside = (
            (m > r.magnitudes[-1])
            | (z < r.z[0])
            | (z > r.z[-1])
            | (wave < r.wavelength[0])
            | (wave > r.wavelength[-1])
        )
        bright = (m < r.magnitudes[0]) & ~outside
        effective = np.maximum(m, r.magnitudes[0])
        raw, scaled = np.zeros(m.shape), np.zeros(m.shape)
        values = np.full(m.shape, 1e20)
        inside = ~outside
        if np.any(inside):
            raw[inside] = real_array(
                self._interpolator(
                    np.column_stack(
                        (
                            effective[inside],
                            np.full(np.count_nonzero(inside), z),
                            np.full(np.count_nonzero(inside), wave),
                        )
                    )
                ),
                "interpolated SNR",
            )
            if np.any(raw < 0):
                raise ValueError("negative interpolated SNR")
            try:
                with np.errstate(
                    over="raise", under="raise", invalid="raise", divide="raise"
                ):
                    scaled[inside] = (
                        raw[inside]
                        * np.sqrt(pixel)
                        * np.sqrt(np.float64(count) / r.provenance["exposure_count"])
                    )
                    values[inside] = 1 / np.maximum(scaled[inside], 1e-10) ** 2
            except FloatingPointError as error:
                raise ValueError("scaled SNR/variance not representable") from error
        result = _record(
            values,
            raw,
            dict(z_source=z, wavelength=wave, magnitudes=m),
            dict(z_source=z, wavelength=wave, magnitudes=effective, evaluated=inside),
            dict(
                snr="legacy_floor_clamp",
                snr_floor=1e-10,
                sentinel_variance=1e20,
                pixel_width_angstrom=pixel,
                exposure_count=count,
                scaled_snr=scaled,
            ),
            dict(
                bright_clamp=bright,
                out_of_range=outside,
                snr_floor=inside & (scaled < 1e-10),
            ),
            r.provenance,
        )
        return result


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
    """Sample once into plain provenance suitable for ForestInput."""
    if scalar(z_source, "z_source") <= geometry.z_eval:
        raise ValueError("z_source must exceed z_eval")
    wave = LYA_REST_ANGSTROM * (1 + geometry.z_eval)
    pixel = pixel_width_angstrom_to_velocity(
        pixel_width_angstrom, lambda_obs_angstrom=wave
    )
    if not np.isclose(
        pixel, response.pixel_width_velocity, rtol=8 * np.finfo(float).eps, atol=0
    ):
        raise ValueError("inconsistent pixel widths")
    d = density.sample(z_source, magnitudes)
    s = snr.sample(
        z_source=z_source,
        magnitudes=magnitudes,
        wavelength=wave,
        pixel_width_angstrom=pixel_width_angstrom,
        exposure_count=exposure_count,
        exposure_time=exposure_time,
    )
    return dict(
        rho=density_per_velocity(d["values"], z_source=z_source),
        variance=s["values"],
        provenance=plain(
            dict(compatibility=True, density=d["provenance"], snr=s["provenance"])
        ),
    )
