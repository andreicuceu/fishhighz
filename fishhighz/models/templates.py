"""Vega-format signed templates; optional libraries are preparation-only.

The source convention follows Vega bin/make_template.py: K is h_template/Mpc,
PK is full power and PKSB is smooth power in (Mpc/h_template)^3. No scientific
code or cosmological assets are copied from Vega. No decomposition is regenerated.
"""

from dataclasses import dataclass
from hashlib import sha256
from io import BytesIO
from pathlib import Path
from types import MappingProxyType
from typing import ClassVar, Mapping

import numpy as np

from .._arrays import integer, readonly, real_array, scalar
from ..kernels.templates import _evaluate

_ROUNDOFF = 64 * np.finfo(np.float64).eps


def _positive(value, name):
    result = scalar(value, name)
    if result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _redshift(value, name):
    result = scalar(value, name)
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _agrees(a, b):
    return abs(a - b) <= _ROUNDOFF * max(1.0, abs(a), abs(b))


def _metadata(values):
    output = {}
    for key, value in ({} if values is None else values).items():
        if not isinstance(key, str):
            raise ValueError("metadata keys must be strings")
        if isinstance(value, np.generic):
            value = value.item()
        if value is not None and not isinstance(value, (str, bool, int, float)):
            raise ValueError(f"metadata {key}: require a plain scalar or string")
        if isinstance(value, float) and not np.isfinite(value):
            raise ValueError(f"metadata {key}: require finite values")
        output[key] = value
    return MappingProxyType(output)


@dataclass(frozen=True, eq=False)
class PowerTemplate:
    """Prepared record returned by prepare_template/load_template.

    Arrays are owned read-only float64. ``components`` has (k,component) shape;
    ``coefficients`` has (4,interval,component) shape, descending cubic order in
    ln(k). Source samples, header metadata and h conventions remain auditable.
    Use the validated factories rather than constructing this record directly.
    """

    component_names: ClassVar[tuple[str, str]] = ("smooth", "wiggle")
    k_file: np.ndarray
    pk_file: np.ndarray
    pksb_file: np.ndarray
    k: np.ndarray
    full: np.ndarray
    components: np.ndarray
    log_k: np.ndarray
    coefficients: np.ndarray
    z_ref: float
    h_template: float
    h_fid: float
    metadata: Mapping
    source_path: str | None = None
    source_sha256: str | None = None

    @property
    def domain(self):
        """Closed converted domain in h_fid/Mpc, not forecast scale cuts."""
        return float(self.k[0]), float(self.k[-1])

    def evaluate(self, k, *, derivative=0, z=None, growth=None):
        """Return owned (query,2) smooth/wiggle powers or dP/dk.

        Parameters
        ----------
        k : array_like
            Nonempty 1D positive finite queries in h_fid/Mpc, in any order.
            The closed converted domain is enforced before logarithms.
        derivative : {0, 1}
            0 returns power in (Mpc/h_fid)^3; 1 returns its derivative with
            respect to k (not ln(k)), in (Mpc/h_fid)^4.
        z : float, optional
            Nonnegative redshift. Defaults to the exact reference redshift.
        growth : float, optional
            Explicit positive power factor G=[D(z)/D(z_ref)]^2, applied once.
            Required away from z_ref. At z_ref require unity within
            64*eps64*max(1,abs(G)); no other growth inference is performed.
        """
        derivative = integer(derivative, "derivative")
        if derivative not in (0, 1):
            raise ValueError("derivative must be 0 or 1")
        query = real_array(k, "template query k")
        if query.ndim != 1 or not len(query) or np.any(query <= 0):
            raise ValueError("template query k must be nonempty positive 1D")
        requested = (float(query.min()), float(query.max()))
        if requested[0] < self.k[0] or requested[1] > self.k[-1]:
            raise ValueError(
                f"requested k range {requested} outside template domain {self.domain} "
                "in h_fid/Mpc; provide sufficient template coverage, not extrapolation"
            )
        z = self.z_ref if z is None else _redshift(z, "z")
        if growth is None:
            if z != self.z_ref:
                raise ValueError(
                    "nonreference redshift requires explicit growth power factor G"
                )
            growth = 1.0
        else:
            growth = _positive(growth, "growth power factor G")
            if z == self.z_ref and not _agrees(growth, 1.0):
                raise ValueError("growth power factor G must equal 1 at z_ref")
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            output = (
                _evaluate(self.log_k, self.coefficients, query, derivative) * growth
            )
        if not np.all(np.isfinite(output)):
            raise ValueError(
                f"nonfinite template evaluation: derivative={derivative}, z={z}, "
                f"G={growth}, requested k range {requested}"
            )
        return np.array(output, dtype=np.float64, order="C", copy=True)


def prepare_template(k, pk, pksb, *, z_ref, h_template, h_fid, metadata=None):
    """Prepare a signed not-a-knot cubic spline in ln(k_fid) from source arrays.

    Parameters
    ----------
    k, pk, pksb : array_like
        Matching finite real 1D source samples; K strictly positive/increasing,
        at least four knots. PK is full power and PKSB is smooth power.
        Source units are h_template/Mpc and (Mpc/h_template)^3.
    z_ref, h_template, h_fid : float
        Explicit reference redshift (>=0) and positive source/fiducial h values.
        k_fid=k*h_template/h_fid; P_fid=P*(h_fid/h_template)^3.
    metadata : mapping, optional
        Copied plain scalar/string provenance only, never a source of growth.

    Returns
    -------
    PowerTemplate
        Owned source/converted arrays and coefficients. Component order is
        (smooth,wiggle), with wiggle=PK-PKSB. Neither amplitude is logged.
        SciPy is needed here only; repeated evaluation requires just NumPy.
    """
    z_ref = _redshift(z_ref, "z_ref")
    h_template, h_fid = _positive(h_template, "h_template"), _positive(h_fid, "h_fid")
    metadata = _metadata(metadata)
    k, pk, pksb = (
        real_array(value, name)
        for value, name in ((k, "K"), (pk, "PK"), (pksb, "PKSB"))
    )
    if k.ndim != 1 or len(k) < 4 or pk.shape != k.shape or pksb.shape != k.shape:
        raise ValueError(
            "K/PK/PKSB require matching scalar 1D samples and at least four knots"
        )
    if np.any(k <= 0) or np.any(k[1:] <= k[:-1]):
        raise ValueError(
            "K must be positive and strictly increasing; no sorting or merging"
        )
    with np.errstate(over="ignore", invalid="ignore", divide="ignore", under="ignore"):
        k_factor = np.float64(h_template) / h_fid
        p_factor = (np.float64(h_fid) / h_template) ** 3
        converted_k = k * k_factor
        full = pk * p_factor
        source_components = np.column_stack((pksb, pk - pksb))
        components = source_components * p_factor
    if (
        not np.isfinite(k_factor)
        or k_factor <= 0
        or not np.isfinite(p_factor)
        or p_factor <= 0
        or not np.all(np.isfinite(converted_k))
        or np.any(converted_k <= 0)
        or np.any(converted_k[1:] <= converted_k[:-1])
        or not np.all(np.isfinite(full))
        or not np.all(np.isfinite(components))
        or np.any((pk != 0) & (full == 0))
        or np.any((source_components != 0) & (components == 0))
    ):
        raise ValueError(
            "nonfinite or unrepresentable template decomposition/h-unit conversion"
        )
    log_k = np.log(converted_k)
    if np.any(log_k[1:] <= log_k[:-1]):
        raise ValueError("template knots are not distinct in float64 ln(k_fid)")
    try:
        from scipy.interpolate import CubicSpline
    except ImportError as error:
        raise ImportError(
            "template preparation requires SciPy; install 'fishhighz[templates]'"
        ) from error
    try:
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            spline = CubicSpline(
                log_k, components, axis=0, bc_type="not-a-knot", extrapolate=False
            )
    except (ValueError, np.linalg.LinAlgError) as error:
        raise ValueError(
            f"template spline coefficient preparation failed: {error}"
        ) from error
    if not np.all(np.isfinite(spline.c)):
        raise ValueError("nonfinite template spline coefficients")
    return PowerTemplate(
        *(
            readonly(value, np.float64)
            for value in (k, pk, pksb, converted_k, full, components, log_k, spline.c)
        ),
        z_ref,
        h_template,
        h_fid,
        metadata,
    )


def _required_metadata(header, key, explicit, *, divisor=1.0, positive=False):
    validate = _positive if positive else _redshift
    original = None
    if key in header:
        original = validate(header[key], f"PK header {key}") / divisor
        original = validate(original, f"PK header {key}/{divisor}")
    supplied = None if explicit is None else validate(explicit, f"explicit {key}")
    if original is None and supplied is None:
        raise ValueError(
            f"missing PK header {key}; supply {'h_template' if positive else 'z_ref'} explicitly"
        )
    if (
        original is not None
        and supplied is not None
        and not _agrees(original, supplied)
    ):
        raise ValueError(
            f"contradictory PK header {key} and explicit value: {original} vs {supplied}"
        )
    return original if original is not None else supplied


def _check_unit(name, unit):
    if unit is None:
        return
    # Deliberately bounded spelling table: no physical-Mpc or general conversions.
    normalized = str(unit).replace(" ", "").replace("**", "^")
    allowed = (
        {"h/Mpc", "hMpc^-1", "hMpc-1"}
        if name == "K"
        else {
            "(Mpc/h)^3",
            "(Mpc/h)3",
            "Mpc^3/h^3",
            "Mpc3/h3",
            "Mpc^3h^-3",
            "Mpc3h-3",
        }
    )
    if normalized not in allowed:
        raise ValueError(
            f"{name}: unsupported/inconsistent units {unit!r}; require "
            f"{'h/Mpc' if name == 'K' else '(Mpc/h)^3'} in the Vega h_template convention"
        )


def load_template(path, *, h_fid, z_ref=None, h_template=None):
    """Load the unique Vega PK binary table from a local file snapshot.

    Hash and parse the same bytes; retain resolved path/SHA-256 and plain header
    values after closing the file. Table ZREF/H0 supply z_ref/h_template=H0/100.
    Explicit inputs fill missing metadata; contradictions fail at
    64*eps64*max(1,abs(header_value),abs(explicit_value)). Table values take
    precedence when consistent. Primary-header values are not substitutes.
    Missing TUNIT implies the documented Vega h-scaled convention. Recognized
    spellings are documented in README; physical Mpc units are rejected.
    Astropy is needed only for this operation, SciPy only for preparation.
    """
    h_fid = _positive(h_fid, "h_fid")
    try:
        from astropy.io import fits
    except ImportError as error:
        raise ImportError(
            "FITS template loading requires Astropy; install 'fishhighz[templates]'"
        ) from error
    source = Path(path).expanduser().resolve(strict=True)
    content = source.read_bytes()
    digest = sha256(content).hexdigest()
    try:
        with fits.open(BytesIO(content), memmap=False) as hdus:
            candidates = [hdu for hdu in hdus if hdu.name.upper() == "PK"]
            if len(candidates) != 1 or not isinstance(candidates[0], fits.BinTableHDU):
                raise ValueError(
                    "require a unique named PK binary table (missing, ambiguous or wrong type)"
                )
            table = candidates[0]
            header = table.header
            names = [name.upper() for name in table.columns.names]
            values = []
            for name in ("K", "PK", "PKSB"):
                if names.count(name) != 1:
                    raise ValueError(f"PK table requires exactly one {name} column")
                column = table.columns[names.index(name)]
                _check_unit(name, column.unit)
                values.append(real_array(table.data[column.name], name))
            z_ref = _required_metadata(header, "ZREF", z_ref)
            h_template = _required_metadata(
                header, "H0", h_template, divisor=100.0, positive=True
            )
            # FITS header objects, comments and history objects never reach kernels.
            metadata = {
                key: None
                if isinstance(header[key], fits.card.Undefined)
                else header[key]
                for key in header
                if key not in ("", "COMMENT", "HISTORY")
            }
        prepared = prepare_template(
            *values, z_ref=z_ref, h_template=h_template, h_fid=h_fid, metadata=metadata
        )
    except (ValueError, OSError) as error:
        raise ValueError(f"template {source}: {error}") from error
    # The immutable provenance is attached only after validation/preparation.
    from dataclasses import replace

    return replace(prepared, source_path=str(source), source_sha256=digest)
