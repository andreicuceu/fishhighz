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
    """Import optional interpolation and smoothing routines.

    Returns
    -------
    spline : type
        scipy.interpolate.RectBivariateSpline.
    interpolator : type
        scipy.interpolate.RegularGridInterpolator.
    gaussian_filter : callable
        scipy.ndimage.gaussian_filter1d.

    Raises
    ------
    ImportError
        If the optional survey dependencies are unavailable.
    """
    try:
        from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator
        from scipy.ndimage import gaussian_filter1d
    except ImportError as error:
        raise ImportError(
            "raw survey readers require SciPy; install fishhighz[survey]"
        ) from error
    return RectBivariateSpline, RegularGridInterpolator, gaussian_filter1d


def _read(path):
    """Read and hash a numeric survey table.

    Parameters
    ----------
    path : path-like
        UTF-8 whitespace-delimited table.

    Returns
    -------
    path : str
        Resolved input path.
    digest : str
        SHA256 of the original bytes.
    text : str
        Decoded original contents, including headers.
    data : ndarray
        Finite numeric table, shape (n_row, n_column), in source units.

    Raises
    ------
    ValueError
        If the table is malformed or contains invalid numeric data.
    FileNotFoundError
        If the input path does not exist.
    """
    path = Path(path).resolve(strict=True)
    content = path.read_bytes()
    try:
        data = real_array(np.loadtxt(io.BytesIO(content), ndmin=2), str(path))
    except Exception as error:
        raise ValueError(f"{path}: malformed numeric table: {error}") from error
    return str(path), hashlib.sha256(content).hexdigest(), content.decode("utf-8"), data


def _axis(axis, name, minimum=1):
    """Validate an ordered interpolation axis.

    Parameters
    ----------
    axis : array_like
        Real coordinates, shape (n_node,), in the named coordinate units.
    name : str
        Coordinate label for errors.
    minimum : int, optional
        Minimum number of nodes; default 1.

    Returns
    -------
    coordinates : ndarray
        Immutable float64 coordinates, shape (n_node,).

    Raises
    ------
    ValueError
        If the axis is nonfinite, unordered, duplicated or too short.
    """
    axis = real_array(axis, name)
    if axis.ndim != 1 or len(axis) < minimum or np.any(np.diff(axis) <= 0):
        raise ValueError(f"{name}: require ordered unique axis with >= {minimum} nodes")
    return _immutable(axis)


def _query(axis, values, name, context):
    """Require coordinates to lie within a closed interpolation domain.

    Parameters
    ----------
    axis : ndarray
        Ordered domain nodes, shape (n_node,).
    values : array_like
        Query coordinates, in the same units as axis.
    name : str
        Coordinate label for errors.
    context : str
        Reader identity for error messages.

    Returns
    -------
    coordinates : ndarray
        Validated real query values, preserving input shape.

    Raises
    ------
    ValueError
        If any query is nonfinite or outside the closed domain.
    """
    values = real_array(values, name)
    if np.any((values < axis[0]) | (values > axis[-1])):
        raise ValueError(
            f"{context}: {name}={values.tolist()} outside closed domain [{axis[0]}, {axis[-1]}]"
        )
    return values


DENSITY_INTERPOLATIONS = ("piecewise_constant", "spline")


def _cell_edges(centres, widths, name, context):
    """Construct lyaforecast cell edges from tabulated centres and widths.

    Parameters
    ----------
    centres : ndarray
        Ordered cell centres, shape (n_cell,).
    widths : array_like
        Positive cell widths, scalar or shape (n_cell,), in centre units.
    name : str
        Coordinate label for errors.
    context : str
        Reader identity for error messages.

    Returns
    -------
    edges : ndarray
        Immutable cell edges, shape (n_cell + 1,).
    tiling_residual : float
        Maximum absolute gap or overlap between the nominal upper edge
        centre + width/2 of each cell and the lower edge of the next cell, in
        centre units. It does not exceed the roundoff tolerance
        64*eps64*max(1, max|centres|).

    Raises
    ------
    ValueError
        If the edges are not strictly increasing, or if the cells do not tile
        the axis, i.e. the tiling residual exceeds the roundoff tolerance.

    Notes
    -----
    As in lyaforecast/tracer.py, edges are the lower edges centre - width/2 and
    the upper edge of the last cell; each cell extends to the next lower edge.
    Cells whose nominal upper edges do not meet the next lower edge (for
    example centres rounded to two decimals, or widths that differ from the
    centre spacing) would be silently stretched or shrunk by this construction,
    changing the integral of the density over each cell relative to the
    tabulated count.  Such tilings are therefore rejected rather than
    recorded: the tolerance is the same absolute float64 coordinate-roundoff
    allowance used for the uniform-axis checks, so regular tables and
    explicit, contiguous, nonuniform widths pass.
    """
    widths = np.broadcast_to(np.asarray(widths, dtype=float), centres.shape)
    lower, upper = centres - widths / 2, centres + widths / 2
    edges = np.append(lower, upper[-1])
    if np.any(np.diff(edges) <= 0):
        raise ValueError(
            f"{context}: {name} cell edges centre - width/2 are not strictly increasing"
        )
    residual = float(np.max(np.abs(lower[1:] - upper[:-1]), initial=0.0))

    # A cell extends to the next lower edge, so any residual above roundoff
    # would change the integrated density of the cell relative to the table.
    tolerance = 64 * np.finfo(float).eps * max(1, np.max(np.abs(centres)))
    if residual > tolerance:
        raise ValueError(
            f"{context}: {name} cells do not tile the axis (maximum gap or "
            f"overlap between centre + width/2 and the next lower edge is "
            f"{residual:.3g}, above the roundoff tolerance {tolerance:.3g}); "
            f"piecewise-constant interpolation requires contiguous cells. "
            f"Supply a regular table or contiguous explicit redshift_widths, "
            f"or use interpolation='spline'"
        )
    return _immutable(edges), residual


@dataclass(frozen=True, eq=False)
class CellHistogram2D:
    """Piecewise-constant density on contiguous (redshift, magnitude) cells.

    Parameters
    ----------
    z_edges : ndarray
        Ordered redshift cell edges, shape (n_z + 1,).
    magnitude_edges : ndarray
        Ordered magnitude cell edges, shape (n_magnitude + 1,).
    values : ndarray
        Cell densities, shape (n_z, n_magnitude).
    fill_value : float
        Value returned outside the outer cell edges.

    Notes
    -----
    Cells are closed below and open above, so a query on an interior edge takes
    the upper cell, as in lyaforecast's Histogram2DInterpolator.
    """

    z_edges: np.ndarray
    magnitude_edges: np.ndarray
    values: np.ndarray
    fill_value: float = 0.0

    def __post_init__(self):
        """Pad cell values with the exterior fill value.

        Returns
        -------
        None
            No value is returned.
        """
        object.__setattr__(
            self,
            "_padded",
            _immutable(
                np.pad(
                    np.asarray(self.values, dtype=float),
                    1,
                    mode="constant",
                    constant_values=self.fill_value,
                )
            ),
        )

    def __call__(self, z, magnitudes):
        """Evaluate cell densities at broadcast coordinates.

        Parameters
        ----------
        z : array_like
            Redshift queries.
        magnitudes : array_like
            Magnitude queries, broadcastable against z.

        Returns
        -------
        density : ndarray
            Cell densities with the broadcast shape; fill_value outside the
            outer edges.
        """
        z, magnitudes = np.broadcast_arrays(
            np.asarray(z, dtype=float), np.asarray(magnitudes, dtype=float)
        )
        z_index = np.searchsorted(self.z_edges, z, side="right")
        magnitude_index = np.searchsorted(
            self.magnitude_edges, magnitudes, side="right"
        )
        # The upper outer edge belongs to the last cell (closed outer domain).
        z_index[z == self.z_edges[-1]] = len(self.z_edges) - 1
        magnitude_index[magnitudes == self.magnitude_edges[-1]] = (
            len(self.magnitude_edges) - 1
        )
        return self._padded[z_index, magnitude_index]


@dataclass(frozen=True, init=False, eq=False)
class DensityReader:
    """Density adapter for explicit cell counts per square degree.

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
    interpolation : {'piecewise_constant', 'spline'}, optional
        'piecewise_constant' (default) holds each cell density constant on
        centre +/- width/2 cells, so magnitude and redshift integrals reproduce the
        tabulated counts. It requires the cells to tile each axis: the maximum gap
        or overlap between centre + width/2 and the next lower edge must not exceed
        the coordinate roundoff tolerance below, otherwise ValueError is raised
        (regular tables and contiguous explicit nonuniform widths are accepted).
        'spline' selects the legacy quadratic RectBivariateSpline (kx=ky=2, s=0)
        through the cell centres and has no tiling requirement.

    Notes
    -----
    Uniform-axis tolerance is 64*eps64*max(1,max(abs(axis))), an absolute
    float64 coordinate-roundoff allowance; the same tolerance bounds the cell
    tiling residual under piecewise-constant interpolation. Magnitudes must
    remain uniform; irregular redshifts require explicit widths or
    legacy_first_spacing, and under piecewise-constant interpolation those widths
    must make the cells contiguous (for example explicit widths equal to the
    nodal spacing; legacy_first_spacing on irregular nodes is rejected).
    """

    z: np.ndarray
    magnitudes: np.ndarray
    raw_counts: np.ndarray
    redshift_widths: np.ndarray
    density: np.ndarray
    provenance: object
    interpolation: str
    z_edges: np.ndarray
    magnitude_edges: np.ndarray
    _spline: object
    _histogram: object
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
        interpolation="piecewise_constant",
    ):
        """Read a rectangular source-count table and prepare its density interpolant.

        Parameters
        ----------
        path : path-like
            Three columns: redshift, magnitude and cell counts per deg^2.
        semantics : str
            Required cell_count_per_deg2 declaration.
        target_density : float or None
            Positive target count per deg^2; None skips renormalization.
        z_norm_min : float or None
            Strict lower redshift threshold for normalization; None uses all rows.
        magnitude_bounds : tuple or None, optional
            Inclusive raw magnitude limits, allowing None endpoints; default None
            selects all magnitudes.
        redshift_widths : array_like or None, optional
            Positive physical cell widths in redshift, shape (n_z,); default None.
        width_policy : str or None, optional
            uniform or legacy_first_spacing when explicit widths are absent. None
            requires uniform spacing unless explicit widths are supplied.
        label : str, optional
            Population label for errors and provenance; default density.
        interpolation : str, optional
            piecewise_constant (default) or the legacy quadratic spline.

        Returns
        -------
        None
            No value is returned.

        Raises
        ------
        ValueError
            If table shape, grid, normalization support, width declarations or
            interpolation choice are invalid, or piecewise-constant cell edges are
            not strictly increasing or do not tile the redshift or magnitude axis
            within the coordinate roundoff tolerance.
        ImportError
            If SciPy is unavailable for the spline interpolation.

        Notes
        -----
        Reads the table and creates the selected interpolant. Counts are normalized over the selected raw cells before division by redshift and magnitude cell widths. Explicit physical widths and the legacy first-spacing convention remain distinct.
        """
        if interpolation not in DENSITY_INTERPOLATIONS:
            raise ValueError(
                f"interpolation must be one of {DENSITY_INTERPOLATIONS}, "
                f"got {interpolation!r}"
            )
        spline = _scipy()[0] if interpolation == "spline" else None
        if semantics != "cell_count_per_deg2":
            raise ValueError("require explicit cell_count_per_deg2 semantics")
        caller = validate_label(label, "density label")
        path, digest, _, data = _read(path)
        context = f"{caller} ({path})"
        if data.shape[1] != 3 or np.any(data[:, 2] < 0):
            raise ValueError(
                f"{context}: require three columns and nonnegative cell counts"
            )

        # Reconstruct the rectangular grid before assigning physical cell widths.
        redshift_grid, redshift_indices = np.unique(data[:, 0], return_inverse=True)
        magnitude_grid, magnitude_indices = np.unique(data[:, 1], return_inverse=True)
        redshift_grid, magnitude_grid = (
            _axis(redshift_grid, "redshift", 3),
            _axis(magnitude_grid, "magnitude", 3),
        )
        if redshift_grid[0] < 0:
            raise ValueError(f"{context}: negative redshift")
        flat = redshift_indices * len(magnitude_grid) + magnitude_indices
        if len(np.unique(flat)) != len(flat) or len(flat) != len(redshift_grid) * len(
            magnitude_grid
        ):
            raise ValueError(f"{context}: duplicate or missing grid cells")
        redshift_spacing, magnitude_spacing = (
            np.diff(redshift_grid),
            np.diff(magnitude_grid),
        )
        tolerance = 64 * np.finfo(float).eps * max(1, np.max(np.abs(magnitude_grid)))
        if np.any(np.abs(magnitude_spacing - magnitude_spacing[0]) > tolerance):
            raise ValueError(f"{context}: nonuniform magnitude grid spacing")
        magnitude_spacing = float(magnitude_spacing[0])
        if redshift_widths is not None:
            if width_policy is not None:
                raise ValueError(
                    f"{context}: explicit widths conflict with width_policy"
                )
            widths = real_array(redshift_widths, "redshift_widths")
            if widths.shape != redshift_grid.shape or np.any(widths <= 0):
                raise ValueError(
                    f"{context}: redshift_widths must be positive (n_z,) in sorted z order"
                )
            policy = "explicit"
        else:
            policy = "uniform" if width_policy is None else width_policy
            if policy not in ("uniform", "legacy_first_spacing"):
                raise ValueError(f"{context}: unknown width_policy {policy!r}")
            tolerance = 64 * np.finfo(float).eps * max(1, np.max(np.abs(redshift_grid)))
            if policy == "uniform" and np.any(
                np.abs(redshift_spacing - redshift_spacing[0]) > tolerance
            ):
                raise ValueError(
                    f"{context}: nonuniform redshift nodes require redshift_widths or explicit width_policy='legacy_first_spacing'"
                )
            widths = np.full(len(redshift_grid), redshift_spacing[0])
        widths = _immutable(widths)

        # Raw-table selection determines the count used for renormalization.
        raw_counts = np.empty((len(redshift_grid), len(magnitude_grid)))
        raw_counts[redshift_indices, magnitude_indices] = data[:, 2]
        magnitude_mask = np.ones(len(magnitude_grid), dtype=bool)
        bounds = None
        if magnitude_bounds is not None:
            if len(magnitude_bounds) != 2:
                raise ValueError("magnitude_bounds must have two endpoints")
            magnitude_min, magnitude_max = (
                None if x is None else scalar(x, "magnitude bound")
                for x in magnitude_bounds
            )
            if (
                magnitude_min is not None
                and magnitude_max is not None
                and magnitude_min > magnitude_max
            ):
                raise ValueError("reversed magnitude bounds")
            magnitude_mask = (
                magnitude_grid >= (-np.inf if magnitude_min is None else magnitude_min)
            ) & (magnitude_grid <= (np.inf if magnitude_max is None else magnitude_max))
            bounds = (magnitude_min, magnitude_max)
        threshold = None if z_norm_min is None else scalar(z_norm_min, "z_norm_min")
        redshift_mask = (
            np.ones(len(redshift_grid), dtype=bool)
            if threshold is None
            else redshift_grid > threshold
        )
        target = (
            None
            if target_density is None
            else _positive(target_density, "target_density")
        )

        # Normalize counts first, then divide by the declared cell measures.
        selected_counts = raw_counts * magnitude_mask
        try:
            with np.errstate(
                over="raise", under="raise", invalid="raise", divide="raise"
            ):
                total_raw_count = float(raw_counts.sum())
                normalization_count = float(selected_counts[redshift_mask].sum())
                if target is not None and normalization_count <= 0:
                    raise ValueError(f"{context}: zero target-normalization support")
                scale = (
                    1.0 if target is None else np.float64(target) / normalization_count
                )
                density = (
                    selected_counts * scale / (widths[:, None] * magnitude_spacing)
                )
        except FloatingPointError as error:
            raise ValueError(
                f"{context}: density normalization not representable"
            ) from error
        if not np.all(np.isfinite(density)):
            raise ValueError(f"{context}: nonfinite prepared density")
        if interpolation == "piecewise_constant":
            z_edges, z_tiling_residual = _cell_edges(
                redshift_grid, widths, "redshift", context
            )
            magnitude_edges, magnitude_tiling_residual = _cell_edges(
                magnitude_grid, magnitude_spacing, "magnitude", context
            )
            histogram = CellHistogram2D(z_edges, magnitude_edges, _immutable(density))
            spline_interpolant = None
            z_domain = (z_edges[0], z_edges[-1])
            magnitude_domain = (magnitude_edges[0], magnitude_edges[-1])
            interpolation_label = "piecewise_constant_cells"
        else:
            z_edges = magnitude_edges = histogram = None
            z_tiling_residual = magnitude_tiling_residual = None
            spline_interpolant = spline(
                redshift_grid,
                magnitude_grid,
                density,
                kx=2,
                ky=2,
                s=0,
                bbox=[
                    redshift_grid[0],
                    redshift_grid[-1],
                    magnitude_grid[0],
                    magnitude_grid[-1],
                ],
            )
            z_domain = (redshift_grid[0], redshift_grid[-1])
            magnitude_domain = (magnitude_grid[0], magnitude_grid[-1])
            interpolation_label = "RectBivariateSpline kx=2 ky=2 s=0"
        provenance = freeze(
            dict(
                path=path,
                sha256=digest,
                label=caller,
                semantics=semantics,
                original_measure=total_raw_count,
                selected_measure=normalization_count,
                magnitude_bounds=bounds,
                z_norm_min=threshold,
                target_density=target,
                scale=float(scale),
                dz=None if policy == "explicit" else float(widths[0]),
                dm=magnitude_spacing,
                width_policy=policy,
                redshift_widths=widths,
                redshift_axis=redshift_grid,
                magnitude_axis=magnitude_grid,
                width_order="reconstructed sorted redshift axis",
                units="deg^-2 redshift^-1 mag^-1",
                interpolation=interpolation_label,
                z_domain=z_domain,
                magnitude_domain=magnitude_domain,
                z_edges=z_edges,
                magnitude_edges=magnitude_edges,
                z_cell_tiling_residual=z_tiling_residual,
                magnitude_cell_tiling_residual=magnitude_tiling_residual,
            )
        )
        for name, value in dict(
            z=redshift_grid,
            magnitudes=magnitude_grid,
            raw_counts=_immutable(raw_counts),
            redshift_widths=widths,
            density=_immutable(density),
            provenance=provenance,
            interpolation=interpolation,
            z_edges=z_edges,
            magnitude_edges=magnitude_edges,
            _context=context,
            _spline=spline_interpolant,
            _histogram=histogram,
        ).items():
            object.__setattr__(self, name, value)

    def query(self, z, magnitudes):
        """Evaluate density at one redshift in requested magnitude order.

        Parameters
        ----------
        z : float
            Dimensionless redshift inside the tabulated domain: the closed cell-edge
            range for piecewise_constant, or the closed centre range for spline.
        magnitudes : array_like
            Magnitude nodes, shape (n_magnitude,), in requested order, inside the
            corresponding magnitude domain.

        Returns
        -------
        density : ndarray
            Immutable dN/(dz dm dOmega) in deg^-2 redshift^-1 mag^-1, shape
            (n_magnitude,).

        Raises
        ------
        ValueError
            If queries are invalid or interpolated densities are negative or
            nonfinite.
        """
        z = scalar(z, "redshift")
        histogram = self._histogram
        _query(
            self.z if histogram is None else self.z_edges, z, "redshift", self._context
        )
        magnitude_grid = _query(
            self.magnitudes if histogram is None else self.magnitude_edges,
            magnitudes,
            "magnitude",
            self._context,
        )
        if magnitude_grid.ndim != 1 or not len(magnitude_grid):
            raise ValueError("magnitudes must be nonempty 1D")
        result = (
            self._spline.ev(np.full(magnitude_grid.shape, z), magnitude_grid)
            if histogram is None
            else histogram(z, magnitude_grid)
        )
        if np.any(result < 0) or not np.all(np.isfinite(result)):
            raise ValueError(
                f"{self._context}: negative/nonfinite interpolated density at z={z}, magnitudes={magnitude_grid.tolist()}"
            )
        return _immutable(result)

    def local_galaxy_density(self, geometry, magnitudes, quadrature):
        """Convert the local density at z_eval into a comoving galaxy density.

        Parameters
        ----------
        geometry : BinGeometry
            Conversion geometry at the foreground evaluation redshift.
        magnitudes : array_like
            Magnitude nodes, shape (n_magnitude,), in requested order.
        quadrature : array_like
            Positive magnitude integration weights in mag, shape (n_magnitude,).

        Returns
        -------
        number_density : float
            Comoving density in (h_fid/Mpc)^3.

        Notes
        -----
        Uses the explicit local-redshift approximation, with no bin-volume averaging.
        """
        return local_galaxy_density(
            self.query(geometry.z_eval, magnitudes), quadrature, geometry
        )


def _header(text, path):
    """Parse SNR table metadata and source redshift columns.

    Parameters
    ----------
    text : str
        Table text containing comment headers.
    path : str
        Input identifier for errors.

    Returns
    -------
    band : str
        Photometric band.
    magnitude : float
        Source magnitude.
    exposure_time : float
        Per-exposure duration in seconds.
    exposure_count : float
        Number of reference exposures.
    redshifts : ndarray
        Immutable source-redshift nodes, shape (n_source_redshift,).

    Raises
    ------
    ValueError
        If required headers are missing, duplicated or malformed.
    """
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
        """Read SNR tables and prepare the magnitude/redshift/wavelength interpolator.

        Parameters
        ----------
        paths : sequence of path-like
            Explicit nonempty list of SNR tables, one per magnitude.
        smoothing : str
            Explicit legacy smoothing or none.
        label : str, optional
            Population label for errors and provenance; default SNR.

        Returns
        -------
        None
            No value is returned.

        Raises
        ------
        ValueError
            If tables, headers, grids or smoothing choices are inconsistent.
        ImportError
            If SciPy is unavailable.

        Notes
        -----
        Reads and hashes every table. Legacy smoothing uses a Gaussian width of 10 wavelength sample indices, reflect boundaries and truncation at four standard deviations. Interpolation is linear.
        """
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
            band, magnitude, exposure_time, exposure_count, redshift_grid = _header(
                text, path
            )
            wavelength_grid = _axis(data[:, 0], "wavelength")
            if (
                wavelength_grid[0] <= 0
                or data.shape[1] != len(redshift_grid) + 1
                or np.any(data[:, 1:] < 0)
            ):
                raise ValueError(
                    f"{path}: SNR column count, wavelength or negative raw SNR"
                )
            if rows:
                first = rows[0]
                if (
                    (band, exposure_time, exposure_count) != first[3:6]
                    or not np.array_equal(redshift_grid, first[6])
                    or not np.array_equal(wavelength_grid, first[7])
                ):
                    raise ValueError(f"{path}: inconsistent SNR metadata/grids")
            rows.append(
                (
                    magnitude,
                    path,
                    digest,
                    band,
                    exposure_time,
                    exposure_count,
                    redshift_grid,
                    wavelength_grid,
                    data[:, 1:].T,
                )
            )

        # Match the interpolator axis order: magnitude, source redshift, wavelength.
        rows.sort(key=lambda r: r[0])
        magnitude_grid = _axis([r[0] for r in rows], "SNR magnitude")
        band, exposure_time, exposure_count, redshift_grid, wavelength_grid = rows[0][
            3:8
        ]
        raw_snr = np.stack([r[8] for r in rows])
        smoothed_snr = (
            raw_snr.copy()
            if smoothing == "none"
            else gaussian(raw_snr, 10, axis=2, mode="reflect", truncate=4)
        )
        if not np.all(np.isfinite(smoothed_snr)) or np.any(smoothed_snr < 0):
            raise ValueError("smoothed SNR must be finite nonnegative")
        provenance = freeze(
            dict(
                label=caller,
                paths=[r[1] for r in rows],
                sha256=[r[2] for r in rows],
                band=band,
                exposure_time=exposure_time,
                exposure_count=exposure_count,
                smoothing=smoothing,
                sigma_samples=0 if smoothing == "none" else 10,
                mode="reflect",
                truncate=4,
                units="SNR per Angstrom",
                interpolation="linear RegularGridInterpolator",
                magnitudes=magnitude_grid,
                source_redshifts=redshift_grid,
                wavelengths=wavelength_grid,
            )
        )
        for name, value in dict(
            magnitudes=magnitude_grid,
            z=redshift_grid,
            wavelength=wavelength_grid,
            raw_snr=_immutable(raw_snr),
            smoothed_snr=_immutable(smoothed_snr),
            provenance=provenance,
            _context=f"{caller} ({[r[1] for r in rows]})",
            _interpolator=interpolator(
                (magnitude_grid, redshift_grid, wavelength_grid),
                smoothed_snr,
                method="linear",
                bounds_error=True,
            ),
        ).items():
            object.__setattr__(self, name, value)

    def query(self, *, z_source, magnitudes, wavelength):
        """Interpolate SNR inside the three closed table domains.

        Parameters
        ----------
        z_source : float
            Dimensionless redshift of the background source.
        magnitudes : array_like
            Magnitude nodes, shape (n_magnitude,), in requested order.
        wavelength : float
            Observed wavelength in angstrom.

        Returns
        -------
        snr : ndarray
            Immutable SNR per reference angstrom sampling and exposure set, shape
            (n_magnitude,).

        Raises
        ------
        ValueError
            If coordinates are outside the domain or SNR is invalid.
        """
        source_redshift = scalar(z_source, "z_source")
        observed_wavelength = scalar(wavelength, "wavelength")
        _query(self.z, source_redshift, "source redshift", self._context)
        _query(self.wavelength, observed_wavelength, "wavelength", self._context)
        magnitude_grid = _query(self.magnitudes, magnitudes, "magnitude", self._context)
        if magnitude_grid.ndim != 1 or not len(magnitude_grid):
            raise ValueError("magnitudes must be nonempty 1D")
        values = self._interpolator(
            np.column_stack(
                (
                    magnitude_grid,
                    np.full(len(magnitude_grid), source_redshift),
                    np.full(len(magnitude_grid), observed_wavelength),
                )
            )
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
        """Convert table SNR to dimensionless pixel delta-flux variance.

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
        variance : ndarray
            Immutable dimensionless variance, shape (n_magnitude,).

        Raises
        ------
        ValueError
            If SNR is nonpositive, exposure metadata disagree or variance is not
            representable.

        Notes
        -----
        Scales inverse squared SNR by pixel width and the ratio of exposure counts; the duration per exposure is not rescaled.
        """
        pixel = _positive(pixel_width_angstrom, "pixel_width_angstrom")
        n_exposures = _positive(exposure_count, "exposure_count")
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
                    * (np.float64(n_exposures) / self.provenance["exposure_count"])
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
    sample : mappingproxy
        Immutable rho in deg^-2 (km/s)^-1 mag^-1 and dimensionless variance
        arrays, both shape (n_magnitude,), with reader provenance.

    Raises
    ------
    ValueError
        If source/foreground order or wavelength/velocity pixel widths disagree.

    Notes
    -----
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
