"""Survey source-density and forest S/N tables, adapted from lyaforecast (GPLv3).

Provenance: lyaforecast/tracer.py (_setup_dndzdm_lya/_tracer) and
spectrograph.py (_setup_desi_spectro, get_pixel_rms_noise). No raw assets are
distributed. SciPy is imported only when an S/N table is prepared or
interpolated.

Two layers are provided:

- ``DensityTable`` and ``SNRTable`` read and validate the raw tables, retain
  their explicit normalization conventions and raise on any query outside the
  closed table domains.
- ``SourceDensity`` and ``ForestSNR`` sample those tables with the lyaforecast
  boundary conventions adopted by the native INI path: the source density is
  1e-20 outside the cell edges, and the S/N is clamped at the brightest
  tabulated magnitude, floored at 1e-10 after pixel/exposure scaling, and
  replaced by a 1e20 variance sentinel outside the table. Sample records are
  owned snapshots; no mutable counters or table objects enter forecasts.

Source densities are piecewise constant on contiguous (redshift, magnitude)
cells, as in lyaforecast's Histogram2DInterpolator.
"""

import hashlib
import io
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ._arrays import label as validate_label
from ._arrays import real_array, scalar
from .geometry import LYA_REST_ANGSTROM, _immutable, _positive
from .noise import local_galaxy_density
from .response import pixel_width_angstrom_to_velocity
from .survey import freeze
from .weights import density_per_velocity

# Source density assigned outside the tabulated cell edges, as in lyaforecast.
DENSITY_FLOOR = 1e-20


def _scipy():
    """Import the S/N interpolation and smoothing routines.

    Returns
    -------
    interpolator : type
        scipy.interpolate.RegularGridInterpolator.
    gaussian_filter : callable
        scipy.ndimage.gaussian_filter1d.

    Notes
    -----
    SciPy is a core dependency; the import is deferred so that importing
    FishHighz does not load it.
    """
    from scipy.interpolate import RegularGridInterpolator
    from scipy.ndimage import gaussian_filter1d

    return RegularGridInterpolator, gaussian_filter1d


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
        Table identity for error messages.

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


GRID_CHUNK_POINTS = 1 << 20


def _grid_vector(values, name):
    """Validate a nonempty one-dimensional finite float64 vector.

    Parameters
    ----------
    values : array_like
        Real coordinates, shape (n_value,), in the units of the named quantity.
    name : str
        Quantity name used in validation errors.

    Returns
    -------
    vector : ndarray
        Owned float64 copy, shape (n_value,).

    Raises
    ------
    ValueError
        If the input is not nonempty, one dimensional and finite.
    """
    vector = real_array(values, name)
    if vector.ndim != 1 or not len(vector):
        raise ValueError(f"{name} must be nonempty 1D")
    return vector


def _magnitude_vector(values):
    """Validate a nonempty magnitude vector.

    Parameters
    ----------
    values : array_like
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
    magnitudes = real_array(values, "magnitudes")
    if magnitudes.ndim != 1 or not len(magnitudes):
        raise ValueError("magnitudes must be nonempty 1D")
    return magnitudes


def _paired_points(z_source, wavelength, pixel_width_angstrom):
    """Validate paired source-redshift, wavelength and pixel-width vectors.

    Parameters
    ----------
    z_source : array_like
        Dimensionless source redshifts, shape (n_point,).
    wavelength : array_like
        Observed wavelengths in angstrom, shape (n_point,); strictly positive.
    pixel_width_angstrom : array_like
        Observed pixel widths in angstrom, shape (n_point,); strictly positive.

    Returns
    -------
    z_source, wavelength, pixel_width_angstrom : ndarray
        Owned float64 vectors, each of shape (n_point,).

    Raises
    ------
    ValueError
        If any vector is nonfinite, not one dimensional, of unequal length, or
        if wavelengths or pixel widths are not strictly positive.
    """
    z_source = _grid_vector(z_source, "z_source")
    wavelength = _grid_vector(wavelength, "wavelength")
    pixel_width_angstrom = _grid_vector(pixel_width_angstrom, "pixel_width_angstrom")
    if not len(z_source) == len(wavelength) == len(pixel_width_angstrom):
        raise ValueError(
            "z_source, wavelength and pixel_width_angstrom must have equal length"
        )
    if np.any(wavelength <= 0):
        raise ValueError("wavelength must be positive and representable")
    if np.any(pixel_width_angstrom <= 0):
        raise ValueError("pixel_width_angstrom must be positive and representable")
    return z_source, wavelength, pixel_width_angstrom


def _row_chunks(n_row, n_column):
    """Split rows into slices holding at most GRID_CHUNK_POINTS points.

    Parameters
    ----------
    n_row : int
        Number of rows (paired points).
    n_column : int
        Number of columns (magnitude nodes) per row.

    Yields
    ------
    rows : slice
        Contiguous row selection bounding the interpolation memory.
    """
    step = max(1, GRID_CHUNK_POINTS // max(1, int(n_column)))
    for start in range(0, n_row, step):
        yield slice(start, min(start + step, n_row))


def _axis_interval(axis, values):
    """Locate coordinates in the intervals of a grid axis.

    Parameters
    ----------
    axis : ndarray
        Strictly increasing grid nodes, shape (n_node,), n_node >= 2.
    values : ndarray
        Coordinates inside the closed domain, shape (n,).

    Returns
    -------
    index : ndarray of int, shape (n,)
        Lower node of the interval holding each coordinate; coordinates at the
        upper edge belong to the last interval.
    fraction : ndarray, shape (n,)
        (value - axis[index]) / (axis[index + 1] - axis[index]) in [0, 1], the
        normalized distance used by SciPy's linear RegularGridInterpolator.
    """
    index = np.clip(np.searchsorted(axis, values, side="right") - 1, 0, len(axis) - 2)
    lower = axis[index]
    return index, (values - lower) / (axis[index + 1] - lower)


def snr_grid_separable(interpolator, z_source, wavelength, magnitudes):
    """Interpolate an SNR tensor on paired (z, wavelength) points times magnitudes.

    Parameters
    ----------
    interpolator : scipy.interpolate.RegularGridInterpolator
        Linear interpolator over (magnitude, source redshift, wavelength);
        its grid and values define the table.
    z_source, wavelength : ndarray of shape (n_point,)
        Paired coordinates inside the closed table domains.
    magnitudes : ndarray of shape (n_magnitude,)
        Magnitude coordinates inside the closed table domain.

    Returns
    -------
    snr : ndarray of shape (n_point, n_magnitude)
        Trilinear interpolant at (magnitudes[j], z_source[p], wavelength[p]).

    Notes
    -----
    The trilinear kernel factorizes: the four (z, wavelength) corner weights of
    each point contract the table into one profile over the tabulated
    magnitudes, which the per-magnitude linear weights then interpolate. This
    equals the generic 3-D evaluation up to floating-point rounding (relative
    ~1e-15) at a small fraction of its cost, since the corner search and
    weights are formed per point and per magnitude rather than per element.
    """
    magnitude_axis, z_axis, wavelength_axis = interpolator.grid
    table = np.ascontiguousarray(np.moveaxis(interpolator.values, 0, -1))
    z_index, z_fraction = _axis_interval(z_axis, z_source)
    wavelength_index, wavelength_fraction = _axis_interval(wavelength_axis, wavelength)
    magnitude_index, magnitude_fraction = _axis_interval(magnitude_axis, magnitudes)

    # Bilinear contraction over the (z, wavelength) cell: profile (n_point, n_table_mag).
    z_low, z_high = 1 - z_fraction, z_fraction
    wavelength_low, wavelength_high = 1 - wavelength_fraction, wavelength_fraction
    profile = table[z_index, wavelength_index] * (z_low * wavelength_low)[:, None]
    profile += table[z_index, wavelength_index + 1] * (z_low * wavelength_high)[:, None]
    profile += table[z_index + 1, wavelength_index] * (z_high * wavelength_low)[:, None]
    profile += (
        table[z_index + 1, wavelength_index + 1] * (z_high * wavelength_high)[:, None]
    )

    # Linear interpolation over magnitude.
    snr = np.take(profile, magnitude_index, axis=1)
    snr *= (1 - magnitude_fraction)[None, :]
    snr += np.take(profile, magnitude_index + 1, axis=1) * magnitude_fraction[None, :]
    return snr


def snr_grid_values(interpolator, z_source, wavelength, magnitudes):
    """Interpolate the SNR table on paired points times magnitudes.

    Parameters
    ----------
    interpolator : callable
        Linear interpolator over (magnitude, source redshift, wavelength),
        normally a SciPy RegularGridInterpolator.
    z_source, wavelength : ndarray of shape (n_point,)
        Paired coordinates inside the closed table domains.
    magnitudes : ndarray of shape (n_magnitude,)
        Magnitude coordinates inside the closed table domain.

    Returns
    -------
    snr : ndarray of shape (n_point, n_magnitude)
        Interpolated SNR.

    Notes
    -----
    A RegularGridInterpolator with at least two nodes on every axis uses the
    separable trilinear evaluation (snr_grid_separable). A table with a single
    node on some axis (one magnitude file, say) has no interval to locate, so
    it, like any other callable, is evaluated with one generic call on the
    (n_point, n_magnitude, 3) coordinate array.
    """
    regular_grid_interpolator, _ = _scipy()
    if isinstance(interpolator, regular_grid_interpolator) and all(
        len(axis) >= 2 for axis in interpolator.grid
    ):
        return snr_grid_separable(interpolator, z_source, wavelength, magnitudes)
    points = np.empty((len(z_source), len(magnitudes), 3))
    points[..., 0] = magnitudes[None, :]
    points[..., 1] = z_source[:, None]
    points[..., 2] = wavelength[:, None]
    return interpolator(points)


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
        Table identity for error messages.

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
            f"piecewise-constant cells must be contiguous. Supply a regular "
            f"table or contiguous explicit redshift_widths"
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
class DensityTable:
    """Source-density table of explicit cell counts per square degree.

    Parameters
    ----------
    path : path-like
        Three-column (z, magnitude, cell_count_per_deg2) rectangular table.
    semantics : str
        Must be 'cell_count_per_deg2'; normalized arrays bypass this table.
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
        uniform z. Explicit legacy_first_spacing uses z[1]-z[0] for all rows;
        this is a compatibility convention, not physical width inference.
        Cannot be combined with redshift_widths.
    label : str
        Caller population label, used only for provenance and errors.

    Notes
    -----
    The density is held constant on centre +/- width/2 cells, so magnitude and
    redshift integrals reproduce the tabulated counts. The cells must tile each
    axis: the maximum gap or overlap between centre + width/2 and the next lower
    edge must not exceed the coordinate roundoff tolerance 64*eps64*max(1,
    max(abs(axis))), the same absolute float64 allowance used for the
    uniform-axis checks. Magnitudes must remain uniform; irregular redshifts
    require contiguous explicit widths (for example widths equal to the nodal
    spacing); legacy_first_spacing on irregular nodes is rejected.
    """

    z: np.ndarray
    magnitudes: np.ndarray
    raw_counts: np.ndarray
    redshift_widths: np.ndarray
    density: np.ndarray
    provenance: object
    z_edges: np.ndarray
    magnitude_edges: np.ndarray
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
    ):
        """Read a rectangular source-count table and prepare its cell densities.

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

        Returns
        -------
        None
            No value is returned.

        Raises
        ------
        ValueError
            If table shape, grid, normalization support or width declarations
            are invalid, or the cell edges are not strictly increasing or do not
            tile the redshift or magnitude axis within the coordinate roundoff
            tolerance.

        Notes
        -----
        Counts are normalized over the selected raw cells before division by
        the redshift and magnitude cell widths. Explicit physical widths and the
        legacy first-spacing convention remain distinct.
        """
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

        # Contiguous cells centre +/- width/2 on both axes.
        z_edges, z_tiling_residual = _cell_edges(
            redshift_grid, widths, "redshift", context
        )
        magnitude_edges, magnitude_tiling_residual = _cell_edges(
            magnitude_grid, magnitude_spacing, "magnitude", context
        )
        histogram = CellHistogram2D(z_edges, magnitude_edges, _immutable(density))
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
                interpolation="piecewise_constant_cells",
                z_domain=(z_edges[0], z_edges[-1]),
                magnitude_domain=(magnitude_edges[0], magnitude_edges[-1]),
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
            z_edges=z_edges,
            magnitude_edges=magnitude_edges,
            _context=context,
            _histogram=histogram,
        ).items():
            object.__setattr__(self, name, value)

    def query(self, z, magnitudes):
        """Evaluate the cell density at one redshift in requested magnitude order.

        Parameters
        ----------
        z : float
            Dimensionless redshift inside the closed redshift cell-edge range.
        magnitudes : array_like
            Magnitude nodes, shape (n_magnitude,), in requested order, inside
            the closed magnitude cell-edge range.

        Returns
        -------
        density : ndarray
            Immutable dN/(dz dm dOmega) in deg^-2 redshift^-1 mag^-1, shape
            (n_magnitude,).

        Raises
        ------
        ValueError
            If queries are invalid or outside the cell edges, or densities are
            negative or nonfinite.
        """
        z = scalar(z, "redshift")
        _query(self.z_edges, z, "redshift", self._context)
        magnitude_grid = _query(
            self.magnitude_edges, magnitudes, "magnitude", self._context
        )
        if magnitude_grid.ndim != 1 or not len(magnitude_grid):
            raise ValueError("magnitudes must be nonempty 1D")
        result = self._histogram(z, magnitude_grid)
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
class SNRTable:
    """Header-driven SNR tensor (magnitude, source_redshift, wavelength).

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

        Notes
        -----
        Reads and hashes every table. Legacy smoothing uses a Gaussian width of
        10 wavelength sample indices, reflect boundaries and truncation at four
        standard deviations. Interpolation is linear.
        """
        interpolator, gaussian = _scipy()
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
        Scales inverse squared SNR by pixel width and the ratio of exposure
        counts; the duration per exposure is not rescaled.
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
        Original table provenance.

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
class SourceDensity:
    """Sample a validated DensityTable with the lyaforecast exterior floor.

    Coordinates outside the redshift or magnitude cell edges receive 1e-20, as
    in lyaforecast. In-domain cell densities are never negative, so the
    recorded negative-density policy ('floor_negative') never acts; zero cells
    are retained.

    Parameters
    ----------
    reader : DensityTable
        Validated source-density table. The attribute name is kept because
        provenance and magnitude-partition code read ``density.reader``.
    """

    reader: DensityTable

    def __post_init__(self):
        """Prepare the forest-normalized cell densities.

        Returns
        -------
        None
            No value is returned.

        Raises
        ------
        ValueError
            If the source table is not a validated DensityTable.

        Notes
        -----
        Retains the reference reduction order, including the flat zero-masked
        sum for forest normalization. The source table is unchanged.
        """
        if not isinstance(self.reader, DensityTable):
            raise ValueError("require validated DensityTable")
        # Reference forest normalization sums the entire masked flat array,
        # including zeroed redshift rows. Preserve its arithmetic ordering here;
        # the strict table's accepted selected-row reduction stays unchanged.
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
        object.__setattr__(
            self,
            "_histogram",
            CellHistogram2D(
                reader.z_edges, reader.magnitude_edges, values, DENSITY_FLOOR
            ),
        )

    def sample(self, z, magnitudes):
        """Evaluate the cell density with the exterior floor.

        Parameters
        ----------
        z : float
            Nonnegative dimensionless source redshift; redshifts outside the
            cell edges receive 1e-20.
        magnitudes : array_like
            Magnitude nodes, shape (n_magnitude,), in requested order.

        Returns
        -------
        record : mappingproxy
            Density values in deg^-2 redshift^-1 mag^-1, shape (n_magnitude,),
            with cell values and exterior-floor provenance.

        Raises
        ------
        ValueError
            If the redshift is invalid or negative.
        """
        z, magnitude_grid = scalar(z, "redshift"), _magnitude_vector(magnitudes)
        if z < 0:
            raise ValueError("redshift must be nonnegative")
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
                negative="floor_negative",
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
            If a redshift is negative or the inputs are invalid.
        """
        z = _grid_vector(z, "redshift")
        magnitude_grid = _magnitude_vector(magnitudes)
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
                        negative="floor_negative",
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


@dataclass(frozen=True)
class ForestSNR:
    """Forest pixel variance from an SNRTable with the lyaforecast fallbacks.

    Bright magnitudes are clamped to the brightest tabulated magnitude, the
    scaled SNR is floored at 1e-10, and coordinates outside the table receive
    the exposure-independent variance sentinel 1e20.

    Parameters
    ----------
    reader : SNRTable
        Validated SNR table. The attribute name is kept because provenance and
        magnitude-partition code read ``snr.reader``.
    """

    reader: SNRTable

    def __post_init__(self):
        """Prepare linear interpolation of the validated SNR table.

        Returns
        -------
        None
            No value is returned.

        Raises
        ------
        ValueError
            If the supplied object is not a validated SNRTable.
        """
        if not isinstance(self.reader, SNRTable):
            raise ValueError("require validated SNRTable")
        interpolator, _ = _scipy()
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
        Bright magnitudes clamp to the brightest node. Other out-of-range
        coordinates use variance 1e20 independently of exposure. In-domain SNR
        is scaled before the 1e-10 floor.
        """
        reader = self.reader
        source_redshift, observed_wavelength = (
            scalar(z_source, "z_source"),
            _positive(wavelength, "wavelength"),
        )
        magnitude_grid = _magnitude_vector(magnitudes)
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
        magnitude_grid = _magnitude_vector(magnitudes)
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


def sample_forest_tables(
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
    """Query each strict table once; return immutable normalized rho/variance.

    Parameters
    ----------
    density : DensityTable
        Prepared source-density table (any object exposing ``query``).
    snr : SNRTable
        Prepared signal-to-noise table (any object exposing ``variance``).
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
        arrays, both shape (n_magnitude,), with table provenance.

    Raises
    ------
    ValueError
        If source/foreground order or wavelength/velocity pixel widths disagree.

    Notes
    -----
    Pass these arrays to ForestInput with explicit magnitude quadrature, L_v and
    weighting method. The clustering wavelength uses z_eval; both source tables
    use z_source. No table object survives this boundary.
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


__all__ = [
    "DENSITY_FLOOR",
    "CellHistogram2D",
    "DensityTable",
    "ForestSNR",
    "SNRTable",
    "SourceDensity",
    "sample_forest_tables",
    "snr_grid_separable",
    "snr_grid_values",
]
