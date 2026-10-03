"""Integrated forest-source geometry: quadrature over source redshift and pixel.

A redshift bin observes the Lyman-alpha forest in the wavelength slice
[lambda_min, lambda_max].  Every background source whose forest overlaps that
slice contributes, with its pixels restricted to the slice.  With u = ln(lambda)
and y = ln(1 + z_q), the forest of a source at y occupies u in [y + a, y + b],
where [a, b] = ln[lambda_r_min, lambda_r_max] is the rest-frame forest range.
Sources contribute for y in [u1 - b, u2 - a] (the bin is [u1, u2]) and the
pixels of source y lie in the overlap [max(u1, y + a), min(u2, y + b)].

The measure carried by a quadrature node (pixel p of source y_p, magnitude node
j) is, in deg^-2,

    mu_pj = geom_p * density[y_index_p, j] * quadrature_j,
    geom_p = (1 + z_q,p) * w_y,p * c * w_u,p / L_bin,   L_bin = c (u2 - u1),

where w_y and w_u are the Gauss-Legendre weights in y and u, c is the speed of
light in km/s and density is dN/(dz dm dOmega) in deg^-2 per unit redshift per
magnitude.  The module is NumPy-only and performs no I/O; callers supply the
tabulated densities and variances through the vectorised reader queries.
"""

from dataclasses import dataclass, field
from functools import cached_property

import numpy as np

from ._arrays import immutable_float_array, integer, readonly, real_array
from .geometry import SPEED_LIGHT_KMS
from .magnitude import composite

# Floor of the legacy density adapter (adapters.legacy_compat.DENSITY_FLOOR).
# Duplicated here to keep this module free of adapter and survey imports; a test
# asserts that the two values agree.
LEGACY_DENSITY_FLOOR = 1e-20


def _positive_scalar(value, name):
    """Validate a finite strictly positive real scalar.

    Parameters
    ----------
    value : float
        Scalar in the units of the named quantity.
    name : str
        Quantity name used in validation errors.

    Returns
    -------
    value : float
        Validated scalar.

    Raises
    ------
    ValueError
        If the value is not a finite positive real scalar.
    """
    array = real_array(value, name)
    if array.ndim != 0 or array <= 0:
        raise ValueError(f"{name} must be a positive finite scalar")
    return float(array)


def _merge_breakpoints(lower, upper, candidates):
    """Merge candidate breakpoints into an ordered unique partition.

    Parameters
    ----------
    lower, upper : float
        Fixed partition ends (here ln(1+z_q) window limits).
    candidates : array_like
        Extra breakpoints; those outside the open interval (lower, upper) are
        discarded.

    Returns
    -------
    partition : ndarray of shape (n_boundary,)
        Sorted boundaries starting at lower and ending at upper.

    Notes
    -----
    Adjacent boundaries closer than 64*eps*max(1, |boundary|) are merged, the
    criterion used by magnitude.breakpoints. The fixed ends are never moved:
    interior candidates within the tolerance of either end are dropped.
    """
    eps = 64 * np.finfo(float).eps
    candidates = np.unique(np.asarray(candidates, dtype=np.float64))
    candidates = candidates[(candidates > lower) & (candidates < upper)]
    kept = [lower]
    for point in candidates:
        if point - kept[-1] > eps * max(1.0, abs(point)):
            kept.append(float(point))
    if len(kept) > 1 and upper - kept[-1] <= eps * max(1.0, abs(upper)):
        kept.pop()
    kept.append(upper)
    return np.array(kept, dtype=np.float64)


@dataclass(frozen=True, eq=False)
class IntegrationNodes:
    """Quadrature nodes over source redshift and observed forest pixel.

    Attributes
    ----------
    zq_nodes : ndarray of shape (n_y,)
        Source redshifts z_q = exp(y) - 1 at the y nodes.
    zq_weights : ndarray of shape (n_y,)
        (1 + z_q) * w_y, the quadrature weights in z_q (dz_q = (1+z_q) dy).
    y_index : ndarray of int64, shape (n_p,)
        Index into the y nodes of the source of each pixel node.
    z_q : ndarray of shape (n_p,)
        Source redshift of each pixel node, zq_nodes[y_index].
    lam_obs : ndarray of shape (n_p,)
        Observed wavelength in angstrom of each pixel node.
    velocity_weights : ndarray of shape (n_p,)
        c * w_u in km/s, the pixel-quadrature weights in velocity.
    bin_length_velocity : float
        L_bin = c (u2 - u1) in km/s.
    geom : ndarray of shape (n_p,)
        zq_weights[y_index] * velocity_weights / bin_length_velocity,
        dimensionless; multiplies density and magnitude weights to form the
        measure.
    z_pix : ndarray of shape (n_p,)
        Redshift of the absorbing gas, lam_obs / lambda_Lya - 1.
    info : dict
        Plain metadata: orders, panels, breakpoints, node counts, window.

    Notes
    -----
    Pixel nodes are ordered by y node, then by increasing wavelength within the
    y node; each y node carries lambda_panels * lambda_order pixel nodes.
    """

    zq_nodes: np.ndarray
    zq_weights: np.ndarray
    y_index: np.ndarray
    z_q: np.ndarray
    lam_obs: np.ndarray
    velocity_weights: np.ndarray
    bin_length_velocity: float
    geom: np.ndarray
    z_pix: np.ndarray
    info: dict = field(default_factory=dict)


def integration_nodes(
    *,
    lambda_min,
    lambda_max,
    rest_min,
    rest_max,
    zq_min,
    zq_max,
    zq_breaks,
    zq_order,
    lambda_order,
    lambda_panels,
    lya_rest_angstrom,
    label="",
):
    """Construct Gauss-Legendre nodes over source redshift and forest pixels.

    Parameters
    ----------
    lambda_min, lambda_max : float
        Observed-wavelength limits of the redshift-bin slice in angstrom,
        lambda_Lya * (1 + z_edge).
    rest_min, rest_max : float
        Rest-frame forest limits in angstrom, [lambda_r_min, lambda_r_max].
    zq_min, zq_max : float
        Source-redshift limits (the field's min/max forest source redshift).
    zq_breaks : array_like
        Source redshifts at which the density or S/N tables are non-smooth
        (density cell edges, S/N z_q nodes); those inside the window become
        panel boundaries. May be empty.
    zq_order : int
        Gauss-Legendre order per y panel.
    lambda_order : int
        Gauss-Legendre order per wavelength panel.
    lambda_panels : int
        Number of equal wavelength panels per y node.
    lya_rest_angstrom : float
        Lyman-alpha rest wavelength in angstrom, defining z_pix.
    label : str, optional
        Field and bin identifier used in error messages.

    Returns
    -------
    nodes : IntegrationNodes
        Quadrature nodes and weights; see the class documentation.

    Raises
    ------
    ValueError
        If any argument is invalid, the y window is empty, or no node has a
        positive overlap, i.e. the field has no forest coverage in the slice.

    Notes
    -----
    With y = ln(1+z_q), u = ln(lambda), a = ln(rest_min), b = ln(rest_max),
    the window is y in [u1 - b, u2 - a] intersected with
    [ln(1+zq_min), ln(1+zq_max)], and the overlap for each y is
    [max(u1, y + a), min(u2, y + b)]. The overlap length is piecewise linear
    in y with kinks at y = u1 - a and y = u2 - b, which are panel boundaries
    together with the window ends and ln(1+zq_breaks); the quadrature of
    overlap length is therefore exact, and sum(w_y * overlap) equals
    (b - a) (u2 - u1) when zq_min and zq_max do not truncate the window.
    Panel boundaries closer than 64 eps are merged as in magnitude.breakpoints.
    """
    lambda_min = _positive_scalar(lambda_min, "lambda_min")
    lambda_max = _positive_scalar(lambda_max, "lambda_max")
    rest_min = _positive_scalar(rest_min, "rest_min")
    rest_max = _positive_scalar(rest_max, "rest_max")
    lya_rest_angstrom = _positive_scalar(lya_rest_angstrom, "lya_rest_angstrom")
    zq_min = float(real_array(zq_min, "zq_min"))
    zq_max = float(real_array(zq_max, "zq_max"))
    if lambda_max <= lambda_min:
        raise ValueError("lambda_max must exceed lambda_min")
    if rest_max <= rest_min:
        raise ValueError("rest_max must exceed rest_min")
    if zq_min <= -1 or zq_max <= zq_min:
        raise ValueError("require -1 < zq_min < zq_max")
    zq_order = integer(zq_order, "zq_order", 1)
    lambda_order = integer(lambda_order, "lambda_order", 1)
    lambda_panels = integer(lambda_panels, "lambda_panels", 1)
    zq_breaks = real_array(zq_breaks, "zq_breaks").ravel()
    if np.any(zq_breaks <= -1):
        raise ValueError("zq_breaks must exceed -1")

    # Logarithmic coordinates of the slice and of the rest-frame forest.
    u_low, u_high = np.log(lambda_min), np.log(lambda_max)
    a_rest, b_rest = np.log(rest_min), np.log(rest_max)
    y_low = max(u_low - b_rest, np.log1p(zq_min))
    y_high = min(u_high - a_rest, np.log1p(zq_max))
    prefix = f"{label}: " if label else ""
    window = (
        f"lambda=[{lambda_min}, {lambda_max}], rest=[{rest_min}, {rest_max}], "
        f"zq limits=[{zq_min}, {zq_max}], "
        f"zq window=[{np.expm1(y_low)}, {np.expm1(y_high)}]"
    )
    if not y_high > y_low:
        raise ValueError(f"{prefix}no forest coverage: empty source window ({window})")

    # y panels: window ends, overlap kinks and table non-smooth points.
    candidates = np.r_[u_low - a_rest, u_high - b_rest, np.log1p(zq_breaks)]
    y_partition = _merge_breakpoints(y_low, y_high, candidates)
    y_nodes, y_weights = composite(y_partition, zq_order)

    # Overlap in u for every y node; positive in the open window interior.
    u_start = np.maximum(u_low, y_nodes + a_rest)
    u_stop = np.minimum(u_high, y_nodes + b_rest)
    overlap = u_stop - u_start
    if not np.all(overlap > 0):
        raise ValueError(
            f"{prefix}no forest coverage: no source node overlaps the slice ({window})"
        )

    # u panels: equal sub-intervals of each overlap, one Gauss-Legendre rule each.
    legendre_nodes, legendre_weights = np.polynomial.legendre.leggauss(lambda_order)
    panel_edges = (
        u_start[:, None]
        + overlap[:, None] * np.arange(lambda_panels + 1)[None, :] / lambda_panels
    )
    panel_lower, panel_upper = panel_edges[:, :-1], panel_edges[:, 1:]
    half_width = (panel_upper - panel_lower) / 2
    u_nodes = (
        (panel_lower + panel_upper)[:, :, None] / 2
        + half_width[:, :, None] * legendre_nodes[None, None, :]
    ).reshape(len(y_nodes), -1)
    u_weights = (half_width[:, :, None] * legendre_weights[None, None, :]).reshape(
        len(y_nodes), -1
    )
    n_per_source = u_nodes.shape[1]

    # Flatten pixel nodes y-major and form the measure geometry.
    y_index = np.repeat(np.arange(len(y_nodes), dtype=np.int64), n_per_source)
    zq_nodes = np.expm1(y_nodes)
    zq_weights = np.exp(y_nodes) * y_weights
    lam_obs = np.exp(u_nodes.ravel())
    velocity_weights = SPEED_LIGHT_KMS * u_weights.ravel()
    bin_length_velocity = SPEED_LIGHT_KMS * (u_high - u_low)
    geom = zq_weights[y_index] * velocity_weights / bin_length_velocity
    z_pix = lam_obs / lya_rest_angstrom - 1
    info = dict(
        label=label,
        zq_order=zq_order,
        lambda_order=lambda_order,
        lambda_panels=lambda_panels,
        y_breakpoints=y_partition.tolist(),
        n_y_panels=len(y_partition) - 1,
        n_y=len(y_nodes),
        n_pixel=int(len(y_index)),
        window_y=[float(y_low), float(y_high)],
        window_zq=[float(np.expm1(y_low)), float(np.expm1(y_high))],
        lambda_limits=[lambda_min, lambda_max],
        rest_limits=[rest_min, rest_max],
        zq_limits=[zq_min, zq_max],
        lya_rest_angstrom=lya_rest_angstrom,
        speed_light_kms=SPEED_LIGHT_KMS,
    )
    return IntegrationNodes(
        zq_nodes=readonly(zq_nodes, np.float64),
        zq_weights=readonly(zq_weights, np.float64),
        y_index=readonly(y_index, np.int64),
        z_q=readonly(zq_nodes[y_index], np.float64),
        lam_obs=readonly(lam_obs, np.float64),
        velocity_weights=readonly(velocity_weights, np.float64),
        bin_length_velocity=float(bin_length_velocity),
        geom=readonly(geom, np.float64),
        z_pix=readonly(z_pix, np.float64),
        info=info,
    )


def pixel_width_angstrom(nodes, pixel_width_velocity):
    """Convert the fixed velocity pixel width to angstrom at each pixel node.

    Parameters
    ----------
    nodes : IntegrationNodes
        Pixel nodes supplying lam_obs.
    pixel_width_velocity : float
        Positive pixel width of the bin in km/s.

    Returns
    -------
    width : ndarray of shape (n_p,)
        pixel_width_velocity * lam_obs / c in angstrom.

    Notes
    -----
    The S/N table is per angstrom, so the pixel variance scales as
    1/(SNR_A^2 * pixel_width_A). Holding the velocity width l_pix fixed across
    the bin, the noise power per unit velocity, sigma^2 * l_pix, becomes
    c / (lambda * SNR_A^2), independent of l_pix and equal to the result for a
    fixed width in angstrom. At the central wavelength lambda_c this reduces to
    the single-point value used by the central-source mode.
    """
    pixel_width_velocity = _positive_scalar(pixel_width_velocity, "pixel_width")
    return pixel_width_velocity * nodes.lam_obs / SPEED_LIGHT_KMS


@dataclass(frozen=True, eq=False)
class IntegratedForestSource:
    """Density, variance and measure of one forest field in one redshift bin.

    Attributes
    ----------
    nodes : IntegrationNodes
        Source-redshift and pixel quadrature.
    density : ndarray of shape (n_y, n_m)
        dN/(dz dm dOmega) in deg^-2 redshift^-1 mag^-1 at each y node.
    magnitudes : ndarray of shape (n_m,)
        Magnitude quadrature nodes.
    quadrature : ndarray of shape (n_m,)
        Magnitude quadrature weights in magnitudes.
    variance : ndarray of shape (n_p, n_m)
        Dimensionless pixel noise variance of each (pixel, magnitude) node.
        Shared, not copied, if already backed by immutable bytes.
    pixel_width_velocity : float
        Fixed velocity pixel width of the bin in km/s.
    info : dict
        Plain provenance metadata.

    Raises
    ------
    ValueError
        At construction, if shapes disagree, values are nonfinite, density is
        negative, variance is not strictly positive, or the width is not
        positive.
    """

    nodes: IntegrationNodes
    density: np.ndarray
    magnitudes: np.ndarray
    quadrature: np.ndarray
    variance: np.ndarray
    pixel_width_velocity: float
    info: dict = field(default_factory=dict)

    def __post_init__(self):
        """Validate shapes and values and store read-only float64 arrays.

        Raises
        ------
        ValueError
            If the inputs are inconsistent as described in the class docstring.
        """
        if not isinstance(self.nodes, IntegrationNodes):
            raise ValueError("nodes must be IntegrationNodes")
        density = real_array(self.density, "density")
        magnitudes = real_array(self.magnitudes, "magnitudes")
        quadrature = real_array(self.quadrature, "quadrature")
        # The variance has n_p * n_m elements; arrays that cannot change are shared.
        variance = immutable_float_array(self.variance, "variance")
        n_y, n_p = len(self.nodes.zq_nodes), len(self.nodes.y_index)
        if magnitudes.ndim != 1 or not len(magnitudes):
            raise ValueError("magnitudes must be nonempty 1D")
        n_m = len(magnitudes)
        if quadrature.shape != (n_m,):
            raise ValueError(f"quadrature must have shape ({n_m},)")
        if density.shape != (n_y, n_m):
            raise ValueError(f"density must have shape ({n_y}, {n_m})")
        if variance.shape != (n_p, n_m):
            raise ValueError(f"variance must have shape ({n_p}, {n_m})")
        if np.any(density < 0):
            raise ValueError("density must be nonnegative")
        if variance.size and variance.min() <= 0:
            raise ValueError("variance must be strictly positive")
        width = _positive_scalar(self.pixel_width_velocity, "pixel_width_velocity")
        for name, value in dict(
            density=density,
            magnitudes=magnitudes,
            quadrature=quadrature,
            variance=variance,
        ).items():
            value.flags.writeable = False
            object.__setattr__(self, name, value)
        object.__setattr__(self, "pixel_width_velocity", width)
        object.__setattr__(self, "info", dict(self.info))

    @cached_property
    def measure(self):
        """Quadrature measure of every (pixel, magnitude) node.

        Returns
        -------
        measure : ndarray of shape (n_p, n_m)
            mu_pj = geom_p * density[y_index_p, j] * quadrature_j in deg^-2
            (float64), read-only.

        Notes
        -----
        Evaluated on first access and cached on the (frozen) instance, since
        the array has n_p * n_m elements (about 1e7 at the default orders) and
        every access would otherwise rebuild it. The product is formed in place
        in one allocation with the same operation order as
        ``(geom[:, None] * density[y_index]) * quadrature[None, :]``, so the
        values are bit-identical to that expression.
        """
        measure = np.take(self.density, self.nodes.y_index, axis=0)
        measure *= self.nodes.geom[:, None]
        measure *= self.quadrature[None, :]
        measure.flags.writeable = False
        return measure


def integrated_forest_source(
    nodes,
    density,
    magnitudes,
    quadrature,
    variance,
    pixel_width_velocity,
    label="",
    info=None,
    density_floor=LEGACY_DENSITY_FLOOR,
):
    """Assemble an IntegratedForestSource and reject windows without sources.

    Parameters
    ----------
    nodes : IntegrationNodes
        Source-redshift and pixel quadrature.
    density : array_like of shape (n_y, n_m)
        Density in deg^-2 redshift^-1 mag^-1 at each y node.
    magnitudes : array_like of shape (n_m,)
        Magnitude quadrature nodes.
    quadrature : array_like of shape (n_m,)
        Magnitude quadrature weights.
    variance : array_like of shape (n_p, n_m)
        Dimensionless pixel variance.
    pixel_width_velocity : float
        Fixed velocity pixel width in km/s.
    label : str, optional
        Field and bin identifier used in error messages.
    info : mapping or None, optional
        Extra plain provenance merged into the source info.
    density_floor : float, optional
        Density at or below which a node counts as outside the source table;
        default 1e-20, the legacy density-adapter floor.

    Returns
    -------
    source : IntegratedForestSource
        Validated source.

    Raises
    ------
    ValueError
        If validation fails or every density is at most density_floor (the
        window lies outside the tabulated source density).
    """
    source = IntegratedForestSource(
        nodes=nodes,
        density=density,
        magnitudes=magnitudes,
        quadrature=quadrature,
        variance=variance,
        pixel_width_velocity=pixel_width_velocity,
        info=dict(info or {}),
    )
    if np.all(source.density <= density_floor):
        window = nodes.info.get("window_zq")
        prefix = f"{label}: " if label else ""
        raise ValueError(
            f"{prefix}no forest coverage: source density is at the floor "
            f"({density_floor}) everywhere in the zq window {window}"
        )
    return source


__all__ = [
    "LEGACY_DENSITY_FLOOR",
    "IntegratedForestSource",
    "IntegrationNodes",
    "integrated_forest_source",
    "integration_nodes",
    "pixel_width_angstrom",
]
