"""Geometry, quadrature and validation of the integrated forest-source nodes."""

import numpy as np
import pytest

from fishhighz.adapters.legacy_compat import DENSITY_FLOOR
from fishhighz.forest_integration import (
    LEGACY_DENSITY_FLOOR,
    IntegratedForestSource,
    IntegrationNodes,
    integrated_forest_source,
    integration_nodes,
    pixel_width_angstrom,
)
from fishhighz.geometry import LYA_REST_ANGSTROM, SPEED_LIGHT_KMS

REST_MIN, REST_MAX = 1040.0, 1200.0
Z_EDGES = (2.2, 2.4)


def make_nodes(**overrides):
    """Build nodes for a typical bin with breakpoints strictly inside the window.

    Parameters
    ----------
    **overrides : dict
        Replacements for the default integration_nodes arguments.

    Returns
    -------
    nodes : IntegrationNodes
        Quadrature nodes for the redshift bin Z_EDGES.
    """
    options = dict(
        lambda_min=LYA_REST_ANGSTROM * (1 + Z_EDGES[0]),
        lambda_max=LYA_REST_ANGSTROM * (1 + Z_EDGES[1]),
        rest_min=REST_MIN,
        rest_max=REST_MAX,
        zq_min=1.5,
        zq_max=6.0,
        zq_breaks=[2.5, 2.8, 3.0, 3.3],
        zq_order=4,
        lambda_order=3,
        lambda_panels=2,
        lya_rest_angstrom=LYA_REST_ANGSTROM,
        label="synthetic bin",
    )
    options.update(overrides)
    return integration_nodes(**options)


def overlap_length(nodes):
    """Sum the velocity weights of each y node and convert to ln(lambda) length.

    Parameters
    ----------
    nodes : IntegrationNodes
        Quadrature nodes.

    Returns
    -------
    overlap : ndarray of shape (n_y,)
        Overlap length in u = ln(lambda) per source node.
    """
    return np.bincount(nodes.y_index, weights=nodes.velocity_weights) / SPEED_LIGHT_KMS


def analytic_overlap(y, u1, u2, a, b):
    """Evaluate the overlap length min(u2, y+b) - max(u1, y+a).

    Parameters
    ----------
    y : ndarray
        ln(1 + z_q) values.
    u1, u2 : float
        ln of the bin wavelength limits.
    a, b : float
        ln of the rest-frame forest limits.

    Returns
    -------
    length : ndarray
        Overlap length, clipped at zero.
    """
    return np.clip(np.minimum(u2, y + b) - np.maximum(u1, y + a), 0, None)


def test_overlap_formula_and_identity():
    nodes = make_nodes(zq_min=0.0, zq_max=20.0)
    u1, u2 = np.log(LYA_REST_ANGSTROM * (1 + np.array(Z_EDGES)))
    a, b = np.log(REST_MIN), np.log(REST_MAX)
    y_nodes = np.log1p(nodes.zq_nodes)
    y_weights = nodes.zq_weights / (1 + nodes.zq_nodes)

    # Pixel quadrature reproduces the analytic overlap at every source node.
    np.testing.assert_allclose(
        overlap_length(nodes), analytic_overlap(y_nodes, u1, u2, a, b), rtol=1e-13
    )

    # Exact identity: the window is not truncated by the zq limits.
    total = np.sum(y_weights * overlap_length(nodes))
    assert total == pytest.approx((b - a) * (u2 - u1), rel=1e-13, abs=0)
    assert nodes.info["window_y"] == pytest.approx([u1 - b, u2 - a], rel=1e-14)


def test_breakpoints_sorted_unique_and_contain_kinks():
    nodes = make_nodes(zq_min=0.0, zq_max=20.0)
    breaks = np.array(nodes.info["y_breakpoints"])
    u1, u2 = np.log(LYA_REST_ANGSTROM * (1 + np.array(Z_EDGES)))
    a, b = np.log(REST_MIN), np.log(REST_MAX)
    assert np.all(np.diff(breaks) > 0)
    assert breaks[0] == pytest.approx(u1 - b, rel=1e-14)
    assert breaks[-1] == pytest.approx(u2 - a, rel=1e-14)
    for expected in [u1 - a, u2 - b, *np.log1p([2.5, 2.8, 3.0, 3.3])]:
        inside = (expected > breaks[0]) & (expected < breaks[-1])
        assert (not inside) or np.min(np.abs(breaks - expected)) < 1e-13
    assert nodes.info["n_y"] == nodes.info["n_y_panels"] * 4


def test_near_duplicate_breakpoints_merge():
    reference = make_nodes(zq_breaks=[2.8])
    duplicated = make_nodes(zq_breaks=[2.8, 2.8 * (1 + 1e-15), 2.8])
    assert duplicated.info["y_breakpoints"] == reference.info["y_breakpoints"]

    # Breakpoints at or beyond the window ends never move the ends.
    edge = make_nodes(zq_breaks=[reference.info["window_zq"][1]])
    assert edge.info["window_y"] == reference.info["window_y"]
    assert edge.info["y_breakpoints"][-1] == reference.info["window_y"][1]


def test_truncation_by_zq_limits():
    full = make_nodes(zq_min=0.0, zq_max=20.0)
    truncated = make_nodes(zq_min=2.5, zq_max=2.8)
    y_low, y_high = np.log1p([2.5, 2.8])
    assert truncated.info["window_y"] == pytest.approx([y_low, y_high], rel=1e-14)
    assert truncated.info["window_zq"][1] < full.info["window_zq"][1]
    assert np.all(truncated.zq_nodes > 2.5) and np.all(truncated.zq_nodes < 2.8)

    # Integral of the overlap over the truncated range: piecewise-linear exact.
    u1, u2 = np.log(LYA_REST_ANGSTROM * (1 + np.array(Z_EDGES)))
    a, b = np.log(REST_MIN), np.log(REST_MAX)
    y_fine = np.linspace(y_low, y_high, 400001)
    reference = np.trapezoid(analytic_overlap(y_fine, u1, u2, a, b), y_fine)
    y_weights = truncated.zq_weights / (1 + truncated.zq_nodes)
    assert np.sum(y_weights * overlap_length(truncated)) == pytest.approx(
        reference, rel=1e-9
    )


def test_pixels_inside_slice_and_forest():
    nodes = make_nodes()
    lambda_min = LYA_REST_ANGSTROM * (1 + Z_EDGES[0])
    lambda_max = LYA_REST_ANGSTROM * (1 + Z_EDGES[1])
    rest = nodes.lam_obs / (1 + nodes.z_q)
    tolerance = 1e-12
    assert np.all(nodes.lam_obs >= lambda_min * (1 - tolerance))
    assert np.all(nodes.lam_obs <= lambda_max * (1 + tolerance))
    assert np.all(rest >= REST_MIN * (1 - tolerance))
    assert np.all(rest <= REST_MAX * (1 + tolerance))
    np.testing.assert_array_equal(nodes.z_q, nodes.zq_nodes[nodes.y_index])
    np.testing.assert_allclose(
        nodes.z_pix, nodes.lam_obs / LYA_REST_ANGSTROM - 1, rtol=1e-15
    )
    assert nodes.y_index.dtype == np.int64
    assert nodes.info["n_pixel"] == len(nodes.y_index) == 2 * 3 * nodes.info["n_y"]


def test_geom_and_velocity_weights():
    nodes = make_nodes()
    lambda_limits = LYA_REST_ANGSTROM * (1 + np.array(Z_EDGES))
    expected_length = SPEED_LIGHT_KMS * np.log(lambda_limits[1] / lambda_limits[0])
    # ln(lambda2) - ln(lambda1) cancels to ~1e-12 relative precision.
    assert nodes.bin_length_velocity == pytest.approx(expected_length, rel=1e-11)
    expected_geom = (
        nodes.zq_weights[nodes.y_index]
        * nodes.velocity_weights
        / nodes.bin_length_velocity
    )
    np.testing.assert_array_equal(nodes.geom, expected_geom)
    assert not nodes.geom.flags.writeable

    # Measure normalization: with unit density, the geom sum is the overlap
    # integral over the bin length, in redshift units.
    total = nodes.geom.sum()
    u1, u2 = np.log(lambda_limits)
    y_weights = nodes.zq_weights / (1 + nodes.zq_nodes)
    expected = np.sum(
        (1 + nodes.zq_nodes)
        * y_weights
        * analytic_overlap(
            np.log1p(nodes.zq_nodes), u1, u2, np.log(REST_MIN), np.log(REST_MAX)
        )
        / (u2 - u1)
    )
    assert total == pytest.approx(expected, rel=1e-13)


def test_pixel_width_angstrom_scales_with_wavelength():
    nodes = make_nodes()
    width = pixel_width_angstrom(nodes, 100.0)
    np.testing.assert_allclose(
        width, 100.0 * nodes.lam_obs / SPEED_LIGHT_KMS, rtol=1e-15
    )
    with pytest.raises(ValueError):
        pixel_width_angstrom(nodes, 0.0)


@pytest.mark.parametrize(
    "overrides",
    [
        dict(zq_min=1.0, zq_max=1.2),
        dict(zq_min=10.0, zq_max=11.0),
        dict(rest_min=1040.0, rest_max=1041.0, zq_min=5.0, zq_max=6.0),
    ],
)
def test_empty_window_raises_with_label(overrides):
    with pytest.raises(ValueError, match=r"synthetic bin: no forest coverage"):
        make_nodes(**overrides)


def test_no_node_with_overlap_raises():
    # Window touches the slice at a single point: lower limit equals upper end.
    lambda_min = LYA_REST_ANGSTROM * (1 + Z_EDGES[0])
    zq_edge = lambda_min / REST_MAX - 1
    with pytest.raises(ValueError, match="no forest coverage"):
        make_nodes(zq_min=-0.5, zq_max=zq_edge)


@pytest.mark.parametrize(
    "overrides",
    [
        dict(zq_order=0),
        dict(lambda_order=1.5),
        dict(lambda_panels=0),
        dict(lambda_max=1.0),
        dict(rest_max=1000.0),
        dict(zq_max=1.0),
        dict(zq_breaks=[-2.0]),
        dict(lya_rest_angstrom=-1.0),
    ],
)
def test_invalid_arguments(overrides):
    with pytest.raises(ValueError):
        make_nodes(**overrides)


def make_source(base_nodes, density=None, **overrides):
    """Assemble valid source arrays on given nodes.

    Parameters
    ----------
    base_nodes : IntegrationNodes
        Quadrature nodes.
    density : ndarray or None, optional
        Density of shape (n_y, 3); default constant 2.
    **overrides : dict
        Replacements for the arguments of integrated_forest_source.

    Returns
    -------
    source : IntegratedForestSource
        Validated source.
    """
    n_y, n_p = len(base_nodes.zq_nodes), len(base_nodes.y_index)
    options = dict(
        nodes=base_nodes,
        density=np.full((n_y, 3), 2.0) if density is None else density,
        magnitudes=np.array([20.0, 21.0, 22.0]),
        quadrature=np.array([0.5, 1.0, 0.25]),
        variance=np.full((n_p, 3), 0.3),
        pixel_width_velocity=60.0,
        label="field a, bin 2",
    )
    options.update(overrides)
    return integrated_forest_source(**options)


def test_integrated_source_measure():
    nodes = make_nodes()
    density = np.random.default_rng(1).uniform(0.5, 3, (len(nodes.zq_nodes), 3))
    source = make_source(nodes, density=density)
    expected = (
        nodes.geom[:, None]
        * density[nodes.y_index]
        * np.array([0.5, 1.0, 0.25])[None, :]
    )
    assert source.measure.dtype == np.float64
    assert source.measure.shape == (len(nodes.y_index), 3)
    np.testing.assert_allclose(source.measure, expected, rtol=1e-15, atol=0)
    assert isinstance(source, IntegratedForestSource)
    assert not source.density.flags.writeable


def test_density_all_at_floor_raises_with_label():
    nodes = make_nodes()
    floor = np.full((len(nodes.zq_nodes), 3), LEGACY_DENSITY_FLOOR)
    with pytest.raises(ValueError, match=r"field a, bin 2: no forest coverage"):
        make_source(nodes, density=floor)
    with pytest.raises(ValueError, match="no forest coverage"):
        make_source(nodes, density=np.zeros_like(floor))

    # One populated node suffices; the floor matches the legacy adapter.
    floor[0, 0] = 1.0
    make_source(nodes, density=floor)
    assert LEGACY_DENSITY_FLOOR == DENSITY_FLOOR


@pytest.mark.parametrize(
    "case",
    [
        "density_shape",
        "variance_shape",
        "quadrature_shape",
        "density_negative",
        "variance_zero",
        "variance_nan",
        "density_inf",
        "width",
        "nodes_type",
    ],
)
def test_integrated_source_validation(case):
    nodes = make_nodes()
    n_y, n_p = len(nodes.zq_nodes), len(nodes.y_index)
    overrides = {}
    if case == "density_shape":
        overrides["density"] = np.ones((n_y + 1, 3))
    if case == "variance_shape":
        overrides["variance"] = np.ones((n_p, 2))
    if case == "quadrature_shape":
        overrides["quadrature"] = np.ones(2)
    if case == "density_negative":
        overrides["density"] = np.full((n_y, 3), -1.0)
    if case == "variance_zero":
        overrides["variance"] = np.zeros((n_p, 3))
    if case == "variance_nan":
        overrides["variance"] = np.full((n_p, 3), np.nan)
    if case == "density_inf":
        overrides["density"] = np.full((n_y, 3), np.inf)
    if case == "width":
        overrides["pixel_width_velocity"] = 0.0
    if case == "nodes_type":
        overrides["nodes"] = object()
    with pytest.raises(ValueError):
        make_source(nodes, **overrides)


def test_integration_nodes_is_frozen():
    nodes = make_nodes()
    assert isinstance(nodes, IntegrationNodes)
    with pytest.raises(AttributeError):
        nodes.geom = nodes.geom
