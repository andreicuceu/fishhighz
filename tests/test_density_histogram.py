"""Piecewise-constant dN/dz/dm cells: values, integrals, domains and configuration."""

import numpy as np
import pytest
from test_survey_config import _write_modified_ini

from fishhighz.accuracy import (
    DENSITY_INTERPOLATION_POLICIES,
    NATIVE_REVISION,
    REVISION_DENSITY_INTERPOLATION,
)
from fishhighz.adapters.legacy_compat import LegacyDensity
from fishhighz.adapters.legacy_inputs import CellHistogram2D, DensityReader
from fishhighz.magnitude import breakpoints, composite

REDSHIFTS = np.array([2.1, 2.3, 2.5, 2.7])
MAGNITUDES = np.array([20.05, 20.15, 20.25, 20.35, 20.45])
DZ, DM = 0.2, 0.1


def counts():
    """Return deterministic nonnegative cell counts per deg^2.

    Returns
    -------
    counts : ndarray of shape (4, 5)
        Counts per cell, including an empty cell, in (z, magnitude) order.
    """
    values = np.random.default_rng(7).uniform(
        0.5, 3.0, (len(REDSHIFTS), len(MAGNITUDES))
    )
    values[1, 2] = 0.0
    return values


def write_table(path, cell_counts=None, redshifts=REDSHIFTS):
    """Write a rectangular (z, magnitude, counts) table in shuffled row order.

    Parameters
    ----------
    path : pathlib.Path
        Output table path.
    cell_counts : ndarray or None, optional
        Counts per cell, shape (n_z, n_magnitude); default counts().
    redshifts : ndarray, optional
        Redshift cell centres; default REDSHIFTS.

    Returns
    -------
    path : pathlib.Path
        The written table path.
    """
    cell_counts = counts() if cell_counts is None else cell_counts
    rows = np.array(
        [
            [z, m, cell_counts[i, j]]
            for i, z in enumerate(redshifts)
            for j, m in enumerate(MAGNITUDES)
        ]
    )
    np.random.default_rng(3).shuffle(rows)
    np.savetxt(path, rows)
    return path


def cell_density(reader):
    """Return counts divided by the reader's tabulated cell measures.

    Parameters
    ----------
    reader : DensityReader
        Prepared reader without renormalization.

    Returns
    -------
    density : ndarray of shape (4, 5)
        Expected cell densities in deg^-2 redshift^-1 mag^-1.
    """
    return counts() / (reader.redshift_widths[:, None] * reader.provenance["dm"])


def read(path, **options):
    """Read a synthetic table with piecewise-constant cells by default.

    Parameters
    ----------
    path : pathlib.Path
        Table path.
    **options : dict
        DensityReader options overriding the defaults.

    Returns
    -------
    reader : DensityReader
        Prepared reader.
    """
    settings = dict(
        semantics="cell_count_per_deg2",
        target_density=None,
        z_norm_min=None,
        width_policy="legacy_first_spacing",
    )
    settings.update(options)
    return DensityReader(path, **settings)


def test_default_reader_is_piecewise_constant_with_cell_edges(tmp_path):
    """Default cells reproduce counts/(dz dm) at arbitrary interior points."""
    reader = read(write_table(tmp_path / "density.txt"))
    assert reader.interpolation == "piecewise_constant"
    assert reader.provenance["interpolation"] == "piecewise_constant_cells"
    assert reader._spline is None
    np.testing.assert_allclose(
        reader.z_edges, np.r_[REDSHIFTS - DZ / 2, REDSHIFTS[-1] + DZ / 2], atol=1e-14
    )
    np.testing.assert_allclose(
        reader.magnitude_edges,
        np.r_[MAGNITUDES - DM / 2, MAGNITUDES[-1] + DM / 2],
        atol=1e-13,
    )
    expected = cell_density(reader)
    rng = np.random.default_rng(11)
    for i, z in enumerate(REDSHIFTS):
        offsets = rng.uniform(-0.49, 0.49, len(MAGNITUDES)) * DM
        z_query = z + rng.uniform(-0.49, 0.49) * DZ
        np.testing.assert_allclose(
            reader.query(z_query, MAGNITUDES + offsets), expected[i], rtol=1e-14
        )


def test_cells_are_closed_below_and_outer_domain_is_closed(tmp_path):
    """Interior edges select the upper cell; outer edges remain in the domain."""
    reader = read(write_table(tmp_path / "density.txt"))
    expected = cell_density(reader)
    edges = reader.magnitude_edges
    np.testing.assert_allclose(
        reader.query(reader.z_edges[1], edges[1:-1]), expected[1, 1:], rtol=1e-14
    )
    np.testing.assert_allclose(
        reader.query(reader.z_edges[-1], [edges[0], edges[-1]]),
        expected[-1, [0, -1]],
        rtol=1e-14,
    )


def test_quadrature_reproduces_tabulated_counts(tmp_path):
    """Magnitude and redshift integrals reproduce the normalized table counts."""
    target = 123.0
    reader = read(write_table(tmp_path / "density.txt"), target_density=target)
    density = LegacyDensity(reader, "floor_negative")
    lo, hi = density.reader.magnitude_edges[[0, -1]]
    scaled = counts() * target / counts().sum()
    row_integrals = []
    for z in REDSHIFTS + 0.37 * DZ:
        partition = breakpoints({"tracer": density}, {}, {"tracer": z}, lo, hi)
        np.testing.assert_allclose(partition, density.reader.magnitude_edges, atol=0)
        magnitudes, weights = composite(partition, 3)
        row_integrals.append(np.sum(density.sample(z, magnitudes)["values"] * weights))
    np.testing.assert_allclose(row_integrals, scaled.sum(axis=1) / DZ, rtol=1e-13)
    np.testing.assert_allclose(np.sum(row_integrals) * DZ, target, rtol=1e-13)


def test_partition_merges_cell_edges_with_limits_and_snr_nodes(tmp_path):
    """Partitions contain the integration limits, interior cell edges and SNR nodes."""

    class SNR:
        magnitudes = np.array([20.12, 20.33])

    density = LegacyDensity(read(write_table(tmp_path / "d.txt")), "floor_negative")
    partition = breakpoints({"d": density}, {"s": SNR()}, {"d": 2.3}, 20.07, 20.42)
    np.testing.assert_allclose(
        partition, [20.07, 20.1, 20.12, 20.2, 20.3, 20.33, 20.4, 20.42], atol=1e-13
    )


def test_legacy_density_floors_outside_cell_edges(tmp_path):
    """Exterior redshift and magnitude queries receive the 1e-20 lyaforecast floor."""
    density = LegacyDensity(read(write_table(tmp_path / "d.txt")), "floor_negative")
    expected = cell_density(density.reader)
    z_edges = density.reader.z_edges
    m_edges = density.reader.magnitude_edges
    inside = density.sample(2.3, [m_edges[0] - 0.01, 20.25, m_edges[-1] + 0.01])
    np.testing.assert_array_equal(inside["values"][[0, 2]], [1e-20, 1e-20])
    assert inside["values"][1] == 0.0  # empty in-domain cells are not floored
    assert inside["provenance"]["counts"]["density_floor"] == 2
    for z in (z_edges[0] - 1e-6, z_edges[-1] + 1e-6):
        record = density.sample(z, MAGNITUDES)
        np.testing.assert_array_equal(record["values"], 1e-20)
        assert record["provenance"]["counts"]["redshift_outside_cells"] == 5
    np.testing.assert_allclose(
        density.sample(z_edges[0], MAGNITUDES)["values"], expected[0], rtol=1e-14
    )


def test_strict_reader_rejects_queries_outside_cell_edges(tmp_path):
    """The strict reader keeps closed-domain errors on the cell-edge domain."""
    reader = read(write_table(tmp_path / "d.txt"))
    with pytest.raises(ValueError, match="redshift"):
        reader.query(reader.z_edges[-1] + 1e-9, [20.2])
    with pytest.raises(ValueError, match="magnitude"):
        reader.query(2.3, [reader.magnitude_edges[0] - 1e-9])


def test_spline_option_preserves_legacy_interpolant(tmp_path):
    """The legacy spline remains selectable and unchanged at tabulated nodes."""
    path = write_table(tmp_path / "d.txt")
    spline = read(path, interpolation="spline")
    assert spline._histogram is None
    assert spline.provenance["interpolation"] == "RectBivariateSpline kx=2 ky=2 s=0"
    np.testing.assert_allclose(
        spline.query(REDSHIFTS[2], MAGNITUDES), counts()[2] / (DZ * DM), rtol=1e-10
    )
    with pytest.raises(ValueError, match="outside"):
        spline.query(REDSHIFTS[-1] + DZ / 4, [20.25])
    with pytest.raises(ValueError, match="interpolation"):
        read(path, interpolation="nearest")


def test_explicit_and_rounded_widths_follow_lyaforecast_edges(tmp_path):
    """Cells span consecutive lower edges centre - width/2, as in lyaforecast.

    Notes
    -----
    Explicit nonuniform widths that make the cells contiguous are accepted. A
    tiling residual above the float64 roundoff tolerance (rounded centres under
    the first-spacing width, or explicit widths that leave gaps) would stretch
    the cells and change their integrated density, so it is rejected under
    piecewise-constant interpolation; the spline has no tiling requirement.
    """
    redshifts = np.array([2.1, 2.35, 2.6, 2.8])
    path = write_table(tmp_path / "d.txt", redshifts=redshifts)
    widths = np.array([0.2, 0.3, 0.2, 0.2])
    reader = read(path, redshift_widths=widths, width_policy=None)
    np.testing.assert_allclose(reader.z_edges, [2.0, 2.2, 2.5, 2.7, 2.9], atol=1e-14)
    assert reader.provenance["z_cell_tiling_residual"] < 1e-14
    np.testing.assert_allclose(
        reader.query(2.45, MAGNITUDES), cell_density(reader)[1], rtol=1e-14
    )

    # Rounded centres: the first spacing sets every width and leaves a residual
    # of 0.05 (gap or overlap between neighbouring cells).
    with pytest.raises(ValueError, match=r"redshift cells do not tile.*0\.05"):
        read(path, width_policy="legacy_first_spacing")

    # Contiguous in the first two cells but with a gap after the second.
    with pytest.raises(ValueError, match="contiguous explicit redshift_widths"):
        read(path, redshift_widths=np.array([0.2, 0.3, 0.1, 0.2]), width_policy=None)
    with pytest.raises(ValueError, match="strictly increasing"):
        read(path, redshift_widths=np.array([0.2, 0.9, 0.2, 0.2]), width_policy=None)

    # The spline path is unchanged and accepts the same irregular table.
    assert (
        read(path, width_policy="legacy_first_spacing", interpolation="spline").z_edges
        is None
    )


def test_tiling_tolerance_accepts_roundoff_and_rejects_percent_gaps(tmp_path):
    """Roundoff-sized residuals pass; a regular grid with one displaced node fails."""
    redshifts = REDSHIFTS.copy()
    roundoff = np.array([0.0, 1e-15, -1e-15, 1e-15])
    path = write_table(tmp_path / "d.txt", redshifts=redshifts + roundoff)
    reader = read(path, width_policy="legacy_first_spacing")
    assert 0 < reader.provenance["z_cell_tiling_residual"] < 1e-13

    displaced = redshifts.copy()
    displaced[-1] += 0.02 * DZ
    path = write_table(tmp_path / "e.txt", redshifts=displaced)
    with pytest.raises(ValueError, match="do not tile"):
        read(path, width_policy="legacy_first_spacing")


def test_histogram_matches_lyaforecast_convention():
    """Queries off the upper outer edge follow lyaforecast Histogram2DInterpolator.

    Notes
    -----
    On the exact upper outer edge FishHighz keeps the last cell (closed outer
    domain), whereas lyaforecast returns its fill value; that point is excluded.
    """
    values = counts()
    z_edges = np.r_[REDSHIFTS - DZ / 2, REDSHIFTS[-1] + DZ / 2]
    m_edges = np.r_[MAGNITUDES - DM / 2, MAGNITUDES[-1] + DM / 2]
    histogram = CellHistogram2D(z_edges, m_edges, values, 0.0)
    z = np.array([1.9, 2.0, 2.2, 2.45, 2.79])
    m = np.array([20.0, 20.1, 20.25, 20.49, 20.6])
    padded = np.pad(values, 1)
    reference = padded[
        np.searchsorted(z_edges, z, side="right"),
        np.searchsorted(m_edges, m, side="right"),
    ]
    np.testing.assert_array_equal(histogram(z, m), reference)


def test_ini_default_selects_piecewise_constant_cells(tmp_path):
    """Compact INIs default to cells; explicit spline restores the legacy pair."""
    assert REVISION_DENSITY_INTERPOLATION[NATIVE_REVISION] == "piecewise_constant_cells"
    default = _write_modified_ini(tmp_path, {})
    assert default.input_policies["density_interpolation"] == "piecewise_constant_cells"
    assert (
        default.input_policies["magnitude_partition"]
        == DENSITY_INTERPOLATION_POLICIES["piecewise_constant_cells"][1]
    )
    spline = _write_modified_ini(
        tmp_path,
        {"input policies": {"density_interpolation": "RectBivariateSpline_kx2_ky2_s0"}},
    )
    assert (
        spline.input_policies["magnitude_partition"]
        == "density_knots_support_snr_nodes_negative_roots"
    )


@pytest.mark.parametrize(
    "changes,match",
    [
        ({"density_interpolation": "nearest"}, "density_interpolation"),
        (
            {"magnitude_partition": "density_knots_support_snr_nodes_negative_roots"},
            "inconsistent",
        ),
        (
            {
                "density_interpolation": "RectBivariateSpline_kx2_ky2_s0",
                "magnitude_partition": "density_cell_edges_support_snr_nodes",
            },
            "inconsistent",
        ),
    ],
)
def test_ini_rejects_unknown_or_mismatched_density_policies(tmp_path, changes, match):
    """Unknown interpolations and inconsistent partition pairs are rejected."""
    with pytest.raises(ValueError, match=match):
        _write_modified_ini(tmp_path, {"input policies": changes})
