"""Synthetic survey-table cells, SNR interpolation and normalization oracles."""

import numpy as np
import pytest
from scipy.interpolate import RegularGridInterpolator
from test_forecast import forest_spec

from fishhighz.geometry import LYA_REST_ANGSTROM, SPEED_LIGHT_KMS
from fishhighz.noise import forest_noise
from fishhighz.response import InstrumentResponse, pixel_width_angstrom_to_velocity
from fishhighz.survey_tables import (
    DensityTable,
    ForestSNR,
    SNRTable,
    SourceDensity,
    sample_forest_tables,
    snr_grid_separable,
)
from fishhighz.weights import prepare_forest_weights


def polynomial(z, m):
    """Evaluate the synthetic differential source density.

    Parameters
    ----------
    z : float or ndarray
        Dimensionless redshift.
    m : float or ndarray
        Apparent magnitude.

    Returns
    -------
    density : float or ndarray
        Source density per deg^2 per redshift per magnitude, with broadcast
        input shape.
    """
    return 2 + z * z + 0.1 * (m - 20) ** 2 + 0.3 * z * (m - 20)


def cell_polynomial(z, m):
    """Evaluate the synthetic density at the centre of the containing table cell.

    Parameters
    ----------
    z : float or ndarray
        Dimensionless redshift inside the cell edges [1.75, 3.75].
    m : float or ndarray
        Apparent magnitude inside the cell edges [19.5, 23.5].

    Returns
    -------
    density : float or ndarray
        Piecewise-constant cell density of the table written by density_file,
        per deg^2 per redshift per magnitude.

    Notes
    -----
    The cells of density_file have centres 2, 2.5, 3, 3.5 (width 0.5) in
    redshift and 20, 21, 22, 23 (width 1) in magnitude, closed below.
    """
    z_centre = 2 + 0.5 * np.floor((np.asarray(z) - 1.75) / 0.5)
    m_centre = 20 + np.floor(np.asarray(m) - 19.5)
    return polynomial(z_centre, m_centre)


def density_file(path, *, shuffle=False):
    """Write a regular synthetic redshift/magnitude count table.

    Parameters
    ----------
    path : pathlib.Path
        Path of the temporary test artifact to read or write.
    shuffle : bool, optional
        Whether to permute the generated table rows with a fixed random seed.
        Default is False.

    Returns
    -------
    rows : ndarray of shape (16, 3)
        Redshift, magnitude, and counts per deg^2 per table cell.
    """
    redshift_grid = np.arange(2, 4, 0.5)
    magnitude_grid = np.arange(20, 24.0)
    rows = np.array(
        [[a, b, polynomial(a, b) * 0.5] for a in redshift_grid for b in magnitude_grid]
    )
    if shuffle:
        np.random.default_rng(42).shuffle(rows)
    np.savetxt(path, rows)
    return rows


def reader(path, **kwargs):
    """Read the synthetic source-count table with explicit normalization options.

    Parameters
    ----------
    path : pathlib.Path
        Path of the temporary test artifact to read or write.
    **kwargs : dict
        DensityTable options, including target density per deg^2,
        normalization redshift and magnitude bounds.

    Returns
    -------
    reader : DensityTable
        Piecewise-constant source-density table for the temporary file.
    """
    return DensityTable(
        path,
        semantics="cell_count_per_deg2",
        target_density=kwargs.pop("target_density", None),
        z_norm_min=kwargs.pop("z_norm_min", None),
        **kwargs,
    )


def snr_files(root, function=None):
    """Write four synthetic SNR files in a deliberately unsorted order.

    Parameters
    ----------
    root : pathlib.Path
        Directory containing the synthetic input files.
    function : callable, optional
        Synthetic model or operation evaluated by the helper. Default is None.

    Returns
    -------
    paths : list of pathlib.Path
        Reversed file order for four magnitude samples.
    """
    if function is None:

        def function(m, z, w):
            """Evaluate the linear synthetic pixel-SNR relation.

            Parameters
            ----------
            m : float or ndarray
                Apparent magnitude.
            z : float or ndarray
                Dimensionless redshift.
            w : float or ndarray
                Observed wavelength in Angstrom.

            Returns
            -------
            snr : float or ndarray
                Dimensionless SNR at the supplied magnitude, source redshift, and
                wavelength.
            """
            return 3 + 0.1 * m + 0.2 * z + 0.001 * w

    paths = []
    for name, m in [("z", 20), ("a", 21), ("b", 22), ("c", 23)]:
        path = root / f"{name}.dat"
        source_redshifts = [2.0, 3.0, 4.0]
        wave = np.arange(3500, 5100, 100.0)
        header = f"BAND= r MAG= {m} EXPTIME= 4000 NEXP= 4\nWave " + " ".join(
            f"SN(z={v})" for v in source_redshifts
        )
        np.savetxt(
            path,
            [[w, *[function(m, v, w) for v in source_redshifts]] for w in wave],
            header=header,
        )
        paths.append(path)
    return paths[::-1]


@pytest.mark.parametrize("shuffle", [False, True])
@pytest.mark.parametrize("target", [None, 100.0])
@pytest.mark.parametrize("threshold", [None, 2.5])
def test_density_count_oracle(tmp_path, shuffle, target, threshold):
    """Check density count oracle.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    shuffle : bool
        Whether to permute input ordering, supplied by pytest parametrization.
    target : float or None
        Target quantity or object under examination, supplied by pytest
        parametrization.
    threshold : float or None
        Acceptance threshold, supplied by pytest parametrization.
    """
    path = tmp_path / "density"
    rows = density_file(path, shuffle=shuffle)
    bounds = (21, 23)
    density_reader = reader(
        path, target_density=target, z_norm_min=threshold, magnitude_bounds=bounds
    )
    measure = sum(
        c for z, m, c in rows if 21 <= m <= 23 and (threshold is None or z > threshold)
    )
    np.testing.assert_allclose(
        density_reader.provenance["selected_measure"], measure, rtol=5e-13, atol=0
    )
    scale = 1 if target is None else target / measure
    expected = np.array(
        [
            [
                polynomial(z, m) * scale if 21 <= m <= 23 else 0
                for m in density_reader.magnitudes
            ]
            for z in density_reader.z
        ]
    )
    np.testing.assert_allclose(density_reader.density, expected, rtol=5e-13, atol=0)
    np.testing.assert_allclose(
        density_reader.provenance["original_measure"],
        sum(c for z, m, c in rows),
        rtol=5e-13,
        atol=0,
    )


def test_density_offgrid_and_order(tmp_path):
    """Check density offgrid and order.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    path = tmp_path / "d"
    density_file(path, shuffle=True)
    density_reader = reader(path)
    magnitude_grid = np.array([22.3, 20, 21.1, 23])
    # Off-centre queries return the density of the containing cell.
    for z in [2, 2.73, 3.5]:
        np.testing.assert_allclose(
            density_reader.query(z, magnitude_grid),
            cell_polynomial(z, magnitude_grid),
            rtol=5e-12,
            atol=0,
        )
    for a in [
        density_reader.z,
        density_reader.magnitudes,
        density_reader.raw_counts,
        density_reader.density,
        density_reader.query(2, magnitude_grid),
    ]:
        with pytest.raises(ValueError):
            a.flags.writeable = True


@pytest.mark.parametrize(
    "case",
    [
        "missing",
        "duplicate",
        "spacing",
        "nan",
        "negative",
        "columns",
        "small",
        "zero_support",
        "semantics",
        "outside_z",
        "outside_m",
    ],
)
def test_density_errors(tmp_path, case):
    """Check density errors.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    case : str
        Named forecast or validation case, supplied by pytest parametrization.
    """
    path = tmp_path / "d"
    rows = density_file(path)
    options = {}
    if case == "missing":
        rows = rows[:-1]
    if case == "duplicate":
        rows = np.vstack((rows, rows[0]))
    if case == "spacing":
        rows[rows[:, 1] == 21, 1] = 21.1
    if case == "nan":
        rows[0, 2] = np.nan
    if case == "negative":
        rows[0, 2] = -1
    if case == "columns":
        rows = rows[:, :2]
    if case == "small":
        rows = rows[rows[:, 0] < 3]
    if case == "zero_support":
        options.update(target_density=10, z_norm_min=4)
    np.savetxt(path, rows)
    with pytest.raises(ValueError):
        if case == "semantics":
            DensityTable(path, semantics="dndzdm", target_density=None, z_norm_min=None)
        else:
            density_reader = reader(path, **options)
            # Closed cell-edge domain: redshift [1.75, 3.75], magnitude [19.5, 23.5].
            if case == "outside_z":
                density_reader.query(1.749, [21])
            if case == "outside_m":
                density_reader.query(2, [23.51])


def test_snr_affine_and_variance(tmp_path):
    """Check snr affine and variance.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    snr = SNRTable(snr_files(tmp_path), smoothing="none")
    magnitude_grid = np.array([20, 22.3, 23, 21.1])
    expected = 3 + 0.1 * magnitude_grid + 0.2 * 2.8 + 0.001 * 4177
    np.testing.assert_allclose(
        snr.query(z_source=2.8, magnitudes=magnitude_grid, wavelength=4177),
        expected,
        rtol=5e-13,
        atol=0,
    )
    for width, nexp in [(0.8, 4), (1.6, 4), (0.8, 8)]:
        pixel_variance = snr.variance(
            z_source=2.8,
            magnitudes=magnitude_grid,
            wavelength=4177,
            pixel_width_angstrom=width,
            exposure_count=nexp,
        )
        np.testing.assert_allclose(
            pixel_variance, 1 / (expected**2 * width * nexp / 4), rtol=5e-13, atol=0
        )
    for z in [2, 4]:
        for w in [3500, 5000]:
            assert np.all(snr.query(z_source=z, magnitudes=[20, 23], wavelength=w) > 0)
    assert snr.provenance["paths"][0].endswith("z.dat")


def test_smoothing_reflection_oracle(tmp_path):
    # Independent half-sample symmetric reflection, not scipy or np.pad.
    """Check smoothing reflection oracle.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    paths = snr_files(tmp_path, lambda m, z, w: float(w == 3500))
    snr = SNRTable(paths, smoothing="legacy")
    offsets = np.arange(-40, 41)
    kernel = np.exp(-0.5 * (offsets / 10) ** 2)
    kernel /= sum(kernel)
    wavelength_count = len(snr.wavelength)
    expected = []
    for i in range(wavelength_count):
        value = 0.0
        for offset, weight in zip(offsets, kernel):
            j = (i + offset) % (2 * wavelength_count)
            if j >= wavelength_count:
                j = 2 * wavelength_count - 1 - j
            if j == 0:
                value += weight
        expected.append(value)
    np.testing.assert_allclose(
        snr.smoothed_snr,
        np.broadcast_to(expected, snr.smoothed_snr.shape),
        rtol=5e-13,
        atol=0,
    )
    for p in paths:
        p.unlink()  # Only generated pytest temporary fixtures.
    snr = SNRTable(snr_files(tmp_path, lambda *a: 3), smoothing="legacy")
    np.testing.assert_allclose(snr.smoothed_snr, 3, rtol=5e-13, atol=0)


@pytest.mark.parametrize(
    "case",
    [
        "missing_header",
        "bad_wave",
        "duplicate_mag",
        "exptime",
        "nexp",
        "band",
        "grid",
        "negative",
        "columns",
        "duplicate_wave",
        "unordered_z",
    ],
)
def test_snr_input_errors(tmp_path, case):
    """Check snr input errors.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    case : str
        Named forecast or validation case, supplied by pytest parametrization.
    """
    paths = snr_files(tmp_path)
    snr_path = paths[0]
    text = snr_path.read_text()
    if case == "missing_header":
        text = text.replace("BAND=", "OTHER=")
    if case == "bad_wave":
        text = text.replace("SN(z=2.0)", "SN(oops)")
    if case == "duplicate_mag":
        text = text.replace("MAG= 23", "MAG= 22")
    if case == "exptime":
        text = text.replace("EXPTIME= 4000", "EXPTIME= 3000")
    if case == "nexp":
        text = text.replace("NEXP= 4", "NEXP= 0")
    if case == "band":
        text = text.replace("BAND= r", "BAND= g")
    if case == "unordered_z":
        text = text.replace("SN(z=2.0) SN(z=3.0)", "SN(z=3.0) SN(z=2.0)")
    lines = text.splitlines()
    if case == "grid":
        lines = lines[:-1]
    if case == "negative":
        lines[2] = lines[2].split()[0] + " -1 1 1"
    if case == "columns":
        lines[2] += " 1"
    if case == "duplicate_wave":
        lines[3] = lines[2]
    snr_path.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError):
        SNRTable(paths, smoothing="none")


@pytest.mark.parametrize(
    "options",
    [
        {"z_source": 1.9},
        {"z_source": 4.1},
        {"magnitudes": [19.9]},
        {"magnitudes": [23.1]},
        {"wavelength": 3499},
        {"wavelength": 5001},
        {"pixel_width_angstrom": 0},
        {"exposure_count": 0},
        {"exposure_time": 3000},
    ],
)
def test_snr_query_errors(tmp_path, options):
    """Check snr query errors.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    options : dict
        Optional controls defining this case, supplied by pytest
        parametrization.
    """
    snr = SNRTable(snr_files(tmp_path), smoothing="none")
    kwargs = dict(
        z_source=3,
        magnitudes=[21],
        wavelength=4000,
        pixel_width_angstrom=0.8,
        exposure_count=4,
    )
    kwargs.update(options)
    with pytest.raises(ValueError):
        snr.variance(**kwargs)


def test_zero_snr_fails(tmp_path):
    """Check zero snr fails.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    snr = SNRTable(snr_files(tmp_path, lambda *a: 0), smoothing="none")
    with pytest.raises(ValueError, match="positive SNR"):
        snr.variance(
            z_source=3,
            magnitudes=[21],
            wavelength=4000,
            pixel_width_angstrom=1,
            exposure_count=4,
        )


def test_reader_weight_boundary(tmp_path):
    """Check reader weight boundary.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    spec, _ = forest_spec(auxiliary=False)
    bin_geometry = spec.geometry
    path = tmp_path / "density"
    density_file(path)
    density_reader = reader(path)
    snr = SNRTable(snr_files(tmp_path), smoothing="none")
    wavelength = LYA_REST_ANGSTROM * (1 + bin_geometry.z_eval)
    response = InstrumentResponse(
        pixel_width_angstrom_to_velocity(0.8, lambda_obs_angstrom=wavelength), 10
    )
    magnitude_grid = np.array([20.0, 21.0, 22.0])
    sampled = sample_forest_tables(
        density_reader,
        snr,
        bin_geometry,
        response,
        z_source=3.2,
        magnitudes=magnitude_grid,
        pixel_width_angstrom=0.8,
        exposure_count=4,
    )
    rho = cell_polynomial(3.2, magnitude_grid) * 4.2 / SPEED_LIGHT_KMS
    variance = 1 / (
        (3 + 0.1 * magnitude_grid + 0.2 * 3.2 + 0.001 * wavelength) ** 2 * 0.8
    )
    np.testing.assert_allclose(sampled["rho"], rho, rtol=5e-13, atol=0)
    field = spec.p3d.selection.fields[0]
    options = dict(
        z_source=3.2,
        magnitudes=magnitude_grid,
        quadrature=[0.5, 1, 0.5],
        length_velocity=10000,
        method="early_lyaforecast",
        signal=0.5,
        alias=2,
    )
    reader_weights = prepare_forest_weights(
        field,
        bin_geometry,
        response,
        **options,
        rho=sampled["rho"],
        variance=sampled["variance"],
    )
    direct_weights = prepare_forest_weights(
        field, bin_geometry, response, **options, rho=rho, variance=variance
    )
    np.testing.assert_allclose(
        reader_weights.weights, direct_weights.weights, rtol=5e-13, atol=0
    )
    for w in (reader_weights, direct_weights):
        noise = forest_noise(
            w,
            field,
            bin_geometry,
            response,
            spec.grid.k_flat,
            spec.grid.mu_flat,
            np.ones(12),
        )
        if w is reader_weights:
            expected = noise.total
        else:
            np.testing.assert_allclose(noise.total, expected, rtol=5e-13, atol=0)
    expected_galaxy = (
        sum(cell_polynomial(bin_geometry.z_eval, magnitude_grid) * [0.5, 1, 0.5])
        * (1 + bin_geometry.z_eval)
        / SPEED_LIGHT_KMS
        * bin_geometry.a_v
        / bin_geometry.d_deg**2
    )
    np.testing.assert_allclose(
        density_reader.local_galaxy_density(
            bin_geometry, magnitude_grid, [0.5, 1, 0.5]
        ),
        expected_galaxy,
        rtol=5e-13,
        atol=0,
    )
    with pytest.raises(ValueError, match="inconsistent"):
        sample_forest_tables(
            density_reader,
            snr,
            bin_geometry,
            InstrumentResponse(1, 10),
            z_source=3.2,
            magnitudes=magnitude_grid,
            pixel_width_angstrom=0.8,
            exposure_count=4,
        )


@pytest.mark.parametrize("target", [1e308, 1e-320])
def test_density_range_rejection(tmp_path, target):
    """Check density range rejection.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    target : float
        Target quantity or object under examination, supplied by pytest
        parametrization.
    """
    path = tmp_path / "d"
    density_file(path)
    rows = np.loadtxt(path)
    rows[:, 2] *= 1e-300 if target == 1e308 else 1e300
    np.savetxt(path, rows)
    with pytest.raises(ValueError, match="representable"):
        reader(path, target_density=target)


@pytest.mark.parametrize("value", [1e-300, 1e300])
def test_snr_variance_range_rejection(tmp_path, value):
    """Check snr variance range rejection.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    value : float
        Value at the tested validation boundary, supplied by pytest
        parametrization.
    """
    snr = SNRTable(snr_files(tmp_path, lambda *a: value), smoothing="none")
    with pytest.raises(ValueError, match="representable"):
        snr.variance(
            z_source=3,
            magnitudes=[21],
            wavelength=4000,
            pixel_width_angstrom=1,
            exposure_count=4,
        )


# --- ForestSNR / SourceDensity fallbacks, vectorised queries and validation ---------


@pytest.fixture
def constant_density(tmp_path):
    """Write a quadratic synthetic source-count table and read it.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory supplied by pytest.

    Returns
    -------
    density : DensityTable
        Cell-count table with unnormalized densities, three redshift and three
        magnitude cells.
    """
    path = tmp_path / "density.txt"
    np.savetxt(
        path,
        [[z, m, (z - 3) ** 2 + (m - 21) ** 2] for z in [2, 3, 4] for m in [20, 21, 22]],
    )
    return reader(path)


@pytest.fixture
def constant_snr(tmp_path):
    """Write three constant pixel-S/N tables (S/N = 2) and read them unsmoothed.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory supplied by pytest.

    Returns
    -------
    snr_table : SNRTable
        Table on magnitudes 20-22, redshifts 2-4 and wavelengths 4000-5000 A
        with EXPTIME=1000 s and NEXP=4.
    """
    paths = []
    for magnitude in [20, 21, 22]:
        path = tmp_path / f"snr{magnitude}.dat"
        np.savetxt(
            path,
            [[w, 2, 2, 2] for w in [4000, 4500, 5000]],
            header=(
                f"BAND=r MAG={magnitude} EXPTIME=1000 NEXP=4\n"
                "Wave SN(z=2) SN(z=3) SN(z=4)"
            ),
        )
        paths.append(path)
    return SNRTable(paths, smoothing="none")


@pytest.mark.parametrize("magnitude", [19.9, 20, 21, 22, 22.00001])
@pytest.mark.parametrize("z", [1.99999, 2, 3, 4, 4.00001])
@pytest.mark.parametrize("wavelength", [3999.999, 4000, 4500, 5000, 5000.001])
def test_snr_branches(constant_snr, magnitude, z, wavelength):
    """Check the bright clamp and the 1e20 sentinel at the table edges.

    Parameters
    ----------
    constant_snr : SNRTable
        Constant S/N = 2 table.
    magnitude, z, wavelength : float
        Query coordinates straddling the closed table domains.
    """
    record = ForestSNR(constant_snr).sample(
        z_source=z,
        magnitudes=[magnitude],
        wavelength=wavelength,
        pixel_width_angstrom=2,
        exposure_count=8,
    )

    # S/N = 2 sqrt(2) sqrt(8/4) = 4 inside the table, so variance = 1/16.
    outside = magnitude > 22 or z < 2 or z > 4 or wavelength < 4000 or wavelength > 5000
    assert record["values"][0] == pytest.approx(1e20 if outside else 1 / 16, rel=5e-13)
    counts = record["provenance"]["counts"]
    assert counts["out_of_range"] == int(outside)
    assert counts["bright_clamp"] == int(magnitude < 20 and not outside)
    if outside:
        assert record["values"][0] == 1e20


@pytest.mark.parametrize("pixel,count", [(1e-24, 4), (1, 4), (2, 8)])
def test_post_scale_floor(constant_snr, pixel, count):
    """Check that the 1e-10 S/N floor acts after the pixel and exposure scaling.

    Parameters
    ----------
    constant_snr : SNRTable
        Constant S/N = 2 table.
    pixel : float
        Pixel width in angstrom.
    count : int
        Exposure count.
    """
    record = ForestSNR(constant_snr).sample(
        z_source=3,
        magnitudes=[21],
        wavelength=4500,
        pixel_width_angstrom=pixel,
        exposure_count=count,
    )
    expected = 1 / max(2 * np.sqrt(pixel) * np.sqrt(count / 4), 1e-10) ** 2
    np.testing.assert_allclose(record["values"], [expected], rtol=5e-13, atol=0)


@pytest.mark.parametrize("bad", [True, 1j, "3", np.nan, np.inf])
def test_invalid_queries(constant_density, constant_snr, bad):
    """Check that non-real or nonfinite coordinates are rejected.

    Parameters
    ----------
    constant_density : DensityTable
        Synthetic density table.
    constant_snr : SNRTable
        Constant S/N table.
    bad : bool or complex or str or float
        Invalid redshift or magnitude.
    """
    with pytest.raises(ValueError):
        SourceDensity(constant_density).sample(bad, [21])
    with pytest.raises(ValueError):
        ForestSNR(constant_snr).sample(
            z_source=3,
            magnitudes=[bad],
            wavelength=4500,
            pixel_width_angstrom=1,
            exposure_count=4,
        )


@pytest.mark.parametrize(
    "key,value",
    [
        ("pixel_width_angstrom", 0),
        ("pixel_width_angstrom", 1j),
        ("exposure_count", -1),
        ("exposure_time", 2),
        ("wavelength", 0),
        ("z_source", -1),
    ],
)
def test_invalid_exposure(constant_snr, key, value):
    """Check that invalid pixel, exposure or coordinate settings are rejected.

    Parameters
    ----------
    constant_snr : SNRTable
        Constant S/N table (EXPTIME = 1000 s).
    key : str
        Keyword of ForestSNR.sample set to an invalid value.
    value : int or complex
        The invalid value.
    """
    kwargs = dict(
        z_source=3,
        magnitudes=[30],
        wavelength=4500,
        pixel_width_angstrom=1,
        exposure_count=4,
    )
    kwargs[key] = value
    with pytest.raises(ValueError):
        ForestSNR(constant_snr).sample(**kwargs)


GRID_EXPOSURE_COUNT = 4
GRID_EXPOSURE_TIME = 4000
SNR_MAGNITUDES = np.array([19.0, 20.0, 21.3, 22.9, 23.0, 23.4, 25.0])
SNR_Z_SOURCE = np.array([2.0, 2.37, 3.0, 3.99, 4.0, 1.5, 4.5, 2.8, 3.3, 2.9, 3.1])
SNR_WAVELENGTH = np.array(
    [3500.0, 4177.0, 5000.0, 4100.0, 3999.0, 4000.0]
    + [4200.0, 3400.0, 5100.5, 4600.0, 4650.0]
)
SNR_PIXEL_WIDTH = np.array([0.8, 1.6, 0.4, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 1e-30, 0.8])


@pytest.fixture
def linear_snr(tmp_path):
    """Read the linear synthetic S/N tables (magnitudes 20-23, EXPTIME = 4000 s).

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory supplied by pytest.

    Returns
    -------
    snr_table : SNRTable
        Unsmoothed table of the linear model of snr_files.
    """
    return SNRTable(snr_files(tmp_path), smoothing="none")


def scalar_variance(adapter, **options):
    """Evaluate ForestSNR.sample at every point and stack the variances.

    Parameters
    ----------
    adapter : ForestSNR
        Adapter under test.
    **options : dict
        Extra sample arguments (exposure_count, exposure_time).

    Returns
    -------
    variance : ndarray of shape (n_point, n_magnitude)
        Stacked scalar variances.
    counts : dict
        Fallback counts summed over the points.
    """
    rows, counts = [], dict(bright_clamp=0, out_of_range=0, snr_floor=0)
    for z, wavelength, width in zip(SNR_Z_SOURCE, SNR_WAVELENGTH, SNR_PIXEL_WIDTH):
        record = adapter.sample(
            z_source=z,
            magnitudes=SNR_MAGNITUDES,
            wavelength=wavelength,
            pixel_width_angstrom=width,
            **options,
        )
        rows.append(record["values"])
        for key in counts:
            counts[key] += record["provenance"]["counts"][key]
    return np.array(rows), counts


@pytest.mark.parametrize("exposure_count", [4, 7.5])
def test_snr_grid_matches_scalar(linear_snr, exposure_count):
    """Check that variance_grid equals the scalar samples point by point.

    Parameters
    ----------
    linear_snr : SNRTable
        Linear synthetic S/N table.
    exposure_count : float
        Exposure count of the query.
    """
    adapter = ForestSNR(linear_snr)
    expected, expected_counts = scalar_variance(adapter, exposure_count=exposure_count)
    record = adapter.variance_grid(
        z_source=SNR_Z_SOURCE,
        wavelength=SNR_WAVELENGTH,
        magnitudes=SNR_MAGNITUDES,
        pixel_width_angstrom=SNR_PIXEL_WIDTH,
        exposure_count=exposure_count,
    )
    np.testing.assert_allclose(record["values"], expected, rtol=1e-14, atol=0)
    assert record["values"].shape == (len(SNR_Z_SOURCE), len(SNR_MAGNITUDES))
    assert dict(record["provenance"]["counts"]) == expected_counts

    # Every fallback class is exercised: sentinel, bright clamp and S/N floor.
    assert expected_counts["bright_clamp"] > 0
    assert expected_counts["out_of_range"] > 0
    assert expected_counts["snr_floor"] > 0
    assert np.any(record["values"] == 1e20)


def test_snr_grid_all_outside(linear_snr):
    """Check that a grid entirely outside the table is the sentinel everywhere.

    Parameters
    ----------
    linear_snr : SNRTable
        Linear synthetic S/N table.
    """
    record = ForestSNR(linear_snr).variance_grid(
        z_source=[1.0, 9.0],
        wavelength=[4000.0, 4000.0],
        magnitudes=[21.0, 22.0],
        pixel_width_angstrom=[0.8, 0.8],
        exposure_count=GRID_EXPOSURE_COUNT,
    )
    np.testing.assert_array_equal(record["values"], np.full((2, 2), 1e20))
    assert record["provenance"]["counts"]["out_of_range"] == 4


def test_snr_grid_negative_raw_snr_raises(linear_snr):
    """Check that a negative interpolated S/N raises on the scalar and grid paths.

    Parameters
    ----------
    linear_snr : SNRTable
        Linear synthetic S/N table.
    """
    adapter = ForestSNR(linear_snr)

    class Negative:
        """Interpolator stand-in returning S/N = -1 at every point."""

        def __call__(self, points):
            """Return -1 for each coordinate triple of shape (..., 3)."""
            return -np.ones(points.shape[:-1])

    object.__setattr__(adapter, "_interpolator", Negative())
    with pytest.raises(ValueError, match="negative interpolated SNR"):
        adapter.sample(
            z_source=3.0,
            magnitudes=[21.0],
            wavelength=4000.0,
            pixel_width_angstrom=0.8,
            exposure_count=GRID_EXPOSURE_COUNT,
        )
    with pytest.raises(ValueError, match="negative interpolated SNR"):
        adapter.variance_grid(
            z_source=[3.0],
            magnitudes=[21.0],
            wavelength=[4000.0],
            pixel_width_angstrom=[0.8],
            exposure_count=GRID_EXPOSURE_COUNT,
        )


def test_snr_grid_unrepresentable_variance_raises(linear_snr):
    """Check that underflow of the scaled S/N is trapped on both paths.

    Parameters
    ----------
    linear_snr : SNRTable
        Linear synthetic S/N table.
    """
    adapter = ForestSNR(linear_snr)
    kwargs = dict(magnitudes=[21.0], exposure_count=1e-320, exposure_time=None)
    with pytest.raises(ValueError, match="not representable"):
        adapter.sample(
            z_source=3.0, wavelength=4000.0, pixel_width_angstrom=1e-320, **kwargs
        )
    with pytest.raises(ValueError, match="not representable"):
        adapter.variance_grid(
            z_source=[3.0],
            wavelength=[4000.0],
            pixel_width_angstrom=[1e-320],
            **kwargs,
        )


def test_source_density_grid_matches_scalar(tmp_path):
    """Check that SourceDensity.sample_grid equals the scalar samples per redshift.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory supplied by pytest.
    """
    redshift_centres = [2.1, 2.3, 2.5, 2.7]
    magnitude_centres = [20.05, 20.15, 20.25, 20.35, 20.45]
    counts = np.random.default_rng(5).uniform(0.5, 3.0, (4, 5))
    counts[1, 2] = 0.0
    path = tmp_path / "d.txt"
    np.savetxt(
        path,
        [
            [z, m, counts[i, j]]
            for i, z in enumerate(redshift_centres)
            for j, m in enumerate(magnitude_centres)
        ],
    )
    adapter = SourceDensity(reader(path, target_density=50.0, z_norm_min=2.2))

    # Interior edges, outer edges, cell interiors and exterior points.
    z_grid = np.array([1.9, 2.0, 2.05, 2.2, 2.3, 2.45, 2.6, 2.7999999, 2.8, 2.9])
    magnitude_grid = np.array([19.9, 20.0, 20.1, 20.2, 20.33, 20.5, 20.6])
    scalar = [adapter.sample(z, magnitude_grid) for z in z_grid]
    record = adapter.sample_grid(z_grid, magnitude_grid)
    np.testing.assert_allclose(
        record["values"], [row["values"] for row in scalar], rtol=1e-14, atol=0
    )
    np.testing.assert_allclose(
        record["raw"], [row["raw"] for row in scalar], rtol=1e-14, atol=0
    )
    assert record["values"].shape == (len(z_grid), len(magnitude_grid))
    floor_count = record["provenance"]["counts"]["density_floor"]
    assert floor_count == np.count_nonzero(record["values"] == 1e-20)
    assert floor_count > 0
    with pytest.raises(ValueError, match="nonnegative"):
        adapter.sample_grid([-0.1], magnitude_grid)


def random_snr_table(root, rng):
    """Read non-separable random S/N tables with legacy smoothing.

    Parameters
    ----------
    root : pathlib.Path
        Output directory.
    rng : numpy.random.Generator
        Random number generator.

    Returns
    -------
    snr_table : SNRTable
        Table on magnitudes 20-23, redshifts [2, 2.6, 3, 4], wavelengths
        3500-5000 A.
    """
    paths = []
    wavelength_grid = np.arange(3500.0, 5100.0, 100.0)
    redshifts = [2.0, 2.6, 3.0, 4.0]
    for name, magnitude in [("d", 23), ("a", 20), ("c", 22), ("b", 21)]:
        path = root / f"{name}.dat"
        header = f"BAND= r MAG= {magnitude} EXPTIME= 4000 NEXP= 4\nWave " + " ".join(
            f"SN(z={v})" for v in redshifts
        )
        np.savetxt(
            path,
            np.column_stack(
                [wavelength_grid, rng.uniform(0.01, 4, (len(wavelength_grid), 4))]
            ),
            header=header,
        )
        paths.append(path)
    return SNRTable(paths, smoothing="legacy")


def test_separable_equals_independent_interpolation(tmp_path):
    """Check snr_grid_separable against a 3-D RegularGridInterpolator to 1e-13.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory supplied by pytest.
    """
    table = random_snr_table(tmp_path, np.random.default_rng(6))
    rng = np.random.default_rng(7)
    n_point, n_magnitude = 400, 60
    z_source = rng.uniform(table.z[0], table.z[-1], n_point)
    wavelength = rng.uniform(table.wavelength[0], table.wavelength[-1], n_point)
    magnitudes = rng.uniform(table.magnitudes[0], table.magnitudes[-1], n_magnitude)

    # Table nodes and domain edges, where the interval choice is ambiguous.
    z_source[:4] = table.z[[0, -1, 1, 2]]
    wavelength[:4] = table.wavelength[[0, -1, 3, 5]]
    magnitudes[:4] = table.magnitudes[[0, -1, 1, 2]]

    # Oracle built here on the same smoothed tensor, evaluated element by element.
    oracle = RegularGridInterpolator(
        (table.magnitudes, table.z, table.wavelength),
        table.smoothed_snr,
        method="linear",
        bounds_error=True,
    )
    points = np.empty((n_point, n_magnitude, 3))
    points[..., 0] = magnitudes[None, :]
    points[..., 1] = z_source[:, None]
    points[..., 2] = wavelength[:, None]
    expected = oracle(points)

    result = snr_grid_separable(oracle, z_source, wavelength, magnitudes)
    np.testing.assert_allclose(result, expected, rtol=1e-13, atol=0)


def test_variance_grid_separable_equals_independent_oracle(tmp_path):
    """Check ForestSNR.variance_grid against an independent interpolation oracle.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory supplied by pytest.

    Notes
    -----
    The oracle clamps bright magnitudes, interpolates the smoothed tensor with
    a RegularGridInterpolator, applies the pixel/exposure scaling and the 1e-10
    floor, and sets the sentinel 1e20 outside the table.
    """
    table = random_snr_table(tmp_path, np.random.default_rng(6))
    rng = np.random.default_rng(9)
    n_point = 50
    z_source = rng.uniform(1.9, 4.1, n_point)
    wavelength = rng.uniform(3450.0, 5150.0, n_point)
    pixel_width = rng.uniform(0.2, 2.0, n_point)
    magnitudes = np.array([19.0, 20.0, 21.7, 22.9, 23.0, 24.0])
    exposure_count = 2

    record = ForestSNR(table).variance_grid(
        z_source=z_source,
        wavelength=wavelength,
        magnitudes=magnitudes,
        pixel_width_angstrom=pixel_width,
        exposure_count=exposure_count,
    )

    oracle = RegularGridInterpolator(
        (table.magnitudes, table.z, table.wavelength),
        table.smoothed_snr,
        method="linear",
        bounds_error=True,
    )
    expected = np.full((n_point, len(magnitudes)), 1e20)
    for i in range(n_point):
        point_inside = (
            table.z[0] <= z_source[i] <= table.z[-1]
            and table.wavelength[0] <= wavelength[i] <= table.wavelength[-1]
        )
        if not point_inside:
            continue
        for j, magnitude in enumerate(magnitudes):
            if magnitude > table.magnitudes[-1]:
                continue
            effective_magnitude = max(magnitude, table.magnitudes[0])
            raw_snr = oracle([[effective_magnitude, z_source[i], wavelength[i]]])[0]
            scaled_snr = raw_snr * np.sqrt(pixel_width[i])
            scaled_snr *= np.sqrt(exposure_count / table.provenance["exposure_count"])
            expected[i, j] = 1 / max(scaled_snr, 1e-10) ** 2

    np.testing.assert_allclose(record["values"], expected, rtol=1e-13, atol=0)

    # Sentinels and bright clamps are present in the fixture.
    assert np.any(expected == 1e20)
    assert record["provenance"]["counts"]["bright_clamp"] > 0
