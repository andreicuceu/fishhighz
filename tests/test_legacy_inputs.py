"""Synthetic raw-grid interpolation and independent normalization oracles."""

import numpy as np
import pytest
from test_forecast import forest_spec

from fishhighz.adapters.legacy_inputs import (
    DensityReader,
    SNRReader,
    sample_forest_readers,
)
from fishhighz.geometry import LYA_REST_ANGSTROM, SPEED_LIGHT_KMS
from fishhighz.noise import forest_noise
from fishhighz.response import InstrumentResponse, pixel_width_angstrom_to_velocity
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
        DensityReader options, including target density per deg^2, normalization redshift, magnitude bounds, and interpolation settings.

    Returns
    -------
    reader : DensityReader
        Source-density interpolator for the temporary table.
    """
    return DensityReader(
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
    for z in [2, 2.73, 3.5]:
        np.testing.assert_allclose(
            density_reader.query(z, magnitude_grid),
            polynomial(z, magnitude_grid),
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
        "overshoot",
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
    if case == "overshoot":
        rows[:, 2] = (rows[:, 1] == 20).astype(float)
    np.savetxt(path, rows)
    with pytest.raises(ValueError):
        if case == "semantics":
            DensityReader(
                path, semantics="dndzdm", target_density=None, z_norm_min=None
            )
        else:
            density_reader = reader(path, **options)
            if case == "outside_z":
                density_reader.query(1.999, [21])
            if case == "outside_m":
                density_reader.query(2, [23.01])
            if case == "overshoot":
                density_reader.query(2.7, [21.5])


def test_snr_affine_and_variance(tmp_path):
    """Check snr affine and variance.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    snr = SNRReader(snr_files(tmp_path), smoothing="none")
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
    snr = SNRReader(paths, smoothing="legacy")
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
    snr = SNRReader(snr_files(tmp_path, lambda *a: 3), smoothing="legacy")
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
        SNRReader(paths, smoothing="none")


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
    snr = SNRReader(snr_files(tmp_path), smoothing="none")
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
    snr = SNRReader(snr_files(tmp_path, lambda *a: 0), smoothing="none")
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
    snr = SNRReader(snr_files(tmp_path), smoothing="none")
    wavelength = LYA_REST_ANGSTROM * (1 + bin_geometry.z_eval)
    response = InstrumentResponse(
        pixel_width_angstrom_to_velocity(0.8, lambda_obs_angstrom=wavelength), 10
    )
    magnitude_grid = np.array([20.0, 21.0, 22.0])
    sampled = sample_forest_readers(
        density_reader,
        snr,
        bin_geometry,
        response,
        z_source=3.2,
        magnitudes=magnitude_grid,
        pixel_width_angstrom=0.8,
        exposure_count=4,
    )
    rho = polynomial(3.2, magnitude_grid) * 4.2 / SPEED_LIGHT_KMS
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
        method="legacy",
        iterations=3,
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
        sum(polynomial(bin_geometry.z_eval, magnitude_grid) * [0.5, 1, 0.5])
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
        sample_forest_readers(
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
    snr = SNRReader(snr_files(tmp_path, lambda *a: value), smoothing="none")
    with pytest.raises(ValueError, match="representable"):
        snr.variance(
            z_source=3,
            magnitudes=[21],
            wavelength=4000,
            pixel_width_angstrom=1,
            exposure_count=4,
        )
