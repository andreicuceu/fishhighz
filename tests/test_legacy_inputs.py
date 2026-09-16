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
    return 2 + z * z + 0.1 * (m - 20) ** 2 + 0.3 * z * (m - 20)


def density_file(path, *, shuffle=False):
    z = np.arange(2, 4, 0.5)
    m = np.arange(20, 24.0)
    rows = np.array([[a, b, polynomial(a, b) * 0.5] for a in z for b in m])
    if shuffle:
        np.random.default_rng(42).shuffle(rows)
    np.savetxt(path, rows)
    return rows


def reader(path, **kwargs):
    return DensityReader(
        path,
        semantics="cell_count_per_deg2",
        target_density=kwargs.pop("target_density", None),
        z_norm_min=kwargs.pop("z_norm_min", None),
        **kwargs,
    )


def snr_files(root, function=None):
    if function is None:

        def function(m, z, w):
            return 3 + 0.1 * m + 0.2 * z + 0.001 * w

    paths = []
    for name, m in [("z", 20), ("a", 21), ("b", 22), ("c", 23)]:
        path = root / f"{name}.dat"
        z = [2.0, 3.0, 4.0]
        wave = np.arange(3500, 5100, 100.0)
        header = f"BAND= r MAG= {m} EXPTIME= 4000 NEXP= 4\nWave " + " ".join(
            f"SN(z={v})" for v in z
        )
        np.savetxt(
            path, [[w, *[function(m, v, w) for v in z]] for w in wave], header=header
        )
        paths.append(path)
    return paths[::-1]


@pytest.mark.parametrize("shuffle", [False, True])
@pytest.mark.parametrize("target", [None, 100.0])
@pytest.mark.parametrize("threshold", [None, 2.5])
def test_density_count_oracle(tmp_path, shuffle, target, threshold):
    path = tmp_path / "density"
    rows = density_file(path, shuffle=shuffle)
    bounds = (21, 23)
    d = reader(
        path, target_density=target, z_norm_min=threshold, magnitude_bounds=bounds
    )
    measure = sum(
        c for z, m, c in rows if 21 <= m <= 23 and (threshold is None or z > threshold)
    )
    np.testing.assert_allclose(
        d.provenance["selected_measure"], measure, rtol=5e-13, atol=0
    )
    scale = 1 if target is None else target / measure
    expected = np.array(
        [
            [polynomial(z, m) * scale if 21 <= m <= 23 else 0 for m in d.magnitudes]
            for z in d.z
        ]
    )
    np.testing.assert_allclose(d.density, expected, rtol=5e-13, atol=0)
    np.testing.assert_allclose(
        d.provenance["original_measure"], sum(c for z, m, c in rows), rtol=5e-13, atol=0
    )


def test_density_offgrid_and_order(tmp_path):
    path = tmp_path / "d"
    density_file(path, shuffle=True)
    d = reader(path)
    m = np.array([22.3, 20, 21.1, 23])
    for z in [2, 2.73, 3.5]:
        np.testing.assert_allclose(d.query(z, m), polynomial(z, m), rtol=5e-12, atol=0)
    for a in [d.z, d.magnitudes, d.raw_counts, d.density, d.query(2, m)]:
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
            d = reader(path, **options)
            if case == "outside_z":
                d.query(1.999, [21])
            if case == "outside_m":
                d.query(2, [23.01])
            if case == "overshoot":
                d.query(2.7, [21.5])


def test_snr_affine_and_variance(tmp_path):
    snr = SNRReader(snr_files(tmp_path), smoothing="none")
    m = np.array([20, 22.3, 23, 21.1])
    expected = 3 + 0.1 * m + 0.2 * 2.8 + 0.001 * 4177
    np.testing.assert_allclose(
        snr.query(z_source=2.8, magnitudes=m, wavelength=4177),
        expected,
        rtol=5e-13,
        atol=0,
    )
    for width, nexp in [(0.8, 4), (1.6, 4), (0.8, 8)]:
        v = snr.variance(
            z_source=2.8,
            magnitudes=m,
            wavelength=4177,
            pixel_width_angstrom=width,
            exposure_count=nexp,
        )
        np.testing.assert_allclose(
            v, 1 / (expected**2 * width * nexp / 4), rtol=5e-13, atol=0
        )
    for z in [2, 4]:
        for w in [3500, 5000]:
            assert np.all(snr.query(z_source=z, magnitudes=[20, 23], wavelength=w) > 0)
    assert snr.provenance["paths"][0].endswith("z.dat")


def test_smoothing_reflection_oracle(tmp_path):
    # Independent half-sample symmetric reflection, not scipy or np.pad.
    paths = snr_files(tmp_path, lambda m, z, w: float(w == 3500))
    snr = SNRReader(paths, smoothing="legacy")
    offsets = np.arange(-40, 41)
    kernel = np.exp(-0.5 * (offsets / 10) ** 2)
    kernel /= sum(kernel)
    n = len(snr.wavelength)
    expected = []
    for i in range(n):
        value = 0.0
        for offset, weight in zip(offsets, kernel):
            j = (i + offset) % (2 * n)
            if j >= n:
                j = 2 * n - 1 - j
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
    paths = snr_files(tmp_path)
    p = paths[0]
    text = p.read_text()
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
    p.write_text("\n".join(lines) + "\n")
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
    spec, _ = forest_spec(auxiliary=False)
    g = spec.geometry
    path = tmp_path / "density"
    density_file(path)
    d = reader(path)
    snr = SNRReader(snr_files(tmp_path), smoothing="none")
    wavelength = LYA_REST_ANGSTROM * (1 + g.z_eval)
    response = InstrumentResponse(
        pixel_width_angstrom_to_velocity(0.8, lambda_obs_angstrom=wavelength), 10
    )
    m = np.array([20.0, 21.0, 22.0])
    sampled = sample_forest_readers(
        d,
        snr,
        g,
        response,
        z_source=3.2,
        magnitudes=m,
        pixel_width_angstrom=0.8,
        exposure_count=4,
    )
    rho = polynomial(3.2, m) * 4.2 / SPEED_LIGHT_KMS
    variance = 1 / ((3 + 0.1 * m + 0.2 * 3.2 + 0.001 * wavelength) ** 2 * 0.8)
    np.testing.assert_allclose(sampled["rho"], rho, rtol=5e-13, atol=0)
    field = spec.p3d.selection.fields[0]
    options = dict(
        z_source=3.2,
        magnitudes=m,
        quadrature=[0.5, 1, 0.5],
        length_velocity=10000,
        method="legacy",
        iterations=3,
        signal=0.5,
        alias=2,
    )
    a = prepare_forest_weights(
        field, g, response, **options, rho=sampled["rho"], variance=sampled["variance"]
    )
    b = prepare_forest_weights(
        field, g, response, **options, rho=rho, variance=variance
    )
    np.testing.assert_allclose(a.weights, b.weights, rtol=5e-13, atol=0)
    for w in (a, b):
        n = forest_noise(
            w, field, g, response, spec.grid.k_flat, spec.grid.mu_flat, np.ones(12)
        )
        if w is a:
            expected = n.total
        else:
            np.testing.assert_allclose(n.total, expected, rtol=5e-13, atol=0)
    expected_galaxy = (
        sum(polynomial(g.z_eval, m) * [0.5, 1, 0.5])
        * (1 + g.z_eval)
        / SPEED_LIGHT_KMS
        * g.a_v
        / g.d_deg**2
    )
    np.testing.assert_allclose(
        d.local_galaxy_density(g, m, [0.5, 1, 0.5]), expected_galaxy, rtol=5e-13, atol=0
    )
    with pytest.raises(ValueError, match="inconsistent"):
        sample_forest_readers(
            d,
            snr,
            g,
            InstrumentResponse(1, 10),
            z_source=3.2,
            magnitudes=m,
            pixel_width_angstrom=0.8,
            exposure_count=4,
        )


@pytest.mark.parametrize("target", [1e308, 1e-320])
def test_density_range_rejection(tmp_path, target):
    path = tmp_path / "d"
    density_file(path)
    rows = np.loadtxt(path)
    rows[:, 2] *= 1e-300 if target == 1e308 else 1e300
    np.savetxt(path, rows)
    with pytest.raises(ValueError, match="representable"):
        reader(path, target_density=target)


@pytest.mark.parametrize("value", [1e-300, 1e300])
def test_snr_variance_range_rejection(tmp_path, value):
    snr = SNRReader(snr_files(tmp_path, lambda *a: value), smoothing="none")
    with pytest.raises(ValueError, match="representable"):
        snr.variance(
            z_source=3,
            magnitudes=[21],
            wavelength=4000,
            pixel_width_angstrom=1,
            exposure_count=4,
        )
