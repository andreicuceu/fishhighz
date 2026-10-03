"""Vectorised SNR/density queries against the scalar adapters, point by point."""

import numpy as np
import pytest

from fishhighz.adapters import legacy_inputs
from fishhighz.adapters.legacy_compat import LegacyDensity, LegacySNR
from fishhighz.adapters.legacy_inputs import DensityReader, SNRReader

EXPOSURE_COUNT = 4
EXPOSURE_TIME = 4000

MAGNITUDES = np.array([19.0, 20.0, 21.3, 22.9, 23.0, 23.4, 25.0])
Z_SOURCE = np.array([2.0, 2.37, 3.0, 3.99, 4.0, 1.5, 4.5, 2.8, 3.3, 2.9, 3.1])
WAVELENGTH = np.array(
    [
        3500.0,
        4177.0,
        5000.0,
        4100.0,
        3999.0,
        4000.0,
        4200.0,
        3400.0,
        5100.5,
        4600.0,
        4650.0,
    ]
)
PIXEL_WIDTH = np.array([0.8, 1.6, 0.4, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 1e-30, 0.8])


def snr_table(root):
    """Write four synthetic SNR tables with a linear SNR model.

    Parameters
    ----------
    root : pathlib.Path
        Directory for the generated files.

    Returns
    -------
    paths : list of pathlib.Path
        Table paths in unsorted magnitude order.
    """
    paths = []
    for name, magnitude in [("d", 23), ("a", 20), ("c", 22), ("b", 21)]:
        path = root / f"{name}.dat"
        source_redshifts = [2.0, 3.0, 4.0]
        wave = np.arange(3500, 5100, 100.0)
        header = (
            f"BAND= r MAG= {magnitude} EXPTIME= {EXPOSURE_TIME} NEXP= {EXPOSURE_COUNT}\n"
            "Wave " + " ".join(f"SN(z={v})" for v in source_redshifts)
        )
        np.savetxt(
            path,
            [
                [
                    w,
                    *[
                        3 + 0.1 * magnitude + 0.2 * v + 0.001 * w
                        for v in source_redshifts
                    ],
                ]
                for w in wave
            ],
            header=header,
        )
        paths.append(path)
    return paths


@pytest.fixture
def snr_reader(tmp_path):
    return SNRReader(snr_table(tmp_path), smoothing="none")


def scalar_legacy_variance(adapter, **options):
    """Evaluate the scalar legacy sample for every point and stack the variances.

    Parameters
    ----------
    adapter : LegacySNR
        Adapter under test.
    **options : dict
        Extra sample arguments (exposure_count, exposure_time).

    Returns
    -------
    variance : ndarray of shape (n_point, n_magnitude)
        Stacked scalar variances.
    counts : dict
        Summed fallback counts over the points.
    """
    rows, counts = [], dict(bright_clamp=0, out_of_range=0, snr_floor=0)
    for z, wavelength, width in zip(Z_SOURCE, WAVELENGTH, PIXEL_WIDTH):
        record = adapter.sample(
            z_source=z,
            magnitudes=MAGNITUDES,
            wavelength=wavelength,
            pixel_width_angstrom=width,
            **options,
        )
        rows.append(record["values"])
        for key in counts:
            counts[key] += record["provenance"]["counts"][key]
    return np.array(rows), counts


@pytest.mark.parametrize("exposure_count", [4, 7.5])
def test_legacy_snr_grid_matches_scalar(snr_reader, exposure_count):
    adapter = LegacySNR(snr_reader)
    expected, expected_counts = scalar_legacy_variance(
        adapter, exposure_count=exposure_count
    )
    record = adapter.variance_grid(
        z_source=Z_SOURCE,
        wavelength=WAVELENGTH,
        magnitudes=MAGNITUDES,
        pixel_width_angstrom=PIXEL_WIDTH,
        exposure_count=exposure_count,
    )
    np.testing.assert_allclose(record["values"], expected, rtol=1e-14, atol=0)
    assert record["values"].shape == (len(Z_SOURCE), len(MAGNITUDES))
    assert dict(record["provenance"]["counts"]) == expected_counts

    # Every fallback class is exercised: sentinel, bright clamp and SNR floor.
    assert expected_counts["bright_clamp"] > 0
    assert expected_counts["out_of_range"] > 0
    assert expected_counts["snr_floor"] > 0
    assert np.any(record["values"] == 1e20)


def test_legacy_snr_grid_chunking(snr_reader, monkeypatch):
    adapter = LegacySNR(snr_reader)
    options = dict(
        z_source=Z_SOURCE,
        wavelength=WAVELENGTH,
        magnitudes=MAGNITUDES,
        pixel_width_angstrom=PIXEL_WIDTH,
        exposure_count=EXPOSURE_COUNT,
    )
    reference = adapter.variance_grid(**options)
    monkeypatch.setattr(legacy_inputs, "GRID_CHUNK_POINTS", 5)
    chunked = adapter.variance_grid(**options)
    np.testing.assert_array_equal(chunked["values"], reference["values"])
    assert dict(chunked["provenance"]["counts"]) == dict(
        reference["provenance"]["counts"]
    )


def test_legacy_snr_grid_exposure_time_and_validation(snr_reader):
    adapter = LegacySNR(snr_reader)
    options = dict(
        z_source=Z_SOURCE,
        wavelength=WAVELENGTH,
        magnitudes=MAGNITUDES,
        pixel_width_angstrom=PIXEL_WIDTH,
        exposure_count=EXPOSURE_COUNT,
    )
    adapter.variance_grid(exposure_time=EXPOSURE_TIME, **options)
    with pytest.raises(ValueError, match="EXPTIME"):
        adapter.variance_grid(exposure_time=EXPOSURE_TIME + 1, **options)
    with pytest.raises(ValueError, match="nonnegative"):
        adapter.variance_grid(**{**options, "z_source": Z_SOURCE - 3})
    with pytest.raises(ValueError, match="equal length"):
        adapter.variance_grid(**{**options, "wavelength": WAVELENGTH[:-1]})
    with pytest.raises(ValueError):
        adapter.variance_grid(**{**options, "exposure_count": 0})
    with pytest.raises(ValueError):
        adapter.variance_grid(**{**options, "pixel_width_angstrom": PIXEL_WIDTH * 0})
    with pytest.raises(ValueError):
        adapter.variance_grid(**{**options, "magnitudes": []})


def test_legacy_snr_grid_all_outside(snr_reader):
    adapter = LegacySNR(snr_reader)
    record = adapter.variance_grid(
        z_source=[1.0, 9.0],
        wavelength=[4000.0, 4000.0],
        magnitudes=[21.0, 22.0],
        pixel_width_angstrom=[0.8, 0.8],
        exposure_count=EXPOSURE_COUNT,
    )
    np.testing.assert_array_equal(record["values"], np.full((2, 2), 1e20))
    assert record["provenance"]["counts"]["out_of_range"] == 4


def test_legacy_snr_grid_negative_raw_snr_raises(tmp_path):
    reader = SNRReader(snr_table(tmp_path), smoothing="none")
    adapter = LegacySNR(reader)

    # Force a negative interpolant, as the scalar path test would.
    class Negative:
        def __call__(self, points):
            return -np.ones(points.shape[:-1])

    object.__setattr__(adapter, "_interpolator", Negative())
    with pytest.raises(ValueError, match="negative interpolated SNR"):
        adapter.sample(
            z_source=3.0,
            magnitudes=[21.0],
            wavelength=4000.0,
            pixel_width_angstrom=0.8,
            exposure_count=EXPOSURE_COUNT,
        )
    with pytest.raises(ValueError, match="negative interpolated SNR"):
        adapter.variance_grid(
            z_source=[3.0],
            magnitudes=[21.0],
            wavelength=[4000.0],
            pixel_width_angstrom=[0.8],
            exposure_count=EXPOSURE_COUNT,
        )


def test_legacy_snr_grid_unrepresentable_variance_raises(snr_reader):
    adapter = LegacySNR(snr_reader)

    # Underflow of the scaled SNR is trapped, exactly as in the scalar path.
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


@pytest.mark.parametrize("exposure_count", [4, 6])
def test_strict_snr_grid_matches_scalar(snr_reader, exposure_count):
    inside = (
        (Z_SOURCE >= 2.0)
        & (Z_SOURCE <= 4.0)
        & (WAVELENGTH >= 3500.0)
        & (WAVELENGTH <= 5000.0)
    )
    magnitudes = np.array([20.0, 20.7, 21.3, 22.9, 23.0])
    z, wavelength, width = Z_SOURCE[inside], WAVELENGTH[inside], PIXEL_WIDTH[inside]
    width = np.where(width < 1e-3, 0.5, width)
    expected = np.array(
        [
            snr_reader.variance(
                z_source=zi,
                magnitudes=magnitudes,
                wavelength=wi,
                pixel_width_angstrom=pi,
                exposure_count=exposure_count,
            )
            for zi, wi, pi in zip(z, wavelength, width)
        ]
    )
    result = snr_reader.variance_grid(
        z_source=z,
        magnitudes=magnitudes,
        wavelength=wavelength,
        pixel_width_angstrom=width,
        exposure_count=exposure_count,
        exposure_time=EXPOSURE_TIME,
    )
    assert expected.shape == result.shape == (inside.sum(), len(magnitudes))
    np.testing.assert_allclose(result, expected, rtol=1e-14, atol=0)
    assert not result.flags.writeable


@pytest.mark.parametrize(
    "case", ["z_low", "z_high", "wavelength", "bright", "faint", "exposure_time"]
)
def test_strict_snr_grid_raises_out_of_domain(snr_reader, case):
    options = dict(
        z_source=[2.5, 3.5],
        magnitudes=[20.0, 22.0],
        wavelength=[4000.0, 4500.0],
        pixel_width_angstrom=[0.8, 0.8],
        exposure_count=EXPOSURE_COUNT,
    )
    snr_reader.variance_grid(**options)
    if case == "z_low":
        options["z_source"] = [1.99, 3.5]
    if case == "z_high":
        options["z_source"] = [2.5, 4.01]
    if case == "wavelength":
        options["wavelength"] = [4000.0, 5000.1]
    if case == "bright":
        options["magnitudes"] = [19.9, 22.0]
    if case == "faint":
        options["magnitudes"] = [20.0, 23.1]
    if case == "exposure_time":
        options["exposure_time"] = 1.0
    with pytest.raises(ValueError):
        snr_reader.variance_grid(**options)
    if case != "exposure_time":
        # The scalar path rejects the same point.
        with pytest.raises(ValueError):
            for z, wavelength in zip(options["z_source"], options["wavelength"]):
                snr_reader.variance(
                    z_source=z,
                    magnitudes=options["magnitudes"],
                    wavelength=wavelength,
                    pixel_width_angstrom=0.8,
                    exposure_count=EXPOSURE_COUNT,
                )


def test_strict_snr_grid_rejects_zero_snr(tmp_path):
    paths = snr_table(tmp_path)
    for path in paths:
        text = path.read_text().splitlines()
        rows = [line for line in text if not line.startswith("#")]
        values = rows[0].split()
        values[1:] = ["0.0"] * (len(values) - 1)
        rows[0] = " ".join(values)
        path.write_text("\n".join([*[ln for ln in text if ln.startswith("#")], *rows]))
    reader = SNRReader(paths, smoothing="none")
    with pytest.raises(ValueError, match="strictly positive"):
        reader.variance_grid(
            z_source=[3.0],
            magnitudes=[21.0],
            wavelength=[3500.0],
            pixel_width_angstrom=[0.8],
            exposure_count=EXPOSURE_COUNT,
        )


REDSHIFT_CENTRES = np.array([2.1, 2.3, 2.5, 2.7])
MAGNITUDE_CENTRES = np.array([20.05, 20.15, 20.25, 20.35, 20.45])


def density_table(path):
    """Write a rectangular count table with one empty cell.

    Parameters
    ----------
    path : pathlib.Path
        Output path.

    Returns
    -------
    path : pathlib.Path
        The written path.
    """
    counts = np.random.default_rng(5).uniform(0.5, 3.0, (4, 5))
    counts[1, 2] = 0.0
    rows = [
        [z, m, counts[i, j]]
        for i, z in enumerate(REDSHIFT_CENTRES)
        for j, m in enumerate(MAGNITUDE_CENTRES)
    ]
    np.savetxt(path, rows)
    return path


def density_reader(path, interpolation="piecewise_constant"):
    return DensityReader(
        path,
        semantics="cell_count_per_deg2",
        target_density=50.0,
        z_norm_min=2.2,
        interpolation=interpolation,
    )


# Interior edges, outer edges, cell interiors and exterior points.
GRID_Z = np.array([1.9, 2.0, 2.05, 2.2, 2.3, 2.45, 2.6, 2.7999999, 2.8, 2.9])
GRID_MAGNITUDES = np.array([19.9, 20.0, 20.1, 20.2, 20.33, 20.5, 20.6])


def test_legacy_density_grid_matches_scalar(tmp_path):
    adapter = LegacyDensity(density_reader(density_table(tmp_path / "d.txt")), "reject")
    expected = np.array([adapter.sample(z, GRID_MAGNITUDES)["values"] for z in GRID_Z])
    expected_raw = np.array([adapter.sample(z, GRID_MAGNITUDES)["raw"] for z in GRID_Z])
    record = adapter.sample_grid(GRID_Z, GRID_MAGNITUDES)
    np.testing.assert_allclose(record["values"], expected, rtol=1e-14, atol=0)
    np.testing.assert_allclose(record["raw"], expected_raw, rtol=1e-14, atol=0)
    assert record["values"].shape == (len(GRID_Z), len(GRID_MAGNITUDES))
    counts = record["provenance"]["counts"]
    assert counts["density_floor"] == np.count_nonzero(record["values"] == 1e-20)
    assert counts["density_floor"] > 0
    with pytest.raises(ValueError, match="nonnegative"):
        adapter.sample_grid([-0.1], GRID_MAGNITUDES)


def test_density_reader_query_grid_matches_scalar(tmp_path):
    reader = density_reader(density_table(tmp_path / "d.txt"))
    z = np.array([2.0, 2.2, 2.3, 2.45, 2.8])
    magnitudes = np.array([20.0, 20.1, 20.2, 20.33, 20.5])
    expected = np.array([reader.query(zi, magnitudes) for zi in z])
    result = reader.query_grid(z, magnitudes)
    np.testing.assert_allclose(result, expected, rtol=1e-14, atol=0)
    assert result.shape == (5, 5)
    assert not result.flags.writeable


@pytest.mark.parametrize(
    "z, magnitudes",
    [([1.99], [20.1]), ([2.81], [20.1]), ([2.2], [19.99]), ([2.2], [20.51])],
)
def test_density_reader_query_grid_raises_outside(tmp_path, z, magnitudes):
    reader = density_reader(density_table(tmp_path / "d.txt"))
    with pytest.raises(ValueError, match="outside closed domain"):
        reader.query_grid(z, magnitudes)
    with pytest.raises(ValueError):
        reader.query(z[0], magnitudes)


def test_density_grid_rejects_spline(tmp_path):
    reader = density_reader(density_table(tmp_path / "d.txt"), "spline")
    with pytest.raises(ValueError, match="piecewise_constant"):
        reader.query_grid([2.3], [20.2])
    adapter = LegacyDensity(reader, "floor_negative")
    with pytest.raises(ValueError, match="piecewise_constant"):
        adapter.sample_grid([2.3], [20.2])
