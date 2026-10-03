"""Optimised integrated-mode paths against their reference implementations.

The dedicated 'sum_historical' recurrence for flattened integrated sources
(kernels.integrated_weights, NumPy and Numba backends) must reproduce the
central kernels.full_sum_weights.solve and fixed_weights with identical update
counts and stopping decisions; the separable S/N grid interpolation must equal
the generic RegularGridInterpolator query; the cached measure must equal the
direct expression.
"""

from types import SimpleNamespace

import numpy as np
import pytest
from scipy.interpolate import RegularGridInterpolator
from test_forest_integration import make_nodes, make_source
from test_integrated_weights import integrated, synthetic_source

from fishhighz import _arrays
from fishhighz.adapters import legacy_inputs
from fishhighz.adapters.legacy_compat import LegacySNR
from fishhighz.adapters.legacy_inputs import (
    SNRReader,
    snr_grid_separable,
    snr_grid_values,
)
from fishhighz.kernels import integrated_weights as kernel
from fishhighz.kernels.full_sum_weights import fixed_weights, solve


@pytest.fixture(scope="module")
def numba():
    """Skip the test if Numba is not installed."""
    pytest.importorskip("numba")


@pytest.fixture(params=["numpy", "numba"])
def backend(request):
    """Streaming-pass backend name; Numba is skipped if not installed."""
    if request.param == "numba":
        pytest.importorskip("numba")
    return request.param


def random_source(rng, n_pixel=300, n_magnitude=40, kind="random"):
    """Generate flattened measure and variance arrays.

    Parameters
    ----------
    rng : numpy.random.Generator
        Random number generator.
    n_pixel, n_magnitude : int
        Number of pixel and magnitude nodes.
    kind : {'random', 'srd'}
        'random' draws independent log-normal elements; 'srd' mimics an SRD
        source: the measure is geom_p * density[y_p, j] * quadrature_j with a
        steeply rising faint-end density, and the variance 1 / (snr^2 pixel)
        rises by many decades towards faint magnitudes up to the 1e20 sentinel.

    Returns
    -------
    measure, variance : ndarray of shape (n_pixel * n_magnitude,)
        Flattened C-order measure in deg^-2 and dimensionless variance.
    """
    if kind == "random":
        measure = rng.lognormal(-6, 2, n_pixel * n_magnitude)
        variance = rng.lognormal(0, 2, n_pixel * n_magnitude)
        return measure, variance
    magnitude = np.linspace(19.5, 23.5, n_magnitude)
    geom = rng.uniform(0.5, 1.5, n_pixel) * 1e-3
    density = 10 ** (0.4 * (magnitude - 20))[None, :] * rng.uniform(
        0.8, 1.2, (n_pixel, 1)
    )
    quadrature = np.full(n_magnitude, 4.0 / n_magnitude)
    measure = geom[:, None] * density * quadrature[None, :]
    snr = 10 ** (-0.4 * (magnitude - 20))[None, :] * rng.uniform(3, 6, (n_pixel, 1))
    variance = 1 / np.maximum(snr * np.sqrt(0.8), 1e-10) ** 2
    variance[:, -3:] = 1e20
    return measure.ravel(), variance.ravel()


def reference_solution(measure, variance, pixel, signal, alias, **controls):
    """Run the central full-sum solver on the flattened-measure mapping.

    Parameters
    ----------
    measure, variance : ndarray of shape (n,)
        Flattened arrays.
    pixel, signal, alias : float
        Pixel width, S and B.
    **controls : dict
        Keyword controls of kernels.full_sum_weights.solve.

    Returns
    -------
    result : dict
        Convergence record of the reference solver.
    """
    inputs = SimpleNamespace(
        density=measure,
        quadrature=np.ones_like(measure),
        variance=variance,
        length=1.0,
        pixel=pixel,
        signal=signal,
        p1d=alias,
    )
    with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
        return solve(inputs, "sum_historical", **controls)


def assert_same_record(new, reference, *, metric_atol=1e-12):
    """Compare a dedicated-solver record with the reference record.

    Parameters
    ----------
    new, reference : dict
        Convergence records.
    metric_atol : float, default=1e-12
        Absolute tolerance of the (rounding-limited) change metrics.
    """
    assert set(new) == set(reference)
    for key in (
        "status",
        "reason",
        "updates",
        "state_updates",
        "candidate",
        "forward_residual",
    ):
        assert new[key] == reference[key], key
    np.testing.assert_allclose(
        new["coefficients"], reference["coefficients"], rtol=1e-12, atol=0
    )
    np.testing.assert_allclose(new["weights"], reference["weights"], rtol=1e-12, atol=0)
    for key in ("last_step", "confirmation"):
        if reference[key] is None:
            assert new[key] is None
            continue
        assert set(new[key]) == set(reference[key])
        for metric, value in reference[key].items():
            np.testing.assert_allclose(
                new[key][metric], value, rtol=1e-6, atol=metric_atol
            )


CONTROLS = [
    {},
    {"rtol": 1e-9},
    {"rtol": 1e-2, "min_updates": 1, "stable_steps": 1},
    {"rtol": 1e-6, "min_updates": 2, "stable_steps": 2, "max_updates": 40},
    {"max_updates": 5},
]


@pytest.mark.parametrize("kind", ["random", "srd"])
@pytest.mark.parametrize("controls", CONTROLS, ids=str)
def test_dedicated_recurrence_reproduces_reference(backend, kind, controls):
    """Same status, update count, candidate and record as the central solver."""
    rng = np.random.default_rng(11)
    measure, variance = random_source(rng, kind=kind)
    pixel, signal, alias = 0.8, 3.5, 12.0
    reference = reference_solution(measure, variance, pixel, signal, alias, **controls)
    new = kernel.solve_integrated(
        measure,
        variance,
        pixel=pixel,
        signal=signal,
        alias=alias,
        backend=backend,
        **controls,
    )
    assert reference["status"] in ("converged", "capped")
    assert_same_record(new, reference)
    assert new["weights"].shape == measure.shape


@pytest.mark.parametrize("rtol", [1e-2, 1e-4, 1e-6])
def test_slowly_converging_source_reproduces_reference(backend, rtol, monkeypatch):
    """Many updates, unproved candidates and a cap, with the reference decisions.

    A small intrinsic signal contracts the recurrence slowly (34, 70 and more
    than 96 updates for the three tolerances), so candidates stay alive for
    many steps and some are evaluated exactly because the bound is inconclusive.
    """
    rng = np.random.default_rng(3)
    measure = rng.lognormal(-6, 2, 800) * 230.5879116025487
    variance = rng.lognormal(0, 2, 800)
    options = dict(pixel=29.760654938393564, signal=1.948e-3, alias=0.0154648)
    controls = dict(rtol=rtol)

    undecided = []
    original = kernel._proved_within

    def counting(*args):
        proved = original(*args)
        undecided.append(not proved)
        return proved

    monkeypatch.setattr(kernel, "_proved_within", counting)
    new = kernel.solve_integrated(
        measure, variance, backend=backend, **options, **controls
    )
    reference = reference_solution(measure, variance, **options, **controls)
    assert_same_record(new, reference)
    assert new["status"] == ("capped" if rtol == 1e-6 else "converged")
    assert new["updates"] >= 34
    if new["status"] == "converged":
        assert any(undecided)


def test_candidate_bound_does_not_change_decisions(backend, monkeypatch):
    """Proof by the intermediate steps equals the exact evaluation of every metric."""
    rng = np.random.default_rng(5)
    measure, variance = random_source(rng, kind="srd")
    options = dict(pixel=0.8, signal=3.5, alias=12.0, backend=backend)
    for controls in ({}, {"rtol": 1e-9}, {"rtol": 1e-3, "min_updates": 1}):
        bounded = kernel.solve_integrated(measure, variance, **options, **controls)
        monkeypatch.setattr(kernel, "_proved_within", lambda *args: False)
        exact = kernel.solve_integrated(measure, variance, **options, **controls)
        monkeypatch.undo()
        for key in ("status", "updates", "candidate", "state_updates"):
            assert bounded[key] == exact[key]
        np.testing.assert_allclose(bounded["weights"], exact["weights"], rtol=1e-14)
        if exact["confirmation"] is not None:
            assert bounded["confirmation"] == exact["confirmation"]


@pytest.mark.parametrize("updates", [0, 1, 4, 9])
def test_fixed_count_weights_reproduce_reference(backend, updates):
    """The fixed-count path equals central fixed_weights."""
    measure, variance = random_source(np.random.default_rng(3), kind="srd")
    pixel, signal, alias = 0.8, 3.5, 12.0
    inputs = SimpleNamespace(
        density=measure,
        quadrature=np.ones_like(measure),
        variance=variance,
        length=1.0,
        pixel=pixel,
        signal=signal,
        p1d=alias,
    )
    with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
        expected = fixed_weights(inputs, "sum_historical", updates)
    new = kernel.fixed_weights_integrated(
        measure,
        variance,
        pixel=pixel,
        signal=signal,
        alias=alias,
        updates=updates,
        backend=backend,
    )
    np.testing.assert_allclose(new, expected, rtol=1e-12, atol=0)


def test_ineligible_and_invalid_controls():
    """Nonpositive references are ineligible; bad controls raise like the reference."""
    measure, variance = random_source(np.random.default_rng(0), 4, 3)
    record = kernel.solve_integrated(
        measure, variance, pixel=1.0, signal=-1.0, alias=2.0, backend="numpy"
    )
    assert record["status"] == "ineligible" and record["weights"] is None
    with pytest.raises(ValueError, match="positive integers"):
        kernel.solve_integrated(
            measure, variance, pixel=1.0, signal=1.0, alias=2.0, max_updates=0
        )
    with pytest.raises(ValueError, match="rtol"):
        kernel.solve_integrated(
            measure, variance, pixel=1.0, signal=1.0, alias=2.0, rtol=0.0
        )
    with pytest.raises(ValueError, match="nonnegative integer"):
        kernel.fixed_weights_integrated(
            measure, variance, pixel=1.0, signal=1.0, alias=2.0, updates=-1
        )


def test_amplitude_is_the_array_maximum():
    """max(c/(c+v)) is attained at the smallest variance, bit for bit."""
    rng = np.random.default_rng(2)
    for _ in range(20):
        variance = rng.lognormal(0, 3, 5000)
        coefficient = float(rng.lognormal(0, 3))
        weights = coefficient / (variance + coefficient)
        assert kernel._amplitude(variance.min(), coefficient) == weights.max()


@pytest.mark.usefixtures("numba")
def test_compiled_passes_equal_numpy_passes():
    """Each compiled streaming pass equals its NumPy counterpart."""
    from fishhighz.kernels import _compiled_integrated_weights as compiled

    measure, variance = random_source(np.random.default_rng(8), kind="srd")
    c_new, c_old, r_new, r_old = 7.0, 7.4, 1 / 0.97, 1 / 0.96
    np.testing.assert_allclose(
        compiled.moments_pass(measure, variance, c_new)[:3],
        kernel.numpy_moments_pass(measure, variance, c_new),
        rtol=1e-13,
    )
    compiled_step = compiled.step_pass(measure, variance, c_new, c_old, r_new, r_old)
    numpy_step = kernel.numpy_step_pass(measure, variance, c_new, c_old, r_new, r_old)
    assert compiled_step[-1] == 0
    np.testing.assert_allclose(compiled_step[:3], numpy_step[:3], rtol=1e-13)
    np.testing.assert_allclose(compiled_step[3:5], numpy_step[3:], rtol=1e-13)
    np.testing.assert_allclose(
        compiled.difference_pass(variance, c_new, c_old, r_new, r_old),
        kernel.numpy_difference_pass(variance, c_new, c_old, r_new, r_old),
        rtol=1e-13,
    )
    filled, expected = np.empty_like(variance), np.empty_like(variance)
    compiled.fill_pass(variance, c_new, filled)
    kernel.numpy_fill_pass(variance, c_new, expected)
    np.testing.assert_array_equal(filled, expected)


def test_backend_selection(monkeypatch):
    """The environment variable selects the backend; unknown names are rejected."""
    monkeypatch.setenv("FISHHIGHZ_INTEGRATED_BACKEND", "numpy")
    assert kernel.select_backend().name == "numpy"
    monkeypatch.setenv("FISHHIGHZ_INTEGRATED_BACKEND", "bogus")
    with pytest.raises(ValueError, match="numpy or numba"):
        kernel.select_backend()
    monkeypatch.delenv("FISHHIGHZ_INTEGRATED_BACKEND")
    assert kernel.select_backend("numpy").name == "numpy"
    # Compiled arithmetic cannot trap underflow, so NumPy is used if requested.
    with np.errstate(under="raise"):
        assert kernel.select_backend("numba").name == "numpy"


@pytest.mark.parametrize("backend", ["reference", "numpy", "numba"])
@pytest.mark.parametrize("iterations", [None, 3])
def test_prepare_backends_agree(backend, iterations):
    """Every backend of prepare_integrated_forest_weights gives the same weights."""
    if backend == "numba":
        pytest.importorskip("numba")
    rng = np.random.default_rng(4)
    measure = rng.lognormal(-1, 1, (7, 3))
    variance = rng.lognormal(0, 1, (7, 3))
    source = synthetic_source(measure, variance, np.linspace(2.2, 2.4, 7))
    reference = integrated(source, iterations=iterations, backend="reference")
    prepared = integrated(source, iterations=iterations, backend=backend)
    for name in ("A", "P_pixel", "N1", "N2", "N3", "z_eff"):
        np.testing.assert_allclose(
            getattr(prepared, name), getattr(reference, name), rtol=1e-12, atol=0
        )
    np.testing.assert_allclose(prepared.weights, reference.weights, rtol=1e-12, atol=0)
    assert prepared.convergence["status"] == reference.convergence["status"]
    assert prepared.convergence["updates"] == reference.convergence["updates"]
    assert not prepared.weights.flags.writeable
    assert prepared.weights.shape == measure.shape


def test_prepare_backend_validation(monkeypatch):
    """An unknown backend name is rejected, by argument and by environment."""
    source = synthetic_source(np.ones((2, 3)), np.ones((2, 3)), [2.2, 2.3])
    with pytest.raises(ValueError, match="reference, numpy or numba"):
        integrated(source, iterations=1, backend="bogus")
    monkeypatch.setenv("FISHHIGHZ_INTEGRATED_BACKEND", "bogus")
    with pytest.raises(ValueError, match="reference, numpy or numba"):
        integrated(source, iterations=1)


def test_measure_is_cached_readonly_and_identical():
    """The measure is evaluated once, read-only, bit-equal to the direct product."""
    nodes = make_nodes()
    density = np.random.default_rng(1).uniform(0.5, 3, (len(nodes.zq_nodes), 3))
    source = make_source(nodes, density=density)
    quadrature = np.array([0.5, 1.0, 0.25])
    expected = (nodes.geom[:, None] * density[nodes.y_index]) * quadrature[None, :]
    first = source.measure
    assert source.measure is first
    assert not first.flags.writeable
    np.testing.assert_array_equal(first, expected)


def test_immutable_float_array_shares_only_immutable_memory():
    """Arrays backed by immutable bytes are shared; others are copied."""
    owned = np.arange(6.0).reshape(2, 3)
    assert _arrays.immutable_float_array(owned, "x") is not owned
    frozen = np.frombuffer(owned.tobytes(), dtype=np.float64).reshape(2, 3)
    assert _arrays.immutable_float_array(frozen, "x") is frozen
    readonly_view = owned.copy()
    readonly_view.flags.writeable = False
    assert _arrays.immutable_float_array(readonly_view, "x") is not readonly_view
    bad = np.frombuffer(np.array([1.0, np.nan]).tobytes(), dtype=np.float64)
    with pytest.raises(ValueError, match="finite"):
        _arrays.immutable_float_array(bad, "x")


# --- separable S/N interpolation -------------------------------------------------


def random_snr_tables(root, rng):
    """Write four non-separable random SNR tables.

    Parameters
    ----------
    root : pathlib.Path
        Output directory.
    rng : numpy.random.Generator
        Random number generator.

    Returns
    -------
    paths : list of pathlib.Path
        Table paths.
    """
    paths = []
    wave = np.arange(3500.0, 5100.0, 100.0)
    redshifts = [2.0, 2.6, 3.0, 4.0]
    for name, magnitude in [("d", 23), ("a", 20), ("c", 22), ("b", 21)]:
        path = root / f"{name}.dat"
        header = f"BAND= r MAG= {magnitude} EXPTIME= 4000 NEXP= 4\nWave " + " ".join(
            f"SN(z={v})" for v in redshifts
        )
        np.savetxt(
            path,
            np.column_stack([wave, rng.uniform(0.01, 4, (len(wave), len(redshifts)))]),
            header=header,
        )
        paths.append(path)
    return paths


@pytest.fixture
def random_snr_reader(tmp_path):
    return SNRReader(
        random_snr_tables(tmp_path, np.random.default_rng(6)), smoothing="legacy"
    )


def test_separable_equals_generic_interpolation(random_snr_reader):
    """Separable bilinear x linear equals the 3-D RegularGridInterpolator to 1e-13."""
    reader = random_snr_reader
    rng = np.random.default_rng(7)
    n_point, n_magnitude = 400, 60
    z = rng.uniform(reader.z[0], reader.z[-1], n_point)
    wavelength = rng.uniform(reader.wavelength[0], reader.wavelength[-1], n_point)
    magnitudes = rng.uniform(reader.magnitudes[0], reader.magnitudes[-1], n_magnitude)
    # Table nodes and domain edges, where the interval choice is ambiguous.
    z[:4] = reader.z[[0, -1, 1, 2]]
    wavelength[:4] = reader.wavelength[[0, -1, 3, 5]]
    magnitudes[:4] = reader.magnitudes[[0, -1, 1, 2]]

    points = np.empty((n_point, n_magnitude, 3))
    points[..., 0] = magnitudes[None, :]
    points[..., 1] = z[:, None]
    points[..., 2] = wavelength[:, None]
    expected = reader._interpolator(points)
    result = snr_grid_separable(reader._interpolator, z, wavelength, magnitudes)
    np.testing.assert_allclose(result, expected, rtol=1e-13, atol=0)


def grid_options():
    """Return paired points and magnitudes with every legacy fallback present."""
    rng = np.random.default_rng(9)
    n_point = 50
    z = rng.uniform(1.9, 4.1, n_point)
    wavelength = rng.uniform(3450.0, 5150.0, n_point)
    return dict(
        z_source=z,
        wavelength=wavelength,
        magnitudes=[19.0, 20.0, 21.7, 22.9, 23.0, 24.0],
        pixel_width_angstrom=rng.uniform(0.2, 2.0, n_point),
        exposure_count=2,
    )


@pytest.mark.parametrize("legacy", [True, False])
def test_variance_grid_separable_equals_interpolator(
    legacy, random_snr_reader, monkeypatch
):
    """The grid query agrees between implementations, sentinels and clamps included."""
    reader = random_snr_reader
    options = grid_options()
    if legacy:
        adapter = LegacySNR(reader)
    else:
        # The strict reader rejects out-of-domain points.
        adapter = reader
        z = np.clip(options["z_source"], reader.z[0], reader.z[-1])
        wavelength = np.clip(
            options["wavelength"], reader.wavelength[0], reader.wavelength[-1]
        )
        options = dict(
            options,
            z_source=z,
            wavelength=wavelength,
            magnitudes=[20.0, 20.5, 21.7, 22.9, 23.0],
        )
    monkeypatch.setenv("FISHHIGHZ_SNR_GRID", "interpolator")
    reference = adapter.variance_grid(**options)
    monkeypatch.setenv("FISHHIGHZ_SNR_GRID", "separable")
    separable = adapter.variance_grid(**options)
    values = (
        (reference["values"], separable["values"]) if legacy else (reference, separable)
    )
    np.testing.assert_allclose(values[1], values[0], rtol=1e-13, atol=0)
    if legacy:
        assert reference["provenance"]["counts"] == separable["provenance"]["counts"]
        # Sentinels are present in the fixture.
        assert np.any(values[0] == 1e20)


def test_snr_grid_values_dispatch(random_snr_reader, monkeypatch):
    """Selection by environment; a non-SciPy interpolator is called directly."""
    reader = random_snr_reader
    z = np.array([2.5, 3.5])
    wavelength = np.array([4000.0, 4500.0])
    magnitudes = np.array([20.5, 22.0])
    monkeypatch.setenv("FISHHIGHZ_SNR_GRID", "bogus")
    with pytest.raises(ValueError, match="separable or interpolator"):
        snr_grid_values(reader._interpolator, z, wavelength, magnitudes)
    monkeypatch.delenv("FISHHIGHZ_SNR_GRID")

    class Double:
        """Interpolator test double returning a constant."""

        def __call__(self, points):
            return np.full(points.shape[:-1], 2.0)

    np.testing.assert_array_equal(
        snr_grid_values(Double(), z, wavelength, magnitudes), np.full((2, 2), 2.0)
    )
    assert isinstance(reader._interpolator, RegularGridInterpolator)
    # Two-node axes are the minimum for the separable path; one-node falls back.
    single = RegularGridInterpolator(
        (np.array([20.0, 21.0]), np.array([2.0]), np.array([4000.0, 5000.0])),
        np.ones((2, 1, 2)),
    )
    result = snr_grid_values(
        single, np.array([2.0]), np.array([4500.0]), np.array([20.5])
    )
    np.testing.assert_allclose(result, 1.0)
    assert legacy_inputs.snr_grid_values is snr_grid_values
