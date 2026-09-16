"""Independent physical-width, unequal-bin and dtype-boundary repair oracles."""

import warnings
from dataclasses import replace

import numpy as np
import pytest
from test_forecast import forest_spec, scalar_spec

from fishhighz.adapters.legacy_inputs import DensityReader
from fishhighz.forecast import prepare_bin, run_bin, run_forecast
from fishhighz.geometry import SPEED_LIGHT_KMS, prepare_geometry
from fishhighz.models.external import P3DProvider, PreparedP3D
from fishhighz.survey import freeze
from fishhighz.weights import prepare_forest_weights


def poly(z, m):
    return 2 + z**2 + 0.3 * (m - 20) ** 2 + 0.2 * z * (m - 20)


def irregular_fixture(path, shuffle=False):
    z = np.array([2.0, 2.3, 2.9, 3.7])
    m = np.array([20.0, 20.5, 21.0, 21.5])
    widths = np.array([0.2, 0.4, 0.7, 0.9])
    rows = np.array(
        [[a, b, poly(a, b) * w * 0.5] for a, w in zip(z, widths) for b in m]
    )
    if shuffle:
        np.random.default_rng(4).shuffle(rows)
    np.savetxt(path, rows)
    return z, m, widths, rows


def density(path, **kwargs):
    return DensityReader(
        path,
        semantics="cell_count_per_deg2",
        target_density=kwargs.pop("target_density", None),
        z_norm_min=kwargs.pop("z_norm_min", None),
        **kwargs,
    )


@pytest.mark.parametrize("shuffle", [False, True])
@pytest.mark.parametrize("target", [None, 100.0])
@pytest.mark.parametrize("masked", [False, True])
def test_explicit_density_widths(tmp_path, shuffle, target, masked):
    path = tmp_path / "counts"
    z, m, widths, rows = irregular_fixture(path, shuffle)
    bounds = (20.5, 21.5) if masked else None
    d = density(
        path,
        redshift_widths=widths,
        target_density=target,
        z_norm_min=z[1],
        magnitude_bounds=bounds,
    )
    measure = sum(
        row[2] for row in rows if row[0] > z[1] and (not masked or row[1] >= 20.5)
    )
    scale = 1 if target is None else target / measure
    expected = np.array(
        [[poly(a, b) * scale if not masked or b >= 20.5 else 0 for b in m] for a in z]
    )
    np.testing.assert_allclose(d.density, expected, rtol=5e-13, atol=0)
    total = sum(
        d.density[i, j] * widths[i] * 0.5 for i in range(2, 4) for j in range(4)
    )
    np.testing.assert_allclose(
        total, measure if target is None else target, rtol=5e-13, atol=0
    )
    np.testing.assert_allclose(
        d.provenance["selected_measure"], measure, rtol=5e-13, atol=0
    )
    if not masked:
        for a in [z[0], 2.65, z[-1]]:
            query = np.array([21.3, 20.0, 20.7, 21.5])
            np.testing.assert_allclose(
                d.query(a, query), poly(a, query) * scale, rtol=5e-12, atol=0
            )
    np.testing.assert_array_equal(d.provenance["redshift_axis"], z)
    np.testing.assert_array_equal(d.provenance["redshift_widths"], widths)
    assert d.provenance["width_policy"] == "explicit"
    snapshot = d.redshift_widths.copy()
    widths[:] = 50
    np.testing.assert_array_equal(d.redshift_widths, snapshot)
    with pytest.raises(ValueError):
        d.redshift_widths.flags.writeable = True


def test_legacy_width_oracle_and_uniform_preservation(tmp_path):
    path = tmp_path / "counts"
    z, m, widths, rows = irregular_fixture(path)
    legacy = density(path, width_policy="legacy_first_spacing")
    physical = density(path, redshift_widths=widths)
    expected = np.array(
        [[c / (z[1] - z[0]) / 0.5 for a, b, c in rows if a == zi] for zi in z]
    )
    np.testing.assert_allclose(legacy.density, expected, rtol=5e-13, atol=0)
    assert not np.allclose(legacy.density, physical.density)
    assert legacy.provenance["width_policy"] == "legacy_first_spacing"
    np.testing.assert_array_equal(legacy.redshift_widths, np.full(4, z[1] - z[0]))
    # The old uniform path must give exactly the same values as either policy.
    z = np.arange(2, 4, 0.5)
    np.savetxt(path, [[a, b, poly(a, b) * 0.5 * 0.5] for a in z for b in m])
    uniform = density(path)
    explicit = density(path, redshift_widths=np.full(4, 0.5))
    legacy = density(path, width_policy="legacy_first_spacing")
    np.testing.assert_array_equal(uniform.density, explicit.density)
    np.testing.assert_array_equal(uniform.density, legacy.density)
    assert uniform.provenance["width_policy"] == "uniform"


@pytest.mark.parametrize(
    "widths",
    [
        [0.1, 0.2],
        [[0.1, 0.2, 0.3, 0.4]],
        0.2,
        [0, 0.2, 0.3, 0.4],
        [-1, 0.2, 0.3, 0.4],
        [np.inf, 0.2, 0.3, 0.4],
        [np.nan, 0.2, 0.3, 0.4],
        np.ones(4, dtype=complex) * (1 + 1j),
        np.ones(4, dtype=bool),
        np.array(["1"] * 4),
        np.ones(4, dtype=object),
    ],
)
def test_invalid_widths(tmp_path, widths):
    path = tmp_path / "counts"
    irregular_fixture(path)
    with pytest.raises(ValueError):
        density(path, redshift_widths=widths)


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"width_policy": "uniform"},
        {"width_policy": "unknown"},
        {"redshift_widths": [0.1] * 4, "width_policy": "uniform"},
        {"redshift_widths": [0.1] * 4, "width_policy": "legacy_first_spacing"},
    ],
)
def test_missing_or_conflicting_width_policy(tmp_path, kwargs):
    path = tmp_path / "counts"
    irregular_fixture(path)
    with pytest.raises(ValueError):
        density(path, **kwargs)


@pytest.mark.parametrize(
    "field", ["magnitudes", "quadrature", "rho", "variance", "weights"]
)
@pytest.mark.parametrize(
    "kind,container",
    [
        (kind, container)
        for kind in ("complex", "bool", "string", "object", "nan", "inf")
        for container in ("array", "list")
        if (kind, container) != ("object", "list")
    ],
)
def test_scientific_dtype_rejection(field, kind, container):
    spec, _ = forest_spec(method="supplied")
    source = spec.forests["f"]
    options = dict(source.weight_options)
    original = np.asarray(options[field])
    if kind == "complex":
        value = original.astype(complex) + 1j
    elif kind == "bool":
        value = np.ones(2, dtype=bool)
    elif kind == "string":
        value = original.astype(str)
    elif kind == "object":
        value = original.astype(object)
    else:
        value = original.astype(float)
        value[0] = np.nan if kind == "nan" else np.inf
    if container == "list":
        value = value.tolist()
    options[field] = value
    args = (spec.p3d.selection.fields[0], spec.geometry, spec.responses["f"])
    with pytest.raises(ValueError):
        prepare_forest_weights(*args, **options)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # An imaginary-part-loss warning is a failure.
        with pytest.raises(ValueError):
            source = replace(source, weight_options=options)
            prepare_bin(replace(spec, forests={"f": source}))


@pytest.mark.parametrize("dtype", [np.int32, np.int64, np.float32, np.float64])
@pytest.mark.parametrize("container", ["array", "list"])
def test_valid_scientific_ownership(dtype, container):
    spec, _ = forest_spec(method="supplied")
    source = spec.forests["f"]
    options = dict(source.weight_options)
    values = dict(
        magnitudes=[20, 21],
        quadrature=[1, 2],
        rho=[0, 2],
        variance=[0, 0],
        weights=[0, 1],
    )
    arrays = {k: np.array(v, dtype=dtype) for k, v in values.items()}
    options.update(
        {k: v if container == "array" else v.tolist() for k, v in arrays.items()}
    )
    direct = prepare_forest_weights(
        spec.p3d.selection.fields[0], spec.geometry, spec.responses["f"], **options
    )
    source = replace(source, weight_options=options)
    for value in options.values():
        if isinstance(value, np.ndarray):
            value[:] = 999
    b = prepare_bin(replace(spec, forests={"f": source}))
    for name in values:
        a = getattr(b.weights["f"], name)
        np.testing.assert_array_equal(a, getattr(direct, name))
        assert a.dtype == np.float64
        with pytest.raises(ValueError):
            a.flags.writeable = True
    assert b.weights["f"].P_pixel == 0
    assert np.all(np.isfinite(run_bin(b).result.data_fisher))


@pytest.mark.parametrize(
    "array",
    [
        np.array([True, False]),
        np.array([1, 2], dtype=np.int32),
        np.array([1.5, 2], dtype=np.float32),
        np.array(["a", "bb"]),
        np.array([1 + 2j, 3j]),
    ],
)
def test_metadata_dtype_ownership(array):
    expected = array.copy()
    result = freeze({"data": array})["data"]
    assert result.dtype == expected.dtype
    array[:] = 0
    np.testing.assert_array_equal(result, expected)
    with pytest.raises(ValueError):
        result.flags.writeable = True


@pytest.mark.parametrize("batch", [None, 1, 5, 100])
@pytest.mark.parametrize("gap", [0, 0.15])
@pytest.mark.parametrize("reverse", [False, True])
def test_unequal_bins_independent_oracle(batch, gap, reverse):
    spec = scalar_spec()
    bins = []
    expected = []
    for i, (lo, hi, zeval, area, s) in enumerate(
        [(2, 2.7, 2.2, 10, 1), (2.7 + gap, 2.9 + gap, 2.85 + gap, 30, -1)]
    ):
        g = prepare_geometry(
            lo,
            hi,
            z_eval=zeval,
            area_deg2=area,
            h_fid=0.7,
            z_order=4,
            hubble=lambda z: np.full_like(z, 200),
            transverse_distance=lambda z: 1000 * (1 + z),
        )
        current = scalar_spec(spec.p3d.registry, sign=s, id=str(i))
        b = prepare_bin(replace(current, geometry=g))
        volume = (
            area
            * (np.pi / 180) ** 2
            * 0.7**3
            * SPEED_LIGHT_KMS
            / 200
            * 1e6
            * ((1 + hi) ** 3 - (1 + lo) ** 3)
            / 3
        )
        # Independent analytic k-shell integral, no production q_mode or volume.
        modes = volume * (0.2**3 - 0.02**3) / (6 * np.pi**2)
        np.testing.assert_allclose(g.volume, volume, rtol=5e-13, atol=0)
        expected.append(modes / 18 * np.array([[1, s], [s, 1]]))
        bins.append(b)
    if reverse:
        bins.reverse()
        expected.reverse()
    prior = np.diag([0.0, 4.0])
    result = run_forecast(bins, batch_size=batch, prior_fisher=prior)
    assert result.bin_ids == tuple(b.id for b in bins)
    for run, oracle in zip(result.bins, expected):
        np.testing.assert_allclose(run.result.data_fisher, oracle, rtol=5e-13, atol=0)
        assert run.result.diagnostics.rank == 1
    np.testing.assert_allclose(
        result.combined.data_fisher, sum(expected), rtol=5e-13, atol=0
    )
    np.testing.assert_array_equal(result.combined.prior_fisher, prior)
    assert result.combined.diagnostics.rank == 2


def test_explicit_evaluation_redshift_spies(monkeypatch):
    import fishhighz.forecast as forecast

    response_calls = []
    response = forecast.prepare_response

    def response_spy(*args, **kwargs):
        response_calls.append(kwargs["a_v"])
        return response(*args, **kwargs)

    monkeypatch.setattr(forecast, "prepare_response", response_spy)
    model_calls = []
    p1d_calls = []
    bins = []
    spec, _ = forest_spec(method="supplied")

    def model(t, z, k, mu, p):
        model_calls.append(z)
        return np.full((len(k), len(p)), t[0] * (1 + z))

    def p1d(t, z, k):
        p1d_calls.append(z)
        return np.full_like(k, 1 + z)

    owner = spec.p3d.routes[0].provider
    p3d = PreparedP3D(
        spec.p3d.registry,
        spec.p3d.selection,
        [P3DProvider(owner.label, model, owner.parameters, owner.pairs)],
    )
    for i, (lo, hi, z) in enumerate([(2, 2.7, 2.1), (2.8, 3, 2.97)]):
        g = prepare_geometry(
            lo,
            hi,
            z_eval=z,
            area_deg2=10,
            h_fid=0.7,
            z_order=4,
            hubble=lambda z: np.full_like(z, 200),
            transverse_distance=lambda z: 1000 * (1 + z),
        )
        source = replace(spec.forests["f"], p1d_model=p1d)
        b = prepare_bin(
            replace(spec, id=str(i), p3d=p3d, geometry=g, forests={"f": source})
        )
        bins.append(b)
        np.testing.assert_array_equal(b.power, 2 * (1 + z))
        np.testing.assert_allclose(
            response_calls[-1], 200 / ((1 + z) * 0.7), rtol=5e-13, atol=0
        )
    assert p1d_calls == [2.1, 2.97]
    model_calls.clear()
    run_forecast(bins, batch_size=5)
    assert model_calls == [2.1] * 15 + [2.97] * 15
    assert len(response_calls) == 2 and p1d_calls == [2.1, 2.97]
