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
    """Evaluate the nonseparable synthetic source-density polynomial.

    Parameters
    ----------
    z : float or ndarray
        Dimensionless redshift.
    m : float or ndarray
        Apparent magnitude.

    Returns
    -------
    density : float or ndarray
        Source density per deg^2 per redshift per magnitude.
    """
    return 2 + z**2 + 0.3 * (m - 20) ** 2 + 0.2 * z * (m - 20)


def irregular_fixture(path, shuffle=False):
    """Write source counts on nonuniform redshift cells.

    Parameters
    ----------
    path : pathlib.Path
        Destination of the generated nonuniform source-count table.
    shuffle : bool, optional
        Whether to permute the generated table rows with a fixed random seed.
        Default is False.

    Returns
    -------
    fixture : tuple of ndarray
        Redshift centers (4,), magnitudes (4,), redshift-cell widths (4,), and
        table rows (16, 3).
    """
    redshift_grid = np.array([2.0, 2.3, 2.9, 3.7])
    magnitude_grid = np.array([20.0, 20.5, 21.0, 21.5])
    widths = np.array([0.2, 0.4, 0.7, 0.9])
    rows = np.array(
        [
            [a, b, poly(a, b) * w * 0.5]
            for a, w in zip(redshift_grid, widths)
            for b in magnitude_grid
        ]
    )
    if shuffle:
        np.random.default_rng(4).shuffle(rows)
    np.savetxt(path, rows)
    return redshift_grid, magnitude_grid, widths, rows


def density(path, **kwargs):
    """Read the synthetic count table with explicit cell-width options.

    Parameters
    ----------
    path : pathlib.Path
        Path of the temporary test artifact to read or write.
    **kwargs : dict
        DensityReader options, including explicit redshift cell widths, density normalization, and interpolation settings.

    Returns
    -------
    reader : DensityReader
        Density interpolator with the requested raw-cell semantics.
    """
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
    """Check explicit density widths.

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
    masked : bool
        Whether the fixture contains masked entries, supplied by pytest
        parametrization.
    """
    path = tmp_path / "counts"
    redshift_grid, magnitude_grid, widths, rows = irregular_fixture(path, shuffle)
    bounds = (20.5, 21.5) if masked else None
    density_reader = density(
        path,
        redshift_widths=widths,
        target_density=target,
        z_norm_min=redshift_grid[1],
        magnitude_bounds=bounds,
    )
    measure = sum(
        row[2]
        for row in rows
        if row[0] > redshift_grid[1] and (not masked or row[1] >= 20.5)
    )
    scale = 1 if target is None else target / measure
    expected = np.array(
        [
            [
                poly(a, b) * scale if not masked or b >= 20.5 else 0
                for b in magnitude_grid
            ]
            for a in redshift_grid
        ]
    )
    np.testing.assert_allclose(density_reader.density, expected, rtol=5e-13, atol=0)
    total = sum(
        density_reader.density[i, j] * widths[i] * 0.5
        for i in range(2, 4)
        for j in range(4)
    )
    np.testing.assert_allclose(
        total, measure if target is None else target, rtol=5e-13, atol=0
    )
    np.testing.assert_allclose(
        density_reader.provenance["selected_measure"], measure, rtol=5e-13, atol=0
    )
    if not masked:
        for a in [redshift_grid[0], 2.65, redshift_grid[-1]]:
            query = np.array([21.3, 20.0, 20.7, 21.5])
            np.testing.assert_allclose(
                density_reader.query(a, query),
                poly(a, query) * scale,
                rtol=5e-12,
                atol=0,
            )
    np.testing.assert_array_equal(
        density_reader.provenance["redshift_axis"], redshift_grid
    )
    np.testing.assert_array_equal(density_reader.provenance["redshift_widths"], widths)
    assert density_reader.provenance["width_policy"] == "explicit"
    snapshot = density_reader.redshift_widths.copy()
    widths[:] = 50
    np.testing.assert_array_equal(density_reader.redshift_widths, snapshot)
    with pytest.raises(ValueError):
        density_reader.redshift_widths.flags.writeable = True


def test_legacy_width_oracle_and_uniform_preservation(tmp_path):
    """Check legacy width oracle and uniform preservation.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    path = tmp_path / "counts"
    redshift_grid, magnitude_grid, widths, rows = irregular_fixture(path)
    legacy = density(path, width_policy="legacy_first_spacing")
    physical = density(path, redshift_widths=widths)
    expected = np.array(
        [
            [
                c / (redshift_grid[1] - redshift_grid[0]) / 0.5
                for a, b, c in rows
                if a == zi
            ]
            for zi in redshift_grid
        ]
    )
    np.testing.assert_allclose(legacy.density, expected, rtol=5e-13, atol=0)
    assert not np.allclose(legacy.density, physical.density)
    assert legacy.provenance["width_policy"] == "legacy_first_spacing"
    np.testing.assert_array_equal(
        legacy.redshift_widths, np.full(4, redshift_grid[1] - redshift_grid[0])
    )
    # The old uniform path must give exactly the same values as either policy.
    redshift_grid = np.arange(2, 4, 0.5)
    np.savetxt(
        path,
        [[a, b, poly(a, b) * 0.5 * 0.5] for a in redshift_grid for b in magnitude_grid],
    )
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
    """Check invalid widths.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    widths : float or ndarray or list
        Instrumental, broadening, or table-cell widths, supplied by pytest
        parametrization.
    """
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
    """Check missing or conflicting width policy.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    kwargs : dict
        Keyword arguments selecting the parametrized case, supplied by pytest
        parametrization.
    """
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
    """Check scientific dtype rejection.

    Parameters
    ----------
    field : str
        Forest-weight input array whose dtype is perturbed.
    kind : str
        Invalid numeric dtype or nonfinite-value case: complex, bool, string,
        object, nan, or inf.
    container : str
        Input representation, either array or list.
    """
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
    """Check valid scientific ownership.

    Parameters
    ----------
    dtype : type
        Array dtype under examination, supplied by pytest parametrization.
    container : str
        Input container constructor, supplied by pytest parametrization.
    """
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
    prepared_bin = prepare_bin(replace(spec, forests={"f": source}))
    for name in values:
        weight_array = getattr(prepared_bin.weights["f"], name)
        np.testing.assert_array_equal(weight_array, getattr(direct, name))
        assert weight_array.dtype == np.float64
        with pytest.raises(ValueError):
            weight_array.flags.writeable = True
    assert prepared_bin.weights["f"].P_pixel == 0
    assert np.all(np.isfinite(run_bin(prepared_bin).result.data_fisher))


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
    """Check metadata dtype ownership.

    Parameters
    ----------
    array : ndarray
        Array test input, supplied by pytest parametrization.
    """
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
    """Check unequal bins independent oracle.

    Parameters
    ----------
    batch : int or None
        Fourier-node batch size, supplied by pytest parametrization.
    gap : int or float
        Separation of the two redshift bins, supplied by pytest parametrization.
    reverse : bool
        Whether to reverse input ordering, supplied by pytest parametrization.
    """
    spec = scalar_spec()
    bins = []
    expected = []
    for i, (lo, hi, zeval, area, s) in enumerate(
        [(2, 2.7, 2.2, 10, 1), (2.7 + gap, 2.9 + gap, 2.85 + gap, 30, -1)]
    ):
        bin_geometry = prepare_geometry(
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
        prepared_bin = prepare_bin(replace(current, geometry=bin_geometry))
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
        np.testing.assert_allclose(bin_geometry.volume, volume, rtol=5e-13, atol=0)
        expected.append(modes / 18 * np.array([[1, s], [s, 1]]))
        bins.append(prepared_bin)
    if reverse:
        bins.reverse()
        expected.reverse()
    prior = np.diag([0.0, 4.0])
    result = run_forecast(bins, batch_size=batch, prior_fisher=prior)
    assert result.bin_ids == tuple(prepared_bin.id for prepared_bin in bins)
    for run, oracle in zip(result.bins, expected):
        np.testing.assert_allclose(run.result.data_fisher, oracle, rtol=5e-13, atol=0)
        assert run.result.diagnostics.rank == 1
    np.testing.assert_allclose(
        result.combined.data_fisher, sum(expected), rtol=5e-13, atol=0
    )
    np.testing.assert_array_equal(result.combined.prior_fisher, prior)
    assert result.combined.diagnostics.rank == 2


def test_explicit_evaluation_redshift_spies(monkeypatch):
    """Check explicit evaluation redshift spies.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    import fishhighz.forecast as forecast

    response_calls = []
    response = forecast.prepare_response

    def response_spy(*args, **kwargs):
        """Record the velocity conversion passed to response preparation.

        Parameters
        ----------
        *args : tuple
            Positional arguments forwarded to the original callable or accepted by
            the test callback.
        **kwargs : dict
            Keyword options forwarded to the original callable or inspected by the
            test callback.

        Returns
        -------
        response : ndarray
            Original per-field response amplitude on Fourier nodes.

        Notes
        -----
        Appends to the enclosing test call log so provider dispatch can be checked.
        """
        response_calls.append(kwargs["a_v"])
        return response(*args, **kwargs)

    monkeypatch.setattr(forecast, "prepare_response", response_spy)
    model_calls = []
    p1d_calls = []
    bins = []
    spec, _ = forest_spec(method="supplied")

    def model(t, z, k, mu, p):
        """Evaluate the synthetic spectrum used by the enclosing regression test.

        Parameters
        ----------
        t : ndarray of shape (n_parameters,)
            Local model parameters in the provider binding order.
        z : float or ndarray
            Dimensionless redshift.
        k : ndarray of shape (n_nodes,)
            Comoving wavenumbers in h/Mpc.
        mu : ndarray of shape (n_nodes,)
            Dimensionless line-of-sight direction cosines.
        p : ndarray of int, shape (n_pairs, 2)
            Observed-field indices defining the requested spectra.

        Returns
        -------
        power : ndarray of shape (n_nodes, n_pairs)
            Synthetic intrinsic power in (Mpc/h)^3 before response and noise.

        Notes
        -----
        Appends to the enclosing test call log so provider dispatch can be checked.
        """
        model_calls.append(z)
        return np.full((len(k), len(p)), t[0] * (1 + z))

    def p1d(t, z, k):
        """Evaluate an independent synthetic one-dimensional forest spectrum.

        Parameters
        ----------
        t : ndarray of shape (n_parameters,)
            Local model parameters in the provider binding order.
        z : float or ndarray
            Dimensionless redshift.
        k : ndarray of shape (n_nodes,)
            Line-of-sight velocity wavenumbers in s/km.

        Returns
        -------
        power : ndarray of shape (n_nodes,)
            Intrinsic one-dimensional power in km/s.

        Notes
        -----
        Appends to the enclosing test call log so provider dispatch can be checked.
        """
        p1d_calls.append(z)
        return np.full_like(k, 1 + z)

    owner = spec.p3d.routes[0].provider
    p3d = PreparedP3D(
        spec.p3d.registry,
        spec.p3d.selection,
        [P3DProvider(owner.label, model, owner.parameters, owner.pairs)],
    )
    for i, (lo, hi, z) in enumerate([(2, 2.7, 2.1), (2.8, 3, 2.97)]):
        bin_geometry = prepare_geometry(
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
        prepared_bin = prepare_bin(
            replace(
                spec, id=str(i), p3d=p3d, geometry=bin_geometry, forests={"f": source}
            )
        )
        bins.append(prepared_bin)
        np.testing.assert_array_equal(prepared_bin.power, 2 * (1 + z))
        np.testing.assert_allclose(
            response_calls[-1], 200 / ((1 + z) * 0.7), rtol=5e-13, atol=0
        )
    assert p1d_calls == [2.1, 2.97]
    model_calls.clear()
    run_forecast(bins, batch_size=5)
    assert model_calls == [2.1] * 15 + [2.97] * 15
    assert len(response_calls) == 2 and p1d_calls == [2.1, 2.97]
