"""Independent end-to-end information and fixed-survey ownership oracles."""

from dataclasses import replace

import numpy as np
import pytest

from fishhighz.fields import ObservedField, PairSelection
from fishhighz.forecast import prepare_bin, run_bin, run_forecast
from fishhighz.geometry import prepare_geometry
from fishhighz.grids import gauss_legendre_grid
from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.response import InstrumentResponse
from fishhighz.results import diagonal_prior
from fishhighz.survey import BinSpec, ForestInput


def geometry(lo=2, hi=3, h=0.7):
    """Prepare a synthetic common-volume redshift-bin geometry.

    Parameters
    ----------
    lo : float, optional
        Lower redshift boundary of the synthetic bin. Default is 2.
    hi : float, optional
        Upper redshift boundary of the synthetic bin. Default is 3.
    h : float, optional
        Dimensionless fiducial Hubble parameter H0/(100 km/s/Mpc). Default is
        0.7.

    Returns
    -------
    geometry : BinGeometry
        Distances in Mpc/h, volume in (Mpc/h)^3, and velocity conversion in km/s
        per Mpc/h.
    """
    return prepare_geometry(
        lo,
        hi,
        z_eval=(lo + hi) / 2,
        area_deg2=10,
        h_fid=h,
        z_order=4,
        hubble=lambda z: np.full_like(z, 200),
        transverse_distance=lambda z: 1000 * (1 + z),
    )


def scalar_spec(registry=None, sign=1, lo=2, id="one", nonlinear=False):
    """Construct a one-field amplitude/nuisance forecast bin.

    Parameters
    ----------
    registry : ParameterRegistry or None, optional
        Global parameter definitions; None constructs the helper default.
        Default is None.
    sign : float, optional
        Coefficient of the nuisance parameter in the synthetic spectrum. Default
        is 1.
    lo : float, optional
        Lower redshift boundary of the synthetic bin. Default is 2.
    id : str, optional
        Identifier assigned to the synthetic redshift bin. Default is 'one'.
    nonlinear : bool, optional
        Whether to exponentiate the amplitude parameter. Default is False.

    Returns
    -------
    spec : BinSpec
        Synthetic bin with a constant spectrum and unit fixed noise.
    """
    if registry is None:
        registry = ParameterRegistry(
            [
                Parameter("A", 2, "target", step=0.01),
                Parameter("b", 0, "nuisance", step=0.01),
            ]
        )

    def model(t, z, k, mu, pairs):
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
        pairs : ndarray of int, shape (n_pairs, 2)
            Observed-field indices defining the requested spectra.

        Returns
        -------
        power : ndarray of shape (n_nodes, n_pairs)
            Synthetic intrinsic power in (Mpc/h)^3 before response and noise.
        """
        return np.full(
            (len(k), len(pairs)), (np.exp(t[0]) if nonlinear else t[0]) + sign * t[1]
        )

    fields = [ObservedField("g", "galaxy", "g")]
    selection = PairSelection(fields)
    binding = BoundParameters(registry, ("A", "b"), {"A": "A", "b": "b"})
    p3d = PreparedP3D(
        registry,
        selection,
        [P3DProvider("scalar", model, binding, selection.required_pairs)],
    )
    grid = gauss_legendre_grid([0.02, 0.1, 0.2], k_order=2, mu_order=3, h_fid=0.7)
    return BinSpec(
        id,
        geometry(lo, lo + 1),
        grid,
        p3d,
        {"g": InstrumentResponse(0, 0)},
        full_noise=np.ones((len(grid.k_flat), 1)),
    )


@pytest.mark.parametrize("batch", [None, 1, 5, 100])
def test_joint_scalar_oracle(batch):
    """Check joint scalar oracle.

    Parameters
    ----------
    batch : int or None
        Fourier-node batch size, supplied by pytest parametrization.
    """
    spec = scalar_spec()
    bins = [
        prepare_bin(spec),
        prepare_bin(scalar_spec(spec.p3d.registry, -1, 3, "two")),
    ]
    prior = diagonal_prior(spec.p3d.registry, {"b": 0.5})
    result = run_forecast(bins, prior_fisher=prior, batch_size=batch)
    expected = []
    for b, s, r in zip(bins, [1, -1], result.bins):
        alpha = b.modes.sum() / 18
        bin_fisher = alpha * np.array([[1, s], [s, 1]])
        np.testing.assert_allclose(r.result.data_fisher, bin_fisher, rtol=5e-13, atol=0)
        assert r.result.diagnostics.rank == 1
        np.testing.assert_array_equal(r.result.prior_fisher, 0)
        np.testing.assert_allclose(
            b.factors[:, 0, 0] ** 2, 18 / b.modes, rtol=5e-13, atol=0
        )
        np.testing.assert_array_equal(b.power, 2)
        expected.append(bin_fisher)
    np.testing.assert_allclose(
        result.combined.data_fisher, sum(expected), rtol=5e-13, atol=0
    )
    np.testing.assert_array_equal(result.combined.prior_fisher, prior)
    assert result.combined.diagnostics.rank == 2
    assert result.bin_ids == ("one", "two")
    calls = len(result.bins[0].node_slices) * 5
    assert result.bins[0].calls[0].model == calls


def test_independent_nuisance_joint_null():
    """Check independent nuisance joint null."""
    registry = ParameterRegistry(
        [
            Parameter("A", 2, "target", step=0.01),
            Parameter("b1", 0, "nuisance", step=0.01),
            Parameter("b2", 0, "nuisance", step=0.01),
        ]
    )
    bins = []
    for i, s in enumerate([1, -1]):
        spec = scalar_spec()
        binding = BoundParameters(registry, ("A", "b"), {"A": "A", "b": f"b{i + 1}"})
        provider = P3DProvider(
            "scalar",
            lambda t, z, k, mu, p, s=s: np.full((len(k), len(p)), t[0] + s * t[1]),
            binding,
            [(0, 0)],
        )
        bins.append(
            prepare_bin(
                replace(
                    spec,
                    id=str(i),
                    geometry=geometry(2 + i, 3 + i),
                    p3d=PreparedP3D(registry, spec.p3d.selection, [provider]),
                )
            )
        )
    run = run_forecast(bins)
    assert run.combined.diagnostics.rank == 2
    np.testing.assert_allclose(
        run.combined.data_fisher @ np.array([1, -1, 1]), 0, atol=1e-11, rtol=0
    )


def test_signed_subset_matrix_oracle():
    """Check signed subset matrix oracle."""
    registry = ParameterRegistry([Parameter("A", 1.5, "target", step=0.01)])
    fields = [ObservedField(str(i), "galaxy", str(i)) for i in range(3)]
    selection = PairSelection(fields, [("2", "0"), ("1", "1"), ("0", "0")])
    base = np.array([[4, -1, 0.5], [-1, 3, -0.2], [0.5, -0.2, 2]])

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
        """
        return t[0] * (1 + k[:, None]) * base[p[:, 0], p[:, 1]]

    p3d = PreparedP3D(
        registry,
        selection,
        [
            P3DProvider(
                "matrix",
                model,
                BoundParameters(registry, ("A",), {"A": "A"}),
                selection.required_pairs,
            )
        ],
    )
    grid = scalar_spec().grid
    spec = BinSpec(
        "matrix",
        geometry(),
        grid,
        p3d,
        {f.id: InstrumentResponse(0, 0) for f in fields},
        galaxies={f.id: 1 for f in fields},
        independent_sampling=True,
    )
    prepared_bin = prepare_bin(spec)
    expected = 0
    for k, m in zip(prepared_bin.k, prepared_bin.modes):
        power_derivative = (1 + k) * base
        total_power = 1.5 * power_derivative + np.eye(3)
        pairs = selection.selected_pairs
        covariance_block = np.array(
            [
                [
                    (
                        total_power[a, c] * total_power[d, e]
                        + total_power[a, e] * total_power[d, c]
                    )
                    / m
                    for c, e in pairs
                ]
                for a, d in pairs
            ]
        )
        jac = np.array([power_derivative[a, d] for a, d in pairs])
        expected += jac @ np.linalg.solve(covariance_block, jac)
    np.testing.assert_allclose(
        run_bin(prepared_bin, batch_size=5).result.data_fisher,
        [[expected]],
        rtol=5e-12,
        atol=0,
    )


@pytest.mark.parametrize("batch", [0, -1, True, 1.5, "2"])
def test_bad_batch(batch):
    """Check bad batch.

    Parameters
    ----------
    batch : bool or int or float or str
        Fourier-node batch size, supplied by pytest parametrization.
    """
    with pytest.raises(ValueError, match="batch_size"):
        run_bin(prepare_bin(scalar_spec()), batch_size=batch)


def test_fixed_state_refinement_and_mutation(monkeypatch):
    """Check fixed state refinement and mutation.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    import fishhighz.forecast as module

    spec = scalar_spec(nonlinear=True)
    prepared_bin = prepare_bin(spec)
    arrays = [
        prepared_bin.theta,
        prepared_bin.k,
        prepared_bin.mu,
        prepared_bin.modes,
        prepared_bin.response,
        prepared_bin.products,
        prepared_bin.noise,
        prepared_bin.power,
        prepared_bin.total,
        prepared_bin.factors,
        prepared_bin.p3d.selection.required_pairs,
    ]
    hashes = [a.tobytes() for a in arrays]
    for a in arrays:
        with pytest.raises(ValueError):
            a.flags.writeable = True
    spec.p3d.selection.selected_to_required.flags.writeable = True
    spec.p3d.selection.selected_to_required[:] = 0
    for name in [
        "prepare_response",
        "prepare_noise",
        "prepare_forest_weights",
        "factor_covariance",
        "evaluate_p1d",
    ]:
        monkeypatch.setattr(
            module, name, lambda *a, **k: pytest.fail("fixed work repeated")
        )
    errors = []
    exact = prepared_bin.modes.sum() * np.exp(4) / (2 * (np.exp(2) + 1) ** 2)
    for scale in [1, 0.5, 0.25]:
        run = run_forecast([prepared_bin], step_scale=scale, prior_fisher=np.eye(2))
        errors.append(abs(run.combined.data_fisher[0, 0] - exact))
    assert 0.24 < errors[1] / errors[0] < 0.26
    assert 0.24 < errors[2] / errors[1] < 0.26
    assert hashes == [a.tobytes() for a in arrays]


@pytest.mark.parametrize("case", ["registry", "overlap", "duplicate", "h", "identity"])
def test_bin_structure_failures(case):
    """Check bin structure failures.

    Parameters
    ----------
    case : str
        Named forecast or validation case, supplied by pytest parametrization.
    """
    first = scalar_spec()
    second = scalar_spec(first.p3d.registry, lo=3, id="two")
    if case == "registry":
        second = scalar_spec(lo=3, id="two")
    if case == "overlap":
        second = replace(second, geometry=geometry(2.5, 3.5))
    if case == "duplicate":
        second = replace(second, id="one")
    if case == "h":
        second = replace(
            second,
            geometry=geometry(3, 4, 0.6),
            grid=gauss_legendre_grid([0.02, 0.2], k_order=2, mu_order=3, h_fid=0.6),
            full_noise=np.ones((6, 1)),
        )
    if case == "identity":
        selection = PairSelection([ObservedField("g", "galaxy", "different")])
        second = replace(
            second,
            p3d=PreparedP3D(
                first.p3d.registry, selection, [second.p3d.routes[0].provider]
            ),
        )
    with pytest.raises(ValueError):
        run_forecast([prepare_bin(first), prepare_bin(second)])


@pytest.mark.parametrize(
    "change",
    [
        {"forests": {}},
        {"galaxies": {"g": 1}},
        {"independent_sampling": True},
        {"responses": {}},
        {"full_noise": None},
        {"full_noise": np.ones((1, 1))},
        {"full_noise": -np.ones((12, 1))},
    ],
)
def test_bad_spec(change):
    """Check bad spec.

    Parameters
    ----------
    change : dict
        Input override exercising the specified validation boundary, supplied by
        pytest parametrization.
    """
    with pytest.raises(ValueError):
        prepare_bin(replace(scalar_spec(), **change))


def test_singular_and_bad_model():
    """Check singular and bad model."""
    spec = scalar_spec()
    for value in [0, np.nan, -2]:
        provider = P3DProvider(
            "broken",
            lambda t, z, k, mu, p, v=value: np.full((len(k), len(p)), v),
            spec.p3d.routes[0].provider.parameters,
            [(0, 0)],
        )
        with pytest.raises(ValueError, match="bin one"):
            prepare_bin(
                replace(
                    spec,
                    p3d=PreparedP3D(spec.p3d.registry, spec.p3d.selection, [provider]),
                    full_noise=np.zeros((12, 1)),
                )
            )


def test_error_global_slice_and_schedule():
    """Check error global slice and schedule."""
    spec = scalar_spec()

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

        Raises
        ------
        RuntimeError
            Deliberately raised to exercise the rejection path in the enclosing
            test.
        """
        if t[0] != 2 and np.any(k > 0.1):
            raise RuntimeError("outside domain")
        return np.full((len(k), len(p)), t[0])

    provider = P3DProvider(
        "domain", model, spec.p3d.routes[0].provider.parameters, [(0, 0)]
    )
    prepared_bin = prepare_bin(
        replace(
            spec, p3d=PreparedP3D(spec.p3d.registry, spec.p3d.selection, [provider])
        )
    )
    with pytest.raises(
        ValueError, match=r"global nodes \[6:7\).*parameter.*central"
    ) as error:
        run_bin(prepared_bin, batch_size=1)
    assert error.value.__cause__ is not None
    for steps in [{"A": 1e-30}, {"unknown": 1}]:
        with pytest.raises(ValueError):
            run_bin(prepared_bin, steps=steps)


def forest_spec(method="legacy", auxiliary=True):
    """Construct a forest bin with either computed or supplied weights.

    Parameters
    ----------
    method : str, optional
        Forest-weight prescription passed through to preparation. Default is
        'legacy'.
    auxiliary : bool, optional
        Whether to sample auxiliary P3D/P1D rather than supply their values.
        Default is True.

    Returns
    -------
    fixture : tuple
        BinSpec and the mutable list recording P1D queries.
    """
    spec = scalar_spec()
    fields = [
        ObservedField("f", "forest", "lya", background="qso"),
        ObservedField("unused", "forest", "lya", background="lbg"),
    ]
    selection = PairSelection(fields, [(0, 0)])
    p3d = PreparedP3D(spec.p3d.registry, selection, [spec.p3d.routes[0].provider])
    counts = []

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
        """
        counts.append(k.copy())
        return np.ones_like(k) * 2

    options = dict(
        z_source=4,
        magnitudes=[20, 21],
        quadrature=[1, 1],
        rho=[0.01, 0.02],
        variance=[1, 2],
        length_velocity=10000,
        method=method,
    )
    if method == "supplied":
        options["weights"] = [1, 1]
    else:
        options["iterations"] = 3
        if not auxiliary:
            options.update(signal=0.5, alias=2)
    source = ForestInput(
        options,
        p1d,
        BoundParameters(spec.p3d.registry, (), {}),
        spec.p3d.registry.fiducials,
        auxiliary_coordinates=(24, 0.00035)
        if auxiliary and method == "legacy"
        else None,
    )
    return replace(
        spec,
        p3d=p3d,
        responses={
            "f": InstrumentResponse(30, 10),
            "unused": InstrumentResponse(30, 10),
        },
        forests={"f": source},
        galaxies={},
        full_noise=None,
        independent_sampling=True,
    ), counts


@pytest.mark.parametrize("method", ["legacy", "supplied"])
def test_generated_spies(method):
    """Check generated spies.

    Parameters
    ----------
    method : str
        Weighting or evaluation prescription, supplied by pytest
        parametrization.
    """
    spec, calls = forest_spec(method)
    prepared_bin = prepare_bin(spec)
    assert len(calls) == (2 if method == "legacy" else 1)
    assert set(prepared_bin.weights) == {"f"}
    assert prepared_bin.diagnostics["p1d_calls"]["f"] == len(calls)
    assert prepared_bin.diagnostics["p3d_calls"]["scalar"] == len(calls)
    before = len(calls)
    run_bin(prepared_bin, batch_size=1)
    run_bin(prepared_bin, step_scale=0.5)
    assert len(calls) == before
    if method == "legacy":
        assert (
            prepared_bin.weights["f"].auxiliary.k > spec.grid.k_max
            or prepared_bin.weights["f"].auxiliary.k < spec.grid.k_min
        )


def test_full_noise_bypasses_sources():
    """Check full noise bypasses sources."""
    spec, calls = forest_spec()
    prepared_bin = prepare_bin(
        replace(
            spec,
            forests=None,
            galaxies=None,
            independent_sampling=None,
            full_noise=np.ones((12, 1)),
        )
    )
    run_bin(prepared_bin)
    assert not calls and not prepared_bin.weights


def test_kaiser_external_equivalence_and_ties():
    """Check kaiser external equivalence and ties."""
    from fishhighz.models.kaiser import KaiserModel
    from fishhighz.models.templates import prepare_template

    template_k_grid = np.linspace(0.005, 0.5, 30)
    template = prepare_template(
        template_k_grid,
        np.full_like(template_k_grid, 10),
        np.full_like(template_k_grid, 10),
        z_ref=2.5,
        h_template=0.7,
        h_fid=0.7,
    )
    registry = ParameterRegistry(
        [
            Parameter("bias", 2, "nuisance", step=0.01),
            Parameter("rate", 0.8, "target", step=0.01),
        ]
    )
    fields = [
        ObservedField("a", "galaxy", "same"),
        ObservedField("b", "galaxy", "same"),
    ]
    model = KaiserModel(
        template,
        fields,
        biases={"a": "ba", "b": "bb"},
        betas={},
        widths={"a": (0, 0), "b": (0, 0)},
        f="f",
        local_names=("ba", "bb", "f"),
    )

    def external(t, z, k, mu, p):
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
        """
        factors = t[:2][None, :] + t[2] * mu[:, None] ** 2
        return 10 * factors[:, p[:, 0]] * factors[:, p[:, 1]]

    def jac(t, z, k, mu, p):
        """Evaluate analytic local-parameter derivatives of the synthetic spectrum.

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
        jacobian : ndarray of shape (n_nodes, n_pairs, n_parameters)
            Derivatives of synthetic intrinsic power with respect to local
            parameters; units are power divided by parameter units.
        """
        result = np.zeros((len(k), len(p), 3))
        factors = t[:2][None, :] + t[2] * mu[:, None] ** 2
        for col, (i, j) in enumerate(p):
            result[:, col, i] += 10 * factors[:, j]
            result[:, col, j] += 10 * factors[:, i]
            result[:, col, 2] = 10 * mu**2 * (factors[:, i] + factors[:, j])
        return result

    selection = PairSelection(fields)
    binding = BoundParameters(
        registry, model.local_names, {"ba": "bias", "bb": "bias", "f": "rate"}
    )
    grid = scalar_spec().grid
    matrices = []
    for function, analytic in [(model, None), (external, None), (external, jac)]:
        p3d = PreparedP3D(
            registry,
            selection,
            [
                P3DProvider(
                    "model",
                    function,
                    binding,
                    selection.required_pairs,
                    jacobian=analytic,
                )
            ],
        )
        spec = BinSpec(
            "kaiser",
            geometry(),
            grid,
            p3d,
            {f.id: InstrumentResponse(0, 0) for f in fields},
            galaxies={"a": 1, "b": 1},
            independent_sampling=True,
        )
        prepared_bin = prepare_bin(spec)
        for batch in [None, 1, 5, 100]:
            run = run_bin(prepared_bin, batch_size=batch)
            matrices.append(run.result.data_fisher)
        if function is model:
            for wrong in [
                replace(spec, geometry=geometry(2.1, 3.1)),
                replace(spec, responses={"a": InstrumentResponse(0, 0)}),
            ]:
                with pytest.raises(ValueError):
                    prepare_bin(wrong)
    for f in matrices:
        np.testing.assert_allclose(f, matrices[0], rtol=5e-12, atol=0)


def test_template_perturbation_only_domain_failure():
    """Check template perturbation only domain failure."""
    from fishhighz.models.kaiser import KaiserModel, Scaling
    from fishhighz.models.templates import prepare_template

    grid = gauss_legendre_grid([0.1, 0.2], k_order=1, mu_order=1, h_fid=0.7)
    k_grid = np.linspace(0.1, 0.151, 8)
    template = prepare_template(
        k_grid,
        np.full_like(k_grid, 10),
        np.full_like(k_grid, 9),
        z_ref=2.5,
        h_template=0.7,
        h_fid=0.7,
    )
    fields = [ObservedField("g", "galaxy", "g")]
    registry = ParameterRegistry([Parameter("ap", 1, "target", step=0.2)])
    model = KaiserModel(
        template,
        fields,
        biases={"g": 2},
        betas={},
        widths={"g": (0, 0)},
        f=0.8,
        local_names=("sp", "st", "wp", "wt"),
        smooth=Scaling("ap_at", ap="sp", at="st"),
        wiggle=Scaling("ap_at", ap="wp", at="wt"),
    )
    selection = PairSelection(fields)
    p3d = PreparedP3D(
        registry,
        selection,
        [
            P3DProvider(
                "template",
                model,
                BoundParameters(
                    registry, model.local_names, {n: "ap" for n in model.local_names}
                ),
                [(0, 0)],
            )
        ],
    )
    prepared_bin = prepare_bin(
        BinSpec(
            "template",
            geometry(),
            grid,
            p3d,
            {"g": InstrumentResponse(0, 0)},
            full_noise=[[1]],
        )
    )
    with pytest.raises(ValueError, match="parameter.*central"):
        run_bin(prepared_bin)
    np.testing.assert_array_equal(prepared_bin.k, grid.k_flat)


def test_field_pair_parameter_and_bin_permutations():
    """Check field pair parameter and bin permutations."""
    base = np.array([[4, -1, 0.5], [-1, 3, -0.2], [0.5, -0.2, 2]])
    matrices = []
    for permutation in [(0, 1, 2), (2, 0, 1)]:
        names = ("A", "B") if permutation[0] == 0 else ("B", "A")
        registry = ParameterRegistry(
            [Parameter(n, 1.5 if n == "A" else 0.8, "target", step=0.01) for n in names]
        )
        binding = BoundParameters(registry, ("a", "b"), {"a": "A", "b": "B"})
        fields = [ObservedField(str(i), "galaxy", str(i)) for i in permutation]
        selected = [("0", "0"), ("0", "2"), ("1", "1")]
        if permutation[0] != 0:
            selected = selected[::-1]
        selection = PairSelection(fields, selected)
        reordered = base[np.ix_(permutation, permutation)]

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
            """
            return (t[0] + t[1] * k[:, None]) * reordered[p[:, 0], p[:, 1]]

        p3d = PreparedP3D(
            registry,
            selection,
            [P3DProvider("permuted", model, binding, selection.required_pairs)],
        )
        specs = [
            BinSpec(
                str(i),
                geometry(2 + i, 3 + i),
                scalar_spec().grid,
                p3d,
                {f.id: InstrumentResponse(0, 0) for f in fields},
                galaxies={f.id: 1 for f in fields},
                independent_sampling=True,
            )
            for i in range(2)
        ]
        bins = [prepare_bin(s) for s in specs]
        if permutation[0] != 0:
            bins = bins[::-1]
        result = run_forecast(bins, batch_size=5)
        order = [names.index("A"), names.index("B")]
        matrices.append(result.combined.data_fisher[np.ix_(order, order)])
    np.testing.assert_allclose(*matrices, rtol=5e-12, atol=0)


def test_indefinite_noise_not_hidden_by_signal():
    """Check indefinite noise not hidden by signal."""
    spec = scalar_spec()
    fields = [ObservedField("a", "galaxy", "a"), ObservedField("b", "galaxy", "b")]
    selection = PairSelection(fields)

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
        """
        return np.tile([100, 0, 100], (len(k), 1))

    p3d = PreparedP3D(
        spec.p3d.registry,
        selection,
        [
            P3DProvider(
                "large",
                model,
                spec.p3d.routes[0].provider.parameters,
                selection.required_pairs,
            )
        ],
    )
    with pytest.raises(ValueError, match="known noise"):
        prepare_bin(
            replace(
                spec,
                p3d=p3d,
                responses={f.id: InstrumentResponse(0, 0) for f in fields},
                full_noise=np.tile([1, 2, 1], (12, 1)),
            )
        )


def test_no_file_reads_on_repeated_runs(tmp_path, monkeypatch):
    """Check no file reads on repeated runs.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    from pathlib import Path

    from test_legacy_inputs import density_file, reader, snr_files

    from fishhighz.adapters.legacy_inputs import SNRReader, sample_forest_readers
    from fishhighz.geometry import LYA_REST_ANGSTROM
    from fishhighz.response import pixel_width_angstrom_to_velocity

    spec, _ = forest_spec(auxiliary=False)
    path = tmp_path / "density"
    density_file(path)
    density_reader = reader(path)
    snr = SNRReader(snr_files(tmp_path), smoothing="none")
    response = InstrumentResponse(
        pixel_width_angstrom_to_velocity(
            0.8, lambda_obs_angstrom=LYA_REST_ANGSTROM * (1 + spec.geometry.z_eval)
        ),
        10,
    )
    data = sample_forest_readers(
        density_reader,
        snr,
        spec.geometry,
        response,
        z_source=3.2,
        magnitudes=[20, 21],
        pixel_width_angstrom=0.8,
        exposure_count=4,
    )
    source = spec.forests["f"]
    options = dict(
        source.weight_options, z_source=3.2, rho=data["rho"], variance=data["variance"]
    )
    spec = replace(
        spec,
        forests={
            "f": replace(source, weight_options=options, provenance=data["provenance"])
        },
        responses=dict(spec.responses, f=response),
    )
    prepared_bin = prepare_bin(spec)
    monkeypatch.setattr(
        Path, "read_bytes", lambda *a: pytest.fail("file read during derivative run")
    )
    monkeypatch.setattr(
        np, "loadtxt", lambda *a, **k: pytest.fail("raw read during derivative run")
    )
    run_bin(prepared_bin)
    run_bin(prepared_bin, step_scale=0.5)
