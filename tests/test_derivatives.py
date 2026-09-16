"""Independent derivatives, scheduling, convergence and fixed-Fisher checks."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from fishhighz.derivatives import check_convergence, evaluate_derivatives
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.fisher import fisher_from_factors
from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
from fishhighz.parameters import Parameter, ParameterRegistry


def one(
    model, *, fiducial=0.0, step=0.01, bounds=None, jacobian=None, analytic_ids=None
):
    reg = ParameterRegistry(
        [Parameter("x", fiducial, "target", step=step, bounds=bounds)]
    )
    selection = PairSelection([ObservedField("A", "galaxy", "unused")])
    binding = BoundParameters(reg, ["local"], {"local": "x"})
    prepared = PreparedP3D(
        reg,
        selection,
        [
            P3DProvider(
                "toy",
                model,
                binding,
                [("A", "A")],
                jacobian=jacobian,
                analytic_ids=analytic_ids,
            )
        ],
    )
    return prepared, reg.fiducials, 2.0, np.array([0.1, 0.2]), np.array([0.2, 0.7])


def polynomial(t, z, k, mu, pairs):
    return np.full((len(k), len(pairs)), t[0] ** 2 + 3 * t[0])


def test_analytic_mapping_ties_shared_unused_empty_and_mixed():
    registry = ParameterRegistry(
        [
            Parameter("unused", 8, "target"),
            Parameter("cosmo", -0.4, "target", step=0.001),
            Parameter("nuisance", 0.3, "nuisance", step=0.002),
            Parameter("independent", 0.7, "nuisance", step=0.003),
        ]
    )
    selection = PairSelection(
        [
            ObservedField("a", "galaxy", "same"),
            ObservedField("b", "forest", "same", background="qso"),
        ]
    )
    tied = BoundParameters(
        registry, ["n", "x1", "x2"], {"n": "nuisance", "x1": "cosmo", "x2": "cosmo"}
    )
    other = BoundParameters(registry, ["i", "c"], {"i": "independent", "c": "cosmo"})
    empty = BoundParameters(registry, [], {})
    calls = {"a": [], "j": [], "b": [], "bj": [], "c": []}

    def a(t, z, k, mu, pairs):
        assert t[1] == t[2]  # every call must stay on the equality manifold
        calls["a"].append(t.copy())
        return (t[0] ** 2 + t[1] * t[2] + 2 * t[1]) * (1 + k[:, None])

    def aj(t, z, k, mu, pairs):
        calls["j"].append(t.copy())
        return (1 + k[:, None, None]) * np.array([2 * t[0], t[2] + 2, t[1]])[
            None, None, :
        ]

    def b(t, z, k, mu, pairs):
        calls["b"].append(t.copy())
        return np.full((len(k), 1), t[0] ** 2 + t[1])

    def bj(t, z, k, mu, pairs):
        calls["bj"].append(t.copy())
        return np.tile([2 * t[0], 1], (len(k), 1, 1))

    def constant(t, z, k, mu, pairs):
        calls["c"].append(t.copy())
        assert t.shape == (0,)
        return np.full((len(k), 1), -0.2)

    def prepared(analytic_ids=None, jac=aj):
        return PreparedP3D(
            registry,
            selection,
            [
                P3DProvider(
                    "a", a, tied, [("a", "a")], jacobian=jac, analytic_ids=analytic_ids
                ),
                P3DProvider("b", b, other, [("b", "b")], jacobian=bj),
                P3DProvider("c", constant, empty, [("b", "a")]),
            ],
        )

    k, mu = np.array([0.1, 0.2]), np.array([0.3, 0.6])
    expected = np.zeros((2, 3, 4))
    expected[:, 0, 1] = 1.2 * (1 + k)
    expected[:, 0, 2] = 0.6 * (1 + k)
    expected[:, 2, 1] = 1
    expected[:, 2, 3] = 1.4
    analytic = evaluate_derivatives(prepared(), registry.fiducials, 2, k, mu)
    np.testing.assert_allclose(analytic.jacobian, expected, atol=1e-14)
    assert [(c.model, c.jacobian) for c in analytic.calls] == [(1, 1), (1, 1), (1, 0)]
    numerical = evaluate_derivatives(
        prepared(), registry.fiducials, 2, k, mu, numerical=True
    )
    np.testing.assert_allclose(numerical.jacobian, expected, atol=1e-12)
    assert [(c.model, c.jacobian) for c in numerical.calls] == [(5, 0), (5, 0), (1, 0)]

    def placeholders(*args):
        value = aj(*args)
        value[:, :, 1:] = 12345  # deliberately wrong unused analytic columns
        return value

    mixed = evaluate_derivatives(
        prepared(["nuisance"], placeholders), registry.fiducials, 2, k, mu
    )
    np.testing.assert_allclose(mixed.jacobian, expected, atol=1e-12)
    assert [(c.model, c.jacobian) for c in mixed.calls] == [(3, 1), (1, 1), (1, 0)]
    assert len(calls["a"]) == 9 and len(calls["b"]) == 7 and len(calls["j"]) == 2
    assert len(calls["bj"]) == 2 and len(calls["c"]) == 3
    assert {c.parameter for c in mixed.columns} == {"cosmo", "nuisance", "independent"}
    for ids in (["unused"], ["nuisance", "nuisance"]):
        with pytest.raises(ValueError, match="analytic"):
            prepared(ids)


@pytest.mark.parametrize(
    "x,bounds,method",
    [
        (0, (-1, 1), "central"),
        (-0.7, (-1, 1), "central"),
        (0, (0, 1), "forward"),
        (1, (0, 1), "backward"),
        (0.005, (0, 1), "forward"),
        (0.995, (0, 1), "backward"),
    ],
)
def test_second_order_polynomial_and_bounds(x, bounds, method):
    seen = []

    def model(t, *args):
        seen.append(t[0])
        assert bounds[0] <= t[0] <= bounds[1]
        return polynomial(t, *args)

    arguments = one(model, fiducial=x, bounds=bounds)
    result = evaluate_derivatives(*arguments)
    np.testing.assert_allclose(result.jacobian, 2 * x + 3, rtol=0, atol=2e-13)
    assert result.columns[0].method == method
    assert len(seen) == result.calls[0].model == 3
    assert seen.count(x) == 1


@pytest.mark.parametrize(
    "kwargs",
    [
        {"step": None},
        {"step": 0.2, "bounds": (-0.1, 0.1)},
        {"step": 1e-10, "fiducial": 1e100},
    ],
)
def test_bad_schedule_before_any_dispatch(kwargs):
    seen = []
    arguments = one(lambda *args: seen.append(args), **kwargs)
    with pytest.raises(ValueError, match="x.*step"):
        evaluate_derivatives(*arguments)
    assert seen == []


@pytest.mark.parametrize(
    "kwargs",
    [
        {"steps": {"x": 0}},
        {"steps": {"x": -1}},
        {"steps": {"x": np.nan}},
        {"steps": {"x": True}},
        {"steps": {"unknown": 0.01}},
        {"step_scale": 0},
        {"step_scale": np.inf},
        {"step_scale": True},
        {"step_scale": 1e308, "steps": {"x": 1e308}},
        {"numerical": "yes"},
    ],
)
def test_invalid_overrides(kwargs):
    seen = []
    with pytest.raises(ValueError):
        evaluate_derivatives(*one(lambda *a: seen.append(a)), **kwargs)
    assert seen == []


def test_schedule_is_global_and_analytic_or_unused_needs_no_step():
    args = one(
        polynomial,
        step=None,
        jacobian=lambda t, z, k, mu, pairs: np.full((len(k), len(pairs), 1), 3),
    )
    assert evaluate_derivatives(*args).calls[0].model == 1
    with pytest.raises(ValueError, match="explicit"):
        evaluate_derivatives(*args, numerical=True)
    result = evaluate_derivatives(*args, numerical=True, steps={"x": 0.02})
    np.testing.assert_allclose(result.jacobian, 3, atol=1e-14)
    reg = ParameterRegistry(
        [Parameter("ok", 1, "target", step=0.01), Parameter("bad", 1, "target")]
    )
    selection = PairSelection([ObservedField("a", "galaxy", "x")])
    calls = []
    bound = BoundParameters(reg, ["ok", "bad"], {"ok": "ok", "bad": "bad"})
    p = PreparedP3D(
        reg,
        selection,
        [P3DProvider("preflight", lambda *a: calls.append(a), bound, [(0, 0)])],
    )
    with pytest.raises(ValueError, match="bad.*explicit"):
        evaluate_derivatives(p, reg.fiducials, 0, [0.1], [0.2])
    assert not calls


def test_actual_unequal_offsets_and_tiny_steps():
    # At a power-of-two boundary, upward/downward spacings differ.
    x = 1.0
    h = 0.6 * np.spacing(x)

    def shifted(t, z, k, mu, pairs):
        delta = t[0] - x
        return np.full((len(k), len(pairs)), delta**2 + delta)

    result = evaluate_derivatives(*one(shifted, fiducial=x, step=h))
    stencil = result.columns[0].stencil
    assert abs(stencil.offsets[0]) != abs(stencil.offsets[1])
    np.testing.assert_allclose(result.jacobian, 1, atol=5e-16)
    # Coefficients never form 1/h**2, which would overflow here.
    result = evaluate_derivatives(
        *one(lambda t, z, k, mu, p: np.full((len(k), len(p)), t[0]), step=1e-200)
    )
    np.testing.assert_array_equal(result.jacobian, 1)


def test_finite_output_but_nonfinite_arithmetic_fails():
    def discontinuous(t, z, k, mu, pairs):
        return np.full((len(k), len(pairs)), -1e308 if t[0] < 0 else 1e308)

    with pytest.raises(ValueError, match="nonfinite derivative arithmetic"):
        evaluate_derivatives(*one(discontinuous))


def exponential(t, z, k, mu, pairs):
    return np.full((len(k), len(pairs)), np.exp(t[0]))


def test_second_order_convergence_and_bad_step():
    arguments = one(exponential, fiducial=-0.2, step=0.2)
    errors = []
    for scale in (1, 0.5, 0.25):
        result = evaluate_derivatives(*arguments, step_scale=scale)
        assert result.calls[0].model == 3
        errors.append(abs(result.jacobian[0, 0, 0] - np.exp(-0.2)))
    np.testing.assert_allclose(np.array(errors[:-1]) / errors[1:], 4, rtol=0.002)
    bad = check_convergence(*arguments, atol=1e-10, rtol=1e-5, refinements=2)
    assert not bad.passed and not bad.comparisons[0].passed
    good = check_convergence(
        *arguments, steps={"x": 0.001}, atol=1e-10, rtol=1e-6, refinements=2
    )
    assert good.passed
    assert len(good.comparisons) == 2
    assert all(
        c.absolute_change > 0 and c.relative_change > 0 for c in good.comparisons
    )


def test_convergence_zero_canceling_and_family_change():
    result = check_convergence(
        *one(lambda t, z, k, mu, p: np.ones((len(k), len(p)))), atol=0, rtol=0
    )
    assert result.passed and result.comparisons[0].relative_change == 0
    even = check_convergence(
        *one(lambda t, z, k, mu, p: np.full((len(k), len(p)), t[0] ** 2)),
        atol=0,
        rtol=0,
    )
    assert even.passed and even.comparisons[0].absolute_change == 0
    changed = check_convergence(
        *one(exponential, fiducial=0.075, step=0.1, bounds=(0, 1)), atol=1, rtol=0
    )
    assert changed.comparisons[0].family_changed
    assert changed.comparisons[0].coarse.method == "forward"
    assert changed.comparisons[0].fine.method == "central"


def test_analytic_is_not_a_convergence_test_and_invalid_tolerances():
    arguments = one(
        exponential,
        jacobian=lambda t, z, k, mu, p: np.full((len(k), len(p), 1), np.exp(t[0])),
    )
    with pytest.raises(ValueError, match="no numerical"):
        check_convergence(*arguments, atol=1e-4, rtol=0)
    assert check_convergence(*arguments, numerical=True, atol=1e-4, rtol=0).passed
    for options in (
        {"atol": -1, "rtol": 0},
        {"atol": 0, "rtol": np.nan},
        {"atol": 0, "rtol": 0, "refinements": 0},
    ):
        with pytest.raises(ValueError):
            check_convergence(*arguments, **options)


@pytest.mark.parametrize(
    "output", [np.ones((2, 1)), np.full((2, 1, 1), 1j), np.full((2, 1, 1), np.nan)]
)
def test_invalid_supplied_jacobian(output):
    with pytest.raises(ValueError, match="analytic Jacobian"):
        evaluate_derivatives(*one(polynomial, jacobian=lambda *a: output))


def test_slices_noncontiguous_and_no_stale_cache():
    arguments = list(
        one(lambda t, z, k, mu, p: (np.exp(t[0] * k) + mu)[:, None], fiducial=0.7)
    )
    k = np.linspace(0.01, 0.4, 20)[::2]
    mu = np.linspace(0, 1, 20)[::2]
    k.flags.writeable = mu.flags.writeable = False
    arguments[3:] = [k, mu]
    full = evaluate_derivatives(*arguments)
    pieces = [
        evaluate_derivatives(*arguments[:3], k[s], mu[s])
        for s in (slice(0, 3), slice(3, 7), slice(7, None))
    ]
    np.testing.assert_array_equal(np.concatenate([r.power for r in pieces]), full.power)
    np.testing.assert_array_equal(
        np.concatenate([r.jacobian for r in pieces]), full.jacobian
    )
    factors = np.full((len(k), 1, 1), 2.0)
    np.testing.assert_allclose(
        sum(
            fisher_from_factors(r.jacobian, np.full((len(r.power), 1, 1), 2.0))
            for r in pieces
        ),
        fisher_from_factors(full.jacobian, factors),
        rtol=2e-16,
    )
    arguments[1] = [1.2]
    changed = evaluate_derivatives(*arguments)
    assert not np.array_equal(changed.power, full.power)
    np.testing.assert_allclose(
        changed.jacobian[:, 0, 0], k * np.exp(1.2 * k), rtol=3e-6
    )


def test_end_to_end_example_fixed_covariance(monkeypatch):
    spec = importlib.util.spec_from_file_location(
        "synthetic_example",
        Path(__file__).resolve().parents[1] / "examples/external_forecast.py",
    )
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    covariance_calls, factor_calls, snapshots = [], [], []
    fixed_inputs = []
    grid_factory = example.gauss_legendre_grid
    combine = example.combine_observed_power

    def saved_grid(*args, **kwargs):
        grid = grid_factory(*args, **kwargs)
        for name in (
            "k_flat",
            "mu_flat",
            "q_mode",
            "weights",
            "k",
            "mu",
            "w_k",
            "w_mu",
        ):
            value = getattr(grid, name)
            fixed_inputs.append((value, value.copy()))
        fixed_inputs.append((np.array([grid.k_min, grid.k_max]), np.array([0.05, 0.3])))
        return grid

    def saved_noise(signal, noise):
        fixed_inputs.append((noise, noise.copy()))
        return combine(signal, noise)

    covariance = example.gaussian_covariance
    factor = example.factor_covariance
    fisher = example.fisher_from_factors

    def counted_covariance(total, modes, selection):
        covariance_calls.append((total.copy(), modes.copy(), selection))
        assert selection.selected_to_required.tolist() == [2, 1]
        assert len(selection.required_pairs) == 3
        return covariance(total, modes, selection)

    def counted_factor(value):
        factor_calls.append(value.copy())
        return factor(value)

    def checked_fisher(jac, factors):
        snapshots.append((id(factors), factors.copy()))
        return fisher(jac, factors)

    monkeypatch.setattr(example, "gauss_legendre_grid", saved_grid)
    monkeypatch.setattr(example, "combine_observed_power", saved_noise)
    monkeypatch.setattr(example, "gaussian_covariance", counted_covariance)
    monkeypatch.setattr(example, "factor_covariance", counted_factor)
    monkeypatch.setattr(example, "fisher_from_factors", checked_fisher)
    result = example.run()
    assert result["convergence_passed"]
    for value, snapshot in fixed_inputs:
        np.testing.assert_array_equal(value, snapshot)
    assert len(covariance_calls) == len(factor_calls) == 1
    assert len(snapshots) == 4 and len({s[0] for s in snapshots}) == 1
    for _, value in snapshots:
        np.testing.assert_array_equal(value, snapshots[0][1])
    errors = [r["max_fisher_relative_error"] for r in result["comparisons"]]
    np.testing.assert_allclose(np.array(errors[:-1]) / errors[1:], 4, rtol=1e-4)


@pytest.mark.parametrize("multiple,refinements", [(1.2, 1), (2.4, 2)])
def test_convergence_repeated_actual_stencil_preflight(multiple, refinements):
    calls = []

    def model(*args):
        calls.append("model")
        return exponential(*args)

    def jacobian(t, z, k, mu, pairs):
        calls.append("jacobian")
        return np.full((len(k), len(pairs), 1), np.exp(t[0]))

    h = multiple * np.spacing(1.5)
    arguments = one(model, fiducial=1.5, step=h, jacobian=jacobian)
    with pytest.raises(ValueError, match="no effective refinement") as error:
        check_convergence(
            *arguments, atol=0, rtol=0, refinements=refinements, numerical=True
        )
    message = str(error.value)
    for context in (
        "provider 'toy'",
        "global parameter 'x'",
        f"refinement level {refinements}",
        str(h * 0.5 ** (refinements - 1)),
        str(h * 0.5**refinements),
        "1.4999999999999998",
        "1.5000000000000002",
        "better-resolved starting step",
        "fewer refinements",
    ):
        assert context in message
    assert calls == []
    # Ordinary differentiation is still legal; this guard assesses refinement,
    # not the accuracy of a single rounding-dominated derivative.
    assert evaluate_derivatives(*arguments, numerical=True).calls[0].model == 3


@pytest.mark.parametrize("bad_first", [False, True])
def test_later_repeat_preflights_all_providers_and_mixed_columns(bad_first):
    registry = ParameterRegistry(
        [
            Parameter("valid", 0.5, "target", step=0.01),
            Parameter("analytic", 1.0, "nuisance"),
            Parameter("bad", 1.5, "target", step=2.4 * np.spacing(1.5)),
        ]
    )
    selection = PairSelection(
        [ObservedField("A", "galaxy", "toy"), ObservedField("B", "galaxy", "toy")]
    )
    calls = []

    def model(*args):
        calls.append("model")
        raise AssertionError("study must fail before dispatch")

    def jacobian(*args):
        calls.append("jacobian")
        raise AssertionError("study must fail before dispatch")

    valid = BoundParameters(registry, ["v", "a"], {"v": "valid", "a": "analytic"})
    mixed = BoundParameters(
        registry, ["v", "a", "b"], {"v": "valid", "a": "analytic", "b": "bad"}
    )
    providers = [
        P3DProvider(
            "valid",
            model,
            valid,
            [(0, 0)],
            jacobian=jacobian,
            analytic_ids=["analytic"],
        ),
        P3DProvider(
            "mixed",
            model,
            mixed,
            [(0, 1), (1, 1)],
            jacobian=jacobian,
            analytic_ids=["analytic"],
        ),
    ]
    if bad_first:
        providers.reverse()
    prepared = PreparedP3D(registry, selection, providers)
    with pytest.raises(ValueError, match="provider 'mixed'.*'bad'.*refinement level 2"):
        check_convergence(
            prepared, registry.fiducials, 2, [0.1], [0.5], atol=0, rtol=0, refinements=2
        )
    assert calls == []


def test_convergence_one_point_moves_is_effective_refinement():
    # At x=1 the downward ULP is half the upward ULP. Halving this request
    # moves the lower point only. Approximate equality would reject it wrongly.
    def shifted_linear(t, z, k, mu, pairs):
        return np.full((len(k), len(pairs)), t[0] - 1.0)

    result = check_convergence(
        *one(shifted_linear, fiducial=1.0, step=1.2 * np.spacing(1.0)),
        atol=0,
        rtol=0,
    )
    comparison = result.comparisons[0]
    assert comparison.coarse.points[0] != comparison.fine.points[0]
    assert comparison.coarse.points[1] == comparison.fine.points[1]
    assert result.passed
