"""Batch-boundary numerical parity and bounded accuracy reuse regressions."""

import importlib.util
from collections import Counter, OrderedDict
from types import SimpleNamespace

import numpy as np
import pytest

from fishhighz.covariance import (
    _validate_field_power,
    _validate_field_power_scalar,
    gaussian_covariance,
)
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.fisher import (
    _factor_covariance_scalar,
    factor_covariance,
    fisher_from_factors,
)
from fishhighz.validation.accuracy import DEFAULT, AccuracyRecipe
from fishhighz.validation.numerics import relative
from fishhighz.validation.schema import request
from fishhighz.validation.synthetic import evidence_payload


@pytest.fixture(params=["numpy", "numba"])
def backend(request, monkeypatch):
    """Select a Fisher-contraction backend for the parametrized test.

    Parameters
    ----------
    request : pytest.FixtureRequest
        Pytest request containing the selected fixture parameter.
    monkeypatch : pytest.MonkeyPatch
        Fixture restoring patched attributes and environment variables after the
        test.

    Returns
    -------
    backend : str
        Selected numpy or numba backend; unavailable numba is skipped.
    """
    if request.param == "numba" and importlib.util.find_spec("numba") is None:
        pytest.skip("optional compiler unavailable")
    monkeypatch.setenv("FISHHIGHZ_FISHER_BACKEND", request.param)
    return request.param


@pytest.mark.parametrize("nodes", [1, 255, 256, 257, 513])
@pytest.mark.parametrize("columns", [2, 12])
@pytest.mark.parametrize("selected", [None, [(4, 2), (1, 1), (0, 3)]])
def test_spd_oracles(backend, nodes, columns, selected):
    """Check spd oracles.

    Parameters
    ----------
    backend : str
        Fisher-contraction implementation, supplied by pytest parametrization.
    nodes : int
        Evaluation nodes, supplied by pytest parametrization.
    columns : int
        Selected Jacobian or table columns, supplied by pytest parametrization.
    selected : list or None
        Selected spectrum definitions, supplied by pytest parametrization.
    """
    rng = np.random.default_rng(1204)
    sel = PairSelection(
        [ObservedField(str(i), "galaxy", "m") for i in range(5)], selected
    )
    raw = rng.normal(size=(nodes, 5, 5))
    total = raw @ raw.swapaxes(1, 2) + 5 * np.eye(5)
    i, j = sel.required_pairs.T
    modes = rng.uniform(0.2, 4, nodes)
    covariance = gaussian_covariance(total[:, i, j], modes, sel)
    pairs = sel.selected_pairs
    oracle_c = np.stack(
        [
            np.stack(
                [
                    (total[:, a, m] * total[:, b, n] + total[:, a, n] * total[:, b, m])
                    / modes
                    for m, n in pairs
                ],
                axis=1,
            )
            for a, b in pairs
        ],
        axis=1,
    )
    assert relative(covariance, oracle_c) <= 5e-12
    jac = rng.normal(size=(nodes, len(pairs), columns))
    if columns == 12:
        jac[:, :, 2:] = 0
    saved = covariance.copy(), jac.copy()
    scales = np.geomspace(1e-50, 1e50, len(pairs))
    factors = factor_covariance(covariance)
    reference = _factor_covariance_scalar(covariance)
    assert relative(factors, reference) <= 5e-12
    actual = fisher_from_factors(jac, factors)
    direct = np.einsum("nsi,nsj->ij", jac, np.linalg.solve(covariance, jac))
    assert relative(actual, direct) <= 5e-12
    scaled = fisher_from_factors(
        jac * scales[None, :, None],
        factor_covariance(covariance * scales[None, :, None] * scales[None, None, :]),
    )
    assert relative(scaled, direct) <= 5e-12
    errors = np.sqrt(np.diag(np.linalg.inv(actual[:2, :2])))
    assert relative(errors, np.sqrt(np.diag(np.linalg.inv(direct[:2, :2])))) <= 5e-12
    if columns == 12:
        assert np.count_nonzero(actual[2:]) == 0
    np.testing.assert_array_equal(covariance, saved[0])
    np.testing.assert_array_equal(jac, saved[1])
    assert actual.flags.owndata


def outcome(function, *args):
    """Return a calculation result or its ValueError diagnostic.

    Parameters
    ----------
    function : callable
        Synthetic model or operation evaluated by the helper.
    *args : tuple
        Positional arguments forwarded to the original callable or accepted by
        the test callback.

    Returns
    -------
    outcome : object or str
        Function result on success, or the ValueError message on rejection.
    """
    try:
        return function(*args)
    except ValueError as error:
        return str(error)


@pytest.mark.parametrize("node", [0, 255, 256, 257, 511])
@pytest.mark.parametrize(
    "bad",
    [
        [[-1.0, 0], [0, 1]],
        [[0.0, 0.1], [0.1, 1]],
        [[0.0, 0], [0, 1]],
        [[1.0, 2], [2, 1]],
        [[1.0, 0.1], [0.2, 1]],
        [[np.nan, 0], [0, 1]],
        [[1.0, 1 - 1e-15], [1 - 1e-15, 1]],
        [[1.0, 1 - 1e-12], [1 - 1e-12, 1]],
        [[1.0, 0.1 + 1e-15], [0.1, 1]],
        [[1.0, 0.1 + 1e-12], [0.1, 1]],
        [[np.nextafter(0.0, 1.0), 0], [0, np.nextafter(0.0, 1.0)]],
    ],
)
def test_factor_boundary_diagnostics(node, bad):
    """Check factor boundary diagnostics.

    Parameters
    ----------
    node : int
        Fourier-node index, supplied by pytest parametrization.
    bad : list
        Invalid input exercising the specified rejection path, supplied by
        pytest parametrization.
    """
    covariance = np.tile(np.eye(2), (513, 1, 1))
    covariance[node] = bad
    covariance[512, 0, 0] = -1  # first error must still be in the original cell order
    assert outcome(factor_covariance, covariance) == outcome(
        _factor_covariance_scalar, covariance
    )


@pytest.mark.parametrize("node", [255, 256, 257])
@pytest.mark.parametrize(
    "bad", [[-1, 0, 1], [0, 0.1, 1], [0, 0, 1], [1, 1, 1], [1, 1 + 1e-12, 1], [1, 2, 1]]
)
def test_field_boundary_diagnostics(node, bad):
    """Check field boundary diagnostics.

    Parameters
    ----------
    node : int
        Fourier-node index, supplied by pytest parametrization.
    bad : list
        Invalid input exercising the specified rejection path, supplied by
        pytest parametrization.
    """
    sel = PairSelection([ObservedField(str(i), "galaxy", "m") for i in range(2)])
    packed_power = np.tile([1.0, -0.2, 1.0], (513, 1))
    packed_power[node] = bad
    packed_power[512, 0] = -1
    assert outcome(_validate_field_power, packed_power, sel) == outcome(
        _validate_field_power_scalar, packed_power, sel
    )


@pytest.mark.parametrize("kind", ["triangle", "factor_nan", "jac_nan", "solve", "sum"])
def test_compiled_error_order_and_overflow(backend, monkeypatch, kind):
    """Check compiled error order and overflow.

    Parameters
    ----------
    backend : str
        Fisher-contraction implementation, supplied by pytest parametrization.
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    kind : str
        Tracer, input, or calculation classification for this case, supplied by
        pytest parametrization.
    """
    lower = np.tile(np.eye(2), (258, 1, 1))
    j = np.ones((258, 2, 2))
    if kind == "triangle":
        lower[256, 0, 1] = np.nextafter(0.0, 1.0)
    elif kind == "factor_nan":
        lower[256, 0, 0] = np.nan
    elif kind == "jac_nan":
        j[256, 0, 1] = np.nan
    elif kind == "solve":
        lower[256, 0, 0] = np.nextafter(0.0, 1.0)
    else:
        j[256] = 1e155
    actual = outcome(fisher_from_factors, j, lower)
    monkeypatch.setenv("FISHHIGHZ_FISHER_BACKEND", "numpy")
    assert actual == outcome(fisher_from_factors, j, lower)


@pytest.mark.parametrize("side", ["production", "oracle"])
def test_independent_accuracy_check(monkeypatch, side):
    """Check independent accuracy check.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    side : str
        Implementation whose response is perturbed, supplied by pytest
        parametrization.
    """
    import fishhighz.validation.accuracy as accuracy

    task = request("lya_qso_lbg_lae_15x2pt", 0, "accuracy", kind="synthetic_bao")
    evidence_arrays, evidence_report = evidence_payload(task)
    factors = factor_covariance(evidence_arrays["selected_covariance"])
    prepared = SimpleNamespace(
        p3d=None,
        theta=None,
        k=evidence_arrays["k"],
        mu=evidence_arrays["mu"],
        products=np.ones_like(evidence_arrays["total"]),
        total=evidence_arrays["total"],
        factors=factors,
        power=evidence_arrays["total"],
        response=np.ones((len(evidence_arrays["k"]), 5)),
        noise=np.zeros_like(evidence_arrays["total"]),
    )
    recipe = AccuracyRecipe.__new__(AccuracyRecipe)
    recipe.registry = SimpleNamespace(ids=task["parameters"])
    recipe.selection = SimpleNamespace(selected_to_required=np.arange(15))
    recipe.prepare = lambda *args, **kwargs: (prepared, evidence_report["settings"])
    recipe.z = lambda index: 2.5
    recipe.provenance = {}
    calls = []

    def derivative(*args, **kwargs):
        """Record derivative step sizes and return the fixed synthetic Jacobian.

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
        result : types.SimpleNamespace
            Object containing observed_j as its jacobian attribute.

        Notes
        -----
        Appends to the enclosing test call log so provider dispatch can be checked.
        """
        calls.append(kwargs["step_scale"])
        return SimpleNamespace(jacobian=evidence_arrays["observed_j"])

    monkeypatch.setattr(accuracy, "evaluate_derivatives", derivative)
    actual, _ = recipe.evaluate(task, DEFAULT)
    assert relative(actual["fisher"], evidence_arrays["fisher"]) < 5e-12
    assert calls == [DEFAULT["step"] / 0.001]
    original = accuracy._assemble if side == "production" else accuracy.contract

    def perturb(*args, **kwargs):
        """Perturb the selected side of an independent Fisher comparison.

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
        result : tuple
            Original calculation output with a one-percent Fisher perturbation.
        """
        perturbed_result, unmodified_result = original(*args, **kwargs)
        if side == "production":
            perturbed_result["fisher"] *= 1.01
        else:
            perturbed_result *= 1.01
        return perturbed_result, unmodified_result

    monkeypatch.setattr(
        accuracy, "_assemble" if side == "production" else "contract", perturb
    )
    with pytest.raises(ValueError, match="direct NumPy"):
        recipe.evaluate(task, DEFAULT)


def test_study_reuse_controls_failures_and_order(monkeypatch):
    """Check study reuse controls failures and order.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    import fishhighz.validation.accuracy as accuracy

    caches = []

    class TrackedPayloads(OrderedDict):
        def __init__(self):
            """Initialize the synthetic TrackedPayloads fixture.

            Notes
            -----
            Sets the instance state used by the enclosing test; no scientific calculation is run.
            """
            super().__init__()
            self.evictions = 0
            self.hits = 0
            caches.append(self)

        def __getitem__(self, key):
            """Record a cache hit and retrieve the original value.

            Parameters
            ----------
            key : object
                Cache key whose access is recorded.

            Returns
            -------
            value : object
                Cached payload at the requested key.
            """
            self.hits += 1
            return super().__getitem__(key)

        def popitem(self, last=True):
            """Record eviction and enforce the synthetic three-entry cache bound.

            Parameters
            ----------
            last : bool, optional
                Whether to evict the most recently inserted cache entry. Default is
                True.

            Returns
            -------
            item : tuple
                Evicted key and cached value.
            """
            self.evictions += 1
            value = super().popitem(last=last)
            assert len(self) <= 3
            return value

    monkeypatch.setattr(accuracy, "OrderedDict", TrackedPayloads)

    class SyntheticStudy(AccuracyRecipe):
        def __init__(self):
            """Initialize the synthetic SyntheticStudy fixture.

            Notes
            -----
            Sets the instance state used by the enclosing test; no scientific calculation is run.
            """
            self.selection = SimpleNamespace(
                fields=[SimpleNamespace(kind="forest")],
                selected_pairs=np.array([[0, 0]]),
            )
            self.calls = Counter()
            self._prepared = OrderedDict()

        def evaluate(self, task, controls, **options):
            """Evaluate the synthetic Fisher trial for the supplied numerical controls.

            Parameters
            ----------
            task : dict
                Synthetic forecast request including case, redshift-bin index, profile,
                and pair selection.
            controls : dict
                Numerical quadrature and weighting controls for this synthetic trial.
            **options : dict
                Keyword options forwarded to the original callable or inspected by the
                test callback.

            Returns
            -------
            payload : tuple
                Synthetic numerical arrays and validation report, including the
                requested test modification.

            Raises
            ------
            ValueError
                Deliberately raised to exercise the rejection path in the enclosing
                test.

            Notes
            -----
            Uses small analytic matrices to exercise validation control flow; it does not run a survey forecast.
            """
            key = tuple(sorted(controls.items()))
            self.calls[key] += 1
            if controls["iterations"] == 3:
                raise ValueError("synthetic failed weight trial")
            # All supported controls influence independently computed information.
            amplitude = 1 + sum(float(v) * 1e-12 for v in controls.values())
            fisher = amplitude * np.eye(2)
            return dict(fisher=fisher, pair_fisher=fisher[None]), dict(
                settings=dict(grid=dict(volume=1.0))
            )

    recipe = SyntheticStudy()
    arrays, report = recipe.study({})
    keys = [tuple(sorted(c.items())) for c in report["study_controls"]]
    assert len(keys) == len(set(keys))
    assert set(keys) <= set(recipe.calls)
    assert all(dict(k)["iterations"] != 3 for k in keys)
    for c, fisher in zip(report["study_controls"], arrays["study_fisher"]):
        np.testing.assert_array_equal(
            fisher, (1 + sum(float(v) * 1e-12 for v in c.values())) * np.eye(2)
        )
    assert report["unresolved_controls"] and not report["passed"]
    before = recipe.calls.copy()
    again, repeated = recipe.study({})
    assert recipe.calls == Counter({k: 2 * v for k, v in before.items()})
    np.testing.assert_array_equal(again["study_fisher"], arrays["study_fisher"])
    assert repeated["study_controls"] == report["study_controls"]
    assert len(caches) == 2
    assert all(
        len(cache) <= 3 and cache.evictions > 0 and cache.hits > 0 for cache in caches
    )


def test_numpy_only_import_and_requested_compiler_fallback(tmp_path):
    """Check numpy only import and requested compiler fallback.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    import os
    import subprocess
    import sys

    code = """
import importlib.abc
import sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, *args):
        if fullname.split('.')[0] in {'numba', 'scipy', 'astropy', 'camb', 'lyaforecast'}:
            raise ImportError('optional imports unavailable')
sys.meta_path.insert(0, Block())
import numpy as np
from fishhighz.fisher import fisher_matrix
assert not {'numba', 'scipy', 'astropy'} & set(sys.modules)
assert fisher_matrix(np.ones((1,1,1)),np.ones((1,1,1)))[0,0] == 1
assert not {'numba', 'scipy', 'astropy'} & set(sys.modules)
"""
    # The installed-wheel run executes this same regression outside the checkout.
    for backend in ("numpy", "numba"):
        env = {**os.environ, "FISHHIGHZ_FISHER_BACKEND": backend}
        subprocess.run([sys.executable, "-c", code], env=env, check=True)


@pytest.mark.parametrize("policy", ["raise", "warn"])
def test_explicit_underflow_policy_is_preserved(monkeypatch, policy):
    """Check explicit underflow policy is preserved.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    policy : str
        Input normalization or interpolation policy, supplied by pytest
        parametrization.
    """
    import warnings

    lower = np.ones((1, 1, 1))
    j = np.full((1, 1, 1), 1e-200)
    outcomes = []
    for backend in ("numpy", "numba"):
        monkeypatch.setenv("FISHHIGHZ_FISHER_BACKEND", backend)
        with warnings.catch_warnings(record=True) as caught, np.errstate(under=policy):
            warnings.simplefilter("always")
            try:
                value = fisher_from_factors(j, lower)
                result = value.tolist()
            except FloatingPointError as error:
                result = str(error)
        outcomes.append((result, [str(w.message) for w in caught]))
    assert outcomes[0] == outcomes[1]


def test_unset_backend_selects_compiled_contraction_when_available(monkeypatch):
    """Check unset backend selects compiled contraction when available.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    if importlib.util.find_spec("numba") is None:
        pytest.skip("optional compiler unavailable")
    import fishhighz.kernels._compiled_fisher as compiled

    rng = np.random.default_rng(12)
    raw = rng.normal(size=(9, 3, 3))
    factors = factor_covariance(raw @ raw.swapaxes(1, 2) + 3 * np.eye(3))
    jac = rng.normal(size=(9, 3, 2))
    calls = []
    original = compiled.contract

    def spy(*args):
        """Record compiled-kernel dispatch before invoking the original kernel.

        Parameters
        ----------
        *args : tuple
            Positional arguments forwarded to the original callable or accepted by
            the test callback.

        Returns
        -------
        result : object
            Original contraction result.

        Notes
        -----
        Appends to the enclosing test call log so provider dispatch can be checked.
        """
        calls.append(args)
        return original(*args)

    monkeypatch.setattr(compiled, "contract", spy)
    monkeypatch.delenv("FISHHIGHZ_FISHER_BACKEND", raising=False)
    default = fisher_from_factors(jac, factors)
    assert len(calls) == 1
    monkeypatch.setenv("FISHHIGHZ_FISHER_BACKEND", "numpy")
    reference = fisher_from_factors(jac, factors)
    assert len(calls) == 1
    assert relative(default, reference) <= 5e-15
