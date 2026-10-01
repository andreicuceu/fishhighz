"""Semantic corruption controls with updated hashes/inventories, not stale-byte tests."""

import json

import numpy as np
import pytest

from fishhighz.validation.evidence import (
    canonical,
    check,
    digest,
    execute,
    inspect_legacy,
)
from fishhighz.validation.schema import add_metrics, request, token, validate_payload
from fishhighz.validation.synthetic import evidence_payload


@pytest.mark.parametrize(
    "name,arrays,report",
    [
        (
            "wrong_dimensions",
            dict(fisher=np.array([7.0]), errors=np.ones(3)),
            {"passed": True},
        ),
        (
            "negative_information",
            dict(fisher=-np.eye(2), errors=np.ones(2)),
            {"passed": True},
        ),
        (
            "wrong_errors",
            dict(fisher=np.eye(2), errors=np.array([30.0, 40.0])),
            {"passed": True},
        ),
        (
            "wrong_pair_payload",
            dict(
                fisher=np.eye(2), errors=np.ones(2), selected_pairs=np.array([[99, 99]])
            ),
            {"passed": True},
        ),
        ("missing_pass", dict(fisher=np.eye(2), errors=np.ones(2)), {}),
    ],
)
def test_original_review_probes(tmp_path, name, arrays, report):
    """Check original review probes.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    name : str
        Named quantity or policy under examination, supplied by pytest
        parametrization.
    arrays : dict
        Named numerical evidence arrays, supplied by pytest parametrization.
    report : dict
        Validation report fixture, supplied by pytest parametrization.
    """
    out = tmp_path / name
    manifest = execute(out, suite="quick", worker=lambda task: (arrays, report))
    assert manifest["execution_finished"] and not manifest["complete"]
    with pytest.raises(ValueError):
        check(out)


@pytest.mark.parametrize(
    "kind", ["synthetic_bao", "synthetic_amplitude", "external_amplitude"]
)
@pytest.mark.parametrize("null", [False, True])
def test_valid_information(tmp_path, kind, null):
    """Check valid information.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    kind : str
        Tracer, input, or calculation classification for this case, supplied by
        pytest parametrization.
    null : bool
        Choice of exact null direction or zero-response case, supplied by pytest
        parametrization.
    """
    manifest = execute(
        tmp_path / "bundle",
        suite="quick",
        kind=kind,
        worker=lambda t: evidence_payload(t, null=null),
    )
    assert manifest["complete"]
    assert check(tmp_path / "bundle")["complete"]
    if null and kind == "synthetic_bao":
        with np.load(tmp_path / "bundle/records-000.npz") as a:
            assert a["constrained"].tolist() == [1, 0]
            assert a["errors"][1] == 0 and a["rank"][0] == 1


def mutate(out, operation):
    """Mutate a synthetic payload and refresh its manifest checksums.

    Parameters
    ----------
    out : pathlib.Path
        Directory containing the synthetic evidence bundle.
    operation : callable
        Mutation receiving the numerical-array dictionary and manifest record.
    """
    manifest_path = out / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    record = manifest["records"][0]
    with np.load(out / record["arrays"], allow_pickle=False) as d:
        arrays = {k: d[k] for k in d.files}
    operation(arrays, record)
    np.savez_compressed(out / record["arrays"], **arrays)
    record["sha256"] = digest(out / record["arrays"])
    record["inventory"] = {k: list(v.shape) for k, v in arrays.items()}
    record["effective_hash"] = canonical(record["report"])
    manifest_path.write_text(json.dumps(manifest))


@pytest.mark.parametrize(
    "change",
    [
        "shape",
        "negative",
        "errors",
        "covariance",
        "correlation",
        "pairs",
        "required",
        "parameter_order",
        "bounds",
        "nodes",
        "alignment",
        "pass_absent",
        "pass_integer",
        "pass_false",
        "nonfinite",
        "token",
    ],
)
def test_rehashed_semantic_mutations(tmp_path, change):
    """Check rehashed semantic mutations.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    change : str
        Input override exercising the specified validation boundary, supplied by
        pytest parametrization.
    """
    out = tmp_path / "bundle"
    execute(out, suite="quick", kind="synthetic_bao", worker=evidence_payload)

    def edit(a, r):
        """Apply the selected semantic corruption to arrays or report metadata.

        Parameters
        ----------
        a : dict of ndarray
            Named numerical arrays in the synthetic evidence payload.
        r : dict
            Manifest record containing the report to be modified in place.
        """
        if change == "shape":
            a["fisher"] = np.ones(2)
        if change == "negative":
            a["fisher"] *= -1
        if change == "errors":
            a["errors"] *= 2
        if change == "covariance":
            a["covariance"] *= 2
        if change == "correlation":
            a["correlation"][0, 1] = 0.9
        if change == "pairs":
            a["selected_pairs"][0] = 99
        if change == "required":
            a["required_pairs"] = a["required_pairs"][::-1]
        if change == "parameter_order":
            r["report"]["settings"]["parameters"].reverse()
            r["report"]["effective_hash"] = canonical(r["report"]["settings"])
            a["effective_token"] = token(r["report"]["settings"])
        if change == "bounds":
            a["bin_bounds"] += 0.1
        if change == "nodes":
            a["k"] = a["k"][::-1]
        if change == "alignment":
            a["observed_j"] = np.roll(a["observed_j"], 1, axis=0)
        if change == "pass_absent":
            del r["report"]["passed"]
        if change == "pass_integer":
            r["report"]["passed"] = 1
        if change == "pass_false":
            r["report"]["passed"] = False
        if change == "nonfinite":
            a["total"][0, 0] = np.nan
        if change == "token":
            a["request_token"][0] ^= 1

    mutate(out, edit)
    with pytest.raises(ValueError):
        check(out)


def test_wrong_convergence_not_hidden(tmp_path):
    """Check wrong convergence not hidden.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """

    def worker(t):
        """Construct synthetic evidence for the requested forecast task.

        Parameters
        ----------
        t : dict
            Synthetic forecast request including case, redshift-bin index, profile,
            and pair selection.

        Returns
        -------
        payload : tuple
            Synthetic numerical arrays and validation report, including the
            requested test modification.
        """
        evidence_arrays, evidence_report = evidence_payload(t)
        add_metrics(
            evidence_arrays,
            evidence_report,
            ["k"],
            [[evidence_arrays["fisher"] * 2, evidence_arrays["fisher"]]],
            [[evidence_arrays["pair_fisher"] * 2, evidence_arrays["pair_fisher"]]],
            [[1, 1]],
        )
        evidence_report["passed"] = True
        evidence_report["metrics"][0] = {k: 0.0 for k in evidence_report["metrics"][0]}
        return evidence_arrays, evidence_report

    manifest = execute(
        tmp_path / "b", suite="quick", kind="synthetic_bao", worker=worker
    )
    assert not manifest["complete"]
    with pytest.raises(ValueError):
        check(tmp_path / "b")


def test_exact_78_records_and_swaps(tmp_path):
    """Check exact 78 records and swaps.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    seen = []

    def worker(t):
        """Construct synthetic evidence for the requested forecast task.

        Parameters
        ----------
        t : dict
            Synthetic forecast request including case, redshift-bin index, profile,
            and pair selection.

        Returns
        -------
        payload : tuple
            Synthetic numerical arrays and validation report, including the
            requested test modification.

        Notes
        -----
        Appends to the enclosing test call log so provider dispatch can be checked.
        """
        seen.append((t["case"], t["bin"], t["profile"]))
        return evidence_payload(t)

    out = tmp_path / "full"
    manifest = execute(
        out,
        suite="full",
        profiles=("compatibility", "accuracy"),
        kind="synthetic_bao",
        worker=worker,
    )
    assert len(seen) == 78 and len(set(seen)) == 78
    assert len({x[0] for x in seen}) == 7
    assert check(out)["complete"]
    # Consistent file/hash swapping of plausible same-shaped profiles still fails token binding.
    r0, r1 = manifest["records"][:2]
    for key in ("arrays", "sha256", "inventory"):
        r0[key], r1[key] = r1[key], r0[key]
    (out / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        check(out)


def test_legacy_explicit_limit(tmp_path):
    """Check legacy explicit limit.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    (tmp_path / "manifest.json").write_text(json.dumps(dict(schema=1, complete=True)))
    assert (
        inspect_legacy(tmp_path)["limited"] and not inspect_legacy(tmp_path)["complete"]
    )
    with pytest.raises(ValueError, match="schema 2"):
        check(tmp_path)


def test_roundoff_information_control():
    """Check roundoff information control."""
    task = request("lya_qso_2x2pt", 0, "accuracy", kind="synthetic_bao")
    evidence_arrays, evidence_report = evidence_payload(task)
    evidence_arrays["fisher"][0, 1] *= 1 + 2 * np.finfo(float).eps
    validate_payload(task, evidence_arrays, evidence_report)


def test_literal_backward_derivative(monkeypatch):
    """Check literal backward derivative.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    from fishhighz.validation import numerics

    k_grid = np.linspace(0.01, 0.5, 12)
    direction_cosine = 0.3
    monkeypatch.setattr(numerics, "legacy_peak", lambda model, k: np.array([k**2]))
    j = numerics.legacy_jacobian(
        np.ones((1, len(k_grid))), k_grid, direction_cosine, widths=[[0, 0]]
    )
    expected = np.zeros(len(k_grid))
    expected[1:] = (k_grid[1:] + k_grid[:-1]) * k_grid[1:]
    np.testing.assert_allclose(
        j[:, 0, 0], expected * direction_cosine**2, rtol=5e-13, atol=0
    )
    np.testing.assert_allclose(
        j[:, 0, 1], expected * (1 - direction_cosine**2), rtol=5e-13, atol=0
    )
    assert np.array_equal(j[0], [[0, 0]])


def test_composite_polynomial_and_partition():
    """Check composite polynomial and partition."""
    from fishhighz.validation.accuracy import composite

    magnitude_nodes, quadrature_weights = composite(
        np.array([16.0, 19.0, 20.5, 24.0]), 4
    )
    assert np.all(np.diff(magnitude_nodes) > 0) and np.all(quadrature_weights > 0)
    np.testing.assert_allclose(
        quadrature_weights @ magnitude_nodes**3, (24.0**4 - 16.0**4) / 4, rtol=1e-14
    )


def test_offline_plot_values_and_failed_rows(tmp_path):
    """Check offline plot values and failed rows.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    from fishhighz.validation.plots import difference, tables

    def worker(t):
        """Construct synthetic evidence for the requested forecast task.

        Parameters
        ----------
        t : dict
            Synthetic forecast request including case, redshift-bin index, profile,
            and pair selection.

        Returns
        -------
        payload : tuple
            Synthetic numerical arrays and validation report, including the
            requested test modification.
        """
        evidence_arrays, evidence_report = evidence_payload(t)
        if t["profile"] == "compatibility":
            evidence_arrays["reference_fisher"] = evidence_arrays["fisher"].copy()
            evidence_arrays["reference_pair_fisher"] = evidence_arrays[
                "pair_fisher"
            ].copy()
        else:
            evidence_report["passed"] = False
        return evidence_arrays, evidence_report

    execute(
        tmp_path / "source",
        suite="quick",
        profiles=("compatibility", "accuracy"),
        kind="synthetic_bao",
        worker=worker,
    )
    table = tables(tmp_path / "source", tmp_path / "plots")
    assert len(table["rows"]) == 1
    row = table["rows"][0]
    assert len(row["pairs"]) == 15
    assert not row["profiles"]["accuracy"]["passed"]
    assert row["differences"]["compatibility/reference"] == [0.0, 0.0, 0.0]
    assert difference([2.0, 4.0, -0.5], [1.0, 2.0, -0.4]) == [
        100.0,
        100.0,
        -0.09999999999999998,
    ]
    assert difference([None, 1.0, None], [1.0, 2.0, 0.2]) == [None] * 3


def test_linear_common_node_interpolation():
    """Check linear common node interpolation."""
    from fishhighz.validation.attribution import interpolate

    k_grid = np.linspace(0.01, 0.5, 5)
    mu_grid = np.array([0.05, 0.45, 0.95])
    values = (2 * np.tile(k_grid, 3) + 3 * np.repeat(mu_grid, 5))[:, None]
    tk = np.array([0.02, 0.31, 0.49])
    tm = np.array([0.0, 0.5, 1.0])
    np.testing.assert_allclose(
        interpolate(values, k_grid, mu_grid, tk, tm)[:, 0], 2 * tk + 3 * tm, rtol=1e-14
    )


def test_no_false_pass_with_unresolved_control():
    """Check no false pass with unresolved control."""
    task = request("lya_qso_2x2pt", 0, "accuracy", kind="synthetic_bao")
    evidence_arrays, evidence_report = evidence_payload(task)
    evidence_report["unresolved_controls"] = [
        {"control": "weights", "error": "underflow"}
    ]
    with pytest.raises(ValueError, match="unresolved"):
        validate_payload(task, evidence_arrays, evidence_report)


def test_arithmetic_mean_swap_preserves_fixed_geometry_contract():
    """Check arithmetic mean swap preserves fixed geometry contract."""
    import configparser
    from types import SimpleNamespace

    from fishhighz.models.templates import prepare_template
    from fishhighz.parameters import Parameter, ParameterRegistry
    from fishhighz.validation.accuracy import AccuracyRecipe
    from fishhighz.validation.cases import bins, recipe, selection

    accuracy_recipe = AccuracyRecipe.__new__(AccuracyRecipe)
    accuracy_recipe.case = "lbg_lae_3x2pt"
    accuracy_recipe.selection = selection(accuracy_recipe.case)
    accuracy_recipe.config = configparser.ConfigParser()
    accuracy_recipe.config.read_dict(recipe(accuracy_recipe.case))
    accuracy_recipe.tracers = {
        f.id: {"tracer": f.id} for f in accuracy_recipe.selection.fields
    }
    accuracy_recipe.external = SimpleNamespace(
        bias=SimpleNamespace(_get_density_bias=lambda z, n: 2.0)
    )
    accuracy_recipe.cosmo = SimpleNamespace(sigma8=0.8)
    accuracy_recipe._growth = lambda z: (0.3, 0.95)
    accuracy_recipe.registry = ParameterRegistry(
        [Parameter(n, 1.0, "target", step=0.001) for n in ("ap_0", "at_0")]
    )
    k_grid = np.geomspace(0.001, 1.0, 20)
    accuracy_recipe.template = prepare_template(
        k_grid,
        10 + np.sin(k_grid * 100),
        np.full(len(k_grid), 10.0),
        z_ref=2.4,
        h_template=0.7,
        h_fid=0.7,
    )
    mean = sum(bins(accuracy_recipe.case)[0]) / 2
    provider, settings = accuracy_recipe.model(0, mean_z=mean)
    from fishhighz.models.external import evaluate_p3d

    power = evaluate_p3d(
        provider,
        accuracy_recipe.registry.fiducials,
        accuracy_recipe.z(0),
        np.array([0.1]),
        np.array([0.5]),
    )
    assert power.shape == (1, 3) and settings["z_eval"] == mean
    with pytest.raises(ValueError, match="redshift"):
        evaluate_p3d(
            provider,
            accuracy_recipe.registry.fiducials,
            mean,
            np.array([0.1]),
            np.array([0.5]),
        )


def test_real_orchestrator_dispatch_without_reference_imports(tmp_path, monkeypatch):
    """Check real orchestrator dispatch without reference imports.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    from fishhighz.validation import profiles
    from fishhighz.validation.numerics import change

    source = tmp_path / "source.py"
    source.write_text("synthetic test provenance")
    wheel = tmp_path / "fixture.whl"
    wheel.write_bytes(b"synthetic test identity")
    reference = tmp_path / "reference"
    reference.mkdir()
    (reference / "manifest.json").write_text("{}")
    identity = dict(
        fishhighz=dict(origin="fixture", module_hashes={"fixture.py": "a" * 64}),
        wheel=dict(origin="fixture", modules={"fixture.py": "a" * 64}, sha256="b" * 64),
        reference=dict(
            reference_origin=str(tmp_path),
            sources={str(source): "c" * 64},
            versions={"fixture": "1"},
        ),
        resources={str(source): "c" * 64},
    )
    monkeypatch.setattr(profiles, "verify_inventory", lambda p: None)
    monkeypatch.setattr(profiles, "provenance", lambda *args: identity)
    backgrounds = []
    cases = []
    monkeypatch.setattr(
        profiles,
        "background",
        lambda *args: (backgrounds.append(True) or object(), object()),
    )

    def compatibility(task, *args):
        """Return synthetic compatibility evidence with identical reference Fisher data.

        Parameters
        ----------
        task : dict
            Synthetic forecast request including case, redshift-bin index, profile,
            and pair selection.
        *args : tuple
            Positional arguments forwarded to the original callable or accepted by
            the test callback.

        Returns
        -------
        payload : tuple
            Numerical arrays and report with provenance and comparison metrics.
        """
        evidence_arrays, evidence_report = evidence_payload(task)
        evidence_arrays["reference_fisher"] = evidence_arrays["fisher"].copy()
        evidence_arrays["reference_pair_fisher"] = evidence_arrays["pair_fisher"].copy()
        evidence_report["provenance"] = identity
        evidence_report["comparison"] = change(
            evidence_arrays["fisher"],
            evidence_arrays["fisher"],
            evidence_arrays["pair_fisher"],
            evidence_arrays["pair_fisher"],
            1,
            1,
        )
        return evidence_arrays, evidence_report

    monkeypatch.setattr(profiles, "compatibility", compatibility)

    class FakeRecipe:
        def __init__(self, root, case, *args, **kwargs):
            """Initialize the synthetic FakeRecipe fixture.

            Parameters
            ----------
            root : pathlib.Path
                Directory containing the synthetic input files.
            case : str
                Authoritative forecast-case identifier.
            *args : tuple
                Positional arguments forwarded to the original callable or accepted by
                the test callback.
            **kwargs : dict
                Keyword options forwarded to the original callable or inspected by the
                test callback.

            Notes
            -----
            Sets the instance state used by the enclosing test; no scientific calculation is run.
            """
            cases.append(case)
            self._samples = {}
            self._prepared = {}

        def study(self, task):
            """Return synthetic convergence evidence with explicit provenance.

            Parameters
            ----------
            task : dict
                Synthetic forecast request including case, redshift-bin index, profile,
                and pair selection.

            Returns
            -------
            payload : tuple
                Convergence arrays and report for the supplied task.
            """
            from fishhighz.validation.synthetic import convergence_payload

            evidence_arrays, evidence_report = convergence_payload(task)
            evidence_report["provenance"] = identity
            return evidence_arrays, evidence_report

    monkeypatch.setattr(profiles, "AccuracyRecipe", FakeRecipe)
    manifest = profiles.run(
        tmp_path / "out",
        reference=reference,
        template=source,
        reference_bundle=reference,
        wheel=wheel,
        suite="full",
        sensitivities=False,
        profiles=("compatibility", "accuracy"),
        accuracy_method="legacy",
        recipe_revision=None,
    )
    assert manifest["complete"] and len(manifest["records"]) == 78
    assert len(backgrounds) == 1 and len(cases) == 7


def test_full_assignment_gate_rejects_narrowed_inventory(tmp_path):
    """Check full assignment gate rejects narrowed inventory.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    import runpy
    from pathlib import Path

    from fishhighz.validation.evidence import modern_requests

    gate = runpy.run_path(
        str(Path(__file__).resolve().parents[1] / "scripts/check_desi2_full.py")
    )["full_gate"]
    (tmp_path / "manifest.json").write_text(
        json.dumps(
            dict(
                requested=modern_requests(
                    "full",
                    ["lbg_lae_3x2pt"],
                    [0],
                    profiles=("compatibility", "accuracy"),
                ),
                diagnostics_requested=[],
                complete=True,
            )
        )
    )
    with pytest.raises(ValueError, match="78 primary"):
        gate(tmp_path)


def test_controlled_attribution_chain_endpoints_without_assets(tmp_path, monkeypatch):
    """Check controlled attribution chain endpoints without assets.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    import builtins
    import copy
    import runpy
    from pathlib import Path
    from types import SimpleNamespace

    from fishhighz.models.templates import prepare_template
    from fishhighz.validation.schema import assemble, grid_nodes

    pytest.importorskip("scipy")
    original_import = builtins.__import__

    def no_models(name, *args, **kwargs):
        """Reject scientific model imports during offline attribution tests.

        Parameters
        ----------
        name : str
            Name of the artifact, module, or result under examination.
        *args : tuple
            Positional arguments forwarded to the original callable or accepted by
            the test callback.
        **kwargs : dict
            Keyword options forwarded to the original callable or inspected by the
            test callback.

        Returns
        -------
        module : module
            Imported module when it is outside the prohibited model families.

        Raises
        ------
        AssertionError
            Deliberately raised to exercise the rejection path in the enclosing
            test.
        """
        if name.startswith("fishhighz.models") or name.split(".")[0] in (
            "camb",
            "lyaforecast",
        ):
            raise AssertionError("offline attribution imported model code")
        return original_import(name, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(builtins, "__import__", no_models)
        module = runpy.run_path(
            str(Path(__file__).resolve().parents[1] / "scripts/attribute_desi2.py")
        )
    attribute = module["attribution"]
    case = "lbg_lae_3x2pt"
    task = request(case, 0, "accuracy", kind="synthetic_bao")
    k_grid = np.geomspace(0.001, 2, 100)
    template = prepare_template(
        k_grid,
        10 + 0.1 * np.sin(k_grid * 100),
        np.full(len(k_grid), 10.0),
        z_ref=2.4,
        h_template=0.7,
        h_fid=0.7,
    )
    settings = dict(
        profile="accuracy",
        parameters=task["parameters"],
        bounds=task["bounds"],
        fields=[f["id"] for f in task["fields"]],
        grid=dict(
            kind="gauss_legendre",
            k_intervals=4,
            k_order=4,
            mu_order=4,
            h_fid=0.7,
            volume=1e8,
        ),
        model=dict(
            z_eval=2.373040171714532,
            biases={"lbg": 2.0, "lae": 2.0},
            betas={},
            widths={"lbg": [2.0, 1.0], "lae": [2.0, 1.0]},
            f=0.9,
            G=1.0,
        ),
        geometry=dict(a_v=100.0, d_deg=70.0),
        samples={},
    )

    def payload(t, s):
        """Construct a three-spectrum angular Jacobian on the requested grid.

        Parameters
        ----------
        t : dict
            Synthetic forecast request including case, redshift-bin index, profile,
            and pair selection.
        s : dict
            Synthetic grid, selection, and numerical settings.

        Returns
        -------
        payload : tuple
            Assembled numerical arrays and evidence report.
        """
        k_grid, mu_grid, _ = grid_nodes(s)
        total = np.tile([100.0, 1.0, 100.0], (len(k_grid), 1))
        j = np.broadcast_to(
            np.column_stack((mu_grid**2, 1 - mu_grid**2))[:, None, :],
            (len(k_grid), 3, 2),
        ).copy()
        return assemble(t, total, j, s)

    accuracy, report = payload(task, settings)
    accuracy["noise"] = np.tile([10.0, 0.0, 10.0], (len(accuracy["k"]), 1))
    legacy_task = request(case, 0, "compatibility", kind="synthetic_bao")
    legacy, _ = payload(
        legacy_task,
        {
            **settings,
            "profile": "compatibility",
            "grid": dict(kind="legacy", volume=1e8),
        },
    )
    report["legacy_volume"] = 1e8
    report["final_controls"] = {"step": 0.00025}
    report["legacy_pair_inputs"] = {
        n + "_" + n: dict(
            _z_mean=2.373040171714532,
            _distance_to_velocity=100.0,
            _angle_to_distance=70.0,
        )
        for n in ["lbg", "lae"]
    }
    adapters = {
        n: SimpleNamespace(
            sample=lambda z, m: dict(
                raw=np.full(len(m), 1e5),
                provenance=dict(masks=dict(density_floor=np.zeros(len(m), dtype=bool))),
            )
        )
        for n in ["lbg", "lae"]
    }
    arrays, result = attribute(legacy, accuracy, report, task, template, adapters)
    assert len(result["stages"]) == 8 and result["endpoint_relative"] < 5e-12
    np.testing.assert_allclose(arrays["fisher_0"], legacy["fisher"], rtol=5e-12, atol=0)
    np.testing.assert_allclose(
        arrays["fisher_7"], accuracy["fisher"], rtol=5e-12, atol=0
    )
    path = tmp_path / "row.npz"
    np.savez_compressed(path, **arrays)
    result.update(array=path.name, sha256=digest(path))
    with monkeypatch.context() as patch:
        patch.setattr(builtins, "__import__", no_models)
        module["verify_saved"](tmp_path, {"records": [result]})
    changed = copy.deepcopy(result)
    changed["cross_rule"]["reference_values"][0] *= 2
    with pytest.raises(AssertionError):
        module["verify_saved"](tmp_path, {"records": [changed]})
    changed = copy.deepcopy(result)
    changed["negative_density_floor"]["change"]["error_relative"] = 0.5
    with pytest.raises(ValueError, match="branch difference"):
        module["verify_saved"](tmp_path, {"records": [changed]})
    changed = copy.deepcopy(result)
    changed["stages"][1]["fisher_relative_previous"] = 0.5
    with pytest.raises(ValueError, match="stage difference"):
        module["verify_saved"](tmp_path, {"records": [changed]})
    arrays["errors_6"] *= 2
    np.savez_compressed(path, **arrays)
    result["sha256"] = digest(path)
    with pytest.raises(ValueError, match="derived quantities"):
        module["verify_saved"](tmp_path, {"records": [result]})


def test_external_reports_are_bounded_and_semantically_checked(tmp_path):
    """Check external reports are bounded and semantically checked.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    from fishhighz.validation.evidence import record_report

    diagnostic = request(
        "lya_qso_lbg_lae_15x2pt",
        2,
        "accuracy",
        kind="diagnostic",
        diagnostic_id="fixture",
    )

    def worker(task):
        """Construct synthetic evidence for the requested forecast task.

        Parameters
        ----------
        task : dict
            Synthetic forecast request including case, redshift-bin index, profile,
            and pair selection.

        Returns
        -------
        payload : tuple
            Synthetic numerical arrays and validation report, including the
            requested test modification.
        """
        evidence_arrays, evidence_report = evidence_payload(task)
        if task["kind"] == "diagnostic":
            evidence_report["large_notes"] = "x" * 200000
        return evidence_arrays, evidence_report

    root = tmp_path / "bundle"
    manifest = execute(
        root,
        suite="quick",
        kind="synthetic_bao",
        worker=worker,
        diagnostic_requests=[diagnostic],
    )
    assert manifest["complete"] and (root / "manifest.json").stat().st_size < 60000
    row = manifest["diagnostics"][0]
    assert "report" not in row and "report_file" in row
    report = record_report(root, row)
    assert len(report["large_notes"]) == 200000
    check(root)
    report["passed"] = False
    path = root / row["report_file"]
    path.write_text(json.dumps(report))
    row["report_sha256"] = digest(path)
    row["effective_hash"] = canonical(report)
    (root / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        check(root)
    row["report_file"] = "../outside.json"
    with pytest.raises(ValueError, match="report path"):
        record_report(root, row)


def test_reuse_rejects_scientific_changes_and_preserves_producer(tmp_path):
    """Check reuse rejects scientific changes and preserves producer.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    import copy

    from fishhighz.validation.evidence import modern_requests
    from fishhighz.validation.profiles import completed_cache, reassemble_cached

    identity = dict(
        fishhighz=dict(
            module_hashes={
                "fishhighz/kernels/weights.py": "a" * 64,
                "fishhighz/validation/evidence.py": "b" * 64,
            }
        ),
        reference={"versions": {"fixture": "1"}},
        resources={"fixture": "c" * 64},
        wheel={"path": str(tmp_path / "fixture.whl")},
    )
    (tmp_path / "fixture.whl").write_bytes(b"fixture")

    def worker(task):
        """Construct synthetic evidence for the requested forecast task.

        Parameters
        ----------
        task : dict
            Synthetic forecast request including case, redshift-bin index, profile,
            and pair selection.

        Returns
        -------
        payload : tuple
            Synthetic numerical arrays and validation report, including the
            requested test modification.
        """
        evidence_arrays, evidence_report = evidence_payload(task)
        evidence_report["provenance"] = identity
        return evidence_arrays, evidence_report

    root = tmp_path / "old"
    execute(root, suite="quick", kind="synthetic_bao", worker=worker)
    original = (root / "manifest.json").read_bytes()
    work = modern_requests("quick", kind="synthetic_bao")
    changed = copy.deepcopy(identity)
    changed["fishhighz"]["module_hashes"]["fishhighz/kernels/weights.py"] = "d" * 64
    with pytest.raises(ValueError, match="scientific code changed"):
        completed_cache(root, changed, work)
    current = copy.deepcopy(identity)
    current["fishhighz"]["module_hashes"]["fishhighz/validation/evidence.py"] = "d" * 64
    cache, _ = completed_cache(root, current, work)
    arrays, report = reassemble_cached(next(iter(cache.values())), current)
    validate_payload(work[0], arrays, report)
    assert report["cached_numerical_inputs"]["producer"] == identity
    assert report["provenance"] == current
    assert (root / "manifest.json").read_bytes() == original
