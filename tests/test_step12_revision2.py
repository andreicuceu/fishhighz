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
    out = tmp_path / name
    m = execute(out, suite="quick", worker=lambda task: (arrays, report))
    assert m["execution_finished"] and not m["complete"]
    with pytest.raises(ValueError):
        check(out)


@pytest.mark.parametrize(
    "kind", ["synthetic_bao", "synthetic_amplitude", "external_amplitude"]
)
@pytest.mark.parametrize("null", [False, True])
def test_valid_information(tmp_path, kind, null):
    m = execute(
        tmp_path / "bundle",
        suite="quick",
        kind=kind,
        worker=lambda t: evidence_payload(t, null=null),
    )
    assert m["complete"]
    assert check(tmp_path / "bundle")["complete"]
    if null and kind == "synthetic_bao":
        with np.load(tmp_path / "bundle/records-000.npz") as a:
            assert a["constrained"].tolist() == [1, 0]
            assert a["errors"][1] == 0 and a["rank"][0] == 1


def mutate(out, operation):
    p = out / "manifest.json"
    m = json.loads(p.read_text())
    r = m["records"][0]
    with np.load(out / r["arrays"], allow_pickle=False) as d:
        a = {k: d[k] for k in d.files}
    operation(a, r)
    np.savez_compressed(out / r["arrays"], **a)
    r["sha256"] = digest(out / r["arrays"])
    r["inventory"] = {k: list(v.shape) for k, v in a.items()}
    r["effective_hash"] = canonical(r["report"])
    p.write_text(json.dumps(m))


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
    out = tmp_path / "bundle"
    execute(out, suite="quick", kind="synthetic_bao", worker=evidence_payload)

    def edit(a, r):
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
    def worker(t):
        a, r = evidence_payload(t)
        add_metrics(
            a,
            r,
            ["k"],
            [[a["fisher"] * 2, a["fisher"]]],
            [[a["pair_fisher"] * 2, a["pair_fisher"]]],
            [[1, 1]],
        )
        r["passed"] = True
        r["metrics"][0] = {k: 0.0 for k in r["metrics"][0]}
        return a, r

    m = execute(tmp_path / "b", suite="quick", kind="synthetic_bao", worker=worker)
    assert not m["complete"]
    with pytest.raises(ValueError):
        check(tmp_path / "b")


def test_exact_78_records_and_swaps(tmp_path):
    seen = []

    def worker(t):
        seen.append((t["case"], t["bin"], t["profile"]))
        return evidence_payload(t)

    out = tmp_path / "full"
    m = execute(
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
    r0, r1 = m["records"][:2]
    for key in ("arrays", "sha256", "inventory"):
        r0[key], r1[key] = r1[key], r0[key]
    (out / "manifest.json").write_text(json.dumps(m))
    with pytest.raises(ValueError):
        check(out)


def test_legacy_explicit_limit(tmp_path):
    (tmp_path / "manifest.json").write_text(json.dumps(dict(schema=1, complete=True)))
    assert (
        inspect_legacy(tmp_path)["limited"] and not inspect_legacy(tmp_path)["complete"]
    )
    with pytest.raises(ValueError, match="schema 2"):
        check(tmp_path)


def test_roundoff_information_control():
    t = request("lya_qso_2x2pt", 0, "accuracy", kind="synthetic_bao")
    a, r = evidence_payload(t)
    a["fisher"][0, 1] *= 1 + 2 * np.finfo(float).eps
    validate_payload(t, a, r)


def test_literal_backward_derivative(monkeypatch):
    from fishhighz.validation import numerics

    k = np.linspace(0.01, 0.5, 12)
    mu = 0.3
    monkeypatch.setattr(numerics, "legacy_peak", lambda model, k: np.array([k**2]))
    j = numerics.legacy_jacobian(np.ones((1, len(k))), k, mu, widths=[[0, 0]])
    expected = np.zeros(len(k))
    expected[1:] = (k[1:] + k[:-1]) * k[1:]
    np.testing.assert_allclose(j[:, 0, 0], expected * mu**2, rtol=5e-13, atol=0)
    np.testing.assert_allclose(j[:, 0, 1], expected * (1 - mu**2), rtol=5e-13, atol=0)
    assert np.array_equal(j[0], [[0, 0]])


def test_composite_polynomial_and_partition():
    from fishhighz.validation.accuracy import composite

    x, w = composite(np.array([16.0, 19.0, 20.5, 24.0]), 4)
    assert np.all(np.diff(x) > 0) and np.all(w > 0)
    np.testing.assert_allclose(w @ x**3, (24.0**4 - 16.0**4) / 4, rtol=1e-14)


def test_offline_plot_values_and_failed_rows(tmp_path):
    from fishhighz.validation.plots import difference, tables

    def worker(t):
        a, r = evidence_payload(t)
        if t["profile"] == "compatibility":
            a["reference_fisher"] = a["fisher"].copy()
            a["reference_pair_fisher"] = a["pair_fisher"].copy()
        else:
            r["passed"] = False
        return a, r

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
    from fishhighz.validation.attribution import interpolate

    k = np.linspace(0.01, 0.5, 5)
    mu = np.array([0.05, 0.45, 0.95])
    values = (2 * np.tile(k, 3) + 3 * np.repeat(mu, 5))[:, None]
    tk = np.array([0.02, 0.31, 0.49])
    tm = np.array([0.0, 0.5, 1.0])
    np.testing.assert_allclose(
        interpolate(values, k, mu, tk, tm)[:, 0], 2 * tk + 3 * tm, rtol=1e-14
    )


def test_no_false_pass_with_unresolved_control():
    task = request("lya_qso_2x2pt", 0, "accuracy", kind="synthetic_bao")
    a, r = evidence_payload(task)
    r["unresolved_controls"] = [{"control": "weights", "error": "underflow"}]
    with pytest.raises(ValueError, match="unresolved"):
        validate_payload(task, a, r)


def test_arithmetic_mean_swap_preserves_fixed_geometry_contract():
    import configparser
    from types import SimpleNamespace

    from fishhighz.models.templates import prepare_template
    from fishhighz.parameters import Parameter, ParameterRegistry
    from fishhighz.validation.accuracy import AccuracyRecipe
    from fishhighz.validation.cases import bins, recipe, selection

    r = AccuracyRecipe.__new__(AccuracyRecipe)
    r.case = "lbg_lae_3x2pt"
    r.selection = selection(r.case)
    r.config = configparser.ConfigParser()
    r.config.read_dict(recipe(r.case))
    r.tracers = {f.id: {"tracer": f.id} for f in r.selection.fields}
    r.external = SimpleNamespace(
        bias=SimpleNamespace(_get_density_bias=lambda z, n: 2.0)
    )
    r.cosmo = SimpleNamespace(sigma8=0.8)
    r._growth = lambda z: (0.3, 0.95)
    r.registry = ParameterRegistry(
        [Parameter(n, 1.0, "target", step=0.001) for n in ("ap_0", "at_0")]
    )
    k = np.geomspace(0.001, 1.0, 20)
    r.template = prepare_template(
        k,
        10 + np.sin(k * 100),
        np.full(len(k), 10.0),
        z_ref=2.4,
        h_template=0.7,
        h_fid=0.7,
    )
    mean = sum(bins(r.case)[0]) / 2
    provider, settings = r.model(0, mean_z=mean)
    from fishhighz.models.external import evaluate_p3d

    power = evaluate_p3d(
        provider, r.registry.fiducials, r.z(0), np.array([0.1]), np.array([0.5])
    )
    assert power.shape == (1, 3) and settings["z_eval"] == mean
    with pytest.raises(ValueError, match="redshift"):
        evaluate_p3d(
            provider, r.registry.fiducials, mean, np.array([0.1]), np.array([0.5])
        )


def test_real_orchestrator_dispatch_without_reference_imports(tmp_path, monkeypatch):
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
        a, r = evidence_payload(task)
        a["reference_fisher"] = a["fisher"].copy()
        a["reference_pair_fisher"] = a["pair_fisher"].copy()
        r["provenance"] = identity
        r["comparison"] = change(
            a["fisher"], a["fisher"], a["pair_fisher"], a["pair_fisher"], 1, 1
        )
        return a, r

    monkeypatch.setattr(profiles, "compatibility", compatibility)

    class FakeRecipe:
        def __init__(self, root, case, *args):
            cases.append(case)
            self._samples = {}
            self._prepared = {}

        def study(self, task):
            from fishhighz.validation.synthetic import convergence_payload

            a, r = convergence_payload(task)
            r["provenance"] = identity
            return a, r

    monkeypatch.setattr(profiles, "AccuracyRecipe", FakeRecipe)
    m = profiles.run(
        tmp_path / "out",
        reference=reference,
        template=source,
        reference_bundle=reference,
        wheel=wheel,
        suite="full",
        sensitivities=False,
    )
    assert m["complete"] and len(m["records"]) == 78
    assert len(backgrounds) == 1 and len(cases) == 7


def test_full_assignment_gate_rejects_narrowed_inventory(tmp_path):
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
    k = np.geomspace(0.001, 2, 100)
    template = prepare_template(
        k,
        10 + 0.1 * np.sin(k * 100),
        np.full(len(k), 10.0),
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
        k, mu, _ = grid_nodes(s)
        total = np.tile([100.0, 1.0, 100.0], (len(k), 1))
        j = np.broadcast_to(
            np.column_stack((mu**2, 1 - mu**2))[:, None, :], (len(k), 3, 2)
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
    from fishhighz.validation.evidence import record_report

    diagnostic = request(
        "lya_qso_lbg_lae_15x2pt",
        2,
        "accuracy",
        kind="diagnostic",
        diagnostic_id="fixture",
    )

    def worker(task):
        a, r = evidence_payload(task)
        if task["kind"] == "diagnostic":
            r["large_notes"] = "x" * 200000
        return a, r

    root = tmp_path / "bundle"
    m = execute(
        root,
        suite="quick",
        kind="synthetic_bao",
        worker=worker,
        diagnostic_requests=[diagnostic],
    )
    assert m["complete"] and (root / "manifest.json").stat().st_size < 60000
    row = m["diagnostics"][0]
    assert "report" not in row and "report_file" in row
    report = record_report(root, row)
    assert len(report["large_notes"]) == 200000
    check(root)
    report["passed"] = False
    path = root / row["report_file"]
    path.write_text(json.dumps(report))
    row["report_sha256"] = digest(path)
    row["effective_hash"] = canonical(report)
    (root / "manifest.json").write_text(json.dumps(m))
    with pytest.raises(ValueError):
        check(root)
    row["report_file"] = "../outside.json"
    with pytest.raises(ValueError, match="report path"):
        record_report(root, row)


def test_reuse_rejects_scientific_changes_and_preserves_producer(tmp_path):
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
        a, r = evidence_payload(task)
        r["provenance"] = identity
        return a, r

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
