"""W12 fixed-reference accuracy routing and narrow trial-contract semantics."""

import configparser
import copy
from collections import OrderedDict
from types import SimpleNamespace

import numpy as np
import pytest

from fishhighz.fields import ObservedField, PairSelection
from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.response import InstrumentResponse
from fishhighz.validation import accuracy, trials
from fishhighz.validation.study import study


def profile_fixture(monkeypatch):
    fields = (
        ObservedField("lya(qso)", "forest", "lya", "qso"),
        ObservedField("lya(lbg)", "forest", "lya", "lbg"),
    )
    selection = PairSelection(fields, [(0, 0), (1, 1)])
    registry = ParameterRegistry(
        [Parameter("ap_0", 1, "target"), Parameter("at_0", 1, "target")]
    )
    binding = BoundParameters(registry, ("ap", "at"), {"ap": "ap_0", "at": "at_0"})
    p3d_calls = []

    def model(theta, z, k, mu, pairs):
        p3d_calls.append(1)
        return np.ones((len(k), len(pairs)))

    p3d = PreparedP3D(
        registry,
        selection,
        [P3DProvider("fixture", model, binding, selection.required_pairs)],
    )
    recipe = accuracy.AccuracyRecipe.__new__(accuracy.AccuracyRecipe)
    recipe.case = "lya_qso_lbg_lae_15x2pt"
    recipe.selection = selection
    recipe.registry = registry
    recipe.h = 0.7
    recipe.provenance = {}
    recipe.weight_method = "inverse_variance"
    recipe._prepared = OrderedDict()
    recipe._partitions = {0: np.array([20.0, 22.0])}
    recipe.config = configparser.ConfigParser()
    recipe.config.read_dict({"survey": {"survey_area": "1000"}})
    recipe.cosmo = SimpleNamespace(
        results=SimpleNamespace(
            hubble_parameter=lambda z: np.full_like(z, 230.0),
            comoving_radial_distance=lambda z: np.full_like(z, 5300.0),
        )
    )
    recipe.model = lambda *args, **kwargs: (p3d, {"z_eval": recipe.z(0)})
    responses = {
        fields[0].id: InstrumentResponse(30.0, 8.0),
        fields[1].id: InstrumentResponse(60.0, 16.0),
    }
    recipe.responses = lambda *args, **kwargs: responses
    rows = {
        field.id: dict(
            z_source=3.0,
            density=np.array([1.0, 2.0]),
            variance=np.array([2.0, 5.0]),
            length_velocity=1000.0,
            density_diagnostics={"field": field.id},
            snr_diagnostics={"field": field.id},
        )
        for field in fields
    }
    magnitudes = np.array([20.5, 21.5])
    measure = np.ones(2)
    recipe.samples = lambda *args, **kwargs: (rows, magnitudes, measure)
    captured = []
    monkeypatch.setattr(
        accuracy, "prepare_bin", lambda spec: captured.append(spec) or spec
    )
    return recipe, captured, p3d_calls


def test_actual_accuracy_prepare_uses_per_field_fixed_reference(monkeypatch):
    recipe, captured, p3d_calls = profile_fixture(monkeypatch)
    controls = dict(
        k_intervals=1,
        mu_order=2,
        z_order=2,
        magnitude_order=2,
        step=2.5e-4,
    )
    _, settings = recipe.prepare(0, controls)
    spec = captured[-1]
    assert p3d_calls == []
    assert settings["forest_weighting"]["method"] == "inverse_variance"
    assert settings["forest_weighting"]["reference"] == trials.FIXED_REFERENCE
    assert settings["forest_weighting"]["iterations"] == {
        "applicable": False,
        "status": "inapplicable",
    }
    aliases = []
    for field_id, source in spec.forests.items():
        assert source.weight_options["method"] == "inverse_variance"
        assert source.auxiliary_coordinates is None
        assert "iterations" not in source.weight_options
        assert "signal" not in source.weight_options
        reference = settings["forest_weighting"]["forests"][field_id]["reference"]
        assert reference["B_star"] == source.weight_options["alias"]
        assert reference["B_star"] == pytest.approx(
            reference["intrinsic_p1d"] * reference["response_factor"] ** 2,
            rel=5e-15,
        )
        aliases.append(reference["B_star"])
    assert aliases[0] != aliases[1]
    with pytest.raises(ValueError, match="iterations are inapplicable"):
        recipe.prepare(0, {**controls, "iterations": 0})

    recipe.weight_method = "legacy"
    legacy, legacy_settings = recipe.prepare(0, {**controls, "iterations": 3})
    for source in legacy.forests.values():
        assert source.weight_options["method"] == "legacy"
        assert source.weight_options["iterations"] == 3
        assert source.auxiliary_coordinates == (2.4, 0.00035)
    assert legacy_settings["forest_weighting"]["method"] == "legacy"
    assert legacy is not spec


class FixedStudy:
    weight_method = "inverse_variance"
    selection = SimpleNamespace(
        fields=[SimpleNamespace(kind="forest")], selected_pairs=np.array([[0, 0]])
    )

    def __init__(self):
        self._prepared = {}

    def evaluate(self, task, controls):
        assert "iterations" not in controls
        fisher = np.diag([2.0, 1.0])
        return dict(fisher=fisher, pair_fisher=fisher[None]), dict(
            settings=dict(
                controls=dict(controls),
                grid=dict(
                    volume=1.0,
                    k_intervals=controls["k_intervals"],
                    mu_order=controls["mu_order"],
                    k_order=4,
                ),
                forest_weighting=dict(
                    method="inverse_variance",
                    reference=trials.FIXED_REFERENCE,
                    iterations={"applicable": False, "status": "inapplicable"},
                ),
            )
        )


class LegacyStudy(FixedStudy):
    weight_method = "legacy"

    def evaluate(self, task, controls):
        fisher = np.diag([2.0, 1.0])
        return dict(fisher=fisher, pair_fisher=fisher[None]), dict(
            settings=dict(
                controls=dict(controls),
                grid=dict(
                    volume=1.0,
                    k_intervals=controls["k_intervals"],
                    mu_order=controls["mu_order"],
                    k_order=4,
                ),
                forest_weighting=dict(
                    method="legacy",
                    reference=None,
                    iterations={"applicable": True, "count": controls["iterations"]},
                ),
            )
        )


def test_fixed_controller_version2_replay_and_narrow_negative_controls():
    arrays, report = study(FixedStudy(), {})
    assert report["metric_names"] == [
        "k",
        "mu",
        "magnitude",
        "volume",
        "step",
        "combined",
    ]
    assert report["trial_contract"]["version"] == 2
    assert report["trial_contract"]["method"] == "inverse_variance"
    assert "weights" not in report["actual_levels"]
    assert all("iterations" not in row for row in report["study_controls"])
    assert trials.validate(arrays, report)

    missing = copy.deepcopy(report)
    missing["metric_names"].remove("magnitude")
    with pytest.raises(ValueError, match="ordered refinement families"):
        trials.validate(arrays, missing)

    legacy_arrays, legacy_report = study(LegacyStudy(), {})
    historical = copy.deepcopy(legacy_report)
    historical["settings"].pop("forest_weighting")
    assert trials.validate(legacy_arrays, historical)
    relabelled = copy.deepcopy(legacy_report)
    relabelled["settings"]["forest_weighting"] = report["settings"]["forest_weighting"]
    with pytest.raises(ValueError, match="version-1.*legacy"):
        trials.validate(legacy_arrays, relabelled)


def test_three_weights_is_explicitly_legacy_only(monkeypatch):
    recipe = accuracy.AccuracyRecipe.__new__(accuracy.AccuracyRecipe)
    recipe.weight_method = "inverse_variance"
    with pytest.raises(ValueError, match="legacy cumulative diagnostic"):
        recipe.sensitivity({}, {}, "three_weights")

    recipe.weight_method = "legacy"
    seen = []
    recipe.evaluate = lambda task, controls, **options: (
        {},
        {"passed": True, "settings": {}, "metrics": []},
    )
    monkeypatch.setattr(
        recipe,
        "evaluate",
        lambda task, controls, **options: (
            seen.append((controls, options)) or {},
            {"passed": True},
        ),
    )
    _, report = recipe.sensitivity({}, {"iterations": 12}, "three_weights")
    assert seen[0][0]["iterations"] == 3
    assert report["sensitivity"]["weight_method"] == "legacy"


def test_cached_legacy_primary_cannot_be_relabelled_fixed(tmp_path):
    from fishhighz.validation.evidence import execute, modern_requests
    from fishhighz.validation.profiles import completed_cache
    from fishhighz.validation.synthetic import convergence_payload

    identity = dict(
        fishhighz=dict(
            origin="fixture",
            module_hashes={"fishhighz/fixture.py": "a" * 64},
        ),
        wheel=dict(
            path=str(tmp_path / "fixture.whl"),
            origin="fixture",
            modules={"fishhighz/fixture.py": "a" * 64},
            sha256="b" * 64,
        ),
        reference=dict(
            reference_origin="/fixture",
            sources={"/fixture/source.py": "c" * 64},
            versions={"fixture": "1"},
        ),
        resources={"fixture": "d" * 64},
    )
    (tmp_path / "fixture.whl").write_bytes(b"fixture")

    def worker(task):
        arrays, report = convergence_payload(task)
        report["provenance"] = identity
        return arrays, report

    root = tmp_path / "legacy"
    execute(
        root,
        suite="quick",
        profiles=("accuracy",),
        worker=worker,
    )
    work = modern_requests("quick", profiles=("accuracy",))
    cache, _ = completed_cache(root, identity, work, accuracy_method="legacy")
    assert len(cache) == 1
    with pytest.raises(ValueError, match="weight method differs"):
        completed_cache(root, identity, work, accuracy_method="inverse_variance")
