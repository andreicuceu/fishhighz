"""Explicit two-profile real orchestration with separate sensitivity inventory."""

import importlib.metadata
import importlib.resources
import json
import time
import zipfile
from pathlib import Path

import numpy as np

from .accuracy import AccuracyRecipe, background
from .cases import CASE_IDS, recipe, verify_inventory
from .compatibility import run as compatibility
from .evidence import canonical, digest, execute, modern_requests, record_report
from .profile_definitions import ADAPTIVE, REVISION, identity
from .reference_capture import imported_reference, resolved_resources
from .schema import request, token

POLICIES = (
    "floor_low",
    "floor_high",
    "remove_bright",
    "remove_sentinel",
    "width_minus",
    "width_plus",
    "growth_eds",
    "legacy_resolution",
    "rectangular_magnitude",
    "no_reconstruction",
    "arithmetic_mean",
    "three_weights",
)


def wheel_identity(path):
    """Verify actual imported payload against one exact caller-selected wheel."""
    import fishhighz

    wheel = Path(path).resolve()
    root = importlib.resources.files("fishhighz")
    modules = {}
    with zipfile.ZipFile(wheel) as archive:
        for name in archive.namelist():
            if name.startswith("fishhighz/") and name.endswith(".py"):
                current = root.joinpath(name.split("/", 1)[1]).read_bytes()
                if current != archive.read(name):
                    raise ValueError(
                        f"imported FishHighz differs from exact wheel: {name}"
                    )
                import hashlib

                modules[name] = hashlib.sha256(current).hexdigest()
    return dict(
        path=str(wheel),
        sha256=digest(wheel),
        origin=fishhighz.__file__,
        version=importlib.metadata.version("fishhighz"),
        modules=modules,
    )


def provenance(reference, template, wheel):
    import configparser

    identity = imported_reference(reference)
    resources = {}
    for case in CASE_IDS:
        c = configparser.ConfigParser()
        c.read_dict(recipe(case))
        resources.update(resolved_resources(reference, c))
    resources[str(Path(template).resolve())] = digest(template)
    installed = wheel_identity(wheel)
    return dict(
        fishhighz=dict(
            origin=installed["origin"],
            version=installed["version"],
            module_hashes=installed["modules"],
        ),
        wheel=installed,
        reference=identity,
        resources=resources,
    )


# Only report I/O/controller files may differ when reusing immutable numerical
# inputs after an interrupted run. Schema, numerical helpers and every scientific
# module must be byte-identical, with identical actual reference dependencies.
_IO_MODULES = {
    "fishhighz/validation/evidence.py",
    "fishhighz/validation/profiles.py",
    "fishhighz/validation/plots.py",
}


def completed_cache(root, ident, work, *, accuracy_method=None):
    from .schema import validate_payload

    root = Path(root).resolve()
    path = root / "manifest.json"
    m = json.loads(path.read_text())
    if m.get("schema") not in (2, 3):
        raise ValueError("reuse requires schema-2 numerical inputs")
    requested = {canonical(t) for t in work}
    cache = {}
    inputs = [path]
    for record in m["records"]:
        key = canonical(record["task"])
        if key not in requested or record.get("status") != "completed":
            continue
        report = record_report(root, record)
        if record["task"]["profile"] == "accuracy" and accuracy_method is not None:
            weighting = report.get("settings", {}).get("forest_weighting")
            recorded_method = "legacy" if weighting is None else weighting.get("method")
            if recorded_method != accuracy_method:
                raise ValueError(
                    "cached accuracy weight method differs from the requested profile"
                )
        if record["task"].get("recipe_revision"):
            expected = identity(
                record["task"]["profile"],
                accuracy_method if record["task"]["profile"] == "accuracy" else None,
            )
            if report.get("settings", {}).get("recipe_identity") != expected:
                raise ValueError(
                    "cached recipe revision/mode/stopping/selection differs"
                )
        old = report["provenance"]
        before = old["fishhighz"]["module_hashes"]
        after = ident["fishhighz"]["module_hashes"]
        if set(before) != set(after) or any(
            before[n] != after[n] for n in before if n not in _IO_MODULES
        ):
            raise ValueError(
                "scientific code changed; cached numerical inputs cannot be reused"
            )
        if (
            old["reference"] != ident["reference"]
            or old["resources"] != ident["resources"]
        ):
            raise ValueError("reference interpreter/source/resources changed; no reuse")
        data = (root / record["arrays"]).resolve()
        if (
            data.parent != root
            or digest(data) != record["sha256"]
            or canonical(report) != record["effective_hash"]
        ):
            raise ValueError("cached input hash/path mismatch")
        # Validate one record at a time; no model imports or retained node arrays.
        with np.load(data, allow_pickle=False) as f:
            arrays = {n: f[n] for n in f.files}
        validate_payload(record["task"], arrays, report, require_pass=False)
        cache[key] = dict(path=data, record=record, root=root, report=report)
        inputs.extend([data, Path(old["wheel"]["path"])])
        if "report_file" in record:
            inputs.append(root / record["report_file"])
    return cache, inputs


def reassemble_cached(entry, ident):
    from .numerics import contract, summaries, wick

    record = entry["record"]
    report = entry["report"]
    task = record["task"]
    with np.load(entry["path"], allow_pickle=False) as f:
        arrays = {n: f[n] for n in f.files}
    c = wick(
        arrays["total"],
        arrays["modes"],
        arrays["required_pairs"],
        arrays["selected_pairs"],
        len(task["fields"]),
    )
    fisher, pairs = contract(c, arrays["observed_j"])
    arrays.update(selected_covariance=c, **summaries(fisher, pairs))
    report["cached_numerical_inputs"] = dict(
        path=str(entry["path"]),
        sha256=record["sha256"],
        producer=report["provenance"],
        byte_identical_scientific_modules={
            n: s
            for n, s in ident["fishhighz"]["module_hashes"].items()
            if n not in _IO_MODULES
        },
        interpretation="Immutable powers/Jacobians/convergence inputs reused after report-I/O repair; C/F reassembled with the current exact wheel",
    )
    report["provenance"] = ident
    print(
        "Reused byte-identical scientific inputs; independently reassembled C/F",
        flush=True,
    )
    return arrays, report


def run(
    output,
    *,
    reference,
    template,
    reference_bundle,
    wheel,
    suite="quick",
    cases=None,
    bin_indices=None,
    profiles=("full-compatibility", "fixed-compatibility", "accuracy"),
    sensitivities=False,
    reuse_completed=None,
    accuracy_method="early_lyaforecast",
    compatibility_bundle=None,
    recipe_revision=REVISION,
):
    """Serial full is explicit; primary and diagnostic inventories are disjoint."""
    if accuracy_method not in ("inverse_variance", "legacy", *ADAPTIVE):
        raise ValueError(
            "accuracy_method must be legacy, inverse_variance, "
            "early_lyaforecast or mcdonald"
        )
    if recipe_revision:
        if cases is None:
            cases = ["lya_qso_lbg_lae_15x2pt"]
        if list(cases) != ["lya_qso_lbg_lae_15x2pt"]:
            raise ValueError("revised profiles are bounded to the DESI-2 15x2pt case")
        profiles = tuple(
            "full-compatibility" if p == "compatibility" else p for p in profiles
        )
    work = modern_requests(
        suite, cases, bin_indices, profiles=profiles, recipe_revision=recipe_revision
    )
    verify_inventory(Path(reference) / "examples/desi2")
    ident = provenance(reference, template, wheel)
    primary_rows = {}
    cached = {}
    cache_inputs = []
    if reuse_completed is not None:
        cached, cache_inputs = completed_cache(
            reuse_completed, ident, work, accuracy_method=accuracy_method
        )
    diagnostics = []
    if sensitivities and "accuracy" in profiles:
        diagnostics = [
            request(
                t["case"],
                t["bin"],
                "accuracy",
                kind="diagnostic",
                diagnostic_id=p,
                recipe_revision=recipe_revision,
            )
            for t in work
            if t["profile"] == "accuracy"
            for p in POLICIES
            if accuracy_method == "legacy" or p != "three_weights"
        ]
    c = t = None
    active_recipe = None
    active_case = None
    out = Path(output).resolve()

    def accuracy_recipe(case):
        nonlocal c, t, active_recipe, active_case
        if c is None:
            print(
                "Preparing shared CAMB background at all actual redshifts", flush=True
            )
            c, t = background(reference, template)
        if active_case != case:
            arguments = (reference, case, c, t, ident)
            active_recipe = AccuracyRecipe(
                *arguments,
                weight_method=accuracy_method,
                recipe_revision=recipe_revision,
            )
            active_case = case
        return active_recipe

    def worker(task):
        started = time.monotonic()
        print(
            f"{task['kind']} {task['case']} bin {task['bin']} {task['profile']} {task['diagnostic_id'] or ''}",
            flush=True,
        )
        key = canonical(task)
        if key in cached:
            arrays, report = reassemble_cached(cached.pop(key), ident)
            if task["profile"] == "accuracy":
                index = work.index(task)
                primary_rows[(task["case"], task["bin"])] = dict(
                    controls=report["final_controls"],
                    fisher=arrays["fisher"],
                    pair_fisher=arrays["pair_fisher"],
                    array=f"records-{index:03d}.npz",
                    index=index,
                )
        elif task["kind"] == "diagnostic":
            key = (task["case"], task["bin"])
            info = primary_rows.get(key)
            if info is None:
                raise ValueError(
                    "primary accuracy execution failed; diagnostic unavailable"
                )
            r = accuracy_recipe(task["case"])
            policy = task["diagnostic_id"]
            sampled, _, _ = r.samples(task["bin"], info["controls"]["magnitude_order"])
            if policy == "remove_bright":
                active = any(
                    row.get("snr_diagnostics", {})
                    .get("counts", {})
                    .get("bright_clamp", 0)
                    > 0
                    for row in sampled.values()
                )
            elif policy == "remove_sentinel":
                active = any(
                    row.get("snr_diagnostics", {})
                    .get("counts", {})
                    .get("out_of_range", 0)
                    > 0
                    for row in sampled.values()
                )
            elif policy.startswith("floor_"):
                active = any(
                    any(
                        row["density_diagnostics"]["counts"][k] > 0
                        for k in ("density_floor", "negative_density")
                    )
                    for row in sampled.values()
                )
            else:
                active = True
            if active:
                arrays, report = r.sensitivity(task, info["controls"], policy)
            else:
                with np.load(out / info["array"], allow_pickle=False) as d:
                    arrays = {
                        k: d[k]
                        for k in d.files
                        if not k.startswith(("metric_", "study_"))
                    }
                manifest = json.loads((out / "manifest.json").read_text())
                report = record_report(out, manifest["records"][info["index"]])
                report.pop("unresolved_controls", None)
                report.update(
                    context=task,
                    passed=True,
                    metric_names=[],
                    metrics=[],
                    sensitivity=dict(policy=policy, active=False, zero_mask_proof=True),
                )
                arrays["request_token"] = token(task)
            report["sensitivity"]["active"] = active
            from .numerics import change

            report["sensitivity"]["change"] = change(
                arrays["fisher"],
                info["fisher"],
                arrays["pair_fisher"],
                info["pair_fisher"],
                1,
                1,
            )
        elif task["profile"] in ("full-compatibility", "fixed-compatibility"):
            from .revised_compatibility import run as selected_compatibility

            if compatibility_bundle is None:
                raise ValueError(
                    "new compatibility profiles require captured compatibility_bundle"
                )
            arrays, report = selected_compatibility(task, compatibility_bundle, ident)
        elif task["profile"] == "compatibility":
            arrays, report = compatibility(task, reference_bundle, ident)
        else:
            r = accuracy_recipe(task["case"])
            arrays, report = r.study(task)
            index = work.index(task)
            primary_rows[(task["case"], task["bin"])] = dict(
                controls=report["final_controls"],
                fisher=arrays["fisher"],
                pair_fisher=arrays["pair_fisher"],
                array=f"records-{index:03d}.npz",
                index=index,
            )
            # Avoid retaining large sampled metadata for every bin in the controller.
            r._samples.clear()
            r._prepared.clear()
        print(
            f"Finished in {time.monotonic() - started:.2f}s; scientific passed={report['passed']}",
            flush=True,
        )
        return arrays, report

    inputs = [
        wheel,
        template,
        Path(reference_bundle) / "manifest.json",
        *ident["resources"],
        *ident["reference"]["sources"],
        *cache_inputs,
    ]
    return execute(
        out,
        suite=suite,
        cases=cases,
        bin_indices=bin_indices,
        profiles=profiles,
        worker=worker,
        inputs=inputs,
        diagnostic_requests=diagnostics,
        recipe_revision=recipe_revision,
    )
