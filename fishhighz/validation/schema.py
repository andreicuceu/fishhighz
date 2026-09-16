"""Schema-2 scientific payload validation, shared by writer and offline checker.

Hashes establish byte identity; this module separately establishes dimensions,
request identity, node measures, Wick covariance, Fisher and rank-aware summaries.
It imports no external model, reference package, pickle or plotting dependency.
"""

import hashlib
import json

import numpy as np

from ..adapters.legacy_compat import plain
from ..fisher import fisher_from_factors
from ..grids import gauss_legendre_grid
from .cases import bins, recipe, selection
from .numerics import change, contract, field_matrix, relative, summaries, wick

LIMITS = dict(
    fisher_relative=1e-3,
    error_relative=1e-3,
    pair_error_relative=1e-3,
    volume_relative=1e-6,
)
KINDS = (
    "real_bao",
    "synthetic_bao",
    "synthetic_amplitude",
    "external_amplitude",
    "diagnostic",
)
PROFILES = ("compatibility", "accuracy")


def canonical(value):
    return hashlib.sha256(
        json.dumps(plain(value), sort_keys=True, allow_nan=False).encode()
    ).hexdigest()


def token(value):
    return np.frombuffer(bytes.fromhex(canonical(value)), dtype=np.uint8).copy()


def request(case, index, profile, *, kind="real_bao", diagnostic_id=None):
    """Bind fixed recipe identities, targets and thresholds before worker execution."""
    if profile not in PROFILES or kind not in KINDS:
        raise ValueError("unknown profile/payload kind")
    if (
        isinstance(index, bool)
        or not isinstance(index, int)
        or not 0 <= index < len(bins(case))
    ):
        raise ValueError("invalid bin index")
    if (kind == "diagnostic") != (
        isinstance(diagnostic_id, str) and bool(diagnostic_id)
    ):
        raise ValueError("diagnostic kind requires an explicit diagnostic ID")
    sel = selection(case)
    parameters = ["A"] if kind.endswith("amplitude") else [f"ap_{index}", f"at_{index}"]
    return dict(
        case=case,
        bin=index,
        bounds=list(bins(case)[index]),
        profile=profile,
        kind=kind,
        diagnostic_id=diagnostic_id,
        fields=[
            dict(
                id=f.id, kind=f.kind, physical=f.physical_model, background=f.background
            )
            for f in sel.fields
        ],
        selected_pairs=sel.selected_pairs.tolist(),
        required_pairs=sel.required_pairs.tolist(),
        parameters=parameters,
        original_hash=canonical(recipe(case)),
        thresholds=dict(LIMITS),
        response_ownership="observed J, response already applied exactly once",
    )


def validate_request(task):
    if task != request(
        task["case"],
        task["bin"],
        task["profile"],
        kind=task["kind"],
        diagnostic_id=task["diagnostic_id"],
    ):
        raise ValueError(
            "request does not match exact recipe/field/pair/parameter contract"
        )


def grid_nodes(settings):
    """Reconstruct expected ordered paired nodes and measure from declared controls."""
    g = settings["grid"]
    volume = float(g["volume"])
    if not np.isfinite(volume) or volume <= 0:
        raise ValueError("positive finite volume required")
    if g["kind"] == "gauss_legendre":
        for name in ("k_intervals", "k_order", "mu_order"):
            if (
                isinstance(g[name], bool)
                or not isinstance(g[name], int)
                or g[name] <= 0
            ):
                raise ValueError("invalid quadrature count")
        grid = gauss_legendre_grid(
            np.linspace(0.01, 0.5, g["k_intervals"] + 1),
            k_order=g["k_order"],
            mu_order=g["mu_order"],
            h_fid=g["h_fid"],
        )
        return grid.k_flat, grid.mu_flat, volume * grid.q_mode
    if g["kind"] == "legacy":
        # These seven explicit INIs share the literal legacy grid; no normalized endpoint fix.
        k = np.linspace(0.01, 0.5, 500)
        mu = (np.arange(10) + 0.5) / 10
        nodes = np.tile(k, len(mu))
        angles = np.repeat(mu, len(k))
        return nodes, angles, volume * nodes**2 * (k[1] - k[0]) * 0.1 / (2 * np.pi**2)
    raise ValueError("unknown evidence grid contract")


def assemble(task, total, observed_j, settings, *, extra=None):
    """Construct a scientifically bound payload with FishHighz independent F."""
    return _assemble(task, total, observed_j, settings, extra=extra)


def _assemble(task, total, observed_j, settings, *, extra=None, factors=None):
    """Internal assembly; only owned prepared state may supply reusable factors."""
    validate_request(task)
    k, mu, modes = grid_nodes(settings)
    selected = np.asarray(task["selected_pairs"], dtype=np.int64)
    required = np.asarray(task["required_pairs"], dtype=np.int64)
    c = wick(total, modes, required, selected, len(task["fields"]))
    if factors is None:
        f, single = contract(c, np.asarray(observed_j))
    else:
        f = fisher_from_factors(observed_j, factors)
        single = np.einsum(
            "nsi,nsj,ns->sij",
            observed_j,
            observed_j,
            1 / np.diagonal(c, axis1=1, axis2=2),
        )
    arrays = dict(
        k=k,
        mu=mu,
        modes=modes,
        total=np.asarray(total),
        observed_j=np.asarray(observed_j),
        selected_covariance=c,
        selected_pairs=selected,
        required_pairs=required,
        bin_bounds=np.asarray(task["bounds"]),
        request_token=token(task),
        effective_token=token(settings),
        field_min_eigenvalue=np.linalg.eigvalsh(
            field_matrix(total, required, len(task["fields"]))
        )[:, 0],
        **summaries(f, single),
    )
    arrays.update(extra or {})
    report = dict(
        context=task,
        settings=plain(settings),
        effective_hash=canonical(settings),
        passed=True,
        metric_names=[],
        metrics=[],
    )
    return arrays, report


def metric_values(arrays, names):
    """Recompute convergence values from saved paired Fisher/volume controls."""
    return [
        change(
            *arrays["metric_fisher"][i],
            *arrays["metric_pair_fisher"][i],
            *arrays["metric_volume"][i],
        )
        for i in range(len(names))
    ]


def add_metrics(arrays, report, names, fisher, pairs, volumes):
    arrays.update(
        metric_fisher=np.asarray(fisher),
        metric_pair_fisher=np.asarray(pairs),
        metric_volume=np.asarray(volumes),
    )
    report["metric_names"] = list(names)
    report["metrics"] = metric_values(arrays, names)
    report["passed"] = all(
        all(np.isfinite(m[k]) and m[k] <= LIMITS[k] for k in LIMITS)
        for m in report["metrics"]
    )


def _close(a, b, name, rtol=5e-12):
    if (
        np.shape(a) != np.shape(b)
        or not np.isfinite(relative(a, b))
        or relative(a, b) > rtol
    ):
        raise ValueError(f"{name}: inconsistent numerical content")


def _close_blocks(a, b, name, *, ndim=2):
    """Compare each matrix (or scalar) on its own scale, never a stack norm."""
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or a.ndim < ndim:
        raise ValueError(f"{name}: inconsistent numerical dimensions")
    leading = a.shape[:-ndim] if ndim else a.shape
    for index in np.ndindex(leading):
        _close(a[index], b[index], f"{name} block {index}")


def validate_payload(task, arrays, report, *, require_pass=True, schema=3):
    """Validate semantic meaning even for consistently rehashed corruptions."""
    validate_request(task)
    if not isinstance(report.get("passed"), bool):
        raise ValueError("passed must be an explicit boolean")
    if report.get("context") != task:
        raise ValueError("swapped/missing report context")
    settings = report["settings"]
    if report.get("effective_hash") != canonical(settings):
        raise ValueError("effective settings hash mismatch")
    if (
        settings.get("profile") != task["profile"]
        or settings.get("parameters") != task["parameters"]
    ):
        raise ValueError("wrong profile or named parameter order")
    if settings.get("bounds") != task["bounds"] or settings.get("fields") != [
        f["id"] for f in task["fields"]
    ]:
        raise ValueError("wrong effective bin/field identity")
    npar = len(task["parameters"])
    nsel = len(task["selected_pairs"])
    nreq = len(task["required_pairs"])
    nf = len(task["fields"])
    k, mu, modes = grid_nodes(settings)
    nn = len(k)
    shapes = dict(
        k=(nn,),
        mu=(nn,),
        modes=(nn,),
        total=(nn, nreq),
        observed_j=(nn, nsel, npar),
        selected_covariance=(nn, nsel, nsel),
        selected_pairs=(nsel, 2),
        required_pairs=(nreq, 2),
        bin_bounds=(2,),
        request_token=(32,),
        effective_token=(32,),
        field_min_eigenvalue=(nn,),
        fisher=(npar, npar),
        covariance=(npar, npar),
        errors=(npar,),
        correlation=(npar, npar),
        constrained=(npar,),
        rank=(1,),
        pair_fisher=(nsel, npar, npar),
        pair_covariance=(nsel, npar, npar),
        pair_errors=(nsel, npar),
        pair_correlation=(nsel, npar, npar),
        pair_constrained=(nsel, npar),
        pair_rank=(nsel,),
    )
    if not set(shapes) <= arrays.keys():
        raise ValueError(f"missing scientific payload: {set(shapes) - arrays.keys()}")
    indices = {
        "selected_pairs",
        "required_pairs",
        "request_token",
        "effective_token",
        "constrained",
        "rank",
        "pair_constrained",
        "pair_rank",
    }
    for name, a in arrays.items():
        a = np.asarray(a)
        if a.dtype.kind not in ("iu" if name in indices else "fiu") or not np.all(
            np.isfinite(a)
        ):
            raise ValueError(
                f"{name}: require finite real or declared integer index data"
            )
        if name in shapes and a.shape != shapes[name]:
            raise ValueError(f"{name}: wrong declared scientific dimensions")
    for name, expected in [
        ("request_token", token(task)),
        ("effective_token", token(settings)),
        ("bin_bounds", task["bounds"]),
        ("selected_pairs", task["selected_pairs"]),
        ("required_pairs", task["required_pairs"]),
    ]:
        if not np.array_equal(arrays[name], expected):
            raise ValueError(f"{name}: wrong payload identity/content/order")
    for name, expected in [("k", k), ("mu", mu), ("modes", modes)]:
        _close(arrays[name], expected, name, 5e-13)
    c = wick(
        arrays["total"], modes, arrays["required_pairs"], arrays["selected_pairs"], nf
    )
    _close(arrays["selected_covariance"], c, "Wick covariance")
    eigen = np.linalg.eigvalsh(
        field_matrix(arrays["total"], arrays["required_pairs"], nf)
    )[:, 0]
    _close(arrays["field_min_eigenvalue"], eigen, "field eigenvalues")
    if task["profile"] == "accuracy" and task["kind"] in ("real_bao", "diagnostic"):
        from ..covariance import _validate_field_power

        _validate_field_power(arrays["total"], selection(task["case"]))
    # Direct solve, independently of the factor kernel used by the writer.
    f, single = contract(c, arrays["observed_j"], independent=True)
    expected = summaries(f, single)
    for name, value in expected.items():
        if name.startswith("pair_"):
            _close_blocks(arrays[name], value, name, ndim=value.ndim - 1)
        else:
            _close(arrays[name], value, name)
    # Also inspect the supplied F itself: a valid reconstruction does not license bad sign/rank.
    from .numerics import information

    information(arrays["fisher"])
    for matrix in arrays["pair_fisher"]:
        information(matrix)
    names = report.get("metric_names")
    stored = report.get("metrics")
    if (
        not isinstance(names, list)
        or len(set(names)) != len(names)
        or not isinstance(stored, list)
        or len(stored) != len(names)
    ):
        raise ValueError("missing/duplicate metric inventory")
    if schema >= 3 and task["kind"] == "real_bao" and task["profile"] == "accuracy":
        from .trials import validate

        validate(arrays, report)
    if names:
        for key, shape in [
            ("metric_fisher", (len(names), 2, npar, npar)),
            ("metric_pair_fisher", (len(names), 2, nsel, npar, npar)),
            ("metric_volume", (len(names), 2)),
        ]:
            if key not in arrays or arrays[key].shape != shape:
                raise ValueError("wrong convergence dimensions")
        computed = metric_values(arrays, names)
        for a, b in zip(stored, computed):
            if set(a) != set(LIMITS):
                raise ValueError("wrong metric keys")
            for key in LIMITS:
                if not np.isfinite(a[key]) or not np.isclose(
                    a[key], b[key], rtol=5e-12, atol=1e-15
                ):
                    raise ValueError("false convergence metric")
    if task["kind"] == "real_bao":
        provenance = report.get("provenance", {})
        if not {"fishhighz", "reference", "resources", "wheel"} <= provenance.keys():
            raise ValueError("missing real imported provenance")
        wheel = provenance["wheel"]
        imported = provenance["fishhighz"]
        ref = provenance["reference"]
        if (
            not wheel.get("modules")
            or wheel.get("modules") != imported.get("module_hashes")
            or wheel.get("origin") != imported.get("origin")
            or not ref.get("sources")
            or not ref.get("versions")
            or not provenance["resources"]
        ):
            raise ValueError("inconsistent imported source provenance")
        from pathlib import Path

        if any(Path(p).parent != Path(ref["reference_origin"]) for p in ref["sources"]):
            raise ValueError("inventoried source differs from imported reference")
        for sha in [
            wheel.get("sha256"),
            *wheel["modules"].values(),
            *ref["sources"].values(),
            *provenance["resources"].values(),
        ]:
            if (
                not isinstance(sha, str)
                or len(sha) != 64
                or any(c not in "0123456789abcdef" for c in sha)
            ):
                raise ValueError("invalid source/resource digest")
        if task["profile"] == "accuracy":
            version = report.get("trial_contract", {}).get("version")
            expected = {
                "k",
                "mu",
                "magnitude",
                "volume",
                "step",
                "combined",
            }
            if version == 1:
                expected.add("weights")
            elif version != 2:
                raise ValueError("unknown accuracy trial contract")
            if set(names) != expected:
                raise ValueError("incomplete accuracy convergence inventory")
        if task["profile"] == "compatibility":
            if arrays.get("reference_fisher", np.empty(0)).shape != (
                npar,
                npar,
            ) or arrays.get("reference_pair_fisher", np.empty(0)).shape != (
                nsel,
                npar,
                npar,
            ):
                raise ValueError("missing matched reference information")
            comp = change(
                arrays["fisher"],
                arrays["reference_fisher"],
                arrays["pair_fisher"],
                arrays["reference_pair_fisher"],
                1,
                1,
            )
            comparison = report.get("comparison", {})
            if set(comparison) != set(comp) or any(
                not np.isfinite(comparison[key])
                or not np.isclose(comparison[key], comp[key], rtol=5e-12, atol=1e-15)
                for key in comp
            ):
                raise ValueError("incorrect reference comparison metrics")
            if (
                comp["fisher_relative"] > 5e-12
                or max(comp["error_relative"], comp["pair_error_relative"]) > 1e-6
            ):
                if require_pass:
                    raise ValueError("compatibility comparison failed")
    metrics_pass = all(all(m[k] <= LIMITS[k] for k in LIMITS) for m in stored)
    if require_pass and (report["passed"] is not True or not metrics_pass):
        raise ValueError("scientific validation failed")
    if report["passed"] is True and report.get("unresolved_controls"):
        raise ValueError("unresolved controls disguised as passed")
    if report["passed"] is True and not metrics_pass:
        raise ValueError("failed convergence disguised as passed")
    return True
