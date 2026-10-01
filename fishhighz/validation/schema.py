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
PROFILES = ("compatibility", "accuracy", "full-compatibility", "fixed-compatibility")


def canonical(value):
    """Hash a deterministic JSON representation of scientific metadata.

    Parameters
    ----------
    value : object
        JSON-compatible metadata, allowing NumPy values through plain
        conversion.

    Returns
    -------
    digest : str
        SHA-256 hexadecimal digest of the sorted finite JSON representation.
    """
    return hashlib.sha256(
        json.dumps(plain(value), sort_keys=True, allow_nan=False).encode()
    ).hexdigest()


def token(value):
    """Encode a canonical metadata hash as a numerical evidence token.

    Parameters
    ----------
    value : object
        JSON-compatible scientific metadata.

    Returns
    -------
    token : ndarray of uint8, shape (32,)
        Owned byte array containing the canonical SHA-256 digest.
    """
    return np.frombuffer(bytes.fromhex(canonical(value)), dtype=np.uint8).copy()


def request(
    case, index, profile, *, kind="real_bao", diagnostic_id=None, recipe_revision=None
):
    """Bind fixed recipe identities, targets and thresholds before worker execution.

    Parameters
    ----------
    case : str
        Identifier of one of the seven original DESI-2 validation
        configurations.
    index : int
        Zero-based redshift-bin index.
    profile : str
        Declared accuracy or compatibility profile.
    kind : str
        Scientific record type used to distinguish primary and diagnostic
        evidence. Default is ``'real_bao'``.
    diagnostic_id : str or None
        Nonempty identifier required only for diagnostic records. Default is
        ``None``.
    recipe_revision : str or None
        Declared revision of the physical recipe; None retains historical
        request semantics. Default is ``None``.

    Returns
    -------
    request : dict
        Exact case/bin/profile identity with ordered fields, pairs, parameters
        and thresholds.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
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
    if (
        profile in ("full-compatibility", "fixed-compatibility")
        and recipe_revision is None
    ):
        from .profile_definitions import REVISION

        recipe_revision = REVISION
    if recipe_revision is not None:
        from .profile_definitions import REVISION, forecast_selection

        if recipe_revision != REVISION:
            raise ValueError("unknown recipe revision")
        pair_selection = forecast_selection(case, index)
    else:
        pair_selection = selection(case)
    parameters = ["A"] if kind.endswith("amplitude") else [f"ap_{index}", f"at_{index}"]
    result = dict(
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
            for f in pair_selection.fields
        ],
        selected_pairs=pair_selection.selected_pairs.tolist(),
        required_pairs=pair_selection.required_pairs.tolist(),
        parameters=parameters,
        original_hash=canonical(recipe(case)),
        thresholds=dict(LIMITS),
        response_ownership="observed J, response already applied exactly once",
    )

    if recipe_revision is not None:
        result["recipe_revision"] = recipe_revision
    return result


def validate_request(task):
    """Check a request against the declared scientific recipe.

    Parameters
    ----------
    task : dict
        Declared case, bin, selected field pairs, parameter order and validation
        thresholds.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    if task != request(
        task["case"],
        task["bin"],
        task["profile"],
        kind=task["kind"],
        diagnostic_id=task["diagnostic_id"],
        recipe_revision=task.get("recipe_revision"),
    ):
        raise ValueError(
            "request does not match exact recipe/field/pair/parameter contract"
        )


def grid_nodes(settings):
    """Reconstruct expected ordered paired nodes and measure from declared controls.

    Parameters
    ----------
    settings : dict
        Declared grid, volume, derivative and physical-model settings.

    Returns
    -------
    k_grid : ndarray, shape (n_cell,)
        Paired Fourier wavenumbers in h/Mpc.
    mu_grid : ndarray, shape (n_cell,)
        Paired dimensionless direction cosines.
    modes : ndarray, shape (n_cell,)
        Independent Fourier-mode counts including the declared volume.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    grid_settings = settings["grid"]
    volume = float(grid_settings["volume"])
    if not np.isfinite(volume) or volume <= 0:
        raise ValueError("positive finite volume required")
    if grid_settings["kind"] == "gauss_legendre":
        for name in ("k_intervals", "k_order", "mu_order"):
            if (
                isinstance(grid_settings[name], bool)
                or not isinstance(grid_settings[name], int)
                or grid_settings[name] <= 0
            ):
                raise ValueError("invalid quadrature count")
        grid = gauss_legendre_grid(
            np.linspace(0.01, 0.5, grid_settings["k_intervals"] + 1),
            k_order=grid_settings["k_order"],
            mu_order=grid_settings["mu_order"],
            h_fid=grid_settings["h_fid"],
        )
        return grid.k_flat, grid.mu_flat, volume * grid.q_mode
    if grid_settings["kind"] == "legacy":
        # These seven explicit INIs share the literal legacy grid; no normalized endpoint fix.
        k_grid = np.linspace(0.01, 0.5, 500)
        mu_grid = (np.arange(10) + 0.5) / 10
        nodes = np.tile(k_grid, len(mu_grid))
        angles = np.repeat(mu_grid, len(k_grid))
        return (
            nodes,
            angles,
            volume * nodes**2 * (k_grid[1] - k_grid[0]) * 0.1 / (2 * np.pi**2),
        )
    raise ValueError("unknown evidence grid contract")


def assemble(task, total, observed_j, settings, *, extra=None):
    """Construct a scientifically bound payload with FishHighz independent F.

    Parameters
    ----------
    task : dict
        Declared case, bin, selected field pairs, parameter order and validation
        thresholds.
    total : array_like, shape (n_cell, n_required_pair)
        Signed total field-pair powers, including noise, in (Mpc/h)^3.
    observed_j : ndarray, shape (n_cell, n_selected_pair, n_parameter)
        Observed mean-spectrum derivatives, including field responses, in power
        units per parameter unit.
    settings : dict
        Declared grid, volume, derivative and physical-model settings.
    extra : dict or None
        Additional numerical evidence arrays to append. Default is ``None``.

    Returns
    -------
    arrays : dict of str to ndarray
        Ordered grid, powers, derivatives, Wick covariance and rank-aware Fisher
        summaries.
    report : dict
        Bound request, physical settings and initial validation metadata.
    """
    return _assemble(task, total, observed_j, settings, extra=extra)


def _assemble(task, total, observed_j, settings, *, extra=None, factors=None):
    """Internal assembly; only owned prepared state may supply reusable factors.

    Parameters
    ----------
    task : dict
        Declared case, bin, selected field pairs, parameter order and validation
        thresholds.
    total : array_like, shape (n_cell, n_required_pair)
        Signed total field-pair powers, including noise, in (Mpc/h)^3.
    observed_j : ndarray, shape (n_cell, n_selected_pair, n_parameter)
        Observed mean-spectrum derivatives, including field responses, in power
        units per parameter unit.
    settings : dict
        Declared grid, volume, derivative and physical-model settings.
    extra : dict or None
        Additional numerical evidence arrays to append. Default is ``None``.
    factors : ndarray, shape (n_cell, n_selected_pair, n_selected_pair), or None
        Lower Cholesky factors in (Mpc/h)^3 from the same owned prepared state; None
        computes factors during contraction. Default is ``None``.

    Returns
    -------
    arrays : dict of str to ndarray
        Ordered grid, powers, derivatives, Wick covariance and rank-aware Fisher
        summaries.
    report : dict
        Bound request, physical settings and initial validation metadata.
    """
    validate_request(task)
    k_grid, mu_grid, modes = grid_nodes(settings)
    selected = np.asarray(task["selected_pairs"], dtype=np.int64)
    required = np.asarray(task["required_pairs"], dtype=np.int64)
    covariance = wick(total, modes, required, selected, len(task["fields"]))
    if factors is None:
        joint_fisher, single = contract(covariance, np.asarray(observed_j))
    else:
        joint_fisher = fisher_from_factors(observed_j, factors)
        single = np.einsum(
            "nsi,nsj,ns->sij",
            observed_j,
            observed_j,
            1 / np.diagonal(covariance, axis1=1, axis2=2),
        )

    arrays = dict(
        k=k_grid,
        mu=mu_grid,
        modes=modes,
        total=np.asarray(total),
        observed_j=np.asarray(observed_j),
        selected_covariance=covariance,
        selected_pairs=selected,
        required_pairs=required,
        bin_bounds=np.asarray(task["bounds"]),
        request_token=token(task),
        effective_token=token(settings),
        field_min_eigenvalue=np.linalg.eigvalsh(
            field_matrix(total, required, len(task["fields"]))
        )[:, 0],
        **summaries(joint_fisher, single),
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
    """Recompute convergence values from saved paired Fisher/volume controls.

    Parameters
    ----------
    arrays : dict of str to ndarray
        Numerical evidence arrays; Fourier-cell axes and pair order follow the
        declared task. Powers use (Mpc/h)^3 and volumes use (Mpc/h)^3 unless
        separately labeled.
    names : sequence of str
        Ordered names of the quantities to process.

    Returns
    -------
    metrics : list of dict
        Relative information, error, individual-error and volume changes for
        each named refinement.
    """
    return [
        change(
            *arrays["metric_fisher"][i],
            *arrays["metric_pair_fisher"][i],
            *arrays["metric_volume"][i],
        )
        for i in range(len(names))
    ]


def add_metrics(arrays, report, names, fisher, pairs, volumes):
    """Attach actual refinement operands and recomputed convergence metrics.

    Parameters
    ----------
    arrays : dict of str to ndarray
        Numerical evidence arrays; Fourier-cell axes and pair order follow the
        declared task. Powers use (Mpc/h)^3 and volumes use (Mpc/h)^3 unless
        separately labeled.
    report : dict
        Scientific settings, provenance, array inventory and validation outcomes
        associated with the numerical evidence.
    names : sequence of str
        Ordered names of the quantities to process.
    fisher : array_like, shape (n_metric, 2, n_parameter, n_parameter)
        Lower and upper joint Fisher matrices for each refinement.
    pairs : array_like, shape (n_metric, 2, n_pair, n_parameter, n_parameter)
        Corresponding individual-spectrum Fisher matrices.
    volumes : array_like, shape (n_metric, 2)
        Corresponding volumes in (Mpc/h)^3.

    Notes
    -----
    Updates arrays and report in place, including the numerical pass flag.
    """
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
    """Check numerical evidence against an independent reconstruction.

    Parameters
    ----------
    a : array_like
        Saved quantity or stack of quantities.
    b : array_like
        Reconstructed quantity with the same shape and units.
    name : str
        Quantity or record label used in diagnostics.
    rtol : float
        Dimensionless relative numerical tolerance. Default is ``5e-12``.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    if (
        np.shape(a) != np.shape(b)
        or not np.isfinite(relative(a, b))
        or relative(a, b) > rtol
    ):
        raise ValueError(f"{name}: inconsistent numerical content")


def _close_blocks(a, b, name, *, ndim=2):
    """Compare each matrix (or scalar) on its own scale, never a stack norm.

    Parameters
    ----------
    a : array_like
        Saved quantity or stack of quantities.
    b : array_like
        Reconstructed quantity with the same shape and units.
    name : str
        Quantity or record label used in diagnostics.
    ndim : int
        Number of trailing dimensions in each independently scaled block; zero
        compares scalars. Default is ``2``.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or a.ndim < ndim:
        raise ValueError(f"{name}: inconsistent numerical dimensions")
    leading = a.shape[:-ndim] if ndim else a.shape
    for index in np.ndindex(leading):
        _close(a[index], b[index], f"{name} block {index}")


def validate_payload(task, arrays, report, *, require_pass=True, schema=3):
    """Validate semantic meaning even for consistently rehashed corruptions.

    Parameters
    ----------
    task : dict
        Declared case, bin, selected field pairs, parameter order and validation
        thresholds.
    arrays : dict of str to ndarray
        Numerical evidence arrays; Fourier-cell axes and pair order follow the
        declared task. Powers use (Mpc/h)^3 and volumes use (Mpc/h)^3 unless
        separately labeled.
    report : dict
        Scientific settings, provenance, array inventory and validation outcomes
        associated with the numerical evidence.
    require_pass : bool
        Whether to require numerical qualification as well as internally
        consistent evidence. Default is ``True``.
    schema : int
        Evidence schema version controlling trial-binding requirements. Default
        is ``3``.

    Returns
    -------
    valid : bool
        True when all required identities, numerical reconstructions and
        requested pass conditions agree.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
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
        field["id"] for field in task["fields"]
    ]:
        raise ValueError("wrong effective bin/field identity")

    # Establish physical axis sizes before checking saved numerical arrays.
    n_parameters = len(task["parameters"])
    n_selected_pairs = len(task["selected_pairs"])
    n_required_pairs = len(task["required_pairs"])
    n_fields = len(task["fields"])
    k_grid, mu_grid, modes = grid_nodes(settings)
    n_cells = len(k_grid)
    shapes = dict(
        k=(n_cells,),
        mu=(n_cells,),
        modes=(n_cells,),
        total=(n_cells, n_required_pairs),
        observed_j=(n_cells, n_selected_pairs, n_parameters),
        selected_covariance=(n_cells, n_selected_pairs, n_selected_pairs),
        selected_pairs=(n_selected_pairs, 2),
        required_pairs=(n_required_pairs, 2),
        bin_bounds=(2,),
        request_token=(32,),
        effective_token=(32,),
        field_min_eigenvalue=(n_cells,),
        fisher=(n_parameters, n_parameters),
        covariance=(n_parameters, n_parameters),
        errors=(n_parameters,),
        correlation=(n_parameters, n_parameters),
        constrained=(n_parameters,),
        rank=(1,),
        pair_fisher=(n_selected_pairs, n_parameters, n_parameters),
        pair_covariance=(n_selected_pairs, n_parameters, n_parameters),
        pair_errors=(n_selected_pairs, n_parameters),
        pair_correlation=(n_selected_pairs, n_parameters, n_parameters),
        pair_constrained=(n_selected_pairs, n_parameters),
        pair_rank=(n_selected_pairs,),
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
    for name, expected in [("k", k_grid), ("mu", mu_grid), ("modes", modes)]:
        _close(arrays[name], expected, name, 5e-13)

    # Rebuild Wick covariance from total powers, independently of saved factors.
    covariance = wick(
        arrays["total"],
        modes,
        arrays["required_pairs"],
        arrays["selected_pairs"],
        n_fields,
    )
    _close(arrays["selected_covariance"], covariance, "Wick covariance")
    eigen = np.linalg.eigvalsh(
        field_matrix(arrays["total"], arrays["required_pairs"], n_fields)
    )[:, 0]
    _close(arrays["field_min_eigenvalue"], eigen, "field eigenvalues")
    if task["profile"] == "accuracy" and task["kind"] in ("real_bao", "diagnostic"):
        from ..covariance import _validate_field_power
        from .profile_definitions import forecast_selection

        pair_selection = (
            forecast_selection(task["case"], task["bin"])
            if task.get("recipe_revision")
            else selection(task["case"])
        )
        _validate_field_power(arrays["total"], pair_selection)
    # Direct solve, independently of the factor kernel used by the writer.
    joint_fisher, single = contract(covariance, arrays["observed_j"], independent=True)
    expected = summaries(joint_fisher, single)
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
            ("metric_fisher", (len(names), 2, n_parameters, n_parameters)),
            (
                "metric_pair_fisher",
                (len(names), 2, n_selected_pairs, n_parameters, n_parameters),
            ),
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
                or any(character not in "0123456789abcdef" for character in sha)
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
            if version in (1, 3):
                expected.add("weights")
            elif version not in (2, 3):
                raise ValueError("unknown accuracy trial contract")
            if set(names) != expected:
                raise ValueError("incomplete accuracy convergence inventory")
        if task["profile"] == "compatibility":
            if arrays.get("reference_fisher", np.empty(0)).shape != (
                n_parameters,
                n_parameters,
            ) or arrays.get("reference_pair_fisher", np.empty(0)).shape != (
                n_selected_pairs,
                n_parameters,
                n_parameters,
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
    metrics_pass = all(
        all(m[k_grid] <= LIMITS[k_grid] for k_grid in LIMITS) for m in stored
    )
    if require_pass and (report["passed"] is not True or not metrics_pass):
        raise ValueError("scientific validation failed")
    if report["passed"] is True and report.get("unresolved_controls"):
        raise ValueError("unresolved controls disguised as passed")
    if report["passed"] is True and not metrics_pass:
        raise ValueError("failed convergence disguised as passed")
    return True
