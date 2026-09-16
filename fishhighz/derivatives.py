"""Global finite differences and supplied local Jacobians at fixed coordinates."""

from dataclasses import dataclass

import numpy as np

from ._arrays import integer, scalar
from .kernels.derivatives import _combine_three, _scatter_column
from .models.external import _inputs, _invoke
from .parameters import map_jacobian


@dataclass(frozen=True)
class Stencil:
    """Actual float64 points/offsets, requested step and scaled coefficients."""

    method: str
    step: float
    points: tuple[float, float]
    offsets: tuple[float, float]
    requested_offsets: tuple[float, float]
    weights: tuple[float, float]
    scale: float


@dataclass(frozen=True)
class ColumnDiagnostic:
    """One provider/global strategy; analytic columns have no stencil."""

    provider: str
    parameter: str
    method: str
    stencil: Stencil | None


@dataclass(frozen=True)
class CallCount:
    """Actual successful model and Jacobian invocations for one provider."""

    provider: str
    model: int
    jacobian: int


@dataclass(frozen=True)
class DerivativeResult:
    """Owned required-pair powers/Jacobian and invocation-local diagnostics.

    Gather axis 1 with selection.selected_to_required before Fisher assembly;
    multiply mean derivatives by supplied fixed response products exactly once.
    """

    power: np.ndarray
    jacobian: np.ndarray
    columns: tuple[ColumnDiagnostic, ...]
    calls: tuple[CallCount, ...]


def _stencil(parameter, x, h):
    lo, hi = parameter.bounds or (None, None)
    lower, upper = (-np.inf if lo is None else lo), (np.inf if hi is None else hi)
    for method, multiples in (
        ("central", (-1, 1)),
        ("forward", (1, 2)),
        ("backward", (-1, -2)),
    ):
        with np.errstate(over="ignore", invalid="ignore", under="ignore"):
            requested = np.asarray(multiples, dtype=np.float64) * h
            points = np.float64(x) + requested
            offsets = points - x
        if not np.all(np.isfinite(points)) or not np.all(np.isfinite(offsets)):
            continue
        if not np.all((points >= lower) & (points <= upper)):
            continue
        # Both points must be resolved and follow the intended direction/order.
        if (
            np.any(offsets == 0)
            or points[0] == points[1]
            or not np.all(np.sign(offsets) == np.sign(requested))
        ):
            break
        if (method == "central" and not points[0] < x < points[1]) or (
            method != "central" and not abs(offsets[0]) < abs(offsets[1])
        ):
            break
        scale = float(np.max(np.abs(offsets)))
        u, v = offsets / scale
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            weights = (-v / (u * (u - v)), -u / (v * (v - u)))
        if not np.all(np.isfinite(weights)):
            break
        return Stencil(
            method,
            h,
            tuple(points),
            tuple(offsets),
            tuple(requested),
            tuple(weights),
            scale,
        )
    raise ValueError(
        f"{parameter.id}: cannot resolve/fit stencil at {x}, bounds={parameter.bounds}, "
        f"requested step={h}; supply a smaller step for bounds, or a representable step"
    )


def _schedule(prepared, theta, steps, step_scale, numerical):
    scale = scalar(step_scale, "step_scale")
    if scale <= 0:
        raise ValueError("step_scale must be positive")
    if not isinstance(numerical, (bool, np.bool_)):
        raise ValueError("numerical must be boolean")
    overrides = {} if steps is None else dict(steps)
    for name, value in overrides.items():
        if name not in prepared.registry.ids:
            raise ValueError(f"unknown step override ID {name}")
        overrides[name] = scalar(value, f"step for {name}")
        if overrides[name] <= 0:
            raise ValueError(f"step for {name} must be positive")
    schedules = []
    # Validate the entire request before any potentially expensive external call.
    for route in prepared.routes:
        provider = route.provider
        columns = []
        for index in np.unique(provider.parameters.binding.local_to_global):
            parameter = prepared.registry.parameters[index]
            analytic = parameter.id in provider.analytic_ids and not numerical
            stencil = None
            if not analytic:
                step = overrides.get(parameter.id, parameter.step)
                if step is None:
                    raise ValueError(
                        f"{parameter.id}: explicit numerical step required"
                    )
                with np.errstate(over="ignore", under="ignore"):
                    h = float(np.float64(step) * scale)
                if not np.isfinite(h) or h <= 0:
                    raise ValueError(
                        f"{parameter.id}: scaled step must be positive finite"
                    )
                stencil = _stencil(parameter, theta[index], h)
            columns.append((int(index), stencil))
        schedules.append(tuple(columns))
    return tuple(schedules)


def evaluate_derivatives(
    prepared, theta, z, k, mu, *, steps=None, step_scale=1.0, numerical=False
):
    """Evaluate mean and global derivatives without changing survey inputs.

    Parameters
    ----------
    prepared : PreparedP3D
        Explicit provider routes with registry-aware bindings.
    theta : array_like
        Global evaluation point, within inclusive bounds.
    z : float
        Fixed nonnegative redshift.
    k, mu : array_like
        Paired nodes, including arbitrary nonempty consecutive grid slices.
    steps : mapping, optional
        Positive finite absolute step overrides keyed by global ID.
    step_scale : float
        Positive finite multiplier, never an implicit guessed step.
    numerical : bool
        Explicitly force all dependent columns numerical, including columns
        normally analytic. Useful for analytic-versus-FD verification.

    Returns
    -------
    DerivativeResult
        Powers (node,required), Jacobian (node,required,global), methods, actual
        stencils and call counts. Unused global columns are exactly zero.

    Notes
    -----
    Prefer central, then forward, then backward second-order stencils. Never
    shrink/clip steps. Perturb global entries before gathering tied locals.
    All schedules are checked before dispatch. Each provider is called once at
    the fiducial and twice per unique numerical dependency; one supplied
    Jacobian call serves all analytic columns. No persistent value cache exists.
    """
    theta, z, k, mu = _inputs(prepared, theta, z, k, mu)
    schedules = _schedule(prepared, theta, steps, step_scale, numerical)
    power = np.empty((len(k), len(prepared.selection.required_pairs)), dtype=np.float64)
    jacobian = np.zeros((*power.shape, len(theta)), dtype=np.float64)
    diagnostics, counts = [], []
    for route, schedule in zip(prepared.routes, schedules):
        provider = route.provider
        fiducial = _invoke(route, theta, z, k, mu, "fiducial")
        power[:, route.columns] = fiducial
        model_calls, jacobian_calls = 1, 0
        if any(stencil is None for _, stencil in schedule):
            local = _invoke(route, theta, z, k, mu, "analytic Jacobian", jacobian=True)
            with np.errstate(over="ignore", invalid="ignore"):
                analytic = map_jacobian(
                    local, provider.parameters.binding.local_to_global, len(theta)
                )
            jacobian_calls = 1
        for index, stencil in schedule:
            parameter = prepared.registry.ids[index]
            if stencil is None:
                values = analytic[:, :, index]
                method = "analytic"
            else:
                method = stencil.method
                values_at_points = []
                for point, offset in zip(stencil.points, stencil.offsets):
                    perturbed = theta.copy()
                    perturbed[index] = point
                    values_at_points.append(
                        _invoke(
                            route,
                            perturbed,
                            z,
                            k,
                            mu,
                            f"parameter {parameter!r}, {method}, step={stencil.step}, offset={offset}",
                        )
                    )
                    model_calls += 1
                with np.errstate(
                    over="ignore", invalid="ignore", divide="ignore", under="ignore"
                ):
                    values = _combine_three(
                        fiducial, *values_at_points, stencil.weights, stencil.scale
                    )
            if not np.all(np.isfinite(values)):
                raise ValueError(
                    f"provider {provider.label!r}, z={z}, parameter {parameter!r}, {method}, pairs={route.pairs.tolist()}: nonfinite derivative arithmetic"
                )
            _scatter_column(jacobian, values, route.columns, index)
            diagnostics.append(
                ColumnDiagnostic(provider.label, parameter, method, stencil)
            )
        counts.append(CallCount(provider.label, model_calls, jacobian_calls))
    return DerivativeResult(power, jacobian, tuple(diagnostics), tuple(counts))


@dataclass(frozen=True)
class ConvergenceColumn:
    """Infinity-norm changes for one provider/global column between two steps."""

    provider: str
    parameter: str
    coarse: Stencil
    fine: Stencil
    absolute_change: float
    relative_change: float
    passed: bool
    family_changed: bool


@dataclass(frozen=True)
class ConvergenceResult:
    """Opt-in comparisons; agreement is evidence, not a proof of accuracy."""

    comparisons: tuple[ConvergenceColumn, ...]
    passed: bool


def check_convergence(
    prepared,
    theta,
    z,
    k,
    mu,
    *,
    atol,
    rtol,
    refinements=1,
    steps=None,
    step_scale=1.0,
    numerical=False,
):
    """Compare h,h/2 (and optionally further halvings), retaining two results.

    Only numerical columns are assessed; fail if there are none. Set numerical
    explicitly to verify normally analytic columns. Use infinity norm across
    owned pairs and nodes: change <= atol + rtol*max(norm(coarse),norm(fine)).
    Relative change is zero for two zero columns (otherwise change/reference).
    Report stencil-family changes; never switch strategies based on agreement.
    Preflight rejects consecutive steps with identical actual point sets before
    any provider/Jacobian call, even at later refinements. Choose a better-resolved
    starting step or fewer refinements; moving just one point is sufficient.
    """
    atol, rtol = scalar(atol, "atol"), scalar(rtol, "rtol")
    if atol < 0 or rtol < 0:
        raise ValueError("convergence tolerances must be nonnegative")
    refinements = integer(refinements, "refinements", 1)
    theta, z, k, mu = _inputs(prepared, theta, z, k, mu)
    scale = scalar(step_scale, "step_scale")
    schedule = _schedule(prepared, theta, steps, scale, numerical)
    if not any(stencil is not None for route in schedule for _, stencil in route):
        raise ValueError(
            "no numerical columns to assess; use numerical=True and explicit steps"
        )
    # Preflight refinements too, so invalid schedules do not partially run a study.
    for level in range(1, refinements + 1):
        refined = _schedule(prepared, theta, steps, scale * 0.5**level, numerical)
        for route, previous, current in zip(prepared.routes, schedule, refined):
            for (index, old), (_, new) in zip(previous, current):
                if old is not None and sorted(old.points) == sorted(new.points):
                    points = tuple(float(point) for point in sorted(old.points))
                    raise ValueError(
                        f"provider {route.provider.label!r}, global parameter "
                        f"{prepared.registry.ids[index]!r}, refinement level {level}: "
                        f"requested steps {old.step} and {new.step} repeat actual "
                        f"points {points}; no effective refinement occurred. "
                        "Use a better-resolved starting step or fewer refinements."
                    )
        schedule = refined
    coarse = evaluate_derivatives(
        prepared, theta, z, k, mu, steps=steps, step_scale=scale, numerical=numerical
    )
    routes = {route.provider.label: route for route in prepared.routes}
    comparisons = []
    for level in range(1, refinements + 1):
        fine = evaluate_derivatives(
            prepared,
            theta,
            z,
            k,
            mu,
            steps=steps,
            step_scale=scale * 0.5**level,
            numerical=numerical,
        )
        for old, new in zip(coarse.columns, fine.columns):
            if old.stencil is None:
                continue
            index = prepared.registry.ids.index(old.parameter)
            columns = routes[old.provider].columns
            a, b = coarse.jacobian[:, columns, index], fine.jacobian[:, columns, index]
            with np.errstate(over="ignore", invalid="ignore"):
                change = float(np.max(np.abs(a - b)))
                reference = float(max(np.max(np.abs(a)), np.max(np.abs(b))))
                threshold = atol + rtol * reference
            if not np.isfinite(change) or not np.isfinite(threshold):
                raise ValueError("nonfinite convergence arithmetic")
            relative = change / reference if reference else 0.0
            comparisons.append(
                ConvergenceColumn(
                    old.provider,
                    old.parameter,
                    old.stencil,
                    new.stencil,
                    change,
                    relative,
                    change <= threshold,
                    old.method != new.method,
                )
            )
        coarse = fine
    return ConvergenceResult(tuple(comparisons), all(c.passed for c in comparisons))
