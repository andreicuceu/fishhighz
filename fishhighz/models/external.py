"""Explicit external routes; no inferred physics, units, parameters, or noise."""

from dataclasses import dataclass

import numpy as np

from .._arrays import label, readonly, real_array, scalar
from ..fields import PairSelection
from ..parameters import ParameterBinding, ParameterRegistry, gather_local
from .protocols import validate_p1d, validate_p3d, validate_p3d_jacobian


@dataclass(frozen=True, init=False, eq=False)
class BoundParameters:
    """Create an equality binding and retain its exact registry identity.

    Parameters
    ----------
    registry : ParameterRegistry
        Common global basis. A separately constructed registry is not identical.
    local_names : sequence of str
        Provider argument order, possibly empty or repeatedly bound.
    bindings : mapping
        Explicit local-name to global-ID equality mapping.
    """

    registry: ParameterRegistry
    binding: ParameterBinding

    def __init__(self, registry, local_names, bindings):
        if not isinstance(registry, ParameterRegistry):
            raise ValueError("expected ParameterRegistry")
        object.__setattr__(self, "registry", registry)
        object.__setattr__(
            self, "binding", ParameterBinding(registry, local_names, bindings)
        )


@dataclass(frozen=True, init=False, eq=False)
class P3DProvider:
    """Declare a callable, bound parameters, exact pair ownership and derivatives.

    Parameters
    ----------
    label : str
        Unique diagnostic label.
    model : callable
        ``model(theta_local, z, k, mu, pairs) -> (node,pair)``.
    parameters : BoundParameters
        Explicit registry-aware local binding.
    pairs : sequence of pairs
        Field ID pairs or original field-index pairs; reversed aliases allowed.
    jacobian : callable, optional
        Same arguments, returning every local column ``(node,pair,local)``.
    analytic_ids : sequence of str, optional
        Global columns consumed analytically. None means all dependencies when
        a Jacobian is supplied, otherwise none. Empty means entirely numerical.
    """

    label: str
    model: object
    parameters: BoundParameters
    pairs: tuple
    jacobian: object
    analytic_ids: tuple

    def __init__(
        self, label, model, parameters, pairs, *, jacobian=None, analytic_ids=None
    ):
        _provider(label, model, parameters)
        if jacobian is not None and not callable(jacobian):
            raise ValueError("jacobian must be callable")
        if isinstance(analytic_ids, str):
            raise ValueError("analytic_ids must be a sequence of global IDs")
        pairs = tuple(pairs)
        if any(isinstance(pair, str) for pair in pairs):
            raise ValueError("a route pair must contain two IDs or indices")
        dependencies = tuple(
            parameters.registry.ids[i]
            for i in np.unique(parameters.binding.local_to_global)
        )
        analytic = (
            (dependencies if jacobian is not None else ())
            if analytic_ids is None
            else tuple(analytic_ids)
        )
        if len(set(analytic)) != len(analytic) or any(
            i not in dependencies for i in analytic
        ):
            raise ValueError(
                f"{label}: duplicate, unknown or nondependent analytic IDs"
            )
        if analytic and jacobian is None:
            raise ValueError(f"{label}: analytic selection requires a Jacobian")
        for name, value in (
            ("label", label),
            ("model", model),
            ("parameters", parameters),
            ("pairs", tuple(tuple(p) for p in pairs)),
            ("jacobian", jacobian),
            ("analytic_ids", analytic),
        ):
            object.__setattr__(self, name, value)


def _provider(name, model, parameters):
    label(name, "provider label")
    if not callable(model) or not isinstance(parameters, BoundParameters):
        raise ValueError("provider requires a callable and BoundParameters")


@dataclass(frozen=True, eq=False)
class _Route:
    provider: P3DProvider
    pairs: np.ndarray
    columns: np.ndarray


@dataclass(frozen=True, init=False, eq=False)
class PreparedP3D:
    """Prepare exact closure ownership; unused declared routes are rejected.

    Providers receive their required pairs in required-pair order (increasing
    original field indices). Provider order never changes output column order.
    Preparation stores structural maps only, never model values.
    """

    registry: ParameterRegistry
    selection: PairSelection
    routes: tuple[_Route, ...]

    def __init__(self, registry, selection, providers):
        if not isinstance(registry, ParameterRegistry) or not isinstance(
            selection, PairSelection
        ):
            raise ValueError("expected registry and pair selection")
        required = {tuple(p): i for i, p in enumerate(selection.required_pairs)}
        owners, labels, routes = set(), set(), []
        for provider in providers:
            if not isinstance(provider, P3DProvider):
                raise ValueError("expected P3DProvider")
            if provider.parameters.registry is not registry:
                raise ValueError(
                    f"{provider.label}: binding belongs to a different registry"
                )
            if provider.label in labels:
                raise ValueError("duplicate provider label")
            labels.add(provider.label)
            declared = PairSelection(selection.fields, provider.pairs).selected_pairs
            columns = []
            for pair in declared:
                pair = tuple(pair)
                if pair not in required:
                    raise ValueError(f"{provider.label}: unused declared route {pair}")
                if pair in owners:
                    raise ValueError(f"{provider.label}: overlapping ownership {pair}")
                owners.add(pair)
                columns.append(required[pair])
            columns = readonly(sorted(columns), np.int64)
            routes.append(
                _Route(
                    provider,
                    readonly(selection.required_pairs[columns], np.int64),
                    columns,
                )
            )
        if owners != set(required):
            raise ValueError(f"missing required powers: {set(required) - owners}")
        object.__setattr__(self, "registry", registry)
        object.__setattr__(self, "selection", selection)
        object.__setattr__(self, "routes", tuple(routes))


def _state(registry, theta, z):
    theta = real_array(theta, "global parameters")
    if theta.shape != (len(registry.ids),):
        raise ValueError("global parameter vector has wrong shape")
    for parameter, value in zip(registry.parameters, theta):
        lo, hi = parameter.bounds or (None, None)
        if (lo is not None and value < lo) or (hi is not None and value > hi):
            raise ValueError(f"{parameter.id}: parameter outside inclusive bounds")
    z = scalar(z, "redshift")
    if z < 0:
        raise ValueError("redshift must be nonnegative")
    theta.flags.writeable = False
    return theta, z


def _nodes(value, name):
    array = real_array(value, name)
    if array.ndim != 1 or not len(array):
        raise ValueError(f"{name}: require nonempty 1D coordinates")
    array.flags.writeable = False
    return array


def _inputs(prepared, theta, z, k, mu):
    theta, z = _state(prepared.registry, theta, z)
    k, mu = _nodes(k, "k"), _nodes(mu, "mu")
    if k.shape != mu.shape or np.any(k <= 0) or np.any((mu < 0) | (mu > 1)):
        raise ValueError("require paired equal shapes, k>0 and mu in [0,1]")
    return theta, z, k, mu


def _invoke(route, theta, z, k, mu, context, *, jacobian=False):
    provider = route.provider
    binding = provider.parameters.binding
    # Per-call snapshots protect caller/prepared state even if external code
    # deliberately resets a NumPy writeable flag. No full stencil history.
    local = readonly(gather_local(theta, binding.local_to_global), np.float64)
    args = (
        local,
        z,
        readonly(k, np.float64),
        readonly(mu, np.float64),
        readonly(route.pairs, np.int64),
    )
    try:
        if jacobian:
            return validate_p3d_jacobian(
                provider.jacobian(*args), len(k), len(route.pairs), len(local)
            )
        return validate_p3d(provider.model(*args), len(k), len(route.pairs))
    except Exception as error:
        raise ValueError(
            f"provider {provider.label!r}, z={z}, {context}, "
            f"theta_local={local.tolist()}, pairs={route.pairs.tolist()}: {error}"
        ) from error


def evaluate_p3d(prepared, theta, z, k, mu):
    """Evaluate owned float64 intrinsic powers in exact required-pair order.

    Inputs are read-only snapshots. k (h_fid/Mpc) and mu are paired nodes or
    slices; output units are (Mpc/h_fid)^3. Providers own physics and hidden
    cache invalidation. No transformations, covariance or noise are evaluated.
    """
    theta, z, k, mu = _inputs(prepared, theta, z, k, mu)
    output = np.empty(
        (len(k), len(prepared.selection.required_pairs)), dtype=np.float64
    )
    for route in prepared.routes:
        output[:, route.columns] = _invoke(
            route, theta, z, k, mu, "fiducial evaluation"
        )
    return output


def evaluate_p1d(model, parameters, theta, z, k_parallel_velocity, *, label="p1d"):
    """Explicit independent P1D call: s/km nodes (including zero), km/s power.

    ``parameters`` is its own BoundParameters; ``theta`` is that registry's
    global vector. No default model, P3D integration or unit conversion exists.
    Callable inputs are read-only owned snapshots; outputs are copied at once.
    """
    _provider(label, model, parameters)
    theta, z = _state(parameters.registry, theta, z)
    k = _nodes(k_parallel_velocity, "k_parallel_velocity")
    if np.any(k < 0):
        raise ValueError("velocity wavenumbers must be nonnegative")
    local = readonly(
        gather_local(theta, parameters.binding.local_to_global), np.float64
    )
    try:
        return validate_p1d(model(local, z, k), len(k))
    except Exception as error:
        raise ValueError(
            f"provider {label!r}, z={z}, P1D evaluation, theta_local={local.tolist()}: {error}"
        ) from error
