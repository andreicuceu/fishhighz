"""Explicit equality bindings; fixed settings are outside the free vector."""

from dataclasses import dataclass

import numpy as np

from ._arrays import indices, integer, label, readonly, real_array, scalar


@dataclass(frozen=True)
class Parameter:
    """Free-parameter metadata; bounds are inclusive with optional open ends.

    Parameters
    ----------
    id : str
        Opaque unique identity, without automatic scope parsing.
    fiducial : float
        Finite initial value, including negative or zero values.
    role : {'target', 'nuisance'}
        Result metadata only.
    bounds : tuple, optional
        (lower, upper); either endpoint may be None. Finite endpoints must
        satisfy lower < upper and include the fiducial.
    step : float, optional
        Positive absolute finite-difference step, independent of bound width.
    """

    id: str
    fiducial: float
    role: str
    bounds: tuple[float | None, float | None] | None = None
    step: float | None = None

    def __post_init__(self):
        label(self.id, "parameter ID")
        value = scalar(self.fiducial, "fiducial")
        object.__setattr__(self, "fiducial", value)
        if self.role not in ("target", "nuisance"):
            raise ValueError("role must be target or nuisance")
        if self.bounds is not None:
            if len(self.bounds) != 2:
                raise ValueError("bounds must have two endpoints")
            lo, hi = (None if x is None else scalar(x, "bound") for x in self.bounds)
            if lo is not None and hi is not None and lo >= hi:
                raise ValueError("lower bound must be below upper bound")
            if (lo is not None and value < lo) or (hi is not None and value > hi):
                raise ValueError("fiducial outside bounds")
            object.__setattr__(self, "bounds", (lo, hi))
        if self.step is not None:
            step = scalar(self.step, "step")
            if step <= 0:
                raise ValueError("step must be positive")
            object.__setattr__(self, "step", step)


@dataclass(frozen=True, init=False, eq=False)
class ParameterRegistry:
    """Ordered metadata and owned read-only C-order float64 fiducials."""

    parameters: tuple[Parameter, ...]
    ids: tuple[str, ...]
    fiducials: np.ndarray

    def __init__(self, parameters):
        parameters = tuple(parameters)
        if not parameters or any(not isinstance(p, Parameter) for p in parameters):
            raise ValueError("registry requires at least one Parameter")
        ids = tuple(p.id for p in parameters)
        if len(set(ids)) != len(ids):
            raise ValueError("duplicate parameter ID")
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "ids", ids)
        object.__setattr__(
            self, "fiducials", readonly([p.fiducial for p in parameters], np.float64)
        )


@dataclass(frozen=True, init=False, eq=False)
class ParameterBinding:
    """Bind local_names to an explicit name-to-registry-ID mapping.

    Empty dependencies are supported. Stored indices own read-only int64 data.
    Multiple local names may map to the same global ID; nothing is inferred.
    """

    local_names: tuple[str, ...]
    local_to_global: np.ndarray

    def __init__(self, registry, local_names, bindings):
        names = tuple(local_names)
        for name in names:
            label(name, "local name")
        if len(set(names)) != len(names):
            raise ValueError("duplicate local name")
        if set(bindings) != set(names):
            raise ValueError("bindings must exactly match local names")
        lookup = {name: i for i, name in enumerate(registry.ids)}
        if any(bindings[name] not in lookup for name in names):
            raise ValueError("unknown global parameter")
        object.__setattr__(self, "local_names", names)
        object.__setattr__(
            self,
            "local_to_global",
            readonly([lookup[bindings[n]] for n in names], np.int64),
        )


def gather_local(theta_global, local_to_global):
    """Gather a finite 1D global vector into declared local order (float64)."""
    theta = real_array(theta_global, "theta_global")
    if theta.ndim != 1:
        raise ValueError("global vector must be one-dimensional")
    index = indices(local_to_global, len(theta))
    return theta[index]


def map_jacobian(local_jacobian, local_to_global, n_global):
    """Sum (node, pair, local) derivatives into (node, pair, global).

    Repeated equality bindings accumulate by the chain rule. Unused columns
    are zero. Inputs are validated before the array-only accumulation.
    """
    n_global = integer(n_global, "n_global", 1)
    index = indices(local_to_global, n_global)
    jac = real_array(local_jacobian, "local_jacobian")
    if jac.ndim != 3 or jac.shape[2] != len(index):
        raise ValueError("Jacobian must have shape (n_node, n_pair, n_local)")
    result = np.zeros((*jac.shape[:2], n_global), dtype=np.float64)
    for local, global_ in enumerate(index):
        result[:, :, global_] += jac[:, :, local]
    return result
