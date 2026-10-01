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
        """Validate the parameter fiducial, role, bounds, and difference step.

        Returns
        -------
        None
            Normalize numeric metadata to Python floats in the frozen dataclass.

        Raises
        ------
        ValueError
            If the parameter label, role, bounds, fiducial, or positive absolute
            difference step is invalid. Bounds and steps use the parameter's units.
        """
        label(self.id, "parameter ID")
        value = scalar(self.fiducial, "fiducial")
        object.__setattr__(self, "fiducial", value)
        if self.role not in ("target", "nuisance"):
            raise ValueError("role must be target or nuisance")
        if self.bounds is not None:
            if len(self.bounds) != 2:
                raise ValueError("bounds must have two endpoints")
            lower_bound, upper_bound = (
                None if x is None else scalar(x, "bound") for x in self.bounds
            )
            if (
                lower_bound is not None
                and upper_bound is not None
                and lower_bound >= upper_bound
            ):
                raise ValueError("lower bound must be below upper bound")
            if (lower_bound is not None and value < lower_bound) or (
                upper_bound is not None and value > upper_bound
            ):
                raise ValueError("fiducial outside bounds")
            object.__setattr__(self, "bounds", (lower_bound, upper_bound))
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
        """Store an ordered set of uniquely identified free parameters.

        Parameters
        ----------
        parameters : iterable of Parameter
            Nonempty parameter sequence defining global vector order. Values retain
            each parameter's physical units; no unit conversion is performed.

        Returns
        -------
        None
            Store immutable metadata tuples and a read-only float64 fiducial vector.

        Raises
        ------
        ValueError
            If the sequence is empty, contains non-Parameter entries, or repeats IDs.
        """
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
        """Bind each declared local parameter name to an explicit global ID.

        Parameters
        ----------
        registry : ParameterRegistry
            Global free-parameter definitions.
        local_names : iterable of str
            Unique provider parameter names in local vector order; may be empty.
        bindings : mapping of str to str
            Mapping from every local name to a registry ID. Repeated destination IDs
            express equality constraints; parameter units must already agree.

        Returns
        -------
        None
            Store local names and a read-only int64 local-to-global index array.

        Raises
        ------
        ValueError
            If local names or bindings are incomplete, repeated, or unknown.
        """
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
    """Gather a global parameter vector into the declared local order.

    Parameters
    ----------
    theta_global : array_like of shape (n_global,)
        Finite real global parameter values, in their respective physical units.
    local_to_global : array_like of int, shape (n_local,)
        Global index of each local parameter, including any repeated bindings.

    Returns
    -------
    theta_local : ndarray of shape (n_local,)
        Float64 local values with units inherited from their global parameters.

    Raises
    ------
    ValueError
        If the vector or indices are invalid.
    """
    theta = real_array(theta_global, "theta_global")
    if theta.ndim != 1:
        raise ValueError("global vector must be one-dimensional")
    index = indices(local_to_global, len(theta))
    return theta[index]


def map_jacobian(local_jacobian, local_to_global, n_global):
    """Sum local mean derivatives into their explicitly bound global columns.

    Parameters
    ----------
    local_jacobian : array_like of shape (n_node, n_pair, n_local)
        Finite derivatives in power units per local parameter unit.
    local_to_global : array_like of int, shape (n_local,)
        Global destination indices, allowing repeated equality bindings.
    n_global : int
        Positive number of global free parameters.

    Returns
    -------
    jacobian : ndarray of shape (n_node, n_pair, n_global)
        Float64 derivatives in power units per global parameter unit. Unused
        columns remain zero; repeated bindings accumulate by the chain rule.

    Raises
    ------
    ValueError
        If the Jacobian shape, index mapping, or global size is invalid.
    """
    n_global = integer(n_global, "n_global", 1)
    index = indices(local_to_global, n_global)
    local_derivatives = real_array(local_jacobian, "local_jacobian")
    if local_derivatives.ndim != 3 or local_derivatives.shape[2] != len(index):
        raise ValueError("Jacobian must have shape (n_node, n_pair, n_local)")
    result = np.zeros((*local_derivatives.shape[:2], n_global), dtype=np.float64)
    for local, global_ in enumerate(index):
        result[:, :, global_] += local_derivatives[:, :, local]
    return result
