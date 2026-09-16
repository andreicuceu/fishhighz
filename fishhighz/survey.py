"""Explicit host-side survey records; normalized-array use requires only NumPy."""

from dataclasses import dataclass
from types import MappingProxyType

import numpy as np

from ._arrays import label, real_array
from .geometry import BinGeometry, _immutable
from .grids import IntegrationGrid
from .models.external import BoundParameters, PreparedP3D, _state


def freeze(value):
    """Own input snapshots without casting ahead of scientific validation.

    Non-object array dtypes (including metadata indices/Booleans) are preserved.
    Object arrays are rejected: their elements cannot be frozen as byte-backed
    values. Scientific arrays become float64 only at their existing validator.
    """
    if isinstance(value, np.ndarray):
        if value.dtype.hasobject:
            raise ValueError(
                "object arrays cannot be frozen; use typed arrays or plain metadata"
            )
        return np.frombuffer(value.tobytes(), dtype=value.dtype).reshape(value.shape)
    if hasattr(value, "items"):
        return MappingProxyType({k: freeze(v) for k, v in value.items()})
    if isinstance(value, (tuple, list)):
        return tuple(freeze(v) for v in value)
    if value is None or isinstance(value, (str, int, float, complex, bool, np.number)):
        return value
    raise ValueError("provenance must contain only arrays and plain metadata")


@dataclass(frozen=True)
class ForestInput:
    """Explicit normalized weight kwargs and independent P1D state.

    weight_options are the keyword arguments of prepare_forest_weights. Optional
    auxiliary_coordinates=(k_t_deg,k_p_velocity) selects one callable S/B query
    for method='legacy'; omit it for direct scalar S/B or supplied weights.
    Reader provenance is plain metadata, never a live interpolator.
    """

    weight_options: object
    p1d_model: object
    p1d_parameters: BoundParameters
    theta_p1d: np.ndarray
    auxiliary_coordinates: tuple | None = None
    provenance: object = None

    def __post_init__(self):
        if not callable(self.p1d_model) or not isinstance(
            self.p1d_parameters, BoundParameters
        ):
            raise ValueError("forest requires independent P1D callable and binding")
        theta, _ = _state(self.p1d_parameters.registry, self.theta_p1d, 0)
        object.__setattr__(self, "theta_p1d", _immutable(theta))
        object.__setattr__(self, "weight_options", freeze(self.weight_options))
        object.__setattr__(self, "provenance", freeze(self.provenance))
        if not hasattr(self.weight_options, "keys"):
            raise ValueError("weight_options must be a mapping")
        if self.auxiliary_coordinates is not None:
            coords = tuple(self.auxiliary_coordinates)
            if (
                len(coords) != 2
                or self.weight_options.get("method") != "legacy"
                or any(
                    self.weight_options.get(k) is not None
                    for k in ("signal", "alias", "auxiliary", "weights")
                )
            ):
                raise ValueError("auxiliary coordinates conflict with weight settings")
            object.__setattr__(self, "auxiliary_coordinates", coords)


@dataclass(frozen=True)
class BinSpec:
    """A bin's exact geometry, grid, P3D selection, responses and noise choice.

    forests maps active forest IDs to ForestInput; galaxies maps active galaxy
    IDs to positive comoving n_bar. Full noise replaces both, and independence.
    """

    id: str
    geometry: BinGeometry
    grid: IntegrationGrid
    p3d: PreparedP3D
    responses: object
    forests: object = None
    galaxies: object = None
    independent_sampling: bool | None = None
    full_noise: np.ndarray | None = None

    def __post_init__(self):
        label(self.id, "bin ID")
        if (
            not isinstance(self.geometry, BinGeometry)
            or not isinstance(self.grid, IntegrationGrid)
            or not isinstance(self.p3d, PreparedP3D)
        ):
            raise ValueError("require BinGeometry, IntegrationGrid and PreparedP3D")
        if self.geometry.h_fid != self.grid.h_fid:
            raise ValueError("geometry/grid h_fid mismatch")
        for name in ("responses", "forests", "galaxies"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, MappingProxyType(dict(value)))
        if self.full_noise is not None:
            object.__setattr__(
                self,
                "full_noise",
                _immutable(real_array(self.full_noise, "full noise")),
            )


@dataclass(frozen=True)
class PreparedBin:
    """Fixed owned numerical state; no raw adapters or file handles.

    Full factors require O(n_node*n_selected**2) storage. Callables in p3d must
    be deterministic and not externally mutated; changed models/fiducials need
    fresh preparation. Arbitrary callable internals are not hashed or copied.
    """

    id: str
    geometry: BinGeometry
    p3d: PreparedP3D
    theta: np.ndarray
    k: np.ndarray
    mu: np.ndarray
    modes: np.ndarray
    response: np.ndarray
    products: np.ndarray
    noise: np.ndarray
    power: np.ndarray
    total: np.ndarray
    factors: np.ndarray
    weights: object
    diagnostics: object


def snapshot_p3d(prepared):
    """Own structural arrays while retaining exact registry and callable identities."""
    from copy import copy

    result = copy(prepared)
    selection = copy(prepared.selection)
    for name, value in vars(selection).items():
        if isinstance(value, np.ndarray):
            object.__setattr__(
                selection,
                name,
                np.frombuffer(value.tobytes(), dtype=value.dtype).reshape(value.shape),
            )
    object.__setattr__(result, "selection", selection)
    routes = []
    for route in prepared.routes:
        new = copy(route)
        provider = copy(route.provider)
        parameters = copy(provider.parameters)
        binding = copy(parameters.binding)
        index = binding.local_to_global
        object.__setattr__(
            binding,
            "local_to_global",
            np.frombuffer(index.tobytes(), dtype=index.dtype).reshape(index.shape),
        )
        object.__setattr__(parameters, "binding", binding)
        object.__setattr__(provider, "parameters", parameters)
        object.__setattr__(new, "provider", provider)
        for name in ("pairs", "columns"):
            value = getattr(route, name)
            object.__setattr__(
                new,
                name,
                np.frombuffer(value.tobytes(), dtype=value.dtype).reshape(value.shape),
            )
        routes.append(new)
    object.__setattr__(result, "routes", tuple(routes))
    return result
