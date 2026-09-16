"""Explicit intrinsic-P3D bridge; never imports or configures lyaforecast."""

from dataclasses import dataclass

import numpy as np

from .._arrays import integer, real_array, scalar
from ..geometry import _positive
from ..survey import freeze


@dataclass(frozen=True, init=False)
class IntrinsicP3D:
    """Wrap caller-prepared compute_p3d_hmpc(z,k,mu,corr).

    routes maps original integer (i,j) pairs to reference correlation strings.
    k_domain uses source h/Mpc; z_domain is closed. Unequal h is explicitly
    converted with k_source=k_fid*h_fid/h_source and P_fid=P_source*ratio**3.
    The object has zero local parameters and owns no response, noise or P1D.
    """

    provider: object
    routes: object
    n_fields: int
    h_source: float
    h_fid: float
    k_domain: tuple
    z_domain: tuple

    def __init__(
        self, provider, *, routes, n_fields, h_source, h_fid, k_domain, z_domain
    ):
        if not callable(getattr(provider, "compute_p3d_hmpc", None)):
            raise ValueError("provider must expose compute_p3d_hmpc")
        n = integer(n_fields, "n_fields", 1)
        prepared = {}
        for pair, name in routes.items():
            if len(pair) != 2 or not isinstance(name, str) or not name.strip():
                raise ValueError("routes require integer pairs and nonempty labels")
            key = tuple(integer(x, "field index") for x in pair)
            if max(key) >= n:
                raise ValueError("route outside field bounds")
            prepared[key] = name
        if not prepared:
            raise ValueError("routes must not be empty")
        domains = []
        for domain, name, minimum in ((k_domain, "k", 0), (z_domain, "z", -1)):
            a = real_array(domain, name)
            if (
                a.shape != (2,)
                or a[0] <= minimum
                or a[1] <= a[0]
                or (name == "z" and a[0] < 0)
            ):
                raise ValueError("invalid closed domain")
            domains.append(tuple(a))
        for key, value in dict(
            provider=provider,
            routes=freeze(prepared),
            n_fields=n,
            h_source=_positive(h_source, "h_source"),
            h_fid=_positive(h_fid, "h_fid"),
            k_domain=domains[0],
            z_domain=domains[1],
        ).items():
            object.__setattr__(self, key, value)

    def __call__(self, theta_local, z, k, mu, pairs):
        """Validate all queries before any reference call; preserve pair order."""
        if real_array(theta_local, "theta_local").shape != (0,):
            raise ValueError("intrinsic bridge has zero local parameters")
        z = scalar(z, "z")
        k, mu = real_array(k, "k"), real_array(mu, "mu")
        if (
            k.ndim != 1
            or not len(k)
            or mu.shape != k.shape
            or np.any((mu < 0) | (mu > 1))
        ):
            raise ValueError("require paired 1D k/mu with mu in [0,1]")
        p = np.asarray(pairs, dtype=object)
        if p.ndim != 2 or p.shape[1] != 2 or not len(p):
            raise ValueError("require nonempty integer pairs")
        keys = [tuple(integer(v, "field index") for v in row) for row in p]
        if any(max(key) >= self.n_fields or key not in self.routes for key in keys):
            raise ValueError("pair outside bounds or missing explicit route")
        try:
            with np.errstate(
                over="raise", under="raise", invalid="raise", divide="raise"
            ):
                ratio = np.float64(self.h_fid) / self.h_source
                source = k * ratio
                factor = ratio**3
        except FloatingPointError as error:
            raise ValueError("h conversion not representable") from error
        if not self.z_domain[0] <= z <= self.z_domain[1] or np.any(
            (source < self.k_domain[0]) | (source > self.k_domain[1])
        ):
            raise ValueError("query outside declared source k/z domain")
        columns = []
        for key in keys:
            raw = real_array(
                self.provider.compute_p3d_hmpc(
                    z, source.copy(), mu.copy(), self.routes[key]
                ),
                "intrinsic output",
            )
            if raw.shape != k.shape:
                raise ValueError("intrinsic output must match paired nodes")
            with np.errstate(over="raise", under="raise", invalid="raise"):
                columns.append(raw * factor)
        return np.array(np.column_stack(columns), copy=True)
