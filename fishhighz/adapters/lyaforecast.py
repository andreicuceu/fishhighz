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
        """Bind an explicit intrinsic-power provider and unit conventions.

        Parameters
        ----------
        provider : object
            Object exposing compute_p3d_hmpc(z, k, mu, corr).
        routes : mapping
            Integer field-index pairs mapped to reference correlation labels.
        n_fields : int
            Positive number of observed fields.
        h_source : float
            Positive dimensionless reduced Hubble constant used by the source.
        h_fid : float
            Positive dimensionless fiducial reduced Hubble constant.
        k_domain : array_like
            Closed source wavenumber bounds in h_source/Mpc, shape (2,).
        z_domain : array_like
            Closed dimensionless redshift bounds, shape (2,).

        Returns
        -------
        None
            No value is returned.

        Raises
        ------
        ValueError
            If provider, routes, units or declared domains are invalid.
        """
        if not callable(getattr(provider, "compute_p3d_hmpc", None)):
            raise ValueError("provider must expose compute_p3d_hmpc")
        field_count = integer(n_fields, "n_fields", 1)
        prepared = {}
        for pair, name in routes.items():
            if len(pair) != 2 or not isinstance(name, str) or not name.strip():
                raise ValueError("routes require integer pairs and nonempty labels")
            key = tuple(integer(x, "field index") for x in pair)
            if max(key) >= field_count:
                raise ValueError("route outside field bounds")
            prepared[key] = name
        if not prepared:
            raise ValueError("routes must not be empty")
        domains = []
        for domain, name, minimum in ((k_domain, "k", 0), (z_domain, "z", -1)):
            domain_bounds = real_array(domain, name)
            if (
                domain_bounds.shape != (2,)
                or domain_bounds[0] <= minimum
                or domain_bounds[1] <= domain_bounds[0]
                or (name == "z" and domain_bounds[0] < 0)
            ):
                raise ValueError("invalid closed domain")
            domains.append(tuple(domain_bounds))
        for key, value in dict(
            provider=provider,
            routes=freeze(prepared),
            n_fields=field_count,
            h_source=_positive(h_source, "h_source"),
            h_fid=_positive(h_fid, "h_fid"),
            k_domain=domains[0],
            z_domain=domains[1],
        ).items():
            object.__setattr__(self, key, value)

    def __call__(self, theta_local, z, k, mu, pairs):
        """Evaluate intrinsic power after validating all requested coordinates.

        Parameters
        ----------
        theta_local : array_like
            Empty parameter vector, shape (0,).
        z : float
            Dimensionless evaluation redshift.
        k : array_like
            Paired observed wavenumbers in h_fid/Mpc, shape (n_node,).
        mu : array_like
            Paired line-of-sight cosines in [0, 1], shape (n_node,).
        pairs : array_like
            Integer observed-field pairs, shape (n_pair, 2).

        Returns
        -------
        power : ndarray
            Owned intrinsic power in (Mpc/h_fid)^3, shape (n_node, n_pair),
            preserving pair order.

        Raises
        ------
        ValueError
            If parameters, coordinates, pair routes or output shapes are invalid.

        Notes
        -----
        Converts k_source=k*h_fid/h_source and P_fid=P_source*(h_fid/h_source)^3. No response, noise or P1D is added.
        """
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
        pair_array = np.asarray(pairs, dtype=object)
        if pair_array.ndim != 2 or pair_array.shape[1] != 2 or not len(pair_array):
            raise ValueError("require nonempty integer pairs")
        keys = [tuple(integer(v, "field index") for v in row) for row in pair_array]
        if any(max(key) >= self.n_fields or key not in self.routes for key in keys):
            raise ValueError("pair outside bounds or missing explicit route")

        # Source and forecast h conventions rescale both k and power volume.
        try:
            with np.errstate(
                over="raise", under="raise", invalid="raise", divide="raise"
            ):
                hubble_ratio = np.float64(self.h_fid) / self.h_source
                source_wavenumbers = k * hubble_ratio
                power_conversion = hubble_ratio**3
        except FloatingPointError as error:
            raise ValueError("h conversion not representable") from error
        if not self.z_domain[0] <= z <= self.z_domain[1] or np.any(
            (source_wavenumbers < self.k_domain[0])
            | (source_wavenumbers > self.k_domain[1])
        ):
            raise ValueError("query outside declared source k/z domain")

        columns = []
        for key in keys:
            source_power = real_array(
                self.provider.compute_p3d_hmpc(
                    z, source_wavenumbers.copy(), mu.copy(), self.routes[key]
                ),
                "intrinsic output",
            )
            if source_power.shape != k.shape:
                raise ValueError("intrinsic output must match paired nodes")
            with np.errstate(over="raise", under="raise", invalid="raise"):
                columns.append(source_power * power_conversion)
        return np.array(np.column_stack(columns), copy=True)
