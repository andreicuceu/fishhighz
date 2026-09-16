"""Thin callable contracts: no inheritance, autodiff, or FishHighz import required.

Evaluate a scalar redshift z and a local vector in binding order. P3D uses paired
1D k/mu evaluation points (not necessarily tensor axes) and prepared field-index
pairs (n_pair,2). On an IntegrationGrid, nodes follow C-order (n_k,n_mu).
P3D is real clustering power in (Mpc/h_fid)^3 at fixed fiducial k in h_fid/Mpc,
before instrument response and known sampling noise. Signed cross powers are
valid. Providers/adapters own cosmological/AP physics; it must not be reapplied
implicitly. Outputs already containing response/noise need explicit adapter
handling to avoid double counting; there is no automatic interpretation here.

P1D independently uses k_parallel_velocity in s/km and power in km/s. Custom
P3D does not require custom P1D, imply it, or supply a comoving conversion.
Validators check outputs only, not domains, noise positivity, or PSD science.
"""

from typing import Protocol

import numpy as np
from numpy.typing import ArrayLike

from .._arrays import integer, real_array


class P3D(Protocol):
    """Return real power (n_node,n_pair) in requested pair order."""

    def __call__(
        self,
        theta_local: np.ndarray,
        z: float,
        k: np.ndarray,
        mu: np.ndarray,
        pairs: np.ndarray,
    ) -> ArrayLike: ...


class P1D(Protocol):
    """Return power_1d (n_node,) at the requested velocity wavenumbers."""

    def __call__(
        self, theta_local: np.ndarray, z: float, k_parallel_velocity: np.ndarray
    ) -> ArrayLike: ...


class P3DJacobian(Protocol):
    """Return signed derivatives (n_node,n_pair,n_local) in binding order."""

    def __call__(
        self,
        theta_local: np.ndarray,
        z: float,
        k: np.ndarray,
        mu: np.ndarray,
        pairs: np.ndarray,
    ) -> ArrayLike: ...


def _output(value, shape):
    output = real_array(value, "provider output")
    if output.shape != shape:
        raise ValueError(f"expected output shape {shape}, got {output.shape}")
    return output


def validate_p3d(value, n_node, n_pair):
    """Copy real finite power into exact (node,pair), C-order float64."""
    return _output(value, (integer(n_node, "n_node", 1), integer(n_pair, "n_pair", 1)))


def validate_p1d(value, n_node):
    """Copy real finite P1D into exact (node,), C-order float64."""
    return _output(value, (integer(n_node, "n_node", 1),))


def validate_p3d_jacobian(value, n_node, n_pair, n_local):
    """Copy signed finite derivatives into exact (node,pair,local) float64."""
    return _output(
        value,
        (
            integer(n_node, "n_node", 1),
            integer(n_pair, "n_pair", 1),
            integer(n_local, "n_local"),
        ),
    )
