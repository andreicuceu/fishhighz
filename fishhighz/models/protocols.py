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
    ) -> ArrayLike:
        """Evaluate intrinsic three-dimensional clustering power.

        Parameters
        ----------
        theta_local : ndarray of shape (n_local,)
            Local parameter values in binding order, in their declared physical
            units.
        z : float
            Scalar dimensionless evaluation redshift.
        k : ndarray of shape (n_node,)
            Positive observed wavenumbers in h_fid/Mpc.
        mu : ndarray of shape (n_node,)
            Paired dimensionless direction cosines.
        pairs : ndarray of int, shape (n_pair, 2)
            Original field indices, in requested spectrum order.

        Returns
        -------
        power : array_like of shape (n_node, n_pair)
            Real intrinsic clustering power in (Mpc/h_fid)^3; signed cross powers
            are allowed.

        Notes
        -----
        Providers own cosmological and dilation physics. Instrumental response
        and known sampling noise are handled explicitly by the consumer.
        """
        ...


class P1D(Protocol):
    """Return power_1d (n_node,) at the requested velocity wavenumbers."""

    def __call__(
        self, theta_local: np.ndarray, z: float, k_parallel_velocity: np.ndarray
    ) -> ArrayLike:
        """Evaluate intrinsic one-dimensional forest power independently of P3D.

        Parameters
        ----------
        theta_local : ndarray of shape (n_local,)
            Local parameter values in binding order, in their declared physical
            units.
        z : float
            Scalar dimensionless evaluation redshift.
        k_parallel_velocity : ndarray of shape (n_node,)
            Line-of-sight wavenumbers in s/km.

        Returns
        -------
        power : array_like of shape (n_node,)
            Intrinsic one-dimensional forest power in km/s.

        Notes
        -----
        No relation to an external three-dimensional model or comoving conversion
        is inferred from this callable contract.
        """
        ...


class P3DJacobian(Protocol):
    """Return signed derivatives (n_node,n_pair,n_local) in binding order."""

    def __call__(
        self,
        theta_local: np.ndarray,
        z: float,
        k: np.ndarray,
        mu: np.ndarray,
        pairs: np.ndarray,
    ) -> ArrayLike:
        """Evaluate derivatives of intrinsic power with respect to local parameters.

        Parameters
        ----------
        theta_local : ndarray of shape (n_local,)
            Local parameter values in binding order, in their declared physical
            units.
        z : float
            Scalar dimensionless evaluation redshift.
        k : ndarray of shape (n_node,)
            Positive observed wavenumbers in h_fid/Mpc.
        mu : ndarray of shape (n_node,)
            Paired dimensionless direction cosines.
        pairs : ndarray of int, shape (n_pair, 2)
            Original field indices, in requested spectrum order.

        Returns
        -------
        jacobian : array_like of shape (n_node, n_pair, n_local)
            Signed local derivatives in (Mpc/h_fid)^3 per parameter unit.

        Notes
        -----
        Local columns follow binding order. Equality-bound local derivatives are
        summed into global columns by the consumer.
        """
        ...


def _output(value, shape):
    """Copy a provider result and enforce its exact shape.

    Parameters
    ----------
    value : array_like
        Real finite provider values in the relevant model units.
    shape : tuple of int
        Required result shape; broadcasting is not permitted.

    Returns
    -------
    output : ndarray
        Owned C-contiguous float64 copy with the requested shape and unchanged
        units.

    Raises
    ------
    ValueError
        If dtype, finiteness, or exact shape validation fails.
    """
    output = real_array(value, "provider output")
    if output.shape != shape:
        raise ValueError(f"expected output shape {shape}, got {output.shape}")
    return output


def validate_p3d(value, n_node, n_pair):
    """Validate and copy intrinsic three-dimensional power without changing units.

    Parameters
    ----------
    value : array_like of shape (n_node, n_pair)
        Finite real intrinsic three-dimensional power in (Mpc/h_fid)^3.
    n_node : int
        Positive number of Fourier nodes.
    n_pair : int
        Positive number of requested spectra.

    Returns
    -------
    output : ndarray of shape (n_node, n_pair)
        Owned C-contiguous float64 values in (Mpc/h_fid)^3.

    Raises
    ------
    ValueError
        If dimensions, dtype, values, or exact output shape are invalid.
    """
    return _output(value, (integer(n_node, "n_node", 1), integer(n_pair, "n_pair", 1)))


def validate_p1d(value, n_node):
    """Validate and copy intrinsic one-dimensional power without changing units.

    Parameters
    ----------
    value : array_like of shape (n_node,)
        Finite real intrinsic one-dimensional power in km/s.
    n_node : int
        Positive number of Fourier nodes.

    Returns
    -------
    output : ndarray of shape (n_node,)
        Owned C-contiguous float64 values in km/s.

    Raises
    ------
    ValueError
        If dimensions, dtype, values, or exact output shape are invalid.
    """
    return _output(value, (integer(n_node, "n_node", 1),))


def validate_p3d_jacobian(value, n_node, n_pair, n_local):
    """Validate and copy local power derivatives without changing units.

    Parameters
    ----------
    value : array_like of shape (n_node, n_pair, n_local)
        Finite real local power derivatives in (Mpc/h_fid)^3 per local parameter
        unit.
    n_node : int
        Positive number of Fourier nodes.
    n_pair : int
        Positive number of requested spectra.
    n_local : int
        Nonnegative number of local free parameters.

    Returns
    -------
    output : ndarray of shape (n_node, n_pair, n_local)
        Owned C-contiguous float64 values in (Mpc/h_fid)^3 per local parameter
        unit.

    Raises
    ------
    ValueError
        If dimensions, dtype, values, or exact output shape are invalid.
    """
    return _output(
        value,
        (
            integer(n_node, "n_node", 1),
            integer(n_pair, "n_pair", 1),
            integer(n_local, "n_local"),
        ),
    )
