"""Fixed tensor quadrature in fiducial coordinates for real, even spectra.

k is in h_fid/Mpc; P3D and volume are in (Mpc/h_fid)^3. Adapters own unit
conversion and keep h_fid fixed even when physical h varies. Model-domain
padding is independent of this integration grid and adds no forecast modes.
Initial Fisher calculations will hold covariance and survey weights fiducial,
differentiating only the mean. No Fisher calculation is implemented here.
"""

from dataclasses import dataclass

import numpy as np

from ._arrays import integer, readonly, real_array, scalar


def _axis(values, name):
    """Validate a nonempty, strictly increasing quadrature axis.

    Parameters
    ----------
    values : array_like of shape (n_node,)
        Real coordinate values in the units of the named axis.
    name : str
        Axis name used in validation errors.

    Returns
    -------
    array : ndarray of shape (n_node,)
        Owned float64 axis with unchanged units and ordering.

    Raises
    ------
    ValueError
        If the axis is empty, nonfinite, not one-dimensional, or not increasing.
    """
    array = real_array(values, name)
    if array.ndim != 1 or not array.size or np.any(np.diff(array) <= 0):
        raise ValueError(f"{name} must be nonempty, 1D, and strictly increasing")
    return array


@dataclass(frozen=True, init=False, eq=False)
class IntegrationGrid:
    """Custom tensor axes with positive dk/dmu weights and fixed cuts.

    Parameters
    ----------
    k, w_k, mu, w_mu : array_like
        Ordered 1D coordinates and matching positive integration weights.
        Endpoints are allowed. Weight sums must match k_max-k_min and 1 with
        relative tolerance 1e-12 and zero absolute tolerance. No sorting,
        clipping, or normalization is performed. Scientific convergence still
        needs testing even for a positive, correctly normalized rule.
    k_min, k_max : float
        Positive fixed observed-coordinate cuts enclosing all k nodes.
    h_fid : float
        Positive reference h metadata; coordinates are not rescaled.

    Notes
    -----
    All stored arrays own read-only C-contiguous float64 data. Flattening is
    C-order (n_k, n_mu), mu fastest: node = i_k*n_mu + i_mu. Product weights
    integrate dk*dmu. q_mode counts both conjugate hemispheres on mu in [0,1];
    N_modes = V_fid*q_mode, with no additional half factor.
    """

    k: np.ndarray
    w_k: np.ndarray
    mu: np.ndarray
    w_mu: np.ndarray
    k_min: float
    k_max: float
    h_fid: float
    n_k: int
    n_mu: int
    k_flat: np.ndarray
    mu_flat: np.ndarray
    w_k_flat: np.ndarray
    w_mu_flat: np.ndarray
    weights: np.ndarray
    q_mode: np.ndarray

    def __init__(self, k, w_k, mu, w_mu, *, k_min, k_max, h_fid):
        """Prepare fixed Fourier coordinates and positive tensor quadrature.

        Parameters
        ----------
        k, w_k : array_like of shape (n_k,)
            Increasing wavenumber nodes and positive dk weights, both in h_fid/Mpc.
        mu, w_mu : array_like of shape (n_mu,)
            Increasing direction cosines on [0, 1] and positive dimensionless
            dmu weights.
        k_min, k_max : float
            Positive fixed observed-coordinate cuts in h_fid/Mpc.
        h_fid : float
            Positive dimensionless reference Hubble parameter.

        Returns
        -------
        None
            Store owned read-only axes, flattened coordinates and quadrature weights.
            Flattening uses C order with mu fastest. The stored q_mode has units
            (h_fid/Mpc)^3; multiplication by the fiducial volume gives mode counts.

        Raises
        ------
        ValueError
            If coordinates, bounds, or weights are invalid, or mode weights cannot
            be represented as positive float64 values.

        Notes
        -----
        Weight sums must match the declared integration intervals to relative
        precision 1e-12. No sorting or renormalization is performed.
        """
        k_min, k_max, h_fid = (
            scalar(x, name)
            for x, name in ((k_min, "k_min"), (k_max, "k_max"), (h_fid, "h_fid"))
        )
        if not 0 < k_min < k_max or h_fid <= 0:
            raise ValueError("require 0 < k_min < k_max and h_fid > 0")
        k, mu = _axis(k, "k"), _axis(mu, "mu")
        if np.any((k < k_min) | (k > k_max)) or np.any((mu < 0) | (mu > 1)):
            raise ValueError("coordinates outside integration cuts")
        w_k, w_mu = real_array(w_k, "w_k"), real_array(w_mu, "w_mu")
        for axis, weight, interval in ((k, w_k, k_max - k_min), (mu, w_mu, 1)):
            if weight.shape != axis.shape or np.any(weight <= 0):
                raise ValueError("weights must match their axis and be positive")
            if not np.isclose(weight.sum(), interval, rtol=1e-12, atol=0):
                raise ValueError("weights do not integrate the declared interval")
        k_flat = np.repeat(k, mu.size)
        mu_flat = np.tile(mu, k.size)
        k_weights_flat = np.repeat(w_k, mu.size)
        mu_weights_flat = np.tile(w_mu, k.size)
        with np.errstate(over="ignore", invalid="ignore", under="ignore"):
            weights = k_weights_flat * mu_weights_flat
            q_mode = k_flat**2 * weights / (2 * np.pi**2)
        if not np.all(np.isfinite(q_mode)) or np.any(q_mode <= 0):
            raise ValueError("mode weights cannot be represented as positive float64")
        for name, value in (("k_min", k_min), ("k_max", k_max), ("h_fid", h_fid)):
            object.__setattr__(self, name, value)
        object.__setattr__(self, "n_k", k.size)
        object.__setattr__(self, "n_mu", mu.size)
        for name, value in (
            ("k", k),
            ("mu", mu),
            ("w_k", w_k),
            ("w_mu", w_mu),
            ("k_flat", k_flat),
            ("mu_flat", mu_flat),
            ("w_k_flat", k_weights_flat),
            ("w_mu_flat", mu_weights_flat),
            ("weights", weights),
            ("q_mode", q_mode),
        ):
            object.__setattr__(self, name, readonly(value, np.float64))


def gauss_legendre_grid(k_edges, *, k_order, mu_order, h_fid):
    """Construct Gauss–Legendre nodes within each k bin and over mu in [0, 1].

    Parameters
    ----------
    k_edges : array_like of shape (n_bin + 1,)
        Positive, strictly increasing observed wavenumber edges in h_fid/Mpc.
    k_order : int
        Positive quadrature order within each wavenumber bin.
    mu_order : int
        Positive angular quadrature order over [0, 1].
    h_fid : float
        Positive dimensionless reference Hubble parameter.

    Returns
    -------
    grid : IntegrationGrid
        Fixed tensor grid with n_bin*k_order radial and mu_order angular nodes.

    Raises
    ------
    ValueError
        If edges or orders are invalid, or interior nodes cannot be represented.

    Notes
    -----
    Nodes are evaluation points, not bin-averaged bandpowers. The chosen orders
    do not by themselves establish numerical convergence.
    """
    edges = _axis(k_edges, "k_edges")
    if edges.size < 2 or edges[0] <= 0:
        raise ValueError("at least two positive k edges are required")
    k_order = integer(k_order, "k_order", 1)
    mu_order = integer(mu_order, "mu_order", 1)
    k_legendre_nodes, k_legendre_weights = np.polynomial.legendre.leggauss(k_order)
    mu_legendre_nodes, mu_legendre_weights = np.polynomial.legendre.leggauss(mu_order)
    half_width = np.diff(edges) / 2
    k_grid = (edges[:-1, None] + half_width[:, None] * (k_legendre_nodes + 1)).ravel()
    per_bin = k_grid.reshape(-1, k_order)
    if np.any(per_bin <= edges[:-1, None]) or np.any(per_bin >= edges[1:, None]):
        raise ValueError("k bins cannot represent interior quadrature nodes")
    weights = (half_width[:, None] * k_legendre_weights).ravel()
    return IntegrationGrid(
        k_grid,
        weights,
        (mu_legendre_nodes + 1) / 2,
        mu_legendre_weights / 2,
        k_min=edges[0],
        k_max=edges[-1],
        h_fid=h_fid,
    )
