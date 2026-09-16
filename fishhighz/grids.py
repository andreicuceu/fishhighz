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
        wk = np.repeat(w_k, mu.size)
        wm = np.tile(w_mu, k.size)
        with np.errstate(over="ignore", invalid="ignore", under="ignore"):
            weights = wk * wm
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
            ("w_k_flat", wk),
            ("w_mu_flat", wm),
            ("weights", weights),
            ("q_mode", q_mode),
        ):
            object.__setattr__(self, name, readonly(value, np.float64))


def gauss_legendre_grid(k_edges, *, k_order, mu_order, h_fid):
    """Construct interior Gauss–Legendre nodes in each k bin and mu in [0,1].

    Orders are explicit positive integers, with no universal convergence claim.
    Nodes are evaluation points, not bin-averaged bandpowers. The first/last
    positive, strictly increasing edge sets the fixed k cuts.
    """
    edges = _axis(k_edges, "k_edges")
    if edges.size < 2 or edges[0] <= 0:
        raise ValueError("at least two positive k edges are required")
    k_order = integer(k_order, "k_order", 1)
    mu_order = integer(mu_order, "mu_order", 1)
    xk, wk = np.polynomial.legendre.leggauss(k_order)
    xm, wm = np.polynomial.legendre.leggauss(mu_order)
    half_width = np.diff(edges) / 2
    k = (edges[:-1, None] + half_width[:, None] * (xk + 1)).ravel()
    per_bin = k.reshape(-1, k_order)
    if np.any(per_bin <= edges[:-1, None]) or np.any(per_bin >= edges[1:, None]):
        raise ValueError("k bins cannot represent interior quadrature nodes")
    weights = (half_width[:, None] * wk).ravel()
    return IntegrationGrid(
        k,
        weights,
        (xm + 1) / 2,
        wm / 2,
        k_min=edges[0],
        k_max=edges[-1],
        h_fid=h_fid,
    )
