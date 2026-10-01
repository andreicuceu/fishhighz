"""Array-only Gaussian covariance accumulation, without boundary validation."""

import numpy as np


def _gaussian_covariance_kernel(total_power, mode_counts, im, jn, in_, jm, out):
    """Fill the Gaussian covariance triangle and mirror it into the output.

    Parameters
    ----------
    total_power : ndarray of shape (n_node, n_required_pair)
        Validated float64 observed signal plus noise in (Mpc/h_fid)^3.
    mode_counts : ndarray of shape (n_node,)
        Positive dimensionless Fourier mode counts.
    im, jn, in_, jm : ndarray of int64, shape (n_selected, n_selected)
        Required-pair lookup tables for the two Gaussian covariance products.
    out : ndarray of shape (n_node, n_selected, n_selected)
        Preallocated float64 covariance destination; must not alias the inputs.

    Returns
    -------
    None
        Overwrite out with covariance in (Mpc/h_fid)^6.

    Notes
    -----
    The caller owns shape, physical, and arithmetic validation. Only one
    additional node-vector is allocated; signed cross products are retained.
    """
    work = np.empty(total_power.shape[0], dtype=np.float64)
    for a in range(im.shape[0]):
        for b in range(a, im.shape[1]):
            target = out[:, a, b]
            np.multiply(total_power[:, im[a, b]], total_power[:, jn[a, b]], out=target)
            np.multiply(total_power[:, in_[a, b]], total_power[:, jm[a, b]], out=work)
            np.add(target, work, out=target)
            np.divide(target, mode_counts, out=target)
            if a != b:
                out[:, b, a] = target


def _gaussian_variance_kernel(total_power, mode_counts, im, jn, in_, jm, out):
    """Fill the selected-spectrum variances with the Gaussian diagonal terms.

    Parameters
    ----------
    total_power : ndarray of shape (n_node, n_required_pair)
        Validated float64 observed signal plus noise in (Mpc/h_fid)^3.
    mode_counts : ndarray of shape (n_node,)
        Positive dimensionless Fourier mode counts.
    im, jn, in_, jm : ndarray of int64, shape (n_selected, n_selected)
        Required-pair lookup tables for the two Gaussian covariance products.
    out : ndarray of shape (n_node, n_selected)
        Preallocated float64 variance destination; must not alias the inputs.

    Returns
    -------
    None
        Overwrite out with variances in (Mpc/h_fid)^6.

    Notes
    -----
    The arithmetic matches the covariance kernel on each diagonal, preserving
    the variance of an independently prepared one-spectrum calculation.
    """
    work = np.empty(total_power.shape[0], dtype=np.float64)
    for a in range(im.shape[0]):
        target = out[:, a]
        np.multiply(total_power[:, im[a, a]], total_power[:, jn[a, a]], out=target)
        np.multiply(total_power[:, in_[a, a]], total_power[:, jm[a, a]], out=work)
        np.add(target, work, out=target)
        np.divide(target, mode_counts, out=target)
