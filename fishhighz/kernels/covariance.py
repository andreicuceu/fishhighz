"""Array-only Gaussian covariance accumulation, without boundary validation."""

import numpy as np


def _gaussian_covariance_kernel(total_power, mode_counts, im, jn, in_, jm, out):
    """Fill one triangle and mirror into out; return None.

    Inputs are validated float64 power (node, required), counts (node,), and
    int64 lookup tables (selected, selected). out is preallocated float64
    (node, selected, selected), does not alias inputs, and is overwritten.
    Only one additional node-vector is allocated. No physical or shape checks
    occur here; the public wrapper owns validation and arithmetic diagnostics.
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
    """Fill out (node, selected) with the diagonal of the Gaussian covariance.

    The arithmetic is exactly that of _gaussian_covariance_kernel for a == b, so
    each column equals the variance of an independently prepared one-spectrum
    covariance built from the same total power and mode counts.
    """
    work = np.empty(total_power.shape[0], dtype=np.float64)
    for a in range(im.shape[0]):
        target = out[:, a]
        np.multiply(total_power[:, im[a, a]], total_power[:, jn[a, a]], out=target)
        np.multiply(total_power[:, in_[a, a]], total_power[:, jm[a, a]], out=work)
        np.add(target, work, out=target)
        np.divide(target, mode_counts, out=target)
