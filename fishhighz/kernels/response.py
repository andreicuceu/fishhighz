"""Array-only instrumental transfer; widths and queries are validated by host."""

import numpy as np


def _transfer(q, pixel_width, gaussian_sigma):
    """Evaluate pixel integration and Gaussian instrumental smoothing.

    Parameters
    ----------
    q : ndarray of shape (n_node,)
        Line-of-sight wavenumbers in (km/s)^-1.
    pixel_width, gaussian_sigma : ndarray of shape (n_field,)
        Pixel widths and Gaussian standard deviations in km/s.

    Returns
    -------
    transfer : ndarray of shape (n_node, n_field)
        Signed, dimensionless amplitude response; square it for auto power.
    pixel_phase, gaussian_phase : ndarray of shape (n_node, n_field)
        Dimensionless q*pixel_width/2 and q*gaussian_sigma.

    Notes
    -----
    NumPy's sinc includes a factor of pi in its argument convention.
    """
    pixel_phase = q[:, None] * (pixel_width[None, :] / 2)
    gaussian_phase = q[:, None] * gaussian_sigma[None, :]
    return (
        np.sinc(pixel_phase / np.pi) * np.exp(-0.5 * gaussian_phase**2),
        pixel_phase,
        gaussian_phase,
    )
