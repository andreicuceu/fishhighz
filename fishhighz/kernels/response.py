"""Array-only instrumental transfer; widths and queries are validated by host."""

import numpy as np


def _transfer(q, pixel_width, gaussian_sigma):
    x = q[:, None] * (pixel_width[None, :] / 2)
    gaussian = q[:, None] * gaussian_sigma[None, :]
    return np.sinc(x / np.pi) * np.exp(-0.5 * gaussian**2), x, gaussian
