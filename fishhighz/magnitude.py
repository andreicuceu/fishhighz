"""Magnitude partitions and composite Gauss--Legendre quadrature.

These functions preserve the adopted accuracy algorithm exactly.  They accept
the production reader/adapters and do not depend on validation or legacy
forecast packages.
"""

import numpy as np


def composite(partition, order):
    """Return ordered Gauss--Legendre nodes and weights on ``partition``."""

    x, w = np.polynomial.legendre.leggauss(int(order))
    lo, hi = np.asarray(partition[:-1]), np.asarray(partition[1:])
    return ((lo[:, None] + hi[:, None]) / 2 + (hi - lo)[:, None] * x / 2).ravel(), (
        (hi - lo)[:, None] * w / 2
    ).ravel()


def breakpoints(densities, snrs, z_queries, lo, hi):
    """Build the adopted density/support/SNR magnitude partition.

    Quadratic roots are located only to partition the existing density spline;
    interpolation and negative-value policy remain owned by the reader adapter.
    """

    points = [lo, hi]
    for name, density in densities.items():
        reader = getattr(density, "reader", density)
        magnitude_axis = getattr(reader, "magnitudes", None)
        spline = getattr(density, "_spline", None)
        if magnitude_axis is None:
            continue
        spline_knots = spline.get_knots()[1] if spline is not None else magnitude_axis
        knots = np.unique(np.r_[magnitude_axis, spline_knots, lo, hi])
        knots = knots[(knots >= lo) & (knots <= hi)]
        points.extend(knots)
        for a, b in zip(knots[:-1], knots[1:]):
            if a < magnitude_axis[0] or b > magnitude_axis[-1] or spline is None:
                continue
            y = spline.ev(np.full(3, z_queries[name]), [a, (a + b) / 2, b])
            c = y[0]
            aa = 2 * (y[2] - 2 * y[1] + y[0])
            bb = y[2] - y[0] - aa
            roots = np.roots([aa, bb, c]) if aa != 0 else ([-c / bb] if bb != 0 else [])
            for root in roots:
                if np.isreal(root) and 1e-12 < float(np.real(root)) < 1 - 1e-12:
                    points.append(a + (b - a) * float(np.real(root)))
    for snr in snrs.values():
        reader = getattr(snr, "reader", snr)
        if hasattr(reader, "magnitudes"):
            points.extend(reader.magnitudes)
    p = np.unique(np.asarray(points))
    p = p[(p >= lo) & (p <= hi)]
    return p[
        np.r_[True, np.diff(p) > 64 * np.finfo(float).eps * np.maximum(1, abs(p[1:]))]
    ]


__all__ = ["breakpoints", "composite"]
