"""Magnitude partitions and composite Gauss--Legendre quadrature.

These functions preserve the adopted accuracy algorithm exactly.  They accept
the production reader/adapters and do not depend on validation or legacy
forecast packages.
"""

import numpy as np


def composite(partition, order):
    """Construct magnitude quadrature on an explicit interval partition.

    Parameters
    ----------
    partition : sequence of float, shape (n_interval + 1,)
        Ordered magnitude boundaries.
    order : int
        Gauss–Legendre order within each interval.

    Returns
    -------
    magnitudes : ndarray of shape (n_interval*order,)
        Interior magnitude nodes in interval order.
    weights : ndarray of shape (n_interval*order,)
        Corresponding integration weights, in magnitudes.

    Notes
    -----
    The caller supplies the partition. Nodes and weights retain the original
    interval-by-interval quadrature arithmetic.
    """

    legendre_nodes, legendre_weights = np.polynomial.legendre.leggauss(int(order))
    lower_magnitudes, upper_magnitudes = (
        np.asarray(partition[:-1]),
        np.asarray(partition[1:]),
    )
    return (
        (lower_magnitudes[:, None] + upper_magnitudes[:, None]) / 2
        + (upper_magnitudes - lower_magnitudes)[:, None] * legendre_nodes / 2
    ).ravel(), (
        (upper_magnitudes - lower_magnitudes)[:, None] * legendre_weights / 2
    ).ravel()


def breakpoints(densities, snrs, z_queries, lo, hi):
    """Partition magnitude integration at density, support, and SNR boundaries.

    Parameters
    ----------
    densities : mapping
        Named density readers or adapters, optionally exposing a quadratic
        spline.
    snrs : mapping
        SNR readers or adapters exposing tabulated magnitude nodes.
    z_queries : mapping of str to float
        Dimensionless redshift query for each named density reader.
    lo, hi : float
        Lower and upper magnitude integration limits.

    Returns
    -------
    partition : ndarray of shape (n_boundary,)
        Sorted magnitude boundaries after merging numerically indistinguishable
        nodes.

    Notes
    -----
    Spline roots partition the existing interpolation; they do not change
    the reader's density interpolation or negative-value policy. Interior roots
    exclude a relative 1e-12 endpoint neighbourhood. Adjacent boundaries are
    merged with the existing 64*eps magnitude criterion.
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
            density_samples = spline.ev(
                np.full(3, z_queries[name]), [a, (a + b) / 2, b]
            )
            constant_coefficient = density_samples[0]
            quadratic_coefficient = 2 * (
                density_samples[2] - 2 * density_samples[1] + density_samples[0]
            )
            linear_coefficient = (
                density_samples[2] - density_samples[0] - quadratic_coefficient
            )
            roots = (
                np.roots(
                    [quadratic_coefficient, linear_coefficient, constant_coefficient]
                )
                if quadratic_coefficient != 0
                else (
                    [-constant_coefficient / linear_coefficient]
                    if linear_coefficient != 0
                    else []
                )
            )
            for root in roots:
                if np.isreal(root) and 1e-12 < float(np.real(root)) < 1 - 1e-12:
                    points.append(a + (b - a) * float(np.real(root)))
    for snr in snrs.values():
        reader = getattr(snr, "reader", snr)
        if hasattr(reader, "magnitudes"):
            points.extend(reader.magnitudes)
    magnitude_partition = np.unique(np.asarray(points))
    magnitude_partition = magnitude_partition[
        (magnitude_partition >= lo) & (magnitude_partition <= hi)
    ]
    return magnitude_partition[
        np.r_[
            True,
            np.diff(magnitude_partition)
            > 64 * np.finfo(float).eps * np.maximum(1, abs(magnitude_partition[1:])),
        ]
    ]


__all__ = ["breakpoints", "composite"]
