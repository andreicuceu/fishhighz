"""Magnitude partitions and composite Gauss--Legendre quadrature.

These functions preserve the adopted accuracy algorithm exactly.  They accept
the production survey tables and samplers and do not depend on validation or
legacy forecast packages.
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


def breakpoints(densities, snrs, lo, hi):
    """Partition magnitude integration at density, support, and SNR boundaries.

    Parameters
    ----------
    densities : mapping
        Named density tables or samplers, optionally exposing a
        piecewise-constant cell histogram (``_histogram``).
    snrs : mapping
        SNR tables or samplers exposing tabulated magnitude nodes.
    lo, hi : float
        Lower and upper magnitude integration limits.

    Returns
    -------
    partition : ndarray of shape (n_boundary,)
        Sorted magnitude boundaries after merging numerically indistinguishable
        nodes.

    Notes
    -----
    Piecewise-constant densities contribute their magnitude cell edges only, so
    Gauss–Legendre integration of the density alone is exact within each
    interval. A density without cell edges contributes its tabulated magnitude
    nodes. Adjacent boundaries are merged with the existing 64*eps magnitude
    criterion.
    """

    points = [lo, hi]
    for density in densities.values():
        reader = getattr(density, "reader", density)
        magnitude_axis = getattr(reader, "magnitudes", None)
        histogram = getattr(density, "_histogram", None)
        if histogram is not None:
            # Constant cells: the density is discontinuous only at cell edges.
            edges = np.asarray(histogram.magnitude_edges)
            points.extend(edges[(edges >= lo) & (edges <= hi)])
            continue
        if magnitude_axis is None:
            continue
        # Without cell edges the tabulated magnitude nodes partition the range.
        knots = np.unique(np.r_[magnitude_axis, lo, hi])
        points.extend(knots[(knots >= lo) & (knots <= hi)])
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
