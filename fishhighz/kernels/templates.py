"""Array-only lookup and Horner evaluation for prepared log-k cubic splines."""

import numpy as np


def _intervals(log_knots, log_query):
    """Locate log-wavenumber spline intervals including the upper endpoint.

    Parameters
    ----------
    log_knots : ndarray of shape (n_knot,)
        Increasing logarithms of knot wavenumbers expressed in h/Mpc.
    log_query : ndarray of shape (n_node,)
        Logarithms of query wavenumbers in the same convention.

    Returns
    -------
    indices : ndarray of int, shape (n_node,)
        Spline interval indices; the closed upper endpoint belongs to the
        final interval. The caller validates the physical wavenumber domain.
    """
    return np.minimum(
        np.searchsorted(log_knots, log_query, side="right") - 1,
        len(log_knots) - 2,
    )


def _evaluate(log_knots, coefficients, k_query, derivative):
    """Evaluate a prepared cubic power spline or its wavenumber derivative.

    Parameters
    ----------
    log_knots : ndarray of shape (n_knot,)
        Increasing logarithms of knot wavenumbers expressed in h/Mpc.
    coefficients : ndarray of shape (4, n_interval, n_component)
        Cubic coefficients in descending powers of the log-wavenumber offset.
        Coefficients carry the template power units, (Mpc/h)^3.
    k_query : ndarray of shape (n_node,)
        Query wavenumbers in h/Mpc, within the validated template domain.
    derivative : int
        Return power for 0 or dP/dk otherwise; callers restrict this to 0 or 1.

    Returns
    -------
    values : ndarray of shape (n_node, n_component)
        Power in (Mpc/h)^3 or dP/dk in (Mpc/h)^4.

    Notes
    -----
    Horner evaluation avoids a query-by-knot matrix. For the derivative,
    division by k converts dP/dln(k) to dP/dk.
    """
    log_query = np.log(k_query)
    index = _intervals(log_knots, log_query)
    log_k_offset = (log_query - log_knots[index])[:, None]
    if derivative == 0:
        return (
            (coefficients[0, index] * log_k_offset + coefficients[1, index])
            * log_k_offset
            + coefficients[2, index]
        ) * log_k_offset + coefficients[3, index]
    return (
        (3 * coefficients[0, index] * log_k_offset + 2 * coefficients[1, index])
        * log_k_offset
        + coefficients[2, index]
    ) / k_query[:, None]
