"""Array-only lookup and Horner evaluation for prepared log-k cubic splines."""

import numpy as np


def _intervals(log_knots, log_query):
    """Locate intervals; the closed upper endpoint belongs to the last interval.

    The caller validates the physical k domain before taking logarithms.
    """
    return np.minimum(
        np.searchsorted(log_knots, log_query, side="right") - 1,
        len(log_knots) - 2,
    )


def _evaluate(log_knots, coefficients, k_query, derivative):
    """Evaluate (query,component) power or dP/dk, with no query-by-knot matrix.

    Coefficients have shape (4,interval,component), descending powers of
    dx=ln(k)-log_knots[interval]. derivative is the integer 0 or 1.
    """
    log_query = np.log(k_query)
    index = _intervals(log_knots, log_query)
    dx = (log_query - log_knots[index])[:, None]
    if derivative == 0:
        return (
            (coefficients[0, index] * dx + coefficients[1, index]) * dx
            + coefficients[2, index]
        ) * dx + coefficients[3, index]
    return (
        (3 * coefficients[0, index] * dx + 2 * coefficients[1, index]) * dx
        + coefficients[2, index]
    ) / k_query[:, None]
