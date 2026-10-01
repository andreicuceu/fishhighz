"""Array-only derivative combinations; dispatch and coefficients live outside."""


def _combine_three(fiducial, first, second, weights, scale):
    """Combine three model evaluations into a finite-difference derivative.

    Parameters
    ----------
    fiducial, first, second : ndarray of shape (n_node, n_pair)
        Predicted powers at the fiducial and two displaced parameter values,
        all in the same power units.
    weights : ndarray of shape (2,)
        Dimensionless coefficients for the actual scaled parameter offsets.
    scale : float
        Parameter displacement scale, in the differentiated parameter's units.

    Returns
    -------
    derivative : ndarray of shape (n_node, n_pair)
        Power derivative in power units per parameter unit.

    Notes
    -----
    Subtract the fiducial before weighting to cancel constant model terms.
    """
    return (weights[0] * (first - fiducial) + weights[1] * (second - fiducial)) / scale


def _scatter_column(output, values, columns, global_index):
    """Assign one provider's contribution to a global Jacobian column.

    Parameters
    ----------
    output : ndarray of shape (n_node, n_pair, n_global_parameter)
        Destination Jacobian, modified in place. Entries carry power units
        divided by the corresponding global parameter unit.
    values : ndarray of shape (n_node, n_provider_pair)
        Derivatives for this provider's exclusively owned pairs.
    columns : ndarray of int, shape (n_provider_pair,)
        Destination pair indices.
    global_index : int
        Destination global parameter index.

    Returns
    -------
    None
        The selected output entries are overwritten.
    """
    output[:, columns, global_index] = values
