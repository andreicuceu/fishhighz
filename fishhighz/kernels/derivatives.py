"""Array-only derivative combinations; dispatch and coefficients live outside."""


def _combine_three(fiducial, first, second, weights, scale):
    """Combine scaled slopes using actual offsets; cancel constant terms first."""
    return (weights[0] * (first - fiducial) + weights[1] * (second - fiducial)) / scale


def _scatter_column(output, values, columns, global_index):
    """Assign one provider's exclusively owned pair/global column."""
    output[:, columns, global_index] = values
