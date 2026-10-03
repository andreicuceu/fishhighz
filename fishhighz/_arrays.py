"""Internal boundary helpers; scientific arrays use owned C-order float64."""

import numpy as np


def real_array(value, name):
    """Copy finite real numeric values into a C-contiguous float64 array.

    Parameters
    ----------
    value : array_like
        Real integer or floating-point data of any shape; units are preserved.
    name : str
        Quantity name used in validation errors.

    Returns
    -------
    array : ndarray
        Owned float64 copy with the input shape and units.

    Raises
    ------
    ValueError
        If values are Boolean, complex, nonnumeric, or nonfinite.
    """
    array = np.asarray(value)
    if array.dtype.kind not in "iuf":
        raise ValueError(f"{name} must contain real numeric values")
    array = np.array(array, dtype=np.float64, order="C", copy=True)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def immutable_float_array(value, name):
    """Return a finite float64 array, sharing memory that cannot change.

    Parameters
    ----------
    value : array_like
        Real numeric data of any shape; units are preserved.
    name : str
        Quantity name used in validation errors.

    Returns
    -------
    array : ndarray
        The input itself if it is a C-contiguous float64 array backed by
        immutable bytes (as built by ``freeze`` and ``geometry._immutable``,
        whose writeability cannot be re-enabled), after checking that it is
        finite; otherwise an owned copy as ``real_array`` returns. Sharing
        avoids copying the 1e7-element arrays of an integrated forest source.

    Raises
    ------
    ValueError
        If values are Boolean, complex, nonnumeric, or nonfinite.
    """
    base = value
    while isinstance(base, np.ndarray) and not base.flags.writeable:
        base = base.base
    if (
        isinstance(value, np.ndarray)
        and isinstance(base, bytes)
        and value.dtype == np.float64
        and value.flags.c_contiguous
    ):
        if not np.all(np.isfinite(value)):
            raise ValueError(f"{name} must be finite")
        return value
    return real_array(value, name)


def readonly(value, dtype):
    """Copy an array and mark the owned C-contiguous data read-only.

    Parameters
    ----------
    value : array_like
        Data to copy, with arbitrary shape and units.
    dtype : numpy.dtype or dtype-like
        Requested destination dtype.

    Returns
    -------
    array : ndarray
        Read-only copy with the input shape and units and requested dtype.
    """
    array = np.array(value, dtype=dtype, order="C", copy=True)
    array.flags.writeable = False
    return array


def integer(value, name, minimum=0):
    """Validate an integer count or index against its lower bound.

    Parameters
    ----------
    value : int or numpy.integer
        Dimensionless value to validate; Boolean values are rejected.
    name : str
        Quantity name used in validation errors.
    minimum : int, default=0
        Inclusive lower bound.

    Returns
    -------
    value : int
        Validated value converted to a Python integer.

    Raises
    ------
    ValueError
        If the value is not an integer or lies below minimum.
    """
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer")
    if value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return int(value)


def indices(value, size):
    # Check elements before NumPy can coerce mixed bool/integer inputs.
    """Validate a one-dimensional sequence of array indices.

    Parameters
    ----------
    value : array_like of int, shape (n_index,)
        Indices to validate without coercing Boolean entries into integers.
    size : int
        Exclusive upper bound for each index.

    Returns
    -------
    indices : ndarray of int64, shape (n_index,)
        Validated nonnegative indices in their original order.

    Raises
    ------
    ValueError
        If the input is not one-dimensional or contains invalid indices.
    """
    raw = np.asarray(value, dtype=object)
    if raw.ndim != 1:
        raise ValueError("indices must be one-dimensional")
    result = [integer(item, "index") for item in raw]
    if any(item >= size for item in result):
        raise ValueError("index out of range")
    return np.array(result, dtype=np.int64)


def label(value, name):
    """Validate a nonempty string label without normalizing its contents.

    Parameters
    ----------
    value : str
        Label to validate; surrounding whitespace is retained.
    name : str
        Quantity name used in validation errors.

    Returns
    -------
    value : str
        Original label.

    Raises
    ------
    ValueError
        If the value is not a string or contains only whitespace.
    """
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def scalar(value, name):
    """Validate a finite real scalar and return a Python float.

    Parameters
    ----------
    value : float or scalar array_like
        Real numeric scalar; its physical units are unchanged.
    name : str
        Quantity name used in validation errors.

    Returns
    -------
    value : float
        Validated scalar in the input units.

    Raises
    ------
    ValueError
        If the input is nonscalar, nonfinite, or not real numeric data.
    """
    array = real_array(value, name)
    if array.ndim != 0:
        raise ValueError(f"{name} must be scalar")
    return float(array)
