"""Internal boundary helpers; scientific arrays use owned C-order float64."""

import numpy as np


def real_array(value, name):
    array = np.asarray(value)
    if array.dtype.kind not in "iuf":
        raise ValueError(f"{name} must contain real numeric values")
    array = np.array(array, dtype=np.float64, order="C", copy=True)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def readonly(value, dtype):
    array = np.array(value, dtype=dtype, order="C", copy=True)
    array.flags.writeable = False
    return array


def integer(value, name, minimum=0):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer")
    if value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return int(value)


def indices(value, size):
    # Check elements before NumPy can coerce mixed bool/integer inputs.
    raw = np.asarray(value, dtype=object)
    if raw.ndim != 1:
        raise ValueError("indices must be one-dimensional")
    result = [integer(item, "index") for item in raw]
    if any(item >= size for item in result):
        raise ValueError("index out of range")
    return np.array(result, dtype=np.int64)


def label(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def scalar(value, name):
    array = real_array(value, name)
    if array.ndim != 0:
        raise ValueError(f"{name} must be scalar")
    return float(array)
