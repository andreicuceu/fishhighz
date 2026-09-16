"""Explicit scope sharing and independent chain-rule checks."""

import numpy as np
import pytest

from fishhighz.parameters import (
    Parameter,
    ParameterBinding,
    ParameterRegistry,
    gather_local,
    map_jacobian,
)


def registry():
    return ParameterRegistry(
        [
            Parameter("cosmo_(f)", 0, "target"),
            Parameter("bin_1:b", -2, "nuisance", (-3, 0), 4),
            Parameter("bin_2:b", -1, "nuisance"),
        ]
    )


def test_scopes_order_and_chain_rule():
    reg = registry()
    one = ParameterBinding(
        reg, ["b", "f1", "f2"], {"b": "bin_1:b", "f1": "cosmo_(f)", "f2": "cosmo_(f)"}
    )
    two = ParameterBinding(reg, ["f", "b"], {"f": "cosmo_(f)", "b": "bin_2:b"})
    np.testing.assert_array_equal(
        gather_local(reg.fiducials, one.local_to_global), [-2, 0, 0]
    )
    np.testing.assert_array_equal(two.local_to_global, [0, 2])
    local_jac = np.array([[[7.0, 2.0, 3.0]]])
    mapped = map_jacobian(local_jac, one.local_to_global, 3)
    np.testing.assert_array_equal(mapped, [[[5, 7, 0]]])

    def analytic(global_theta):
        b, x, y = gather_local(global_theta, one.local_to_global)
        return 7 * b + 2 * x + 3 * y

    eps = 1e-5
    finite = [
        (analytic(reg.fiducials + eps * unit) - analytic(reg.fiducials - eps * unit))
        / (2 * eps)
        for unit in np.eye(3)
    ]
    np.testing.assert_allclose(mapped[0, 0], finite, atol=1e-9)
    # Per-field and per-pair scopes are opaque IDs, with no automatic sharing.
    scoped = ParameterRegistry(
        [Parameter(x, 1, "nuisance") for x in ("field_A", "pair(A,B)")]
    )
    binding = ParameterBinding(
        scoped, ["pair", "field"], {"pair": "pair(A,B)", "field": "field_A"}
    )
    np.testing.assert_array_equal(binding.local_to_global, [1, 0])


def test_empty_and_ownership():
    values = [Parameter("x", -1, "target")]
    reg = ParameterRegistry(values)
    values.clear()
    names, ties = ["a"], {"a": "x"}
    binding = ParameterBinding(reg, names, ties)
    names[0], ties["a"] = "b", "missing"
    assert binding.local_names == ("a",)
    assert binding.local_to_global.tolist() == [0]
    for array in (reg.fiducials, binding.local_to_global):
        assert array.flags.owndata and array.flags.c_contiguous
        with pytest.raises(ValueError):
            array[0] = 4
    empty = ParameterBinding(reg, [], {})
    assert gather_local(reg.fiducials, empty.local_to_global).shape == (0,)
    np.testing.assert_array_equal(
        map_jacobian(np.empty((2, 3, 0)), [], 1), np.zeros((2, 3, 1))
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(id=""),
        dict(fiducial=np.nan),
        dict(fiducial=np.inf),
        dict(fiducial=True),
        dict(role="fixed"),
        dict(bounds=(2, 1)),
        dict(bounds=(0, 0)),
        dict(bounds=(1, 2)),
        dict(bounds=(None, np.inf)),
        dict(step=0),
        dict(step=-1),
        dict(step=np.nan),
        dict(step=True),
    ],
)
def test_parameter_errors(kwargs):
    args = dict(id="x", fiducial=0, role="target")
    args.update(kwargs)
    with pytest.raises(ValueError):
        Parameter(**args)


@pytest.mark.parametrize(
    "names,ties",
    [
        (["x", "x"], {"x": "cosmo_(f)"}),
        (["x"], {}),
        ([], {"x": "cosmo_(f)"}),
        ([""], {"": "cosmo_(f)"}),
        (["x"], {"x": "absent"}),
    ],
)
def test_binding_errors(names, ties):
    with pytest.raises(ValueError):
        ParameterBinding(registry(), names, ties)


@pytest.mark.parametrize("idx", [[True], [0, False], [0.0], [-1], [3], [[0]], ["0"]])
def test_index_errors(idx):
    with pytest.raises(ValueError):
        gather_local([0, 1, 2], idx)
    with pytest.raises(ValueError):
        map_jacobian(np.zeros((2, 3, 1)), idx, 3)


def test_shape_registry_errors_and_open_bounds():
    with pytest.raises(ValueError):
        ParameterRegistry([])
    p = Parameter("x", 0, "target", (None, 1))
    with pytest.raises(ValueError):
        ParameterRegistry([p, p])
    for theta in ([[1, 2]], [np.nan], [True]):
        with pytest.raises(ValueError):
            gather_local(theta, [])
    for jac in (np.zeros((2, 3)), np.zeros((2, 3, 2))):
        with pytest.raises(ValueError):
            map_jacobian(jac, [0], 1)
    with pytest.raises(ValueError):
        map_jacobian(np.zeros((2, 3, 0)), [], 0)
