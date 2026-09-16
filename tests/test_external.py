"""Provider boundary tests using plain functions and callable objects."""

import numpy as np
import pytest

from fishhighz.derivatives import evaluate_derivatives
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.models.external import (
    BoundParameters,
    P3DProvider,
    PreparedP3D,
    evaluate_p1d,
    evaluate_p3d,
)
from fishhighz.parameters import Parameter, ParameterRegistry


def setup_model(model=None, **kwargs):
    registry = ParameterRegistry([Parameter("x", 1.0, "target", step=0.01)])
    selection = PairSelection(
        [
            ObservedField("unused", "galaxy", "x"),
            ObservedField("A", "galaxy", "x"),
            ObservedField("B", "forest", "x", background="qso"),
        ],
        [("B", "B"), ("B", "A")],
    )
    bound = BoundParameters(registry, ["x"], {"x": "x"})
    if model is None:

        def model(t, z, k, mu, pairs):
            return (
                t[0]
                * (1 + k[:, None])
                * np.where(pairs[:, 0] == pairs[:, 1], 2.0, -1.0)
            )

    provider = P3DProvider(
        "all", model, bound, [("A", "A"), ("A", "B"), ("B", "B")], **kwargs
    )
    return PreparedP3D(registry, selection, [provider])


def args(prepared):
    return (
        prepared,
        prepared.registry.fiducials,
        2.0,
        np.array([0.2, 0.1]),
        np.array([1.0, 0.0]),
    )


def test_routes_signed_and_original_indices():
    prepared = setup_model()
    whole = evaluate_p3d(*args(prepared))
    model = prepared.routes[0].provider.model
    bound = prepared.routes[0].provider.parameters
    seen = []

    class Plain:
        def __call__(self, t, z, k, mu, pairs):
            seen.append(pairs.copy())
            return model(t, z, k, mu, pairs)

    split = PreparedP3D(
        prepared.registry,
        prepared.selection,
        [
            P3DProvider("b", Plain(), bound, [("B", "B")]),
            P3DProvider("a", Plain(), bound, [("B", "A"), ("A", "A")]),
        ],
    )
    np.testing.assert_array_equal(evaluate_p3d(*args(split)), whole)
    assert [x.tolist() for x in seen] == [[[2, 2]], [[1, 1], [1, 2]]]
    assert np.all(whole[:, 1] < 0)
    reordered = PreparedP3D(
        prepared.registry,
        PairSelection(prepared.selection.fields, [("A", "B"), ("B", "B")]),
        [r.provider for r in split.routes][::-1],
    )
    np.testing.assert_array_equal(evaluate_p3d(*args(reordered)), whole)


@pytest.mark.parametrize(
    "pairs,match",
    [
        (["AB"], "route pair"),
        ([("A", "A")], "missing"),
        ([("A", "B"), ("B", "A")], "duplicate"),
        ([("X", "A")], "unknown"),
        ([("unused", "unused")], "unused"),
    ],
)
def test_bad_ownership(pairs, match):
    p = setup_model()
    spec = p.routes[0].provider
    with pytest.raises(ValueError, match=match):
        PreparedP3D(
            p.registry,
            p.selection,
            [P3DProvider("bad", spec.model, spec.parameters, pairs)],
        )


def test_overlap_and_registry_identity():
    p = setup_model()
    spec = p.routes[0].provider
    with pytest.raises(ValueError, match="overlapping"):
        PreparedP3D(
            p.registry,
            p.selection,
            [spec, P3DProvider("other", spec.model, spec.parameters, [("A", "A")])],
        )
    with pytest.raises(ValueError, match="label"):
        PreparedP3D(p.registry, p.selection, [spec, spec])
    with pytest.raises(ValueError, match="different registry"):
        PreparedP3D(ParameterRegistry(p.registry.parameters), p.selection, [spec])


@pytest.mark.parametrize("analytic", [["missing"], ["x", "x"], ["x"]])
def test_bad_analytic_ids(analytic):
    with pytest.raises(ValueError, match="analytic"):
        setup_model(analytic_ids=analytic)


@pytest.mark.parametrize("value", [np.nan, 1j, True, object(), "1"])
def test_bad_output_types(value):
    p = setup_model(lambda t, z, k, mu, pairs: np.full((len(k), len(pairs)), value))
    with pytest.raises(ValueError, match="provider 'all'.*fiducial.*pairs") as failure:
        evaluate_derivatives(*args(p))
    assert isinstance(failure.value.__cause__, ValueError)


@pytest.mark.parametrize(
    "change",
    [
        {1: [1, 2]},
        {1: [True]},
        {1: [np.inf]},
        {2: -1},
        {2: 1j},
        {3: [0, 1]},
        {3: [-1, 1]},
        {3: []},
        {3: [[1, 2]]},
        {4: [0]},
        {4: [0, 1.1]},
        {4: [-0.1, 1]},
    ],
)
def test_input_validation_before_dispatch(change):
    calls = []
    p = setup_model(lambda *a: calls.append(a))
    values = list(args(p))
    for index, value in change.items():
        values[index] = value
    with pytest.raises(ValueError):
        evaluate_derivatives(*values)
    assert not calls


def test_outputs_shapes_and_domain_failure():
    p = setup_model(lambda *a: np.ones((2, 1)))
    with pytest.raises(ValueError, match="shape"):
        evaluate_p3d(*args(p))

    def domain(t, z, k, mu, pairs):
        if t[0] > 1:
            raise RuntimeError("external model domain")
        return np.ones((len(k), len(pairs)))

    p = setup_model(domain)
    with pytest.raises(
        ValueError, match="all.*parameter 'x'.*central.*pairs"
    ) as failure:
        evaluate_derivatives(*args(p))
    assert isinstance(failure.value.__cause__, RuntimeError)


@pytest.mark.parametrize("slot", [0, 2, 3, 4])
def test_readonly_inputs(slot):
    def model(*values):
        values[slot].flat[0] = 9
        return np.ones((2, 3))

    p = setup_model(model)
    with pytest.raises(ValueError, match="read-only"):
        evaluate_p3d(*args(p))
    assert p.registry.fiducials.tolist() == [1.0]
    assert p.routes[0].pairs.tolist() == [[1, 1], [1, 2], [2, 2]]


def test_hostile_flag_reset_cannot_change_prepared_or_caller():
    def model(t, z, k, mu, pairs):
        for value in (t, k, mu, pairs):
            value.flags.writeable = True
            value.flat[0] = 0
        return np.ones((len(k), len(pairs)))

    p = setup_model(model)
    values = args(p)
    snapshots = [a.copy() for a in values[1:] if isinstance(a, np.ndarray)]
    evaluate_p3d(*values)
    for a, b in zip([a for a in values[1:] if isinstance(a, np.ndarray)], snapshots):
        np.testing.assert_array_equal(a, b)
    assert p.routes[0].pairs[0, 0] == 1


def test_reused_buffer_copied_immediately_and_output_owned():
    buffer = np.empty((2, 3))

    def model(t, z, k, mu, pairs):
        buffer[:] = t[0] ** 2
        return buffer

    p = setup_model(model)
    result = evaluate_derivatives(*args(p))
    buffer[:] = -99
    np.testing.assert_allclose(result.jacobian, 2, atol=1e-13)
    np.testing.assert_array_equal(result.power, 1)
    assert result.power.flags.owndata and result.jacobian.flags.owndata


def test_p1d_independent_reordered_binding_and_units():
    reg = ParameterRegistry(
        [Parameter("unused", 99, "target"), Parameter("v", 4, "nuisance")]
    )
    bound = BoundParameters(reg, ["amplitude"], {"amplitude": "v"})
    k = np.array([0.01, 0, 0.001])
    calls = []

    def velocity(t, z, q):
        calls.append((t.copy(), q.copy()))
        return t[0] / (1 + 100 * q)

    output = evaluate_p1d(velocity, bound, reg.fiducials, 2, k)
    np.testing.assert_allclose(output, [2, 4, 4 / 1.1])
    np.testing.assert_array_equal(calls[0][1], k)
    p = setup_model()
    evaluate_derivatives(*args(p))
    assert len(calls) == 1
    for bad in (
        lambda *a: [1],
        lambda *a: [True] * 3,
        lambda *a: [1j] * 3,
        lambda *a: [np.nan] * 3,
    ):
        with pytest.raises(ValueError, match="P1D evaluation"):
            evaluate_p1d(bad, bound, reg.fiducials, 2, k)
    with pytest.raises(ValueError, match="nonnegative"):
        evaluate_p1d(velocity, bound, reg.fiducials, 2, [-1])
    with pytest.raises(ValueError, match="callable"):
        evaluate_p1d(None, bound, reg.fiducials, 2, k)
