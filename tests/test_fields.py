"""Field identities and independently named covariance products."""

import numpy as np
import pytest

from fishhighz.fields import ObservedField, PairSelection


def fields():
    return [
        ObservedField(x, "galaxy", "shared") for x in ("A_(q)", "B_b", "C", "D", "E")
    ]


def test_dependencies_and_opaque_identity():
    fs = fields()
    ab = PairSelection(fs, [("B_b", "A_(q)")])
    np.testing.assert_array_equal(ab.selected_pairs, [[0, 1]])
    np.testing.assert_array_equal(ab.required_pairs, [[0, 0], [0, 1], [1, 1]])
    autos = PairSelection(fs, [(0, 0), (1, 1)])
    np.testing.assert_array_equal(autos.required_pairs, ab.required_pairs)
    all_pairs = PairSelection(fs)
    assert all_pairs.selected_pairs.shape == (15, 2)
    assert all_pairs.available_pairs.tolist() == [
        [0, 0],
        [0, 1],
        [0, 2],
        [0, 3],
        [0, 4],
        [1, 1],
        [1, 2],
        [1, 3],
        [1, 4],
        [2, 2],
        [2, 3],
        [2, 4],
        [3, 3],
        [3, 4],
        [4, 4],
    ]
    forests = PairSelection(
        [
            ObservedField("lya(qso)", "forest", "lya", "qso"),
            ObservedField("lya_lbg", "forest", "lya", "lbg"),
        ]
    )
    assert len(forests.fields) == 2
    assert len(forests.available_pairs) == 3


def test_named_lookup_products_and_permutation():
    fs = fields()
    chosen = [("C", "B_b"), ("A_(q)", "A_(q)"), ("A_(q)", "B_b")]
    prepared = PairSelection(fs, chosen)
    names = [f.id for f in fs]
    required = [frozenset((names[i], names[j])) for i, j in prepared.required_pairs]
    # Explicit products for BC, AA, AB in that order, each table row-major.
    expected = [
        [
            ("B_b", "B_b"),
            ("B_b", "A_(q)"),
            ("B_b", "A_(q)"),
            ("A_(q)", "B_b"),
            ("A_(q)", "A_(q)"),
            ("A_(q)", "A_(q)"),
            ("A_(q)", "B_b"),
            ("A_(q)", "A_(q)"),
            ("A_(q)", "A_(q)"),
        ],
        [
            ("C", "C"),
            ("C", "A_(q)"),
            ("C", "B_b"),
            ("A_(q)", "C"),
            ("A_(q)", "A_(q)"),
            ("A_(q)", "B_b"),
            ("B_b", "C"),
            ("B_b", "A_(q)"),
            ("B_b", "B_b"),
        ],
        [
            ("B_b", "C"),
            ("B_b", "A_(q)"),
            ("B_b", "B_b"),
            ("A_(q)", "C"),
            ("A_(q)", "A_(q)"),
            ("A_(q)", "B_b"),
            ("A_(q)", "C"),
            ("A_(q)", "A_(q)"),
            ("A_(q)", "B_b"),
        ],
        [
            ("C", "B_b"),
            ("C", "A_(q)"),
            ("C", "A_(q)"),
            ("A_(q)", "B_b"),
            ("A_(q)", "A_(q)"),
            ("A_(q)", "A_(q)"),
            ("B_b", "B_b"),
            ("B_b", "A_(q)"),
            ("B_b", "A_(q)"),
        ],
    ]
    for table, products in zip(
        (prepared.im, prepared.jn, prepared.in_, prepared.jm), expected
    ):
        assert table.shape == (3, 3)
        assert [required[x] for x in table.ravel()] == [frozenset(p) for p in products]
    np.testing.assert_array_equal(
        prepared.required_pairs[prepared.selected_to_required], prepared.selected_pairs
    )


@pytest.mark.parametrize(
    "selected",
    [
        [],
        [(0, 1), (1, 0)],
        [(0, 5)],
        [(-1, 0)],
        [(True, 0)],
        [(0.0, 1)],
        [("missing", "B_b")],
        ["A_B"],
        [(0,)],
        [("B_b", 0)],
    ],
)
def test_invalid_pairs(selected):
    with pytest.raises(ValueError):
        PairSelection(fields(), selected)


@pytest.mark.parametrize(
    "args",
    [
        ("", "galaxy", "p"),
        ("a", "bad", "p"),
        ("a", "forest", "p"),
        ("a", "forest", "p", ""),
        ("a", "galaxy", ""),
    ],
)
def test_invalid_fields(args):
    with pytest.raises(ValueError):
        ObservedField(*args)


def test_duplicate_and_ownership():
    fs = fields()
    with pytest.raises(ValueError):
        PairSelection([fs[0], fs[0]])
    selected = np.array([[0, 1]])
    prepared = PairSelection(fs, selected)
    selected[:] = 3
    fs.clear()
    assert len(prepared.fields) == 5
    assert prepared.selected_pairs.tolist() == [[0, 1]]
    for name in (
        "available_pairs",
        "selected_pairs",
        "required_pairs",
        "selected_to_required",
        "im",
        "jn",
        "in_",
        "jm",
    ):
        array = getattr(prepared, name)
        assert (
            array.dtype == np.int64 and array.flags.owndata and array.flags.c_contiguous
        )
        with pytest.raises(ValueError):
            array.flat[0] = 9
