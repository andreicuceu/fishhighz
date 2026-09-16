"""Observed identities and covariance dependencies, without covariance values.

All stored mappings own read-only C-contiguous int64 arrays. Providers must
supply every required spectrum; missing spectra are never zero-filled or inferred.
"""

from dataclasses import dataclass, field

import numpy as np

from ._arrays import integer, label, readonly


@dataclass(frozen=True)
class ObservedField:
    """Identify one sample; shared physical labels do not imply parameter ties."""

    id: str
    kind: str
    physical_model: str
    background: str | None = None

    def __post_init__(self):
        label(self.id, "field ID")
        label(self.physical_model, "physical model")
        if self.kind not in ("galaxy", "forest"):
            raise ValueError("kind must be galaxy or forest")
        if self.kind == "forest" or self.background is not None:
            label(self.background, "background")


@dataclass(frozen=True, init=False, eq=False)
class PairSelection:
    """Prepare ordered fields and selected ID or integer pairs.

    ``im, jn, in_, jm`` index required_pairs for selected A=(i,j), B=(m,n).
    ``selected_to_required`` gathers means in the requested selection order.
    Omitting selected requests all canonical pairs.
    """

    fields: tuple[ObservedField, ...]
    available_pairs: np.ndarray
    selected_pairs: np.ndarray
    required_pairs: np.ndarray
    selected_to_required: np.ndarray
    im: np.ndarray = field(repr=False)
    jn: np.ndarray = field(repr=False)
    in_: np.ndarray = field(repr=False)
    jm: np.ndarray = field(repr=False)

    def __init__(self, fields, selected=None):
        fields = tuple(fields)
        if not fields or any(not isinstance(f, ObservedField) for f in fields):
            raise ValueError("fields must be nonempty ObservedField records")
        ids = {f.id: i for i, f in enumerate(fields)}
        if len(ids) != len(fields):
            raise ValueError("duplicate field ID")
        available = [(i, j) for i in range(len(fields)) for j in range(i, len(fields))]
        pairs = []
        for pair in available if selected is None else selected:
            if isinstance(pair, str) or len(pair) != 2:
                raise ValueError("a pair must contain two IDs or two indices")
            if all(isinstance(x, str) for x in pair):
                if any(x not in ids for x in pair):
                    raise ValueError("unknown field ID")
                pair = tuple(ids[x] for x in pair)
            else:
                pair = tuple(integer(x, "field index") for x in pair)
            if any(x >= len(fields) for x in pair):
                raise ValueError("field index out of range")
            canonical = tuple(sorted(pair))
            if canonical in pairs:
                raise ValueError("duplicate selected pair")
            pairs.append(canonical)
        if not pairs:
            raise ValueError("empty pair selection")
        products = [[], [], [], []]
        required = set(pairs)
        for i, j in pairs:
            for m, n in pairs:
                for table, product in zip(products, ((i, m), (j, n), (i, n), (j, m))):
                    canonical = tuple(sorted(product))
                    table.append(canonical)
                    required.add(canonical)
        required = [p for p in available if p in required]
        lookup = {p: i for i, p in enumerate(required)}
        object.__setattr__(self, "fields", fields)
        for name, data in (
            ("available_pairs", available),
            ("selected_pairs", pairs),
            ("required_pairs", required),
            ("selected_to_required", [lookup[p] for p in pairs]),
        ):
            object.__setattr__(self, name, readonly(data, np.int64))
        for name, table in zip(("im", "jn", "in_", "jm"), products):
            array = np.array([lookup[p] for p in table]).reshape(len(pairs), len(pairs))
            object.__setattr__(self, name, readonly(array, np.int64))
