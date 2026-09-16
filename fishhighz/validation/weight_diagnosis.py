"""Saved-sample diagnosis of the cumulative update; no alternate forecast policy."""

from decimal import Decimal, localcontext

import numpy as np

from ..kernels.weights import _integrals, _iterate


def first_underflow(masses, weights, variance):
    """Locate the first guarded multiply and report its exact decimal operands."""
    value = masses
    for name, factor in (
        ("r*w", weights),
        ("r*w*w", weights),
        ("r*w*w*variance", variance),
    ):
        try:
            with np.errstate(under="raise", over="raise", invalid="raise"):
                product = value * factor
        except FloatingPointError:
            for i, (a, b) in enumerate(zip(value, factor)):
                try:
                    with np.errstate(under="raise", over="raise", invalid="raise"):
                        np.multiply(a, b)
                except FloatingPointError as error:
                    with localcontext() as context:
                        context.prec = 80
                        exact = Decimal.from_float(float(a)) * Decimal.from_float(
                            float(b)
                        )
                    return dict(
                        operation=name,
                        index=i,
                        left=float(a),
                        right=float(b),
                        exact_product=str(exact),
                        error=str(error),
                    )
            raise AssertionError("array exception without scalar counterpart")
        value = product
    return None


def diagnose(masses, variance, length, pixel, signal, alias, counts=(3, 6, 12, 24)):
    """Replay the accepted finite updates and preserve unavailable coefficients."""
    masses, variance = np.asarray(masses), np.asarray(variance)
    rows = []
    arrays = {}
    for count in counts:
        try:
            with np.errstate(
                over="raise", invalid="raise", divide="raise", under="ignore"
            ):
                w, changes = _iterate(
                    masses, variance, length, pixel, signal, alias, count
                )
            arrays[f"weights_{count}"] = w
            arrays[f"changes_{count}"] = changes
            failure = first_underflow(masses, w, variance)
            if failure is not None:
                rows.append(
                    dict(
                        iterations=count, available=False, first_unrepresentable=failure
                    )
                )
                continue
            _, _, _, a, p = _integrals(masses, w, variance, length, pixel)
            rows.append(
                dict(
                    iterations=count,
                    available=True,
                    A=float(a),
                    P_pixel=float(p),
                    max_weight=float(w.max()),
                    min_weight=float(w.min()),
                )
            )
        except (FloatingPointError, ValueError) as error:
            rows.append(dict(iterations=count, available=False, error=str(error)))
    ratio = masses * length / pixel * signal / variance
    return arrays, dict(
        rows=rows,
        linearized_diagonal_max=float(ratio.max()),
        first_cell_linearized_ratio=float(ratio[0]),
        first_cell_positive_fixed_point=max(0.0, 1 - 1 / float(ratio[0])),
        definition="Jacobian diagonal at zero weights: rho*dm*L*S/(pixel*variance)",
    )
