"""W02: exact two-bin identities and independent float64 equation checks."""

import runpy
import sys
from fractions import Fraction as F
from pathlib import Path

import numpy as np


def main():
    """Check only dimensionless synthetic examples; no survey inputs."""
    w = (F(1, 2), F(1, 5))
    v = (F(1), F(4))
    i1 = sum(w)
    i2 = sum(x * x for x in w)
    i3 = sum(y * x * x for x, y in zip(w, v, strict=True))
    assert (i1, i2, i3) == (F(7, 10), F(29, 100), F(41, 100))
    prefix = (w[0] / (w[0] + v[0]), i1 / (i1 + v[1]))
    full = tuple(i1 / (i1 + y) for y in v)
    assert prefix == (F(1, 3), F(7, 47))
    assert full == (F(7, 17), F(7, 47))
    assert (i2 / i1**2, i3 / i1**2) == (F(29, 49), F(41, 49))
    assert (4 * i2 / (2 * i1) ** 2, 4 * i3 / (2 * i1) ** 2) == (F(29, 49), F(41, 49))
    doubled_prefix = (2 * w[0] / (2 * w[0] + 1), 2 * i1 / (2 * i1 + 4))
    doubled_full = tuple(2 * i1 / (2 * i1 + y) for y in v)
    assert doubled_prefix == (F(1, 2), F(7, 27))
    assert doubled_full == (F(7, 12), F(7, 27))

    # Float64 equations evaluated independently of the Fraction arithmetic.
    variance = np.array([1.0, 4.0])
    initial = 1 / (1 + variance)
    for scale, expected_p, expected_f in (
        (1, prefix, full),
        (2, doubled_prefix, doubled_full),
    ):
        weights = scale * initial
        density = np.cumsum(weights)
        np.testing.assert_allclose(
            1 / (1 + variance / density),
            [float(x) for x in expected_p],
            rtol=1e-12,
            atol=1e-14,
        )
        np.testing.assert_allclose(
            1 / (1 + variance / weights.sum()),
            [float(x) for x in expected_f],
            rtol=1e-12,
            atol=1e-14,
        )
        np.testing.assert_allclose(
            [weights @ weights, (weights * variance) @ weights] / weights.sum() ** 2,
            [29 / 49, 41 / 49],
            rtol=1e-12,
            atol=1e-14,
        )

    # Read the live array kernel without importing survey/package orchestration.
    kernel = runpy.run_path(
        str(Path(__file__).resolve().parents[1] / "fishhighz/kernels/weights.py")
    )
    actual, _ = kernel["_iterate"](np.ones(2), variance, 1, 1, 1, 1, 1)
    np.testing.assert_allclose(actual, [1 / 3, 7 / 47], rtol=1e-12, atol=1e-14)
    for scale in (1, 2):
        *_, a, p = kernel["_integrals"](np.ones(2), scale * initial, variance, 1, 1)
        np.testing.assert_allclose([a, p], [29 / 49, 41 / 49], rtol=1e-12, atol=1e-14)

    for r in (F(1, 2), F(1), F(2)):
        for fixed in (F(0), 1 - 1 / r):
            assert r * fixed / (r * fixed + 1) == fixed
            derivative = r / (r * fixed + 1) ** 2
            assert derivative == (r if fixed == 0 else 1 / r)
            x, mass = float(fixed), float(r)
            np.testing.assert_allclose(
                mass * x / (mass * x + 1), x, rtol=1e-12, atol=1e-14
            )
        print(f"r={r}: roots 0, {1 - 1 / r}; slopes {r}, {1 / r}")
    # At r=1, reciprocal weights increase exactly by one each update.
    assert 1 / (w[0] / (1 + w[0])) == 1 / w[0] + 1
    print(f"Python {sys.version.split()[0]}; NumPy {np.__version__}")
    print("PASS: rational identities, float64 equations, live tiny kernel checks")
    print("prefix=(1/3, 7/47); full=(7/17, 7/47); A=29/49; P_pixel=41/49")
    print("doubled: prefix=(1/2, 7/27); full=(7/12, 7/27)")


if __name__ == "__main__":
    main()
