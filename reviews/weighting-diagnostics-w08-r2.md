# W08 revision 2: reject underflow in positive instrumental power

**The W08-R1 range finding is corrected: the explicit inverse-variance method
now rejects inexact underflow in `l_p*v` before forming its denominator.**
The concrete review reproducer raises a contextual `ValueError`; exact-zero
variance still returns unit weights on positive support. All **111 focused tests
pass**, including the unchanged rational, legacy/supplied, noise and synthetic
survey checks. The saved order-32 QSO coefficients are identical to revision 1
and agree with independent scalar 80-digit Decimal arithmetic within
**3.18e-15** relative error. Ready for independent/user review.

## Correction and reproducer

In `fishhighz/weights.py`, only the inverse-variance product is evaluated inside
`np.errstate(under="raise")`. The resulting `instrumental_power` enters the
same addition and division as before. The existing exception handler supplies
field context and raises `ValueError`. No helper, public API, arithmetic
reordering, normalization or range extension was introduced. The multiplication
still acts only on positive-density support; exact zero products remain valid.
The legacy branch and `_iterate`/`_integrals` are unchanged.

Before editing, the review reproducer was confirmed with

```text
tiny = nextafter(0, 1)
l_p = B_star = tiny
q = (1,1), rho = (1/8,1/8), v = (0,1/2), L = 1.
```

Revision 1 accepted weights=(1,1), A=4 and P_pixel=tiny after losing the positive
product tiny/2. The new regression requires a field-context error containing
`not representable` and `underflow`. Its exact-zero control sets v=(0,0),
retaining weights=(1,1), A=4, I3=(0,0) and P_pixel=0. These are the only two
added tests. Ordinary zero-density behavior remains covered by the existing
QSO/LBG tests. The pre-edit result is retained in `r1-reproducer.json`.

README documentation remains applicable and was not changed. No production file
other than `weights.py` needed modification; there is no out-of-scope obstruction.

## Bounded validation

[Evidence directory](../.validation/forest-weight-diagnostics/w08-r2-20260916T005452Z/)
contains the logs, initial source copies, preservation manifest, copied saved
check script and its new results. Commands from the package root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 .venv/bin/python -B -m pytest -q tests/test_inverse_variance_weights.py tests/test_weights.py tests/test_weight_range.py tests/test_weight_models.py tests/test_noise.py tests/test_forecast.py::test_generated_spies
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 30 .venv/bin/python -B .validation/forest-weight-diagnostics/w08-r2-20260916T005452Z/check_saved.py
.venv/bin/ruff check --no-respect-gitignore fishhighz/weights.py tests/test_inverse_variance_weights.py .validation/forest-weight-diagnostics/w08-r2-20260916T005452Z/check_saved.py
.venv/bin/ruff format --check --no-respect-gitignore fishhighz/weights.py tests/test_inverse_variance_weights.py .validation/forest-weight-diagnostics/w08-r2-20260916T005452Z/check_saved.py
git diff --check
```

The single focused pytest invocation passed **111 tests in 2.35 s**. The single
saved-check invocation passed in **0.374 s** internally. Both used the existing
Python 3.13.15 / NumPy 2.5.3 environment on login32, one process and numerical
thread limits of one; neither cap was reached. Ruff lint/format and
`git diff --check` passed. No numerical invocation or lint check failed in this
revision; the pre-edit reproducer deliberately confirmed the known r1 defect.

The rational control still gives weights=(1/2,1/5), A=29/49,
P_pixel=41/49 and Q=10/7 from q=(1/2,2), rho=(2,1/2), v=(1,4) and
L=l_p=B_star=1. It agrees with independently supplied weights and legacy zero
iterations. Existing cumulative counts 0/1/2/3/6, metadata, immutability, context,
range and model-call regressions pass.

The scalar two-mode noise check retains intrinsic P1D=(4,9) at
q=(0.003,0.02) s/km, with the sinc/Gaussian response applied once to aliasing
and no smoothing of pixel noise. B_star only prepares the weights. In the
synthetic survey comparison, inverse-variance and supplied weights agree in
weights, coefficients, noise and the information matrix. Each preparation makes
one P3D and one P1D call, with no auxiliary query. Two central-difference runs
add five mean P3D calls each, giving 12 P3D calls in total while P1D remains at
two. Runtime traps and snapshots establish that weights, coefficients, response
and noise remain fixed during differentiation. These are rerun r1 tests, not
new historical Fisher comparisons or real forecasts.

## Saved-array coefficients and independent arithmetic

The unchanged r1 check script was copied into the new evidence directory and
executed against the revised package. It independently forms rho_v, nu and all
three moments using scalar 80-digit Decimal arithmetic from exact saved binary
inputs, without calling implementation weight/integral helpers for the oracle.
Only the assigned report was loaded:

```text
.validation/step12-r5-20260914T191855Z/profiles-checked/records-001.report.json
SHA-256 9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666
```

Its hash matched before and after execution. Identity checks select accuracy
bin 0, bounds [2.0,2.235], `settings.samples["lya(qso)"]`, physical model lya,
QSO background, and magnitude order 32: 5,184 nodes across 162 intervals.
The calculation uses saved density converted by (1+z_source)/299792.458,
saved quadrature and variance, z_source=2.3830099582078716,
L=44147.09373976799 km/s, l_p=63.328211161339375 km/s and
B_star=16.356748967852234 km/s.

The coefficient-only geometry is synthetic, with z_eval=2.1152848986890427,
constant H=200 and transverse distance=1000, h_fid=0.7, area=1 deg² and
z_order=2. The matching response has pixel width l_p and zero Gaussian sigma.
This is not a reconstruction of historical geometry or instrumental response:
A/P_pixel depend on the supplied arrays, L/l_p and already-smoothed B_star.
No raw reader, interpolation or P3D/P1D provider was called.

| Quantity | Revision 2 | Independent Decimal, rounded |
| --- | ---: | ---: |
| A [deg²] | 0.029419845663300446 | 0.029419845663300400 |
| P_pixel [deg² km/s] | 0.31435492218062416 | 0.31435492218062315 |
| Q [deg² km/s] | 0.7955679523681858 | 0.7955679523681840 |

The three returned coefficients equal the r1 saved result exactly. Maximum
relative discrepancies are 2.85e-16 for weights, 2.16e-15 for prefixes,
3.18e-15 for coefficients and 1.54e-15 against W06's printed reference. All
comparisons require finite values, rtol=5e-12, atol=0 and exact reference-zero
handling. The retained density/SNR policies were neither reapplied nor
physically validated.

## Source identity, preservation and limits

The supplied home path resolves to the requested CFS checkout. HEAD remains
`0d69786a06d5d676564a51fad14a7156951c2228`; source hashes identify this
uncommitted implementation. The initial 184-file manifest covers package
source, tests, reports and planning documents. Only `fishhighz/weights.py` and
`tests/test_inverse_variance_weights.py` changed among those files. The r1
handoff/review, README, planning documents, profile recipes, survey/forecast/noise
code, weighting kernels and `IMPLEMENTATION_STEP.md` remain unchanged.
Historical evidence was not overwritten.

```text
5e9110b63d13eded2cb84dbd0e82b958418ed59602d910bfeb1428a27d67f663  fishhighz/weights.py
a529a18855a66ead5682777bba668eae8525e4737dea8973b90b7bb5784a9ca0  tests/test_inverse_variance_weights.py
3867c9eee6d43f22fca878b16c9f38bafc65a55061355d79828fcb135cc5e66f  fishhighz/kernels/weights.py
56c1a312c44203e7e34cdcc7b542028130da6bd0e58bc59f8b552380a0e92b87  copied check_saved.py
3675cd861fc7909d464d9f7e60f7725f1a9e2f5cb3fa8a9f9849f3a30552c310  IMPLEMENTATION_STEP.md
```

The full hashes and exact Decimal values are in `saved-result.json`; preservation
checks are in `preservation.json`. Ruff explicitly checked the untracked tests
and ignored evidence script, which tracked `git diff --check` alone cannot cover.

This closes the concrete implementation finding, subject to independent review;
it establishes no wider extreme-input domain or new scientific optimum. W07's
BAO improvements remain historical evidence. W08 measures no new survey
improvement and changes no profile result. No real forecast, full suite,
wheel/environment work, performance study, physical input-policy change,
mode-dependent weighting, W01 repair, package Step 13 work, subsequent step,
Slurm action, agent dispatch, commit or push was performed. Stop for
independent/user review.
