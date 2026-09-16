# Scientific conclusion and impact of proposed changes

**W08 revision 2 passes independent scientific review. No scientifically
consequential correction is needed.** The repair closes W08-R1 by rejecting
inexact underflow of a positive `l_p*v` product before it can alter the fixed
inverse-variance weight. Exact zero remains valid, and an exactly representable
subnormal product remains accepted, so the repair does not over-narrow the
declared numerical domain. Ordinary rational, supplied/legacy, noise, survey
and saved-reference results remain unchanged.

No user feedback beyond the request to review revision 2 was supplied. No repair
revision is proposed. This pass is ready for user review and acceptance; it is
not profile adoption, forecast acceptance or permission to advance.

## Exact repair and preserved scope

The exact revision-2 assignment is archived at
`reviews/weighting-diagnostics-w08-instructions-r2.md` with SHA-256
`3a3e8b13...428f`. The 184-file pre-edit manifest independently reconstructs
exactly two changed existing files: `fishhighz/weights.py` and
`tests/test_inverse_variance_weights.py`; none are missing.

The production diff is confined to `weights.py:291–296`. On positive-density
support it now computes `instrumental_power=pixel*variance` inside a local
`np.errstate(under="raise")`, then uses the same
`alias/(alias+instrumental_power)` expression. The surrounding handler converts
`FloatingPointError` to the existing field-contextual `ValueError`. Overflow,
invalid and divide policies are inherited unchanged from the outer context.
The legacy and supplied branches, `_iterate`, `_integrals`, noise, survey,
forecast and README are byte-identical to revision 1.

The only two new tests are the exact W08-R1 reproducer and an exact-zero control.
No public API, estimator, arithmetic normalization, response convention,
provider route, profile recipe or physical input policy changed. Live decisive
hashes match the handoff:

```text
5e9110b63d13eded2cb84dbd0e82b958418ed59602d910bfeb1428a27d67f663  fishhighz/weights.py
a529a18855a66ead5682777bba668eae8525e4737dea8973b90b7bb5784a9ca0  tests/test_inverse_variance_weights.py
3867c9eee6d43f22fca878b16c9f38bafc65a55061355d79828fcb135cc5e66f  fishhighz/kernels/weights.py
3675cd861fc7909d464d9f7e60f7725f1a9e2f5cb3fa8a9f9849f3a30552c310  IMPLEMENTATION_STEP.md
```

## Independent range controls

The independent review check does not import the implementation test or saved
check script. For

```text
tiny = nextafter(0,1)
l_p = B_star = tiny
q = (1,1), rho = (1/8,1/8), v = (0,1/2), L = 1,
```

the positive product `tiny/2` rounds to zero. Revision 2 rejects it with

```text
f: weighting integrals/coefficients are not representable or lack support:
underflow encountered in multiply
```

and retains the originating `FloatingPointError`. This is the requested
field-contextual closure of W08-R1.

Two independent boundary controls pass:

- With v=(0,0), the products are exact zeros and the method returns weights
  (1,1), A=4, I3=0 and P_pixel=0.
- With v=(0,1), the second product equals the smallest positive float exactly.
  It remains accepted with weights (1,1/2), A=4.444444444444445 and
  P_pixel=`tiny`.

Thus the branch rejects the inexact positive underflow without rejecting exact
zeros or exact subnormal arithmetic. Zero-density nodes remain outside the
product calculation and retain canonical zero weights.

## Focused regressions and coefficients

The assigned focused invocation passed **111 tests in 1.51 s**:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 \
  .venv/bin/python -B -m pytest -q \
  tests/test_inverse_variance_weights.py tests/test_weights.py \
  tests/test_weight_range.py tests/test_weight_models.py tests/test_noise.py \
  tests/test_forecast.py::test_generated_spies
```

This reruns the explicit keyword conflicts, immutable metadata, QSO/LBG support,
ordinary and extreme range checks, unchanged legacy counts, two-mode scalar
noise response and tiny survey equivalence/provider-call checks. The unchanged
survey tests verify one P3D and one P1D call per preparation, no auxiliary query,
five mean P3D calls per two-parameter run and no P1D/noise/weight reevaluation
during differentiation. The unchanged source hashes independently support the
fixed-noise boundary.

Exact Fraction arithmetic still gives weights=(1/2,1/5), A=29/49,
P_pixel=41/49 and Q=10/7. The public result differs by at most `1.88e-16`;
independently supplied weights and legacy zero iterations agree exactly.

The saved check verifies the assigned report hash
`9b64dd91...666`, accuracy bin 0, `lya(qso)`, order 32 and 5,184 nodes. An
independent scalar 80-digit Decimal calculation from the exact saved binary
density, quadrature and variance gives:

| Quantity | Revision 2 | Independent Decimal |
| --- | ---: | ---: |
| A [deg^2] | 0.029419845663300446 | 0.029419845663300400 |
| P_pixel [deg^2 km/s] | 0.31435492218062416 | 0.31435492218062315 |
| Q [deg^2 km/s] | 0.7955679523681858 | 0.79556795236818397 |

The maximum coefficient discrepancy is `3.18e-15`; agreement with W06's printed
reference is `1.54e-15`, below rtol=5e-12 with zero absolute tolerance. This is
only a coefficient check in synthetic geometry. It does not replay a forecast
or physically validate inherited density/SNR policies.

## Evidence, execution and disposition

Independent evidence is in
`.validation/forest-weight-diagnostics/w08-review-r2-independent/check.py` and
`result.json`. The successful calculation completed in **0.273 s** under the
30-second cap, one process and numerical thread limits of one:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 30 \
  .venv/bin/python -B \
  .validation/forest-weight-diagnostics/w08-review-r2-independent/check.py
```

The first review-script invocation used `Fraction` objects directly in public
array inputs and failed validation before rational or saved calculations. The
reviewer corrected the fixture to equivalent float inputs; this was a review
fixture error, not an implementation or scientific failure. Final Ruff
lint/format and `git diff --check` pass.

W08-R1 is closed. W07's BAO comparison remains historical evidence; W08 has no
new survey result and makes no multi-mode optimality claim. No production repair
was made during review. No real forecast, broad suite, wheel/environment,
profile/default change, physical input-policy change, package Step 13 work,
subsequent step, Slurm action, agent dispatch, commit or push was performed.
Stop for user review and acceptance.
