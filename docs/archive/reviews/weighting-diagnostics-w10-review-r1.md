# Scientific finding and effect of possible corrections

**W10 revision 1 passes independent scientific review in its exact bounded
scope. No scientifically consequential correction is needed.** The fixed
inverse-variance reference is stable under the assigned 16/32/64 magnitude
refinements for both the bin-0 `lya(lbg)` and bin-5 `lya(qso)` samples. Their
saved order-32/64 individual-spectrum reference errors are likewise stable to
roundoff. Relative to each sample's own order-32 reference, the historical
cumulative t=3 and t=6 matrices give larger marginalized BAO errors in both
directions:

| Sample | t=3 radial / transverse | t=6 radial / transverse |
| --- | ---: | ---: |
| bin 0 `lya(lbg)` | +115.86120224% / +117.51202930% | +250.07503404% / +253.67973352% |
| bin 5 `lya(qso)` | +6.11227109% / +6.73272186% | +8.13892660% / +8.88846339% |

An error in sample or pair identity, native-order anchoring, the positive
measure, fixed B_star/pixel inputs, matrix selection, or the contrast denominator
could have changed this conclusion. A live rerun and a separate scalar/Decimal
reconstruction found no such error. There is no improvement reversal in these
two samples; had one been supported, it would have been a valid scientific
finding rather than an implementation defect. No revision of
`WEIGHTING_DIAGNOSTIC_STEP.md` is warranted.

This result extends the W07 ordering to one LBG-background forest at the same
redshift and one higher-redshift QSO-background forest. It is not a
population-by-redshift survey, a multi-mode optimum, physical validation of the
input policies, profile adoption, or acceptance of the Step 12 forecasts. A
review pass is not user acceptance or permission to prepare W11.

## Source, identity and held inputs

Reviewed the W10 handoff against the live assignment and inspected the submitted
script, public inverse-variance preparation and integral kernel. The diagnostic
binds the five assigned immutable inputs at
`scripts/check_inverse_variance_samples.py:22-45`, selects identity and native
order from report metadata at lines 111-146, reconstructs held physical inputs
at lines 147-201, and checks B_star, measure, support and native arrays at lines
202-257. It calls public `method="inverse_variance"` exactly once per order at
lines 258-270. The public path rejects conflicting iterative inputs and forms
`B_star/(B_star+l_p*v)` at `fishhighz/weights.py:268-299`; the coefficients use
the existing guarded integral products at `fishhighz/kernels/weights.py:23-37`.

All five assignment SHA-256 values match the live diagnosis, bin arrays and
reports. The handoff disclosed an output-only script edit after its original
run. I reran the final live source; the fresh artifact reports unchanged input
and source identities, and its complete `samples` payload is exactly equal to
the original artifact's payload. This closes the only run/source distinction
without changing any numerical result.

Metadata and independent selection establish:

| Sample | Native report order | Field / selected auto-pair index | Parameters | Nodes at 16/32/64 |
| --- | ---: | ---: | --- | --- |
| bin 0 `lya(lbg)` | 32 | 4 / 14 | `ap_0`, `at_0` | 2592/5184/10368 |
| bin 5 `lya(qso)` | 64 | 0 / 0 | `ap_5`, `at_5` | 2640/5280/10560 |

The native magnitudes, quadrature and variance match the corresponding report
exactly. Native masses reconstructed as
`density*(1+z_source)/299792.458*quadrature` also match exactly. Every saved
quadrature cell and mass is strictly positive, variance is finite and
nonnegative, all interval measures reproduce their partition widths within
1.80e-16 relative, and support is retained at the non-native orders.

The held report-derived scalars reproduce the rerun artifact exactly. In
particular, L=44147.09373976799 km/s for both samples; the bin-0 LBG sample has
`l_p=63.328211161339375 km/s` and `B_star=16.356748967852234 km/s`, while the
bin-5 QSO sample has `l_p=45.9777226171801 km/s` and
`B_star=64.49126195731874 km/s`. B_star is positive and exactly constant over
orders 16/32/64 within each sample. It is not response-smoothed again.

The inherited `floor_negative`, `legacy_floor`,
`legacy_spline_extension`, `legacy_floor_clamp` and
`legacy_first_spacing; unknown physical cells` policies are retained. Their
presence is an explicit limitation, not evidence of their physical validity.

## Independent coefficient and refinement arithmetic

The review calculation imports neither the submitted W10 script nor FishHighz
numerical helpers. It selects both samples from the frozen JSON/NPZ identities,
checks their native anchors and positive measures, then uses Python scalar
products with `math.fsum` for

```text
nu_i = B_star/(B_star + l_p*v_i)
I1 = sum(r_i*nu_i)
I2 = sum(r_i*nu_i^2)
I3 = sum(r_i*nu_i^2*v_i)
A = I2/(L*I1^2)
P_pixel = l_p*I3/(L*I1^2)
Q_star = A*B_star + P_pixel = B_star/(L*I1).
```

The order-32 independent coefficients are:

| Sample | A [deg^2] | P_pixel [deg^2 km/s] | Q_star [deg^2 km/s] |
| --- | ---: | ---: | ---: |
| bin 0 `lya(lbg)` | 0.1548959698171755 | 42.03495896736023 | 44.56855346179179 |
| bin 5 `lya(qso)` | 0.27830404960362665 | 2.5608574918779174 | 20.509036858648034 |

Across all six preparations, the maximum independent/public or
independent/historical coefficient residual is 8.16e-15; the Q_star identity
residual is at most 3.19e-16. The largest submitted signed 16/32/64 coefficient
contrast is 1.42e-14 in fractional units. Independent ratios differ from the
submitted ratios by at most 1.44e-14 because the public cumulative sums and
`math.fsum` round differently. Both are far below the 5e-12 numerical budget
and the 0.001 scientific refinement threshold. This is empirical bounded
magnitude stability, not a continuum statement.

## Independent Fisher selection and denominator

Only after the coefficient gate passed, the review selected each row's named
`pair_fisher_key` at the unique auto-pair index. These are historical
individual-spectrum matrices, not blocks of a joint inverse and not new mode
sums. All eight matrices are exactly symmetric and positive definite; their
condition numbers range from 3.22 to 3.39. No prior, regularization or
pseudoinverse is involved.

For each `F=[[a,b],[b,d]]`, an 80-digit Decimal determinant and
`(sqrt(d/det),sqrt(a/det))` agree with an ordinary NumPy solve and the saved
`pair_errors` within 2.20e-16 relative. The reference errors are:

| Sample / order | sigma_parallel | sigma_perp |
| --- | ---: | ---: |
| bin 0 `lya(lbg)` / 32 | 0.7288907443946463 | 0.9889518040663758 |
| bin 0 `lya(lbg)` / 64 | 0.7288907443946513 | 0.9889518040663827 |
| bin 5 `lya(qso)` / 32 | 0.12640403641329562 | 0.17544181456742458 |
| bin 5 `lya(qso)` / 64 | 0.12640403641329434 | 0.17544181456742267 |

The independent 64/32 fractional changes are
`(+6.8543e-15,+6.9603e-15)` for the LBG sample and
`(-1.0101e-14,-1.0916e-14)` for the QSO sample, so both reference-error gates
pass by many orders of magnitude. Every t=3/t=6 percentage uses the matching
sample's full-precision **order-32 reference error** as denominator; the
independent Decimal divisions reproduce the handoff percentages within
4.45e-16 fractional units. The large bin-0 errors are local Gaussian-curvature
constraints and do not establish that a Gaussian likelihood is adequate for
that weak observable.

The historical bin-0 order-64/t=6 joint diagnostic remains unavailable and was
neither needed nor regenerated. Its absence does not imply failure of the
available order-32 individual-spectrum matrix.

## Execution, evidence and disposition

The final live-source diagnostic rerun was:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 \
  .venv/bin/python -B scripts/check_inverse_variance_samples.py
```

It passed in 0.565 s internally and wrote
`.validation/forest-weight-diagnostics/w10-r1-20260916T030416201604Z/summary.json`
(SHA-256 `5d8a62b189b074ae2411b92b62c311a0beeab6a22000d6190bc42569b5d11304`).

The independent calculation is
check.py (local-only path: `../.validation/forest-weight-diagnostics/w10-review-r1-20260916T030501Z/check.py`)
with full-precision output in
result.json (local-only path: `../.validation/forest-weight-diagnostics/w10-review-r1-20260916T030501Z/result.json`).
Its final run passed in 0.648 s internally:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 \
  .venv/bin/python -B \
  .validation/forest-weight-diagnostics/w10-review-r1-20260916T030501Z/check.py
.venv/bin/ruff check scripts/check_inverse_variance_samples.py \
  .validation/forest-weight-diagnostics/w10-review-r1-20260916T030501Z/check.py
.venv/bin/ruff format --check scripts/check_inverse_variance_samples.py \
  .validation/forest-weight-diagnostics/w10-review-r1-20260916T030501Z/check.py
git diff --check
```

Both Ruff checks and `git diff --check` pass. The final independent checker and
result hashes are respectively
`5b71369a9d7f0c323dcfa86565340d6c8cddbce76c3e0e406e147402ce9b2dc1`
and
`05f81d0eaafaeafdc0f1ffae9d69261fed3443975fabdbca4e2692e392773f97`.

No W09 mode sum, W08 range campaign, broad test suite, new model/derivative/
covariance calculation, forecast, profile or production change, assignment
revision, Slurm action, dispatch, commit, push or W11 preparation occurred.
W10 revision 1 is ready for user review. Acceptance and any progression remain
with the user.
