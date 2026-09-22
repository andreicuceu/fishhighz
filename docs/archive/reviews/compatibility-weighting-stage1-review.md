# Scientific review: compatibility weighting stage 1

## Scientific summary and verdict

**PASS.** No scientifically consequential change is needed before stage 2. The
validation path implements the five specified signed recurrences, preserves the
legacy compatibility seed and three-update convention, refreshes moment-based
aliasing from the current iterate, and evaluates each update simultaneously.
Final `A` and `P_pixel` use full-sample moments for every variant. The saved
compatibility control is reproduced in all 54 pair contexts, including negative
density nodes and negative auxiliary cross-signals.

The retained covariance intrinsic, defined as saved total power minus
independently reconstructed baseline forest noise, is sufficient for this
controlled comparison. The saved total is evaluated at the geometric covariance
redshift, so it must not be forced to equal the saved mean evaluated at the
arithmetic redshift. The recorded relative residual range
`0.0003400642--0.0026990375` is consistent with that fixed legacy convention and
does not indicate missing noise.

## Recurrence and coefficient checks

I inspected `fishhighz/validation/compatibility_weights.py`, the stage-1 tests,
the saved replay, the implementation handoff, the weighting findings, and the
relevant live `lyaforecast` source. The following properties agree with stage 1:

- `seed` evaluates `B/lp / (B/lp + v)`, matching the legacy operation order and
  the required `B/(B+lp*v)` seed.
- Each call to `update` computes `J1` and `J2` from one unchanged input vector.
  Prefix variants retain the full cumulative arrays; sum variants use their final
  elements, preserving the same reduction order. Moment aliasing is recomputed
  before each update. No normalization, clipping, density replacement, or
  regularization is introduced.
- `fixed_weights(..., updates=3)` performs exactly three transitions after the
  common seed. `coefficients` then uses the final full-sample `I1`, `I2`, and
  `I3` for all variants.

As an independent exact control, I used two signed cells with
`r=(-1,3)`, `w=(1/2,1/4)`, `v=(1,4)`, `L=2`, `lp=1`, `P=2`, and `B=3`.
The prefix moments are `J1=(-1/2,1/4)` and `J2=(-1/4,-1/16)`. Direct rational
evaluation gives updated weights

| Variant | Independent result |
| --- | --- |
| `prefix_intrinsic` | `(2, 1/5)` |
| `sum_intrinsic` | `(1/2, 1/5)` |
| `prefix_aliasing` | `(-1, 1/17)` |
| `sum_aliasing` | `(1/5, 1/17)` |
| `sum_historical` | `(4/5, 1/2)` |

All five agree exactly with `update`. The single-cell fixed points, homogeneous
subthreshold decay, equal-noise full-sample solution, permutation invariance,
and invalid zero-moment behavior also pass the focused tests.

An independent replay of the literal saved arrays reconstructed all 54
three-update `prefix_intrinsic` weight vectors and all 54 final `(A,P_pixel)`
pairs with zero relative discrepancy. All five variants also produced finite
three-update states for all 54 saved contexts in this bounded check. This is an
algebraic finite-count result on signed inputs; it supplies no convergence or
positive-density optimality claim.

## Noise, response, and covariance routing

The live `lyaforecast` total-power routing applies additive forest noise only
when both field names are identical forests. Distinct forests and forest-galaxy
pairs use the cross-power route without additive noise. Consistently, a
`sum_aliasing` reassembly changed only required-pair columns 0 and 14, the two
forest autos, while the mixed `lya(qso)`--`lya(lbg)` column remained bitwise
unchanged. Across six bins there are therefore 12 auto contexts entering
additive covariance noise; all 54 contexts remain available for the subsequent
iteration study.

I independently evaluated the source-level noise expression on both forest autos
in saved bin 0:

`[A P1D(z_cov,k_parallel) W_pixel^2 W_resolution^2 + P_pixel]
 * angle_to_distance^2 / distance_to_velocity`.

The result agrees with `forest_noise` to `4.09e-17` relative. This verifies that
the response multiplies P1D exactly once at the power level, the pixel-noise term
is not response-filtered, and the angular/velocity conversion yields
`(Mpc/h)^3`. The formula and constants agree with the inspected
`lyaforecast` P1D, spectrograph-response, and forest-auto total-power paths.

The six-bin baseline replay reports relative discrepancies no larger than
`5.92e-20` for total power, `1.25e-19` for selected covariance, `4.67e-17` for
joint Fisher information, and `9.11e-17` for joint marginalized errors. These
near-zero replay metrics confirm preservation of the saved control and Wick/Fisher
routing, but they are cancellation-based for the baseline total: the same
reconstructed noise is subtracted and added when the control coefficients match.
They are not independent validation of the noise formula. The analytic/source
comparison above is the independent evidence needed for that step.

## Checks run

With `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, and `MKL_NUM_THREADS=1`:

```text
.venv/bin/python -m pytest tests/test_compatibility_weights.py -q
18 passed in 0.12 s

.venv/bin/python scripts/check_compatibility_weight_baseline.py \
  .validation/step12-r5-20260914T191855Z/profiles-checked \
  /tmp/compatibility-weighting-stage1-review-baseline.json
passed=true, contexts=54, 3.35 s

.venv/bin/python - <<'PY'
# Independent rational recurrences, saved-input reconstruction,
# source-formula noise evaluation, and changed-column audit described above.
PY
5 exact recurrences; 54 saved contexts exact; noise relative=4.09e-17;
changed total columns=[0,14]; mixed forest cross unchanged
```

The Python processes emitted the known post-run MUNGE socket warning; all checks
completed successfully before that environment message. No forecast, full suite,
production change, or new scientific prescription was run or adopted.
