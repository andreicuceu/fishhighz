# Compatibility weighting stage 1

Implemented 2026-09-17 against stage 1 of the compatibility weighting plan.
This is a validation-only signed-input implementation; production defaults and
strict nonnegative preparation are unchanged. No scientific acceptance or
variant BAO comparison is claimed.

## Changes

- `fishhighz/validation/compatibility_weights.py`: five explicit recurrences,
  a common seed, simultaneous updates, full-sample final moments, and minimal
  saved pair-context preparation/noise/Wick-covariance reassembly. Full sums
  use the final element of the same cumulative reduction. Literal multiplication
  and division order matches the current legacy control. Negative density and
  negative auxiliary cross signals are retained; invalid arithmetic raises.
- `scripts/check_compatibility_weight_baseline.py`: bounded six-bin baseline
  replay from saved arrays, retaining input hashes and all 54 forest-related
  contexts. Each bin has nine distinct contexts; 12 auto preparations in total
  supply the additive forest noise in the covariance.
- `tests/test_compatibility_weights.py`: signed two-cell simultaneous updates,
  homogeneous fixed points and subthreshold decay, equal-noise populations,
  full-sample permutation invariance, exact update counts, invalid arithmetic,
  and an independent analytic P1D/response/unit-conversion noise calculation.

The mixed forest pair has its own preparation, but its total is noise-free:
`lyaforecast/covariance.py:327` routes only identical forest names to the
forest-auto total; distinct names reach `_compute_total_power_cross`. The
`_corr_type` preparation classification does not change this total-power routing.
All contexts remain accessible for the later iteration study.

## Evidence and checks

Saved evidence: `.validation/compatibility-weighting-stage1/baseline.json`.
Source: `.validation/step12-r5-20260914T191855Z/profiles-checked/records-000`,
`002`, `004`, `006`, `008`, `010` (reports and arrays).

All 54 saved three-update weight vectors and coefficient pairs reproduce exactly.
Across six bins the maximum relative discrepancies are:

| Quantity | Maximum |
| --- | ---: |
| total power | 5.92e-20 |
| selected covariance | 1.25e-19 |
| joint Fisher | 4.67e-17 |
| joint marginalized errors | 9.11e-17 |
| individual Fisher/errors | 0 |

Commands from the workspace root, with `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1` for Python checks:

```bash
lib/fishhighz/.venv/bin/python -m pytest lib/fishhighz/tests/test_compatibility_weights.py -q
lib/fishhighz/.venv/bin/python lib/fishhighz/scripts/check_compatibility_weight_baseline.py lib/fishhighz/.validation/step12-r5-20260914T191855Z/profiles-checked lib/fishhighz/.validation/compatibility-weighting-stage1/baseline.json
lib/fishhighz/.venv/bin/ruff check lib/fishhighz/fishhighz/validation/compatibility_weights.py lib/fishhighz/scripts/check_compatibility_weight_baseline.py lib/fishhighz/tests/test_compatibility_weights.py
git -C lib/fishhighz diff --check
```

Results: 18 tests passed in 0.13 s; six-bin replay passed in 3.71 s; Ruff and
whitespace checks passed. The environment emitted MUNGE socket warnings after
Python checks, without changing successful exit status.

## Scientific limits

Reassembly subtracts independently reconstructed baseline forest noise from the
saved total to preserve the intrinsic covariance-redshift spectrum, then adds
noise evaluated with the chosen weights. This is not merely relabelling cached
totals: changed coefficients alter the covariance with the retained P1D response
and angular/velocity conversion. The analytic noise test independently checks
response ownership and the unsmoothed additive pixel term.

The saved mean is evaluated at arithmetic mean redshift, while total powers use
geometric covariance redshift. Direct residual-vs-mean comparisons therefore
show relative differences 0.000340–0.002699 (including noise-free cross powers),
recorded explicitly in the evidence. An initial direct-equality assertion failed
and was removed after verifying this intentional convention in
`forecast_new.py:232–247`; no input or spectrum was adjusted. The residual is not
claimed to be an independently reconstructed intrinsic P3D model. Saved mean and
Jacobian remain untouched.

These are algebraic comparisons on signed compatibility inputs. Positive-density
optimality/convergence claims do not transfer to them. No convergence trajectories,
adaptive stopping, grid refinement, variant BAO forecasts, reference recapture,
Slurm work, full suite, commits, or production-weighting adoption were performed.
Existing modified governance files and untracked findings were preserved.

Ready for the authorized independent stage-1 scientific review.
