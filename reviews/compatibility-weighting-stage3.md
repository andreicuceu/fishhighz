# Compatibility weighting stage 3

Implemented 2026-09-17 following the stage-2 scientific review. An opt-in,
validation-only adaptive prescription now stops the three full-sample recurrences
in positive-signal forest-auto contexts. Every returned converged state passes
comparison with every later saved stage-2 state through update 96 on all three
grids at both tolerances. No production defaults or fixed-count paths changed.

## Prescription and use

`fishhighz/validation/adaptive_weights.py` adds the separate function:

```python
from fishhighz.validation.adaptive_weights import adaptive_weights

result = adaptive_weights(
    inputs, "sum_aliasing", context="lya(qso)_lya(qso)",
    rtol=1e-4, min_updates=3, stable_steps=3, max_updates=96,
)
if result["status"] == "converged":
    weights = result["weights"]
    aliasing, pixel_power = result["coefficients"][-2:]
```

Eligibility explicitly requires one of `lya(qso)_lya(qso)` or
`lya(lbg)_lya(lbg)`, finite positive intrinsic signal and `sum_intrinsic`,
`sum_aliasing` or `sum_historical`. Both prefixes and all 42 other saved contexts
return `ineligible`, without performing updates or silently using three updates.
The caller must supply the truthful context corresponding to its inputs.
Positive signal is necessary but does not guarantee convergence. Signed density
samples, weight amplitude and literal recurrence arithmetic are retained.

Defaults are minimum 3 updates, three successive stable transitions, rtol=1e-4,
and maximum 96 updates. Each transition tests relative amplitude, signed
normalized shape, A, P_pixel and complete vector changes without an absolute
floor. A candidate at t must remain stable at every intervening transition, and
every candidate-to-current comparison must remain within tolerance until 2t.
The state actually reached at 2t is returned. Simultaneous viable candidates are
retained so a failed earlier candidate does not erase a later qualifying one.
The cap is never extended for confirmation.

Statuses are `converged`, `capped`, `arithmetic_failure` and `ineligible`.
Only `converged` permits interpreting the state as satisfying this prescription.
Cap/failure returns the last state with finite weights and coefficients, its
`state_updates`, completed transition count `updates`, and reason. If an update
succeeds but its coefficient calculation fails, `updates` includes that completed
transition while `state_updates` identifies the retained earlier state. If even
the seed coefficients fail, no finite state is available. There is no substitute
or automatic fallback. `last_step` is explicitly the five changes from the
preceding state, not the forward fixed-point residual. `forward_residual=None`
marks the unavailable residual at the returned state: no additional transition
is hidden in the reported count. `confirmation` contains candidate-to-stop
changes, and `candidate` supplies its count.

The existing `fixed_weights(inputs, variant, updates=3)` remains available for
all five variants, with exactly three transitions after the common seed. The
existing covariance reassembly and public preparation APIs are unchanged.

## Saved-trajectory validation

`scripts/check_adaptive_weights.py` reconstructs inputs only from the saved
stage-2 arrays and original report scalars. No raw sampling/model evaluation is
repeated. It checks 36 auto/variant combinations on three grids at two tolerances:
216 eligible runs, with 192 converged and 24 capped. The caps are precisely
LBG bins 0–1 for intrinsic and moment-aliasing full sums, on every grid at both
tolerances. Historical sums converge in all 12 autos on every grid.

For each eligible result, weights and all five coefficients are bitwise equal
to the saved trajectory at the actual returned count. Confirmed counts match
the reviewed stage-2 classification. Every later state through 96 passes all
five relative comparisons against the returned stopped state. Expected capped
states remain equal to the saved state at 96. The script also verifies 1,404
ineligible variant/context/grid/tolerance combinations and bitwise agreement
of exactly-three-update weights for all 810 stage-2 trajectories.

| Variant | Converged per tolerance, across three grids | Updates at 1e-3 | Updates at 1e-4 | Maximum stopped-to-later discrepancy over five measures, 1e-3 / 1e-4 |
| --- | ---: | --- | --- | --- |
| sum_intrinsic | 30/36 | 10–26 | 14–40 | 1.05e-5 / 8.54e-8 |
| sum_aliasing | 30/36 | 10–26 | 12–34 | 3.98e-7 / 6.46e-9 |
| sum_historical | 36/36 | 8–24 | 10–30 | 1.06e-7 / 2.00e-9 |

Evidence is `.validation/compatibility-weighting-stage3/bin-{0..5}.json`,
including source-report and stage-2 JSON/NPZ SHA256 digests, status, actual count,
coefficients, last-step metrics, and maxima against every remaining saved state.
`summary.json` aggregates these six files. Stage-2 evidence is preserved.

## Checks

From `lib/fishhighz`, with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1` for Python commands:

```text
.venv/bin/python -m pytest tests/test_adaptive_weights.py tests/test_weight_convergence.py tests/test_compatibility_weights.py -q
.venv/bin/python scripts/check_adaptive_weights.py .validation/compatibility-weighting-stage2 .validation/compatibility-weighting-stage3 --bin N
.venv/bin/ruff check fishhighz/validation/adaptive_weights.py tests/test_adaptive_weights.py scripts/check_adaptive_weights.py
.venv/bin/ruff format --check fishhighz/validation/adaptive_weights.py tests/test_adaptive_weights.py scripts/check_adaptive_weights.py
git diff --check
```

30 tests passed in 0.18 s. Each saved-input bin check passed in less than one
second. Ruff, formatting and whitespace checks passed. Focused tests include an
analytic exact fixed point with actual stopping at six updates, homogeneous
subthreshold amplitude decay despite stable coefficients, a cap too short for
doubling, excluded prefixes/positive-signal mixed forest/nonpositive signals,
all five exactly-three paths, and signed arithmetic failure retaining the seed.
Known post-exit MUNGE warnings did not alter successful Python exit status.

## Scope

This is evidence on the tested finite grids through 96 updates, not a global
convergence theorem for arbitrary signed measures or a positive-density
optimality claim. Adaptive convergence does not resolve the historical sum's
213-to-425-node spacing sensitivity, 0.194% in A and 0.117% in P_pixel.
No BAO uncertainty criterion has been established; the 0.1% stopping-tolerance
criterion belongs to stage 4. No BAO forecast, full suite, new environment,
Slurm action, commit, push or production weighting adoption was performed.

Ready for the authorized independent stage-3 scientific review. Stop for review.
