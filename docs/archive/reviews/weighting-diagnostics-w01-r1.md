# Forest weighting diagnostic W01, revision 1 handoff

W01 is implemented and ready for independent/user review. The exact saved
three-update legacy result is reproduced. The first prescribed comparison,
six absolute updates from the same initialization, changes both forest-noise
coefficients beyond the 1e-3 diagnostic trigger. The bounded calculation
therefore stops at six updates and projects only the three- and six-update
states for the saved `lya(qso)` auto-spectrum in bin 0.

This result establishes finite-update sensitivity before any profile changes
for this one saved example. It does not invalidate a deliberately specified
three-update estimator, establish an optimal update count, describe other
populations or bins, or establish a continuum limit.

## Saved input identity and baseline gate

The diagnostic resolves the compatibility/bin-0 record by its manifest task,
not by record number alone. It verifies the schema-3 manifest, report and NPZ
hashes; upstream reference manifest, metadata and bin-array hashes; and the
recorded live lyaforecast weighting, covariance, P1D and Fisher source hashes.
The producer metadata and checked report contain identical selected pair inputs.

The retained sample has 107 endpoint-inclusive magnitude nodes on
16.1--26.75, rectangular width 0.100471698113207, nine negative density
interpolants and no zero density entries. No density conversion, floor, mask,
support removal or weight normalization is applied. The source fingerprint is
`047e2d535356e51e365d9c7d5fdd74cc3135c3dee61ed9faa0d21468ce564baf`.

At three updates the saved weights, A and P_pixel reproduce exactly in float64:

| t | A [deg2] | P_pixel [deg2 km/s] | sigma(ap) | sigma(at) |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 0.029449686440878566 | 0.3140935074916962 | unavailable | unavailable |
| 3 | 0.022118395206768268 | 0.5480269062497443 | 0.031927636457863354 | 0.030637148707328768 |
| 6 | 0.022308973469159235 | 0.5544833486755244 | 0.032107755826537035 | 0.030873464804936995 |
| 12 | not attempted after identification | not attempted | not attempted | not attempted |
| 24 | not attempted after identification | not attempted | not attempted | not attempted |

Independent scalar Decimal recurrences at 80 and 160 digits reproduce t=3 and
t=6. The largest float64/Decimal coefficient discrepancy is 6.67e-16 relative;
the recorded 80-to-160-digit coefficient refinements are zero at the retained
Decimal scale. No arithmetic failure or signed-recurrence singularity occurs
through the stopping point. Signed weights and prefix sums are retained in the
NPZ rather than interpreted as a positive-domain recurrence.

## Identified coefficient and forecast sensitivity

From t=3 to t=6, A changes by +0.008616278921205422 and P_pixel by
+0.011781250796540021 relative to their t=3 values. The absolute changes are
+0.00019057826239096687 deg2 and +0.00645644242578014 deg2 km/s. Either
coefficient independently exceeds the 1e-3 trigger.

The one-spectrum projection retains all 5000 saved Fourier nodes, modes,
observed Jacobian, physical total auto-power, geometry and response. Only the
coefficient-dependent noise is replaced. The two Fisher matrices are

```text
t=3: [[1259.1277304320095, 616.7085322519753],
      [ 616.7085322519753, 1367.4348344974080]]
t=6: [[1245.2640293363638, 608.8569927406515],
      [ 608.8569927406515, 1346.8232463264135]]
```

The Fisher Frobenius change is 0.013251541612760317. The marginalized ap and at
errors increase by 0.005641487709602044 and 0.007713384162009174,
respectively, so the single-spectrum forecast also crosses the existing 1e-3
diagnostic budget. Both Fisher matrices are positive definite and rank two.
The t=3 Fisher agrees with the saved `reference_pair_fisher` and `pair_fisher`
at 3.33e-15 and 3.55e-15 elementwise-relative discrepancy. Direct vectorized
summation and independent scalar accumulation/analytic 2x2 inversion agree
within the required tolerance.

## Implementation and focused checks

The standalone script separates array recurrence, Decimal control, scalar
Fisher projection, source verification, evidence writing and evidence replay.
It records every prescribed attempt, including the explicit t=12/24 early-stop
states. The focused tests cover one-cell d below/equal/above one, unequal
positive and signed small arrays, signed cancellation, scale-invariant
coefficient diagnosis, Fisher normalization/rank, field/pair permutations,
source and attempt mutations, exclusive output paths, and both discrepancy and
failed-baseline stopping.

- `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest -q tests/test_legacy_forest_iterations.py`: **23 passed in 0.31 s**; process wall time 0.80 s.
- `.venv/bin/ruff check scripts/diagnose_legacy_forest_iterations.py tests/test_legacy_forest_iterations.py`: **passed**.
- `.venv/bin/ruff format --check scripts/diagnose_legacy_forest_iterations.py tests/test_legacy_forest_iterations.py`: **2 files already formatted**.
- The explicit `--check-only` command reloaded the NPZ with `allow_pickle=False`, reconstructed its recurrence/Fisher content, and revalidated the authoritative source fingerprint: **passed**.

Canonical evidence is in
`.validation/forest-weight-diagnostics/w01-r1-20260915T203132Z/`:
`summary.json` SHA-256
`defdc7331e298c27c8eaafe2c5969d24083416cc1c78044cbabad38342336dd3`,
`arrays.npz` SHA-256
`4218feb3bcf10fe125db84549d9194c20f48e0365402763e01f6ee5c150cf6af`,
and the exact instruction snapshot SHA-256
`99fb517f5fabd67c48179bab437ab25e5d92b10a9a8abb139c1056a4a1fb0110`.
The final numerical portion took 0.108 s with Python 3.13.15, NumPy 2.5.3 and
all three thread limits set to one. A pre-final development bundle was moved
unchanged to `/tmp/fishhighz-w01-r1-pre-final-20260915T202847Z` after tightening
the failed-baseline early stop; it is not canonical evidence.

## Exact changed-file boundary and exclusions

W01 adds only:

- `scripts/diagnose_legacy_forest_iterations.py`;
- `tests/test_legacy_forest_iterations.py`;
- `reviews/weighting-diagnostics-w01-r1.md`;
- `.validation/forest-weight-diagnostics/w01-r1-20260915T203132Z/` containing
  `summary.json`, `arrays.npz`, `WEIGHTING_DIAGNOSTIC_STEP.md` and `checks.json`.

`checks.json` records hashes for the untouched production weighting/noise files,
`IMPLEMENTATION_STEP.md`, README, AGENTS.md, both weighting-diagnostic planning
documents and previous Step 12/13 handoffs. Their pre-existing modified or
untracked worktree states are preserved. No package source, old test/report,
roadmap or active step file was edited.

No NewForecast/controller, raw reader, interpolation, CAMB/background, external
P3D, other population/bin/spectrum, joint covariance, full suite, wheel,
installation, Slurm action, agent dispatch, commit or push was run. No plot was
needed because the finite comparison table identifies the effect directly.
Stop here for independent/user review; W02 and production policy remain
unauthorized.
