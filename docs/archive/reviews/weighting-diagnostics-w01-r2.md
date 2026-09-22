# Forest weighting diagnostic W01, revision 2 handoff

W01 revision 2 is implemented and ready for independent/user review. The
revision-1 scientific result is unchanged: for the saved `lya(qso)` auto-power
in bin 0, the literal signed recurrence reproduces the captured three-update
state exactly and first exceeds the coefficient trigger at six updates. Counts
12 and 24 were not evaluated. This repairs evidence integrity only; it neither
selects a weighting prescription nor authorizes W02 or a production change.

Canonical evidence is
`.validation/forest-weight-diagnostics/w01-r2-20260915T212017Z/`. Its final
`manifest.json` authenticates the numerical summary/arrays, exact revision-2
instructions, immutable diagnostic-script snapshot, final check inventory and
this handoff. The manifest is external to those artifacts, avoiding a circular
hash. The instruction SHA-256 is
`3ec253993f752fafc144f6df1f6cd11c55a0dde5ae703d0986147b48f351b481`;
the authenticated scientific-input fingerprint is
`047e2d535356e51e365d9c7d5fdd74cc3135c3dee61ed9faa0d21468ce564baf`.

## Scientific result

Held fixed were the 107-point rectangular magnitude grid, signed saved density,
variance, S/B auxiliaries, forest/pixel scales, all 5000 Fourier nodes, response
inputs, total-power baseline, mode counts and observed Jacobian. Only the
absolute recurrence count changed from three to six.

| updates | A [deg2] | P_pixel [deg2 km/s] | sigma(ap) | sigma(at) |
| ---: | ---: | ---: | ---: | ---: |
| 3 | 0.022118395206768268 | 0.5480269062497443 | 0.031927636457863354 | 0.030637148707328768 |
| 6 | 0.022308973469159235 | 0.5544833486755244 | 0.032107755826537035 | 0.030873464804936995 |

Relative to t=3, t=6 changes A by `0.008616278921205422` and
P_pixel by `0.011781250796540021`. The one-spectrum Fisher norm changes by
`0.013251541612760317`; the marginalized ap and at errors change by
`0.005641487709602044` and `0.007713384162009174`. Thus the terminal verdicts
remain `reproduced_legacy`, `coefficient_sensitive`, and `forecast_sensitive`.
This establishes finite-update sensitivity before profile changes for this one
saved example only.

The baseline independently matches `reference_pair_fisher` and `pair_fisher`
with maximum elementwise relative discrepancies `3.331e-15` and `3.553e-15`.
Both t=3 and t=6 Fisher matrices have rank two and both ap/at directions are
constrained. Scalar Fisher accumulation agrees exactly at recorded precision;
analytic 2x2 error discrepancies are at most `2.220e-16`. The 80- and
160-digit Decimal recurrences agree with float64 within the `5e-12` replay
tolerance.

## W01-R1 evidence binding

The final `--check-only` route now resolves the authoritative source with
`load_source` and requires the submitted source identity, full provenance and
scientific-input fingerprint to match. It reconstructs every evaluated
recurrence state, Decimal control, coefficient comparison, first terminal
event, verdict, sign inventory, table row, projected-count choice, P1D,
response, total power, node covariance, Fisher matrix, rank/constrained mask,
errors, scalar-oracle discrepancies and forecast comparison. Every named NPZ
array is checked against this reconstruction.

The finalized bundle contains a byte-identical copy of the executed script.
The external final manifest binds that script, `summary.json`, `arrays.npz`,
the instruction snapshot, `checks.json` and this handoff. Full-route regression
tests alter each reviewed claim or artifact class independently: terminal and
verdict, trigger operands, table coefficients/errors, forecast comparison,
projected errors/total/covariance/P1D/response, Decimal coefficients, script
identity, source-report path, source-report hash and manifest artifact hash.
All altered bundles reject.

## W01-R2 zero-safe baseline gate

A and P_pixel now use separate zero-safe scalar checks. A zero reference accepts
only exact zero; nonzero references retain the per-quantity `5e-12` relative
tolerance. Tests cover exact-zero matches and zero/nonzero mismatches in both
directions independently for A and P_pixel. Existing elementwise zero handling
for the weight array is retained.

## W01-R3 bounded failed-baseline evidence

Initialization failure, three-update float64 failure, Decimal/float64
disagreement and finite captured-baseline mismatch now return serializable
terminal results. Completed float states and both Decimal attempts are retained;
later counts are marked `not_attempted_after_failure`, projection arrays are
empty, and `reproduced_legacy` is absent. Each route is exercised through CLI
write, finalization and the same full `--check-only` path. A control in which
float64 reports a singularity while both Decimal evaluations remain finite is
classified as `arithmetic_failure`; only a Decimal-supported zero denominator
is classified as `singular_signed_recurrence`.

## Checks and execution

Commands were run from `lib/fishhighz` on the Perlmutter login node with the
package `.venv`; all numerical invocations used one OpenMP, OpenBLAS and MKL
thread.

- `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest -q tests/test_legacy_forest_iterations.py`: **48 passed in 1.37 s**.
- `.venv/bin/ruff check scripts/diagnose_legacy_forest_iterations.py tests/test_legacy_forest_iterations.py`: **passed**.
- `.venv/bin/ruff format --check scripts/diagnose_legacy_forest_iterations.py tests/test_legacy_forest_iterations.py`: **2 files already formatted**.
- The bounded real-source command recorded in `summary.json`: **passed in 0.1024 s**, with the t=6 early stop.
- Finalization and the authenticated real-source `--check-only` replay: recorded in `checks.json`; both must pass before this handoff is complete.

Python was 3.13.15 and NumPy was 2.5.3. No NewForecast/controller, external
P3D, reader interpolation, additional population/bin, joint 15x2pt covariance,
full package suite, installation, wheel, Slurm operation, commit, push or agent
dispatch was used.

## Changed and preserved files

Permitted implementation changes are exactly:

- `scripts/diagnose_legacy_forest_iterations.py` — SHA-256
  `239ee57d07ff485bc2ac77907ba104b5ed694cf7eb56e1b7b550c3f0847e4ed9`;
- `tests/test_legacy_forest_iterations.py` — SHA-256
  `430d53bacdeb2c8105e349c8f01cfb63ba33c45db5a4bd9f6e70ec4cb6a8fffc`;
- this new revision-2 handoff;
- the unique canonical evidence directory named above.

Production weighting files retain their pre-repair hashes:
`fishhighz/kernels/weights.py` `3867c9eee6d43f22fca878b16c9f38bafc65a55061355d79828fcb135cc5e66f`,
`fishhighz/weights.py` `fdfcef33eed0955729f8b4b0086d25c4360f5bda4ed999850dfd9e0f8658aabf`,
and `fishhighz/noise.py` `9cb580b1cb52860d3b65572bb18dfb707795b30a9646ff92778b930b429ef7b6`.
`IMPLEMENTATION_STEP.md`, README and the Step 12/13 handoffs also retain the
hashes recorded in `checks.json`.

Revision-1 artifacts remain byte-for-byte unchanged. In particular, the r1
handoff hash is
`cdd3401cb0caf21d9da2e46ee43480f06690ffba515288a9fb34aa7c2e5a6ce9`,
the r1 review hash is
`501b1dafc8f4e38917970846328033eb410f6f31d6d0e1d1c230ee2165de5149`,
and the four r1 evidence hashes are retained in the final check inventory.
The source manifest, report, NPZ and producer/reference inputs are rehashed by
every authoritative replay.

## Scope limit and stopping point

W01 revision 2 repairs all three review findings while preserving the supported
finite-update result. It does not establish prevalence across bins or
populations, an optimal finite count, or an iteration/continuum limit. The
smallest next action is independent review of this repair. No further diagnostic
or scientific comparison has been started.
