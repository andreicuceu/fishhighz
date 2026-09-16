# Step 13 revision 1 handoff

## Outcome

Step 13 revision 1 is implemented and ready for independent/user review. The
result is **no suitable finite grid-independent limit** for the existing
cumulative forest-weight recurrence, for both `lya(qso)` and `lya(lbg)` in all
six saved 15x2pt bins. Every fixed saved grid has a finite asymptotic shape, but
that shape concentrates onto shrinking quadrature measure and
`P_pixel` grows approximately linearly with magnitude order. The full
derivation and all 12 diagnoses are in
[step-13-weight-limit-r1.md](step-13-weight-limit-r1.md).

This completes the authorized bounded investigation only. It does not adopt a
new production prescription, reassemble a real forecast, accept the Step 12
accuracy profile, start Step 14, commit, or push.

## Implementation

- `fishhighz/validation/weight_limit.py`: validation-only log-amplitude/shape
  recurrence, scale-cancelled coefficient evaluation, concentration diagnostics,
  discrete spectrum and exact fixed-grid eigenvector limit, analytic controls,
  and hash-bound evidence validation.
- `scripts/diagnose_weight_limit.py`: explicit initialize/one-batch/finalize
  controller using saved arrays only; 60 immutable population/order batches,
  decision table, and three SVG figures.
- `tests/test_weight_limit.py`: 60 mathematical, range, corruption, and evidence
  tests. Coverage includes direct and Decimal recurrences, 0/1/3/6/12/24
  updates, independent rescalings, scalar (d<1,d=1,d>1), the continuum
  counterexample, repeated/near-degenerate spectra, conserved-measure
  subdivision, underflow cases, exact-zero branches, invalid inputs, and false
  convergence/evidence mutations.

No production weighting, noise, covariance, Fisher, model, response, reader,
or forecast file was changed by Step 13.

## Scientific evidence

Final evidence directory:
`.validation/step13-r1-20260914T235633Z/`.

Input identity:

- `weight-diagnosis/diagnosis.json` SHA256
  `069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008`;
- its source manifest SHA256
  `021a73d00644a5a365daea1d7579902888652fd1875d5345df378291969c678a`;
- `profiles-checked/manifest.json` SHA256
  `70edc87acf3c64229a144793c63528632a55b5e1763e7993d163083c71684474`;
- all six diagnosis NPZ hashes and all 12 profile NPZ/report hashes were
  rechecked before initialization and before every batch/finalization phase.

The source relationship was checked directly: all 12 order-64 forest inputs
match their independently saved accuracy-profile samples at relative tolerance
`5e-12`. The recovered pixel width is the unchanged observed-frame 0.8 Å
converted with `LYA_REST_ANGSTROM * (1 + z_eval)`; forest length, auxiliary
signal, alias, density, quadrature, and variance come from the hash-bound saved
reports. No reader, CAMB, P3D/P1D, or forecast call was made.

The controller completed 60/60 batches and 840/840 checkpoints through 1024
updates. Maximum numerical trajectory time was 0.812 s, below the 30 s batch
target. Historical 3/6/12/24 results comprise 144 available and 96 unavailable
coefficient attempts; all available pairs were reproduced to maximum relative
discrepancy `8.95e-14` in `A` and `1.22e-13` in `P_pixel`. The unchanged
production integrals reject 56 of 200 saved-weight states, whereas the exact
diagnostic retains their amplitude, shape, and coefficient provenance.

For the 12 fixed-grid asymptotic sequences, order 16/32/64 gives
`P_pixel` order exponents 0.9615--0.9798 and effective-measure exponents from
-0.9079 to -0.2671. `P_pixel` rises by more than 80% on both successive
refinements for every population, while the dominant physical magnitude moves
by at most 0.00234 and stabilizes. Thus all 12 entries in `decision-table.csv`
have `candidate_finite_limit=False`. The exact fixed-grid spectral radii span
`2.99338e-05`--`0.407248`, with unique dominant eigenvalues; this establishes
amplitude decay at fixed mesh, not noise convergence under refinement.

Artifacts include `summary.json`, the 60 hash-bound JSON/NPZ batch pairs,
`decision-table.csv`, `A-concentration.svg`, `P-pixel-concentration.svg`, and
`effective-measure.svg`. Final summary SHA256 is
`94f069593115e521e0787789195037708fccc9506371a5df954c3894288c6aee`.

## Validation

All work ran on the Perlmutter login node with one-thread settings and the NumPy
backend.

| Check | Result |
|:---|:---|
| Focused Step 13 tests | 60 passed |
| `scripts/check.sh` | 1303 passed, 25 skipped; Ruff passed; 161 files formatted |
| Installed affected tests | 121 passed |
| Installed examples | all six passed outside the checkout |
| Wheel/source/install identity | all 53 Python modules byte-identical; METADATA/WHEEL identical |

The ordinary-suite log contains non-fatal MUNGE socket warnings from login-node
Slurm discovery after pytest; the command exited zero and completed Ruff. No
Slurm job or allocation was requested.

The exact wheel is
`.validation/step13-r1-20260914T235633Z/dist/fishhighz-0.1.0.dev0-py3-none-any.whl`,
SHA256 `01956a930d288cb25b0f4b8d9eb3bedb246453a085ad301ecb11b91e440f3d93`.
It was built with Python 3.13.15, installed without dependency resolution into
an isolated target, and exercised under `python -I` from
`/tmp/fishhighz-step13-r1-tests`. NumPy was 2.5.3 and
`FISHHIGHZ_FISHER_BACKEND=numpy`. The build required approved network access for
the isolated declared setuptools dependency; no sibling environment changed.

Principal logs are `quick.log`, `build-approved.log`, `installed-tests.log`,
and `installed-identity.log`; `installed-identity.json` records every packaged
module digest and the runtime identity.

## Preserved status and review boundary

The Step 12 verdict remains unchanged: six compatibility bins pass, zero
accuracy bins pass, 60/72 historical per-spectrum diagnostics complete, and 12
are unavailable. Step 13's negative coefficient-limit result resolves the
bounded mathematical investigation after review; it does not make those
accuracy forecasts acceptable. A future real forecast requires both a
user-adopted forest-weight/noise prescription and explicit authorization.

The scientific choice now required is whether to formulate a different forest
estimator or to define and justify a finite-iteration/discretization convention.
This handoff makes neither choice and stops for review.
