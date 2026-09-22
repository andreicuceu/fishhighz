# Step 12 revision 5: trial binding and 15x2pt weighting diagnosis

Revision-5 implementation is ready for independent/user review. **R2 is repaired;
R3 remains unresolved. Step 12 is not scientifically accepted.** All six 15x2pt
bins and all 15 individually selected spectra are retained. No other real case
was run, no NewForecast reference was recaptured, and no scientific prescription,
input floor, support, arithmetic guard or tolerance was changed.

Evidence is in step12-r5-20260914T191855Z (local-only path: `../.validation/step12-r5-20260914T191855Z/`).
The previous handoff is preserved byte-for-byte as [step-12-r2.md](step-12-r2.md).
The archived revision-3 assignment, revision-4 reports and evidence remain intact.
`preservation.json` verifies all 78 historical primary NPZ/report hashes, the
128-file initial dirty/untracked snapshot, and protected source files. All 34
production modules outside `validation/`, including the reviewed matrix and
compiled kernels, remain unchanged. The NumPy-only performance target disclosed
in revision 4 remains unmet; it was not waived or re-benchmarked here.

## R2: convergence operands now refer to actual trials

New bundles use schema 3. `validation/trials.py` binds the primary combined F,
15 pair F matrices and volume to the final trial. Every metric names its actual
lower/upper trial IDs, whose controls must match the isolated refinement or the
combined lower controls. The ordered successful arrays and explicit failed
outcomes share one inventory. Required families and bounded levels are declared
before checking. Missing, reordered, duplicate, unrelated or failed trials cannot
establish convergence. Distinct controls yielding identical matrices remain valid.

The bounded controller is separated into `validation/study.py`, preserving the
revision-4 three-payload cache. The offline checker replays that controller using
only the saved summaries and errors; it cannot evaluate a model, fill in missing
trials or discard an attempted failure. Separate direct Wick/NumPy oracles still
verify C/F independently of reused Cholesky factors. The writer and checker use
this contract even when a caller consistently updates file hashes and inventories.

Preparation-cache metadata could retain another finite-difference step label;
`AccuracyRecipe.evaluate` now records the actual requested controls. Numerical
Jacobian evaluation already used the requested step. Current primary matrices
agree with their historical counterparts at roundoff precision.

`historical-derived/` contains an exclusive derivation of the six original
accuracy trial bindings and six compatibility records, plus the original 72
15x2pt diagnostic outcomes. Historical failed controls are explicitly identified
as inferred by exact controller replay, matching every original error and all
saved successful trials. They are not presented as newly observed attempts.
Historical numerical failures remain failures; old bundles were not rewritten.
Schema-2 inspection cannot establish the new scoped acceptance contract.

All six copied historical accuracy records reject the original repeated-operand
mutation. In `historical-duplicate-probe.json`, both the writer and a consistently
rehashed offline bundle reject the copied failed 15x2pt bin with
`metric/trial operand: inconsistent numerical content`. The five original R1
probes remain covered. New tests include wrong/repeated trial selection,
missing/reordered trials, controls, primary/pair F, volume, combined lower
controls, suppressed failed attempts, finite/null/roundoff controls and identical
results at distinct valid settings.

## Completed numerical coverage

Both profiles finished in all six bins: **12/12 primary records, 180 individual-spectrum results**. Compatibility passes 6/6; accuracy converges in 0/6. All **72 diagnostic requests** are retained: **60 complete, 12 unavailable**. Unavailable variants: floor_low: bins 0, 1, 2, 3, 4, 5; remove_bright: bins 0, 1, 2, 3, 4; width_plus: bins 2. `coverage.json` records each error. The schema-3 scoped gate rejects the scientific pass, as required; execution completion is not acceptance.

The preserved actual reference and compatibility curves reproduce the literal
legacy estimator from upstream total powers and observed Jacobians. Accuracy
retains CAMB growth, full wiggle-mapping derivatives, physical FWHM resolution,
per-field noise and the fixed-domain magnitude partition. The signed PK-PKSB
wiggle, smooth PKSB, fiducial geometry, response/noise ownership and independent
per-bin ap/at targets are unchanged. The earlier template/background comparison
remains linked in the historical handoff; no template was retuned.

All primary final settings use 128 k intervals with order 4, mu and volume order
32, and ap/at steps 2.5e-4. Six cumulative updates are retained in every bin;
these are explicitly unconverged. Magnitude orders are 32 in bins 0–1 and 64 in
bins 2–5. The reports store every attempted level, actual outcome and metric
reference. k, mu, volume and derivative-step checks pass. Magnitude convergence
fails for selected-spectrum errors in every bin, even where joint errors satisfy
the budget. Weighting fails for both joint and selected-spectrum results.

| Bin | z interval | sigma ap, accuracy | sigma at, accuracy | Delta sigma ap (%) | Delta sigma at (%) | Delta correlation |
|---|---|---:|---:|---:|---:|---:|
| 0 | 2.000–2.235 | 0.0206376 | 0.0147899 | -3.3069 | -0.6100 | +0.01676 |
| 1 | 2.235–2.470 | 0.0199853 | 0.0146776 | -2.3615 | +0.3768 | +0.01730 |
| 2 | 2.470–2.705 | 0.0179195 | 0.0126499 | -0.4295 | +2.2466 | +0.01740 |
| 3 | 2.705–2.940 | 0.0149231 | 0.0098053 | -0.6017 | +1.4549 | +0.01657 |
| 4 | 2.940–3.175 | 0.0162115 | 0.0108724 | -0.8798 | +0.7005 | +0.01620 |
| 5 | 3.175–3.410 | 0.0236808 | 0.0179877 | +2.3648 | +3.3178 | +0.01814 |

Percent differences above are accuracy/reference; accuracy/compatibility agrees
at the shown precision because compatibility reproduces reference. Positive
values mean larger uncertainties. The complete comparison CSV (local-only path: `../.validation/step12-r5-20260914T191855Z/comparison.csv`)
contains all six joint and 90 selected-pair rows, all three profiles, both
FishHighz/reference comparisons and accuracy/compatibility differences.
Pair convergence (local-only path: `../.validation/step12-r5-20260914T191855Z/pair-convergence.csv`)
contains all 630 bin/pair/control comparisons. Full per-bin controls and raw
policy counts/contributions are in `audit.json` and the source reports. Bin
centres and widths are used for plotting; actual geometric evaluation redshifts
are separately recorded. Independent-bin combinations retain 12 distinct target
parameters and zero information between different bins.

## R3: fixed-weight quadrature versus the cumulative update

`diagnose_desi2_weights.py` uses saved model/background settings and the retained
raw readers without solving CAMB. It first isolates forest magnitude integration
using the original initial weight function w0(m), fixed as a function of m at all
orders 4, 8, 16, 32 and 64 on the same breakpoint partition. Galaxy noise, Fourier
nodes, modes and J stay fixed. It then evaluates the existing cumulative rule at
3, 6, 12 and 24 updates; this is a diagnostic, not a replacement primary survey.

The fixed-weight order-32/order-64 maximum selected-spectrum error change is
1.066e-14. Forty coupled magnitude/weight trial matrices independently match the
saved study results: maximum relative joint-F discrepancy 5.29e-15. All six
primary noise/F reconstructions agree at roundoff precision. Thus the observed
magnitude sensitivity is associated with the discretized cumulative update,
rather than unresolved integration of a fixed weight function in this control.
This does not establish quadrature convergence for arbitrary supplied weights.

Write r_i=rho_i q_i, alpha=L/pixel and a_i=variance_i/S. The update is

    N_i(t) = alpha sum_{j<=i} r_j w_j(t)
    w_i(t+1) = N_i(t)/(N_i(t)+a_i).

For the first cell a positive fixed point exists only if alpha*r_1>a_1;
then w_1*=1-a_1/(alpha*r_1). Otherwise its only nonnegative fixed point is zero.
The linearization at zero is lower triangular,
M_ij=alpha*r_j/a_i for j<=i. Every tested field/bin/order has spectral radius
max_i M_ii below one; the largest across the full diagnostic is 0.40725.
Since w(t+1)<=M w(t), these discrete unnormalized weights tend to zero.
The first-cell factors for floored samples can be of order 1e-24–1e-22.
Small synthetic tests verify the first-cell fixed point and subdivision effect.

This explains the very small intermediate products without blaming the
well-conditioned final two-parameter Fisher inversion. At bin 0, magnitude order
16 and 12 updates, F_L first fails at r*w*w in the first magnitude cell:
about 4.9881e-297 times 1.6280e-269, or 8.1204e-566. The diagnostic saves the
operands and an 80-significant-digit Decimal evaluation. Other 24-update trials
encounter divide overflow or zero denominators during the accepted update.
No failed coefficient is substituted into a forecast and no guard is relaxed.

Vanishing unnormalized weights alone do **not** exclude a finite limiting noise:
A and P_pixel are unchanged by a common rescaling of all weights. This revision
has not established such a normalized asymptotic limit, its quadrature
independence, or joint optimality. The actual allowed finite updates still change
joint errors by up to 5.017%, Fisher norm by 9.314%, and individual-spectrum errors
by 38.50%; they do not meet the 1e-3 information/error budget.

Consequential choices for the user, none implemented here:

- Retain an explicitly specified finite-iteration estimator and reconsider what
  numerical accuracy is required of that convention; this would change the
  present acceptance requirement, not demonstrate its convergence.
- Authorize investigation of a scale-normalized/asymptotic formulation of the
  same rule, establishing equivalence and convergence of A/P_pixel before any
  implementation or changed arithmetic policy.
- Specify and independently justify a different weighting estimator. That is a
  scientific prescription change, outside this repair.

Raw support/floor changes are not an automatic remedy. Their diagnostic effects
and failures remain explicit. No calibrated coverage, physical raw cell widths,
source-overlap noise, within-bin evolution or new P1D physics was inferred.

## Comparisons, sensitivities and figures

Each of the six bins has a saved eight-stage endpoint-reproducing supplied-array comparison, plus separate legacy cross-width and density-floor branches. Largest sequential uncertainty changes by bin:

| Bin | Largest stage | ap change (%) | at change (%) |
|---|---|---:|---:|
| 0 | template model with polynomial peak | -4.183 | -1.766 |
| 1 | template model with polynomial peak | -4.640 | -2.180 |
| 2 | template model with polynomial peak | -4.476 | -2.200 |
| 3 | template model with polynomial peak | -3.562 | -1.925 |
| 4 | template model with polynomial peak | -3.053 | -1.763 |
| 5 | per-field accuracy noise | +5.744 | +4.357 |

Maximum endpoint residual is 4.46e-15; the measured common-node Jacobian interpolation residual ranges from 0.0042709 to 0.0043365. This residual is retained rather than assigned to physics. The ordered swaps are not independent additive effects. Growth/redshift, resolution, reconstruction, finite weighting and input-policy one-setting diagnostics complement the chain. `attribution-summary.json` records every stage and its signed change. These comparisons explain finite endpoints within the chosen ordering, not a converged or physically complete forecast.

Retained-input variations include floor levels 1e-22/1e-20/1e-18, bright/sentinel support removal, and artificial +/-10% raw-cell widths at fixed normalized counts. `sensitivity-summary.json` retains signed effects and inactive-policy flags; failures are separately listed in `coverage.json`. These are local model/input diagnostics, not calibrated systematic uncertainties. No original first-spacing cell widths were promoted to measured physical edges.

Main PNG (local-only path: `../.validation/step12-r5-20260914T191855Z/plots/lya_qso_lbg_lae_15x2pt.png`)
and PDF (local-only path: `../.validation/step12-r5-20260914T191855Z/plots/lya_qso_lbg_lae_15x2pt.pdf`),
all 15 individual spectra (local-only path: `../.validation/step12-r5-20260914T191855Z/plots/lya_qso_lbg_lae_15x2pt-pairs.pdf`),
summary (local-only path: `../.validation/step12-r5-20260914T191855Z/plots/all-cases.png`),
controlled attribution (local-only path: `../.validation/step12-r5-20260914T191855Z/attribution/lya_qso_lbg_lae_15x2pt.png`),
sensitivity (local-only path: `../.validation/step12-r5-20260914T191855Z/sensitivities/lya_qso_lbg_lae_15x2pt.png`),
fixed/coupled magnitude convergence (local-only path: `../.validation/step12-r5-20260914T191855Z/weight-figures/magnitude-convergence.png`)
and weight decay (local-only path: `../.validation/step12-r5-20260914T191855Z/weight-figures/weight-decay.png`)
are generated from saved numerical tables. PDF counterparts and filename/hash
manifests accompany the figures. The inherited `all-cases` filename contains
only the six 15x2pt bins. An initial plotting pass emitted empty diagnostic placeholders named for excluded cases; they are preserved in `plot-scope-failure/`. The corrected command renders only the checked case inventory. No numerical curves from excluded cases were regenerated. The inherited grouped diagnostic panels are retained only in `plots/superseded-diagnostics/`; the linked sensitivity figure is the failure-preserving canonical result.

## Validation, identity, execution history and limitations

- `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PATH="$PWD/.venv/bin:$PATH" scripts/check.sh`:
  **1183 passed, 25 optional-compiler skips**, 115.28 s; Ruff lint/format passed.
  All 1188 baseline cases remain, with 20 new cases. Final Ruff and
  `git diff --check` also pass.
- Installed compiled evidence/performance/profile regressions: 373 passes in the
  initial selection, plus the added scoped-inventory test outside the checkout.
  The optional-dependency subprocess test passes under default NumPy. It fails
  if SciPy is deliberately blocked while `FISHHIGHZ_FISHER_BACKEND=numba` is
  forced; the same failure was reproduced with the untouched revision-4 wheel.
  This pre-existing compiled/blocker combination is disclosed, not repaired or
  counted as a compiled pass. Logs retain both failures and the NumPy control.
- All six examples pass outside the checkout with `-I`. All 52 source/wheel/
  installed Python modules and installed METADATA/WHEEL bytes match the exact
  wheel. Wheel SHA256:
  `acb7ae177d2d20bdc270c088c940bb80c6b42fa05c8c26a911e8a438f76e7e69`.
- Fresh isolated interpreter: Python 3.13.15, NumPy 2.3.5, SciPy 1.15.3, Numba
  0.66.0, llvmlite 0.48.0, Astropy 8.0.1 and CAMB 2.0.1. The read-only sibling
  lyaforecast source is explicitly imported and inventoried. No sibling
  installation or historical baseline environment was changed.
- NumPy/compiled comparison on the installed wheel has exactly zero Fisher
  difference. Recorded times for its small fixture are 0.00589 s NumPy,
  3.686 s cold compilation/contraction and 0.000190 s warm compiled contraction.
  These are a numerical/backend check, not a new performance benchmark claim.
  Actual backend and one-thread settings are recorded in `execution.json`.
- The first non-isolated build lacked setuptools; sandboxed build/download
  attempts failed DNS. Approved isolated dependency installation and wheel build
  succeeded. These failed logs are preserved. One exact wheel was built and
  installed in the new environment; scientific runs used that payload.
- An initial launch omitted explicit compiled-backend selection and was stopped
  during background setup, after one compatibility result and before an accuracy
  result. The corrected run completed all 12 primaries and 30 diagnostics, then
  exited 143/SIGTERM. `resume.py` verified and reused those exact-wheel inputs,
  reassembled primary C/F and retained all recorded diagnostic outcomes; only
  remaining diagnostic requests were evaluated. The partial bundles remain
  untouched. Logs, per-record times and resumed provenance distinguish actual
  model work from reassembly and I/O. No case was rerun merely to benchmark it.

The final bundle is `profiles-checked/`. `bind_desi2_sensitivities.py` derives it exclusively from `profiles-final/`, preserving all NPZ bytes and rebinding diagnostic comparison metrics to the final reassembled primary matrices. This fixes a 4.95e-15 discrepancy in a nominally zero effect without changing any tolerance. Both the original comparison reports and failed sensitivity-plot check log are preserved.

The independent final five-field Wick/direct-solve oracle checks all 12 records and 180 pair matrices. Maximum relative C discrepancy is 0, F discrepancy 1.84e-15, pair-F discrepancy 1.1e-15, and error discrepancy 1.52e-15. All six final accuracy trial bindings pass. The offline plot check verifies 240 series over six bins and 12 profiles; it checks values and reference ratios against the actual captured reference. `postprocess.json` stores exact commands, exit codes and timings for tables, audit, attribution, plotting, sensitivity, independent-bin checks and the expected failing scoped gate. Main figures and representative selected-pair pages were visually inspected; `visual-inspection.json` records the checked files.

The changed-file manifest (`changed-files.json`), exact commands (`postprocess.json`), import/resource/wheel identities (`installed-identity.json` and `final-identity.json`),
initial snapshot, logs and failure history are retained with the evidence.
The implementation stops for independent/user review. The cumulative weighting
choice and scientific acceptance remain with the user. No Step 13, Slurm action,
agent dispatch, commit or push occurred.
