# Step 12 revision 2 independent review

Date: 2026-09-14. Outcome: **changes required; scientific acceptance not achieved**.
The implementation handoff correctly reports its failed accuracy convergence.
The original five R1 invalid-payload probes now reject, but a further convergence
evidence defect remains (R2). All 39 compatibility bins reproduce lyaforecast;
only the five galaxy-only accuracy bins converge. No acceptance or Step 13
progression is implied. [Revision 3](../notes/IMPLEMENTATION_STEP.md) specifies the
same step's bounded repairs and the unresolved scientific decision.

## R2 — P2: convergence operands are not bound to the recorded trials

Locations: `fishhighz/validation/schema.py:156` (`metric_values`), `:309`
(metric validation), and `:381` (accuracy metric inventory);
`fishhighz/validation/accuracy.py:703` (trial/control persistence).

The checker recomputes changes between `metric_fisher`, `metric_pair_fisher`
and `metric_volume` operands. It does not bind these operands to `study_fisher`,
`study_pair_fisher`, `study_volume`, `study_controls`, `actual_levels`,
`final_controls` or `combined_lower_controls`. These trial records are not
required by the schema. A worker choosing the wrong operands can therefore
claim convergence while preserving contradictory trial evidence.

The independent probe copied accuracy bin 0 of `lya_lbg_lae_3x2pt`. Its original
weight refinement changes Fisher norm by 0.5231349801 and combined errors by
0.2006780025. Replacing the first operand of each metric with its second,
recomputing the resulting zero metrics, clearing the failure list and setting
`passed=True` makes both `execute()` and `check()` accept the copied bundle with
`complete=True`. The original nonconverged study matrices and their controls
remain unchanged. Hashes and inventories are updated by the normal writer.
The copied bundle deliberately contains one declared case/bin; this is a
per-record and general-checker defect, not a demonstrated bypass of the separate
78-record inventory check.

Evidence: probe result (local-only path: `../.validation/step12-review-r2/convergence-probe.json`),
reproducer (local-only path: `../.validation/step12-review-r2/independent.py`), and the exclusive
`false-convergence/` bundle beside them. Historical bundles were untouched.

Required repair: require a trial/control inventory and explicit metric-to-trial
references, bind all combined and selected-pair Fisher and volume operands,
and bind the final primary result to the final trial. Validate distinct actual
refinement controls and their outcomes, not merely the names of seven metrics.
Retain valid constant-result controls: distinct refinements may legitimately
produce identical information. Add consistently rehashed mutation tests for
duplicate/wrong trial selection, missing/reordered trials, inconsistent final
controls/results and omitted failed trials.

The independent review separately checked these bindings in **all 39 untouched
accuracy records**; their primary and metric matrices match the recorded trials.
R2 does not invalidate the numerical results in the handoff. It prevents relying
on the current checker as the promised convergence acceptance test.

## R3 — P2: the requested converged forest accuracy comparison is unresolved

Locations: `fishhighz/validation/accuracy.py:544` (bounded study),
`fishhighz/kernels/weights.py:10` (existing cumulative iteration), and the
handoff's convergence and diagnostic tables in [step-12.md](step-12.md).

All 34 forest-containing bins fail weighting stabilization at the existing
1e-3 Fisher/error tolerance; 23 also fail magnitude convergence. Final-two
changes reach 52.31% in Fisher norm, 20.07% in combined errors, and 38.50% in
individual-spectrum errors for weighting. k, mu, volume and derivative-step
checks pass in all 39 bins. The result is 44/78 passing primaries: 39
compatibility plus five galaxy-only accuracy records.

Of 468 diagnostic requests, 422 complete and 46 are unavailable due to guarded
underflow: 26 floor-low, eight rectangular-magnitude, ten bright-support-removal
and two width-plus variants. The full offline acceptance command correctly
exits 1 with `partial or scientifically failed validation`. These are documented
scientific limitations, not falsely reported successes or evidence of an error
in the final covariance/Fisher contraction.

The authorized bounded refinement study has already been attempted. Do not lift
arithmetic guards, relax tolerances, adopt a different weighting prescription,
change input floors/support, or launch additional real full runs to force a pass.
The next scientific choice belongs to the user: the current finite-iteration
prescription has not demonstrated the requested accuracy. Revision 3 permits
bounded diagnosis using stored samples and self-contained numerical examples;
any changed prescription or acceptance criterion requires a decision before
implementation. Preserve these curves as explicitly unconverged comparisons.

## The seven scientific comparisons

Each compares marginalized wiggle-only `ap` and `at` errors and their correlation
in independent original redshift bins, with fixed covariance/weights, no priors
or nuisance parameters, and full covariance among the selected spectra. Both
combined and individually selected-spectrum results are retained.

Write F_Q for Ly-alpha measured in QSO spectra and F_L for Ly-alpha measured in
LBG spectra. The first three cases have five bins over 2.26 <= z <= 3.41; the
others have six bins over 2.0 <= z <= 3.41. Entries below give signed ranges of
100*(sigma_accuracy/sigma_legacy - 1) over the original bins. Positive means
larger forecast uncertainty. Each case name links its original handoff plot.

| Case | Selected spectra | Delta sigma_ap (%) | Delta sigma_at (%) | Accuracy converged |
| --- | --- | ---: | ---: | ---: |
| lbg_lae_3x2pt (local-only path: `../.validation/step12-r2-implementation/plots-final/lbg_lae_3x2pt.png`) | LBG auto, LAE auto, LBG-LAE | -1.114 to -0.186 | -0.612 to +0.413 | 5/5 |
| lya_lbg_lae_3x2pt (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_lbg_lae_3x2pt.png`) | F_L auto, F_L-LBG, F_L-LAE | +7.908 to +15.752 | +18.490 to +23.945 | 0/5 |
| lya_lbg_lae_6x2pt (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_lbg_lae_6x2pt.png`) | All six spectra of F_L, LBG, LAE | +0.216 to +4.855 | +1.209 to +6.470 | 0/5 |
| lya_qso_2x2pt (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_2x2pt.png`) | F_Q auto and F_Q-QSO; QSO auto is covariance-only | -4.474 to -1.205 | -3.056 to +0.158 | 0/6 |
| lya_qso_lbg_lae_15x2pt (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_lbg_lae_15x2pt.png`) | All 15 spectra of F_Q, F_L, QSO, LBG, LAE | -3.307 to +2.365 | -0.610 to +3.318 | 0/6 |
| lya_qso_lbg_lae_4x2pt (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_lbg_lae_4x2pt.png`) | F_Q auto and its crosses with QSO, LBG, LAE | -7.248 to -3.901 | -3.503 to -0.769 | 0/6 |
| lya_qso_lbg_lae_8x2pt (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_lbg_lae_8x2pt.png`) | Both forest autos and each forest's three galaxy/QSO crosses; excludes F_Q-F_L | -6.897 to +6.285 | -2.420 to +12.563 | 0/6 |

The compatibility path uses captured upstream total powers and observed
Jacobians, then independently forms Wick C and FishHighz F. It reproduces the
legacy endpoint sums, polynomial peak extraction/backward derivative, pair
noise, redshift and reconstruction conventions. Its maximum relative Fisher
difference from legacy is 1.9565e-15; combined-error difference is 1.4433e-15.
This establishes estimator reproduction with matched inputs.

The accuracy profile instead uses the accepted KaiserModel and production
forecast preparation: signed PK-PKSB wiggles, full wiggle mapping derivatives,
CAMB growth, integrated volume, physical FWHM resolution and per-field noise.
Correlations become less negative in all bins, by 0.0141–0.0228. The full 15x2pt
errors differ by at most about 3.3%, but that agreement does not establish weight
convergence. For the forest-only LBG 3x2pt selection, the sequential attribution
identifies per-field noise as the largest effect (17.5–29.5% error increases at
that stage), partly offset by later model/derivative changes. Other cases often
have model/peak conventions as their largest sequential contribution; the
galaxy-only case is dominated by the full wiggle derivative stage. Attribution
is conditional on the finite controls and replacement order. The recorded
Jacobian interpolation residual, about 0.43–0.45%, is explicitly separate from
a physical attribution.

Useful supplementary figures: all cases (local-only path: `../.validation/step12-r2-implementation/plots-final/all-cases.png`),
LBG-forest attribution (local-only path: `../.validation/step12-r2-implementation/attribution-final/lya_lbg_lae_3x2pt.png`),
and 15x2pt sensitivity (local-only path: `../.validation/step12-r2-implementation/sensitivities-final/lya_qso_lbg_lae_15x2pt.png`).
The implementation handoff links all seven main PDFs and selected-spectrum PDFs.

## New verification and preservation

New evidence is exclusively in `../.validation/step12-review-r2/`:

- All **139** implementation-manifest files match the live checkout before
  review-document edits; all **48** package modules match the exact streamed
  wheel; all **75** input hashes match. The seven live INIs match the recipes.
  The separate reference audit verifies 108 metadata/configuration/result/array/
  source/resource files across all seven captured cases and 39 bins.
  Wheel SHA256: `be9015e0a0f945357018f92079a76d57581ea48d5a8f94862abec7ea63d421dc`.
- `PATH="$PWD/.venv/bin:$PATH" scripts/check.sh`: **1059 passed in 54.93 s**;
  Ruff lint/format pass, 134 files. This includes the original five R1 probes.
- `PYTHONPATH=. .venv/bin/python .validation/step12-review-r2/independent.py`:
  independently unpacked field matrices and looped over Wick pair indices in
  batches, solved C directly, and inverted the 2x2 Fisher analytically.
  All **78** primary combined and **468** individual-spectrum profile/bin
  matrices agree. Maximum relative discrepancy is zero for C, **3.161e-15**
  for combined F, **1.487e-15** for pair F, and **1.999e-15** for errors.
  All 39 accuracy primary/trial/metric bindings pass the independent check.
- Existing isolated `wheel-env-streamed`, from `/tmp` with `-I`: rerun
  `runtime/probe.py`, `step12_probe.py` and `revision2_probe.py` pass exact
  source/wheel/installed module and METADATA/WHEEL identity, five examples,
  numerical and optional-import probes, seven synthetic selections and all
  78 synthetic schema-2 profile records. A separate `-I` run of the copied
  `desi2_synthetic.py` also passes, completing all six examples. These are reruns of the existing
  installation, not a fresh build/install or a new real external-model run.
- `scripts/check_desi2_plots.py`: **720 plotted series, 78 profiles, 37 files**
  pass, including reference and selected-pair values and ratios. The first
  invocation lacked optional Matplotlib in `.venv`; rerunning with the existing
  read-only reference interpreter passed. Both logs remain. All seven main
  PNGs and the overview were visually inspected, plus representative attribution
  and sensitivity figures. Bins, units, overlapping compatibility curves and
  unresolved-accuracy labels are visible. Sensitivity gaps remain explicit.
- `scripts/attribute_desi2.py:verify_saved`: one saved eight-stage chain and two
  branches per case, **seven cases**, pass matrix/summary/difference checks
  without evaluating a model. The remaining attribution chains retain their
  historical implementation verification; they were not rerun in this review.
- `scripts/check_desi2_full.py` on the saved full bundle: expected exit 1 after
  exact 78/468 inventory checks; scientific completion remains false.
  `git diff --check` passes.

All numerical commands set OMP_NUM_THREADS, OPENBLAS_NUM_THREADS and
MKL_NUM_THREADS to 1. MUNGE/font-cache messages accompanied successful checks;
exit status and assertions determine the results. The seven fresh actual legacy
forecasts, 422 successful diagnostics and real external-P3D demonstration are
**historical implementation evidence**. This review did not run CAMB,
NewForecast, a real full forecast suite, or a new installation.

Production source, tests, implementation handoff, historical evidence and sibling
repositories are preserved. Only review artifacts and planning/governance
documents are added/updated. No agent was dispatched; no Slurm action, commit or
push occurred. The user's scientific choice and acceptance remain pending.
