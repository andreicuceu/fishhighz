# Deferred FishHighz scientific tests

Updated: 2026-09-21.

This file records studies that may inform future scientific use but are not
required for the research-ready Python API milestone. It accompanies the
[implementation roadmap](../../FISHHIGHZ_IMPLEMENTATION_PLAN.md). Listing a
study does not dispatch a calculation or change an adopted forecast recipe.
Keep its status, scientific question, proposed comparison and eventual evidence
together when revisiting it.

The D identifiers are separate from implementation steps S1–S5 and all closed
historical assignments. Step 13, W12 and the old W-series remain closed.

| ID | Study | Status |
| --- | --- | --- |
| D1 | Weighting reference mode: fixed angular/velocity versus fixed comoving coordinates | Planned for later; explicitly deferred by the user on 2026-09-21. |
| D2 | Early weights versus the retained fixed-reference inverse-variance method | Optional proposal, lower priority; no execution decision. |

## D1 — Weighting reference-mode sensitivity

**Question.** How much do individual and joint BAO constraints depend on the
representative mode used to choose the early-lyaforecast forest weights?

The current reference is (k_t_deg, k_p_velocity) = (2.4, 0.00035), in inverse
degrees and s/km. Its comoving coordinates vary with redshift through
k_perp = k_t_deg / d_deg and k_parallel = k_p_velocity * a_v. A fixed comoving
reference selects different physical modes across redshift; expressing the same
physical mode in different units would not constitute this sensitivity test.

**Proposed comparison.** Match the existing reference at a declared pivot
redshift, then hold its comoving k and mu fixed in other bins. Decide the pivot
and exact mode when this study is assigned. Recompute each required forest's
fiducial P/B, converged weights, aliasing and pixel-noise coefficients, and the
resulting individual and joint BAO forecasts. Hold the cosmology, survey inputs,
physical model, response prescription, quadrature and forecast selection fixed.
The full BAO Fourier domain remains unchanged; the reference only selects weights.

An initial comparison in bins 1, 2 and 6 would cover the QSO-only selection,
low-redshift LBG forests and high-redshift endpoint. Bin 1 must contain only
QSO-forest auto, QSO auto and their cross; bins 2 and 6 retain all 15 spectra.
Extend the comparison only if its measured behavior warrants it.

**Evidence to retain.** Both BAO components for every selected individual and
joint result; weight convergence and fiducial reference coordinates; A and
P_pixel; rank and covariance diagnostics. Keep the existing 0.1% numerical
refinement target and 1% BAO-error attribution threshold. A change above 1%
would establish sensitivity to the weighting choice, not a convergence failure
or an automatic reason to replace the default. One alternative reference does
not establish global optimality or general reference-mode independence.

**Why deferred.** The adopted reference is explicit and the baseline already
passes its numerical tests. This study is not required for correct implementation
or documented use of that prescription. No new-baseline D1 forecast has run.

## D2 — Matched early/inverse-variance weighting comparison

**Question.** How strongly do forecasts depend on the early recurrence compared
with the retained fixed-reference inverse-variance prescription, at otherwise
matched accuracy-profile inputs?

Use the existing method implemented during W12 as a control within the current
accuracy model. Keep field responses, P1D, reference convention, survey inputs,
quadrature, physical signal and selected spectra matched. Prepare each method's
own weights and covariance; compare individual and joint BAO errors and the
noise coefficients. Selected informative bins can follow D1's proposed selection,
but this study can run independently of D1.

This is a comparison against the newly adopted baseline, not a reopening of W12
or an inference from its historical forecasts. It tests the weighting prescription
and answers a different question from D1. Exact cases and comparison settings
remain to be specified if the study is requested. The early baseline remains
the accepted weighting method; no new comparison or replacement is implied.

## Recording subsequent studies

Add a D identifier, status, scientific question, bounded proposed comparison,
reason for deferral and evidence links. Distinguish user-agreed future studies
from optional proposals. When a study runs, preserve the original proposal and
record its numerical findings separately from any decision to change a recipe.
