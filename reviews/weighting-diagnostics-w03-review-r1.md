# Scientific conclusion and impact of proposed changes

**W03 revision 1 passes independent scientific review within its assigned
finite-count scope. No scientifically consequential changes are needed.**
Replacing magnitude-prefix feedback by the common full-sample integral lowers
both coefficients at matched counts but increases their absolute fractional
change between three and six updates. The independently reproduced full-sample
changes are -0.9196% for A and +2.0171% for P_pixel. Both exceed the assigned
0.1% diagnostic threshold, so the prescribed stop at six updates is correct.
This does not establish asymptotic convergence, physical admissibility of the
signed measure, a magnitude-continuum limit, or a production prescription.
No repair revision is required.

Reviewed on 2026-09-15: [implementation handoff](weighting-diagnostics-w03-r1.md),
`scripts/compare_forest_weight_integrals.py`, and the
[exact archived W03 r1 assignment](weighting-diagnostics-w03-instructions-r1.md).
No substantive user feedback was supplied; the request retained its feedback
placeholder. **Awaiting user review and acceptance.** This pass does not permit
progression. Package Step 13 and its acceptance status are unchanged.

## Controlled comparison and source assessment

The inspected update and moment functions preserve W01's product order and
cumulative summation. The only branch-dependent operation selects either the
prefix array or its last element for N=(L/lp)I1. Separate weight arrays evolve
from the same seed (B/lp)/(B/lp+v); neither branch borrows the other's updated
state. Coefficient extraction uses full-sample moments in both branches:
A=I2/(L I1²), P_pixel=lp I3/(L I1²). There is no inter-update normalization.

I read the workspace/package guidance, design section 0, main handover/progress
register, diagnostic roadmap, exact assignment, W02 report/scientific review,
and W01 baseline report/review. The relevant W01 `_literal_terms` and
`float_snapshot` operations and legacy initialization/update/integral code agree
with this interpretation. W02's source/literature assessment is reused as
historical context; no literature search or broader source audit was needed.
In particular, this controlled substitution retains the existing fixed signal
without adding aliasing and is not a complete paper-prescription comparison.

The authoritative report is
`.validation/step12-r5-20260914T191855Z/profiles-checked/records-000.report.json`,
SHA-256 `3e5200d3e5b1bab3e7570c61fc02651caa907ae30816036b0b1b323098a87ed6`.
I verified compatibility bin 0, z=[2.0,2.235], and selected exactly
`settings.pair_inputs["lya(qso)_lya(qso)"]`. All six saved input/reference arrays
match the W03 NPZ exactly, are finite, and have 107 elements. The nine negative
density indices are 0, 2, 3, 76, 77, 80, 81, 84, 85; none were altered.

| Held quantity | Value | Units |
| --- | ---: | --- |
| Magnitude endpoints | 16.1, 26.75 | mag |
| Rectangular dm, including endpoints | 0.100471698113207 | mag |
| Forest length L | 44148.20436604321 | km/s |
| Pixel width lp | 63.32980433473595 | km/s |
| Fixed signal S | 1.7189455896517238 | deg² km/s |
| Fixed P1D B | 16.333147249645222 | km/s |

Density is deg⁻² (km/s)⁻¹ mag⁻¹, variance is dimensionless, and the integrated
moments have units deg⁻² (km/s)⁻¹. A is deg² and P_pixel is deg² km/s.
The report uses these units consistently; response already contained in S/B/v
is unchanged.

## Independent numerical checks

I wrote a small review-only scalar calculation, importing neither the
implementation script nor W01's controller. It converts the exact saved binary
floats with `Decimal.from_float`, uses 80-digit arithmetic, computes prefix sums
with direct scalar sums, and evolves each branch separately through t=6.
All four checkpoint weight vectors, moments, coefficients and signed normalized
shapes agree with the submitted arrays at rtol=5e-12, atol=0, with explicit
finite checks and exact-zero comparison. The maximum relative discrepancy over
these checks and baseline comparisons is **1.2721e-15**. No higher precision,
extra checkpoint, or singularity investigation was necessary.

Independent coefficients, rounded here to 15 significant digits:

| Feedback | t | A [deg²] | P_pixel [deg² km/s] |
| --- | ---: | ---: | ---: |
| Prefix | 3 | 0.0221183952067683 | 0.548026906249745 |
| Prefix | 6 | 0.0223089734691592 | 0.554483348675524 |
| Full sample | 3 | 0.0220213796793063 | 0.529854638024920 |
| Full sample | 6 | 0.0218188652536176 | 0.540542582564187 |

The prefix t=3 weights reproduce the source `_w_lya`; coefficients reproduce
its final `_aliasing_weights` and `_effective_noise_power`. Prefix t=3/6 also
reproduce the exact historical numbers specified in the assignment, independently
of W01's verdict labels. Their signed fractional changes are +0.8616278921%
and +1.1781250797%. The full-sample changes are -0.9196264205% and
+2.0171465478%. Thus A reverses its direction of change and both coefficients
have larger absolute finite-count sensitivity. The implementation's matched-count
level differences are consistent with these independently reconstructed values.

The exact rational two-bin control gives prefix (1/3,7/47) and full
(7/17,7/47), including the common last-bin first update. Source inspection verifies
that the implementation checks this same control against its float64 update.
The saved four-checkpoint inventory confirms that t=12 was not attempted.
All checked moments and coefficients are finite and positive. Signed densities
and the negative prefix weights at t=3 remain literal diagnostic inputs; these
results do not establish a physically admissible sampling measure.

Review evidence is in `.validation/forest-weight-diagnostics/w03-review-r1/`:
`check.py` contains the independent calculation and `result.json` its output.
The decisive invocation, from `lib/fishhighz`, was:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -B .validation/forest-weight-diagnostics/w03-review-r1/check.py > .validation/forest-weight-diagnostics/w03-review-r1/result.json
```

It passed with a 30-second alarm and completed in under one second. An initial
invocation failed before execution because the review script was written under
a duplicated relative directory; the script was moved to the path above and
only the newly created empty directories were removed. No scientific calculation
ran in that failed invocation. No broad suite or implementation-controller rerun
was needed.

## Disposition and preserved scope

No proposed scientific correction or additional check remains. The smallest
unanswered question is the full-sample branch's behavior beyond these permitted
counts; this review makes no inference about that behavior and assigns no work.
W01-R4/R5 remain deferred, and no package acceptance state is changed.

The exact assignment was archived byte for byte before updating its status
(SHA-256 `efc3218677af21705a9170e715c4df0e91101fc8a2cadd0fe350c72a6bf1d223`).
The diagnostic status in the current-step file, diagnostic roadmap/register,
design section 0, main handover paragraph and package AGENTS.md now records this
passing review and pending user acceptance. Scientific instructions are unchanged.
The package assignment retains SHA-256
`3675cd861fc7909d464d9f7e60f7725f1a9e2f5cb3fa8a9f9849f3a30552c310`.
The implementation's recorded source hashes matched before the status update.
The saved check now reads the byte-identical instruction archive so it remains
reproducible after that update; its final invocation also passed. Ruff lint and
format checks pass for the review script, and `git diff --check` passes.
Historical reports, implementation code, dirty/untracked work and numerical
artifacts were preserved.

No production edits, new survey/model evaluation, aliasing-signal substitution,
magnitude refinement, Fisher/full forecast, Slurm action, agent dispatch,
commit, push or subsequent-step preparation occurred. Stop for user review.
