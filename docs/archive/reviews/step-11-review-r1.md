# Step 11 revision 1: independent review

Date: 2026-09-13. Outcome: **changes required**; revision 2 instructions prepared.
Reviewed the implementation and `reviews/step-11.md` against the preserved
revision-1 assignment. The user requests nonuniform redshift bins and confirmed
explicit raw redshift-cell widths plus a labeled legacy first-spacing conversion.
No production code, tests, examples, or implementation report were changed by
this review. No acceptance, progression, dispatch, commit or push is implied.

## Required changes

### R1: nonuniform redshift support and contradictory reader requirement

This is a user-requested design correction and an acceptance blocker.
`fishhighz/forecast.py:288` already accepts independently specified unequal bins;
there is no uniform-bin constraint there. Independent tests of [2,2.7] and
[2.7,2.9] recover analytic integrated volumes and the joint Fisher matrix across
four batch sizes (maximum relative Fisher error 2.032e-14).

The restriction is `fishhighz/adapters/legacy_inputs.py:124`, which rejects
nonuniform raw redshift nodes before creating its spline. The LBG/LAE reference
nodes [2.38,2.60,2.83,3.07,3.29] have spacings [0.22,0.23,0.24,0.22]. Three
required density comparisons are blocked, as the implementation report honestly
records. The revision-1 plan itself required uniformity and these comparisons;
the implementer followed that incompatible requirement.

Revision 2 must explicitly support unequal forecast widths, gaps and independent
evaluation redshifts, with regression coverage and an unequal-bin example.
For raw cell counts, accept caller-supplied positive per-redshift-cell widths in
sorted axis order and divide each normalized row by its own width times dm.
Do not infer physical cell boundaries from centres. An explicit labeled
legacy_first_spacing policy must reproduce the existing reference's constant
z[1]-z[0] divisor even for irregular nodes. Never select it automatically.
Raw source cells and forecast bins remain distinct. Keep magnitude-grid scope
unchanged, and document both policies and their provenance.

Closure requires analytic nonuniform polynomial/count-measure tests, malformed
width rejection, uniform regressions, legacy-mode scalar comparisons, all three
currently blocked real-table checks and the known negative LBG query. Preserve
the original reference nodes and evidence; physical widths are not to be invented
for those files merely to force legacy equality.

### R2 (P2): input snapshot bypasses scientific dtype validation

`fishhighz/survey.py:17` calls `geometry._immutable` for any ndarray;
`fishhighz/geometry.py:16` casts it to float64. `ForestInput.__post_init__` invokes
this at `survey.py:51` before `prepare_forest_weights` validates the data.
Consequently complex imaginary parts are discarded with ComplexWarning, while
Boolean, numeric-string and numeric-object arrays become accepted floats.

The independent probe reproduces **20 bypasses**: each of magnitudes,
quadrature, rho, variance and supplied weights with each of those four dtypes.
The direct Step 10 API rejects every case, while ForestInput followed by
prepare_bin returns finite noise. For example rho=[0.01+1j,0.02+1j] is retained
as [0.01,0.02] and used in a forecast. Container-dependent coercion can silently
change scientific input and violates the established real-array boundary.

Preserve types/values until scientific validation, then own accepted float64
arrays. Keep metadata dtypes faithful. Do not suppress the warning or relax
Step 10 validation. Closure requires direct-versus-orchestrated rejection
regressions, list/ndarray coverage, nonfinite rejection and valid integer/float,
zero and immutable ownership controls, including installed-wheel checks.

## Newly executed checks

Evidence is separate in `.validation/step11-review-r1/`:

- `before.json`, `git-before.txt` and governance snapshots preserve review input.
  All **101 entries** in the implementation's after.json matched live files
  before review edits (`identity.json`). The existing dirty/untracked state was
  retained; Git HEAD alone does not represent this implementation.
- `PATH="$PWD/.venv/bin:$PATH" scripts/check.sh`: **723 passed in 31.43s**;
  Ruff lint and format passed (96 files). Exit 0. MUNGE diagnostics in the log
  did not prevent successful completion. Thread limits were one.
- `.venv/bin/python .validation/step11-review-r1/independent.py`: two analytic
  unequal-bin volume checks and four joint Fisher batch checks passed; the
  20 R2 bypasses above were reproduced. Script and exact outputs are retained.
  Test helpers supply fixtures; expected volume and shell-mode integrals are
  independently evaluated and do not reuse production volume/q_mode values.
- Copied/read `reference.py` rerun against the new readers: **34 reference
  hashes match; 17 comparisons pass; three remain blocked** by nonuniform z.
  Saved Step 10 legacy numbers are historical comparison values; current reader
  evaluations are new. This did not execute a legacy forecast.
- Copied/read `probe.py`, helper scripts and all five exact examples rerun
  outside the checkout (`/tmp`, Python `-I`) in the **existing reported wheel
  environment**. All 32 source/wheel/installed modules and METADATA/WHEEL match;
  five examples, scalar Fisher/geometry, blocked-optional normalized multi-bin,
  optional Astropy and Step 10 range checks pass. This is a probe rerun, not a
  new fresh installation. Exact wheel SHA-256:
  `3d7f2b9aac75b0a0689e3212d867158e330f2db06ae4fe0010d6a1b52081951b`.

Commands/logs, JSON numerical values and copied probe scripts are in that review
directory. The wheel remains in `.validation/step11-r1-implementation/dist/`.
Python 3.13.15, NumPy 2.5.3, SciPy 1.18.1 and Astropy 8.0.1 were verified by the
installed probe. No rebuild or real full forecast was requested or run.

## Handoff

The current assignment is revision 2, incorporating R1 and R2 and retaining all
other scope and acceptance requirements. Design, roadmap and package guidance
are synchronized. The implementation report and historical evidence remain
unchanged. Await user-controlled repair dispatch and subsequent independent
review; do not plan Step 12.
