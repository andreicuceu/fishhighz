# Step 02: add quick reference validation and close review findings

Revision: 3 (2026-09-11).
State: revision implemented and independently reviewed; no outstanding findings,
awaiting the user's next-step instruction.
Implementation root: `lib/fishhighz`.

Review outcome: [revision-3 review](reviews/step-02-review-r3.md) passed. The
requirements below are retained as the reviewed assignment, not a request to
repeat implementation or validation. No Step 03 plan has been drafted.

Read [AGENTS.md](AGENTS.md), the [roadmap](../../FISHHIGHZ_IMPLEMENTATION_PLAN.md),
the [design](../../FISHHIGHZ_DESIGN.md), the
[revision-1 review](reviews/step-02-review-r1.md), and `../lyaforecast/AGENTS.md`.
Implement only this revision when the user dispatches it; report and stop for
review. Do not start Step 03.

User decision for this revision: run quick checks by default. Execute the real
full suite only when the user explicitly asks for that run. A general request
to implement, test, review, or finish a step does not authorize full execution.

## Starting point and scope

The revision-1 implementation and [handoff](reviews/step-02.md) were reviewed.
The reviewer verified the five reported implementation hashes and acceptance
manifest hash, reran all 13 portable tests and Ruff checks, and checked the saved
seven-case bundle plus repeat. Independent inspection confirmed that its stored
grids/redshifts agree with configurations and both tool snapshots are intact.
No real forecast was rerun during review.

Preserve the existing working capture and its evidence at
`.validation/baseline/20260911T184946Z-7a6a795d/`. Revise the maintained controller,
worker only as needed, tests, README, and implementation report. Fix the four
review findings below and add a quick validation mode. Reuse existing functions;
avoid a new workflow framework or broad rewrite of the controller.

Do not implement FishHighz numerics, add runtime scientific dependencies, modify
neighboring packages/environments, alter authoritative science inputs, or commit
and push. Preserve the accepted direct package layout and import checks. Leave
planning/governance updates to the planning agent unless the user instructs
otherwise. Retain the revision-1 handoff as `reviews/step-02-r1.md` before updating
`reviews/step-02.md` for this revision.

## 1. Introduce explicit validation levels

Implement this command contract with `--suite quick` as the capture default.
Full execution requires an explicit `--suite full` option and, for agents, the
user's explicit request. Preserve old bundle read/check compatibility, not the
old expensive capture default:

```bash
# Ordinary development checks: synthetic tests and lint only.
./scripts/check.sh

# Quick scientific reference check: exactly one full-resolution 15x2pt run.
python scripts/lyaforecast_baseline.py capture --suite quick \
  --baseline .validation/baseline/<existing-full-bundle> \
  --reference-checkout ../lyaforecast \
  --python ../lyaforecast/.validation/dev-env/bin/python

# Complete reference capture: run only when explicitly requested by the user.
python scripts/lyaforecast_baseline.py capture --suite full \
  --reference-checkout ../lyaforecast \
  --python ../lyaforecast/.validation/dev-env/bin/python

# Offline integrity/result checks; never launch scientific workers.
python scripts/lyaforecast_baseline.py check <bundle>
python scripts/lyaforecast_baseline.py check <bundle> --require-suite full
```

### Quick mode

- Run only `lya_qso_lbg_lae_15x2pt.ini`, exactly once in a fresh process. Do not
  run a second repeat or compute other cases. Preserve its complete scientific
  settings, including CAMB, redshift/magnitude bins, k/mu grids, and all 15 pairs.
- Require an explicit saved full baseline and validate it before launching the
  forecast. Reject quick/incomplete/invalid bundles as baseline inputs.
- Before the expensive run, check compatibility of the relevant reference
  source, original 15x2pt INI, resolved input bytes, Python version, scientific
  library versions, and thread settings. Compare content identities and relevant
  environment fields, allowing location-only differences. Controller changes,
  artifact directory names, and Python executable relocation alone must not make
  an otherwise identical science baseline incompatible. Record the exact fields
  compared; report mismatches and explain when a new full capture would be
  needed. Never fall back to full mode automatically. Complete independent quick
  checks and report any blocked comparison without launching a full capture.
- Compare all fresh results with the saved full bundle's primary 15x2pt result,
  checking keys/order, shapes/dtypes, and numbers at `rtol=1e-10`, `atol=1e-12`.
  Retain fresh results and failure evidence when the comparison fails. Never
  report success merely because numbers are finite.
- Preserve provenance, configs, inputs, outputs, logs, timings, and source/input
  checks for this one case. Do not run seven resolution workers just to prepare
  one case. Lightweight validation of the authoritative inventory is fine.
- Save enough baseline comparison evidence in the new quick bundle to check it
  after relocation without the external baseline directory. Record the original
  full manifest digest and the copied comparison payload's digest. Distinguish
  validating the full baseline at capture time from later checking the embedded
  comparison evidence; do not claim to revalidate uncopied full-bundle artifacts.

### Full mode and explicit coverage

Full mode retains seven primary cases and one independent 15x2pt repeat, each in
a fresh subprocess, with all revision-1 preservation, serialization, numerical,
configuration, and repeat checks. Required primary cases remain:

| INI | Selected pairs |
| --- | ---: |
| `lbg_lae_3x2pt.ini` | 3 |
| `lya_lbg_lae_3x2pt.ini` | 3 |
| `lya_lbg_lae_6x2pt.ini` | 6 |
| `lya_qso_2x2pt.ini` | 2 |
| `lya_qso_lbg_lae_4x2pt.ini` | 4 |
| `lya_qso_lbg_lae_8x2pt.ini` | 8 |
| `lya_qso_lbg_lae_15x2pt.ini` | 15 |

Record suite identity and expected coverage explicitly in new versioned
manifests. Derive permitted inventories from the supported suite, not arbitrary
manifest counts. A valid quick bundle must say quick/one case; it cannot satisfy
`--require-suite full`. A partially completed full capture cannot be relabeled
as a successful quick run. Reject conflicting suite, inventory, role, and repeat
records with useful errors.

Retain read/check support for the existing version-1 full bundle so it remains
usable as a baseline. Apply all checks supported by its recorded evidence;
identify new preservation fields that old captures lack without inventing them
or modifying old artifacts. New captures must include the complete revised
requirements. Keep compatibility handling small and explicit.

### Runtime and cadence

Keep cases serial and single-threaded in this revision. Quick mode is expected
to take about two minutes, based on the recorded 134-second 15x2pt case, while
ordinary tests took 10.63 seconds in review. Measure actual revised timings;
these estimates are not hard performance assertions.

Two concurrent workers could shorten full-suite wall time, but do not reduce
work and add memory/process-management costs. Reduced grids barely address the
observed initialization bottleneck. Do not add parallel execution, change grids,
or cache/share CAMB initialization in this revision.

Run ordinary tests after local edits and quick mode when checking the reference
capture/integration path. Do not execute the real full suite at step completion,
review checkpoints, after failures, or after relevant changes unless the user
explicitly requests it. Do not ask for full execution routinely or treat an
unrequested full run as a prerequisite for a review handoff. Record it as not run
under the user-selected validation policy. A request for one full run does not
authorize subsequent full reruns.

Portable synthetic tests may exercise full-mode orchestration without launching
real forecasts; offline checking of saved full bundles is also part of the quick
workflow. Future steps should reuse accepted reference evidence and run relevant
quick FishHighz tests/comparisons. An unchanged lyaforecast rerun alone is not a
test of newly implemented FishHighz numerics.

## 2. Close the four review findings

### A. Validate tool snapshots and required shared evidence

The reviewer appended a comment to the copied controller snapshot without
updating any checksum; the checker still passed. Validate both required tool
snapshot files against the digests recorded in `tool.json`. Require the shared
artifact inventory rather than trusting whichever keys happen to be present.
Reuse containment validation for paths reached through nested metadata. Missing,
corrupted, or out-of-bundle snapshots must fail with a specific diagnostic.

### B. Use the same worker checks for every case role

The checker accepted a repeat with exit code 17 and a failed stored status.
Factor common case-evidence validation for full primaries, repeats, and quick
runs. Require successful manifest/status exit codes and state, verified round
trips, consistent worker response/status evidence, correct case/role, and all
required artifacts. Ensure a repeat is a distinct captured run with its own
request, paths, and result evidence; aliasing primary artifacts cannot serve as
an independent repeat. Numerical agreement never overrides worker failure.

### C. Validate complete coordinate and identity metadata

The checker accepted an interior k value of -1000 and, separately, a result
redshift shifted by 0.001 relative to its unchanged metadata/configuration.

Check all k and mu values, their one-dimensional shapes, ordering, and agreement
with the actual reference grid formulas. Check result redshifts/edges against
worker metadata and original survey settings, including bin count, and check
fiducial redshift against its configured reference value. Reproduce only the
simple reference coordinate conventions in these standard-library checks; do
not import the scientific package. Support the settings used by the seven
actual INIs, reporting unsupported alternatives explicitly.

Derive tracer order, all pair identities, and selected pairs from the original
configuration using the reference naming rules, rather than checking counts
alone. Preserve forest sample identities and zero-filled unselected results.
Use explicit tight coordinate tolerances allowing floating-point roundoff; do
not require byte equality between independently constructed floating arrays.

### D. Preserve authoritative INIs and source/input inventories

Current before/after checks omit authoritative example INIs and miss newly added
source/SNR files. Capture original INIs before input resolution or forecasts,
hash them before/after, and use those captured originals consistently for
resolution, effective configurations, and repeat comparison.

Compare bounded source and consumed input-directory inventories before and after
capture as well as contents, including added/deleted files and directory entries
relevant to the spectrograph header reader. An added source module or SNR file
must invalidate an unchanged-source/input claim. Scope quick-mode preservation
to its consumed data and reference source, and full mode to all seven cases.
Do not broaden filesystem scans outside the specified reference roots.

## 3. Preserve existing guarantees

- Use the unchanged `NewForecast(path).new_run_forecast()` scientific API.
  Allow only output/input path substitutions in effective INIs; never patch
  scientific methods or alter settings to accelerate runs or force agreement.
- Preserve complete dictionaries, order, dtypes/shapes, exact serialization
  round trips, finite-result checks, positive selected/combined uncertainties,
  bounded correlations, and zero-filled unselected pairs.
- Keep NumPy/CAMB/reference imports in the scientific worker. Preserve the
  lexical virtual-environment interpreter path and verify actual module origins.
- Keep thread limits at one and disable reference bytecode writes. All writes
  stay in FishHighz or explicit new artifact locations; no Slurm actions.
- Create new destinations exclusively, retain failed evidence with nonzero
  status, and keep snapshots/outputs ignored. Never overwrite accepted evidence.
- Keep ordinary pytest independent of scientific dependencies and sibling
  installations; keep check commands offline and usable after relocation.

## 4. Required tests and acceptance evidence

Extend portable tests with synthetic bundles and stub workers. Cover:

1. Full mode launches exactly seven primary forecasts plus one repeat; quick mode
   launches exactly one 15x2pt forecast and no repeat. Count scientific `run`
   actions separately from lightweight probes/resolution, using stub evidence.
   Verify that omitting `--suite` selects quick, and that missing/incompatible
   baseline evidence fails without launching a real or stub full capture.
2. Quick success against a compatible full baseline; numerical and structural
   mismatches; invalid/quick baseline rejection; science/environment mismatch
   rejected before forecasting; location-only compatibility; and missing baseline.
3. Quick/full manifest coverage, `--require-suite full`, and malformed/partial
   captures. Check a relocated quick bundle after making the original baseline
   and reference paths unavailable within the test's temporary directories.
4. Missing/corrupted tool snapshots, missing shared evidence, and nested invalid
   artifact paths; failed repeat exit/status and unverified round trips; primary
   artifacts reused as repeat evidence. Exercise semantic failures with updated
   hashes so tests do not stop at checksum mismatch alone.
5. Interior k/mu mutations, redshift/result/metadata disagreement, wrong bin count,
   tracer/pair identity changes, and valid real-case coordinate conventions.
6. Changed original INI and added/removed source or consumed-directory entries
   during a synthetic capture. Preserve failure evidence and exit nonzero.

Retain the existing regression coverage and run `scripts/check.sh`. Then:

- Recheck the existing full bundle with the revised checker, reporting its legacy
  schema and any unavailable new fields honestly.
- Run quick mode once against that compatible full baseline. Record total wall
  time and comparison outcome, and check the relocated quick bundle offline.
- Exercise full-mode orchestration and revised preservation requirements with
  synthetic tests. Do not run a new real full capture by default. If the user
  explicitly requests one, verify seven primary cases, the separate repeat,
  complete revised preservation evidence, and offline relocation checks.
- Record exact commands, code hashes/Git state, environment, suite identities,
  artifact paths, timings, comparisons, failures, and remaining limitations in
  the updated handoff. Link the revision-1 report and review, and address each
  finding explicitly. Do not claim a measured parallel speedup or FishHighz
  scientific validation.

The revision is ready for review when the required portable tests, offline
checks, and quick validation pass. An unrequested real full run is not an unmet
acceptance requirement; report that it was not run and distinguish historical
full evidence from validation of the new code. If the user requests a full run,
include its result or failure in the handoff. The user reviews and requests
independent review, then decides whether to advance. No agent may start or draft
Step 03 as part of this assignment.
