# Step 12 revision 1 independent review

Date: 2026-09-13. Outcome: **changes required** for R1 below. The user also
requested a larger scientific comparison, now specified in revision 2 of
[the same assignment](../IMPLEMENTATION_STEP.md). No acceptance or progression
to Step 13 is implied.

## Scope and preservation

Reviewed the revision-1 assignment and [implementation handoff](step-12.md),
live adapters, recipes, real/synthetic validation, evidence controller, tests,
and relevant read-only lyaforecast source. HEAD remains
`0d69786a06d5d676564a51fad14a7156951c2228`. All 118 entries in the implementation's
final manifest match live files. All 32 previously accepted Step 11 production
modules match their accepted wheel. Existing dirty/untracked work is preserved.

New review evidence is in `.validation/step12-review-r1/` relative to the package:
original assignment/governance snapshots, `before.json`, `git-before.txt`,
`identity.json`, quick/probe logs, and independent numerical probes. No
production code, tests, examples, historical reports, or reference assets were
edited. This review updates only current planning/governance documents and adds
this report and separate review artifacts.

## R1 — P2: evidence checker accepts invalid scientific payloads

Locations: `fishhighz/validation/evidence.py:103` (worker payload acceptance),
`:155` (optional pass flag), and `:166` (payload validation).

`execute()` and `check()` validate stored arrays against an inventory derived
from those same arrays, without binding their scientific meaning to the declared
case. Only finiteness/dtype, hashes, presence of Fisher/errors, and self-reported
shapes are checked. A missing `report.passed` is accepted. Therefore an ordinary
worker bug can produce an apparently complete, offline-validated bundle.

The review's `independent.py` generated five bundles through the public worker
interface; every bundle incorrectly passed `check()` with `complete=True`:

| Probe | Scientifically invalid content |
| --- | --- |
| wrong_dimensions | Fisher is shape `(1,)` with value 7; errors have length 3. |
| negative_information | Fisher is `-I_2` with positive errors. |
| wrong_errors | Fisher is `I_2`, errors are `[30, 40]`. |
| wrong_pair_payload | Payload selected pair is `[[99, 99]]` for the declared 15x2pt case. |
| missing_pass | Report is `{}`, without an explicit successful scientific status. |

These probes have consistent hashes; they expose semantic validation gaps, not
stale-file detection. Results and failing bundles are preserved alongside the
script. This finding does **not** invalidate the actual saved one-bin Fisher,
which passed the independent reconstruction below.

Required repair: define an explicit validation payload contract distinguishing
synthetic, real, profile, and diagnostic results. Bind array dimensions and
contents to case/field/pair/bin/node/parameter inventories; validate Fisher
symmetry, information sign/rank and derived errors/correlations; require boolean
scientific success and independently evaluate stored acceptance metrics. Reject
incomplete or contradictory evidence even if its hashes agree. Add regression
tests for each original probe and self-consistently rehashed corruptions of
valid new bundles. Preserve historical schema limitations rather than silently
upgrading old artifacts. Revision 2 supplies the detailed closure requirements.

## New verification performed during this review

| Check | Result |
| --- | --- |
| Source identity | All 118 implementation-manifest entries match; 32 accepted prerequisite modules unchanged. |
| Ordinary quick runner | **1017 tests passed in 23.93 s**, Ruff lint/format passed (113 files); `quick.log`. |
| Existing isolated installed-wheel probe rerun | Passed outside checkout using `-I`; all 40 source/wheel/installed modules and package metadata match. Five retained examples and earlier numerical oracles pass. This was not a fresh installation. |
| Step 12 installed probe rerun | Seven synthetic case selections, adapter branches, bridge units/order and optional-dependency blockers pass; `step12-probe.json`. |
| Independent real-array reconstruction | Rebuilt the full field matrix, Wick covariance, observed Jacobian and Fisher via independent NumPy solves from historical `real-02` arrays. Fisher relative Frobenius difference **5.016710052738691e-15**; errors `[0.017415176279183802, 0.012261442480497807]`; `independent.json`. |
| Actual reference density/SNR methods | **104 newly evaluated comparisons** pass, maximum relative difference **3.1516201241062103e-16**; all 34 prior reference source/resource hashes match; `policy.json` and `policy.log`. |
| Semantic evidence probes | Five invalid bundles accepted, confirming R1. |

Quick command: `PATH="$PWD/.venv/bin:$PATH" scripts/check.sh`. Numerical probes
used `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`. Installed probes
used the existing `.validation/step12-r1-implementation/wheel-env/bin/python`
from `/tmp`, with `-I` and copied scripts/examples in the new review directory.
The exact wheel SHA256 is
`b240003ad6497d77fe850107152765750065b027c18fb79f92ccfe157f00847a`.

The bounded raw-reference command used the existing read-only
`../lyaforecast/.validation/dev-env/bin/python`, explicit FishHighz `PYTHONPATH`,
and `scripts/legacy_policy_checks.py --reference ../lyaforecast --saved
.validation/step10-r1-implementation/reference.json --output
.validation/step12-review-r1/policy.json`. It does not run CAMB or NewForecast.

## Historical evidence and the user's expanded comparison

Revision 1 correctly ran one real 15x2pt bin, not seven full real cases. Its
one-bin grid/volume/derivative refinements and actual intrinsic-P3D amplitude
demonstration remain historical implementation evidence. Reconstructing saved
arrays here is a new independent calculation, **not** a newly executed forecast.
The reported differences of about -3.23% and -0.89% in ap/at uncertainties
relative to historical legacy output are unmatched comparisons. They do not
isolate the effects of individual conventions.

The user now requires all seven authoritative cases, every original bin, two
FishHighz settings per case, comprehensive numerical comparisons, plots, and
documented explanations of significant differences. This is new revision-2
scope, not a missing revision-1 acceptance requirement. The user confirmed:

1. A validation-only maximum-compatibility path reproducing legacy model,
   derivative, modes and pair-specific inputs, with independent FishHighz
   covariance/Fisher assembly.
2. A converged existing-model accuracy profile: integrated volumes, refined
   k/mu/magnitude quadrature, full wiggle derivatives, physical FWHM resolution,
   CAMB growth and per-field noise, keeping survey samples and ap/at targets.
3. Retained approved floors/clamps and legacy first-spacing conversion where
   physical input information is missing, with sensitivity studies.

Revision 2 authorizes those full reference and FishHighz runs for the dispatched
implementation assignment. They were not executed during this review because
the requested profiles and repaired evidence checker must first be implemented.
No implementation agent was dispatched, no Slurm action was taken, and nothing
was committed or pushed.
