# Forest weighting diagnostic W01 revision 2 independent review

## Outcome

**Requires revision.** W01-R1, W01-R2 and W01-R3 are substantively repaired.
The authoritative source, recurrence, Decimal controls, forecast projection and
original corruption classes now bind correctly, exact-zero scalar baselines are
handled correctly, and failed-baseline outcomes are serializable. The supported
W01 scientific result is unchanged.

Two residual evidence defects prevent revision 2 from passing. A non-finite JSON
number can bypass the semantic comparison for any expected nonzero float. The
final check inventory and execution record are authenticated as bytes but are not
validated as the required checks; the canonical inventory also claims successful
finalization and final `--check-only` before those operations could have produced
the manifest that authenticates it. These defects do not change the unmodified
numerical result, which was independently reproduced again.

The implementation, revision-2 handoff and evidence remain unchanged. The active
instructions are revised to W01 revision 3 for this bounded finalization repair.
No production source or historical artifact was edited, and W02 remains
unauthorized.

Review evidence is in
`.validation/forest-weight-diagnostics/w01-review-r2-20260915T212925Z/`.
The exact revision-2 instructions are preserved in the implementation bundle at
SHA-256 `3ec253993f752fafc144f6df1f6cd11c55a0dde5ae703d0986147b48f351b481`.

## Repairs that pass independent review

### W01-R1: source, trajectory and forecast binding

The finalized canonical bundle passes its authenticated `--check-only` route.
Independent copied-bundle mutations of the terminal/verdict, authoritative
provenance, projected errors and script snapshot reject after their dependent
hashes are refreshed. The 48 focused tests include all fourteen revision-1
corruption classes and pass.

The validator reconstructs the source fingerprint and full provenance from the
authoritative record, all recurrence states and Decimal values, coefficient
comparisons, first terminal event, table, forecast products, both saved pair
Fisher references, and every named NPZ array. This closes W01-R1 for finite JSON
values; the non-finite bypass below is a new boundary defect.

### W01-R2: exact-zero baseline handling

The new scalar comparator distinguishes exact zero from a nonzero mismatch and
retains the `5e-12` relative tolerance for nonzero values. Focused tests exercise
zero/zero and both zero/nonzero directions independently for A and P_pixel. Code
inspection and the passing controls close W01-R2.

### W01-R3: failed-baseline evidence

Initialization singularity, count-three float64 failure, Decimal/float64
disagreement and captured-baseline mismatch now produce terminal bundles with
partial states, explicit later-stop records, empty projection arrays and no
`reproduced_legacy` verdict. Each synthetic route is written, finalized and
rechecked through the CLI. A separate control correctly classifies a float-only
singularity with finite Decimal results as arithmetic failure. This closes
W01-R3.

## Scientific result independently retained

An independent float64 recurrence, 80/160-digit Decimal calculation and scalar
one-spectrum Fisher sum again reproduced the authoritative report/NPZ hashes and
the revision-2 arrays:

| count | A [deg2] | P_pixel [deg2 km/s] | sigma(ap) | sigma(at) |
| ---: | ---: | ---: | ---: | ---: |
| 3 | 0.022118395206768268 | 0.5480269062497443 | 0.031927636457863354 | 0.03063714870732877 |
| 6 | 0.022308973469159228 | 0.5544833486755246 | 0.032107755826537035 | 0.030873464804937006 |

The independent relative coefficient changes are 0.008616278921204978 and
0.011781250796540466. The t=3 and t=6 Fisher matrices reproduce the handoff to
float64 summation precision. Counts 12 and 24 remain correctly unattempted after
the first confirmed discrepancy.

## Required findings

### W01-R4 — high: NaN bypasses reconstructed finite JSON values

Location: `scripts/diagnose_legacy_forest_iterations.py:1267-1303`; the missing
boundary control is visible in `tests/test_legacy_forest_iterations.py:368-452`.

For an expected float, `_assert_semantic_equal` computes a discrepancy and tests
`discrepancy > tolerance`. If the submitted value is NaN, the discrepancy is
NaN and that comparison is false. Python's default `json.load` accepts the
non-standard `NaN` token. On a copied canonical bundle, replacing the finite
Fisher relative change by NaN and refreshing the summary digest in the manifest
passes the complete authoritative `--check-only` route.

The same bypass applies to nonzero coefficient changes, table values, baseline
metrics and projection checks. The writer uses `allow_nan=False`, so the
unmodified canonical bundle contains no such value; this is a validation defect,
not evidence of corrupted submitted results.

**Closure:** reject `NaN`, `Infinity` and `-Infinity` while parsing every JSON
artifact, and require submitted/reconstructed numerical scalars to be finite
where the schema expects a finite value. Continue representing unavailable JSON
quantities as `null` and preserve the intentional NPZ NaNs used for unavailable
errors. Add refreshed-manifest full-route mutations for the coefficient trigger,
forecast comparison, table and projection checks, plus non-finite manifest/check
metadata.

### W01-R5 — medium: execution and required-check claims are not validated

Locations: `scripts/diagnose_legacy_forest_iterations.py:1373-1404` and
`:1643-1703`, `tests/test_legacy_forest_iterations.py:341-357`, and the canonical
`checks.json` rows quoted by `reviews/weighting-diagnostics-w01-r2.md:88-99`.

`validate_payload` deliberately omits the entire `execution` mapping.
`validate_saved_evidence` accepts any nonempty check list whose rows contain a
command and zero exit code. Two copied canonical bundles therefore pass after
refreshing manifest hashes when either:

- the execution command, duration, Python executable/version, platform, NumPy
  version and all thread limits are replaced by fabricated values; or
- all six check rows are replaced by one fabricated passing row.

The canonical check inventory itself is temporally self-referential. Its
`created_utc` is `2026-09-15T21:20:17Z`, while the manifest that makes
`--check-only` possible was created at `21:22:41Z`. Yet the sealed checks claim
that both finalization and final `--check-only` had already passed. Finalization
requires this checks file as an input and hashes it, so those two rows cannot be
observations completed before the file was sealed. The independent review did
run final `--check-only` successfully; the defect is the saved provenance claim.

**Closure:** define and validate the exact pre-finalization check categories,
their unique commands, successful exit codes and finite nonnegative timings.
Cross-check the numerical-run command and one-thread settings against the summary
execution record, and validate nonempty version/executable fields. Seal only
checks that have actually completed before finalization. Record finalization and
post-manifest `--check-only` as later external observations rather than entries in
their own prerequisite artifact, or use a clearly ordered outer record that does
not claim to authenticate itself. Add full-route rejection tests for omitted,
duplicated and fabricated checks and for altered execution metadata.

## Checks run

- Focused suite: **48 passed in 1.59 s** pytest time; the login-node process
  completed successfully.
- Ruff lint: **passed**.
- Ruff format check: **2 files already formatted**.
- Canonical authenticated `--check-only`: **passed**.
- Independent recurrence, Decimal and Fisher calculation: **passed**.
- Representative original mutations of terminal/verdict, provenance, projected
  errors and script snapshot: **all rejected**.
- New copied-bundle probes: NaN scientific summary, fabricated execution metadata
  and replacement of required checks by one fabricated pass were **all
  incorrectly accepted**.

No NewForecast/controller, model or reader evaluation, additional bin/population,
joint forecast, package-wide suite, wheel, installation, Slurm action, agent,
commit or push was used. The smallest next action is W01 revision 3 only; this
review does not prepare or authorize W02.

