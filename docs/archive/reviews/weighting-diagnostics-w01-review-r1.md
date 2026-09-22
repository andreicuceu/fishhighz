# Forest weighting diagnostic W01 revision 1 independent review

## Outcome

**Requires revision.** The submitted numerical result is independently
reproducible: on the saved `lya(qso)` auto-spectrum in bin 0, changing the
absolute update count from three to six changes both forest-noise coefficients
and both marginalized BAO errors beyond the prescribed `1e-3` diagnostic
trigger. This establishes finite-update sensitivity before profile changes for
this one saved example. It does not select a weighting prescription, establish
prevalence across the survey, or address an iteration or continuum limit.

The current saved bundle is numerically truthful, but its validator does not
bind several scientific claims or provenance fields to reconstructed values.
Two additional failure-path defects violate explicit revision-1 requirements.
The implementation and handoff remain unchanged; the active instructions are
revised to W01 revision 2 for a bounded evidence repair. W02 and production
changes remain unauthorized.

Review evidence is in
`.validation/forest-weight-diagnostics/w01-review-r1-20260915T210025Z/`.
The exact reviewed instructions remain preserved at SHA-256
`99fb517f5fabd67c48179bab437ab25e5d92b10a9a8abb139c1056a4a1fb0110`
in the revision-1 implementation bundle. No production source, implementation
script, test, historical report or saved input was edited during review. No new
forecast, reader/model evaluation, full suite, wheel, installation, Slurm
operation, agent dispatch, commit or push occurred.

## Independently verified numerical result

The reviewer resolved the authoritative compatibility/bin-0 report from its
manifest and independently checked its report and NPZ hashes. A separate
float64 recurrence and 80/160-digit scalar Decimal calculation, neither
importing nor calling the W01 implementation, gave:

| count | A [deg2] | P_pixel [deg2 km/s] | sigma(ap) | sigma(at) |
| ---: | ---: | ---: | ---: | ---: |
| 3 | 0.022118395206768268 | 0.5480269062497443 | 0.031927636457863354 | 0.03063714870732877 |
| 6 | 0.022308973469159228 | 0.5544833486755246 | 0.032107755826537035 | 0.030873464804937006 |

The independent relative changes are 0.008616278921204978 in A and
0.011781250796540466 in P_pixel. The corresponding Fisher matrices reproduce
the handoff to float64 summation precision. The marginalized error changes are
0.005641487709602044 and 0.0077133841620094. All 5000 saved Fourier nodes were
used with the stated scalar covariance; no 15-spectrum covariance was formed.

The original bundle's t=0/3/6 trajectories, moments, projected totals,
covariances, Fisher matrices and errors also match this independent calculation
within the assignment tolerance. The recorded early stop at six updates is
therefore scientifically appropriate. Counts 12 and 24 are correctly reported
as not attempted after identification.

## Required findings

### W01-R1 — high: evidence finalization does not bind the reported conclusion

Locations: `scripts/diagnose_legacy_forest_iterations.py:1072-1189` and
`:1319-1362`; the handoff's validator claim is at
`reviews/weighting-diagnostics-w01-r1.md:73-87`.

`validate_payload` reconstructs the source fingerprint, float64 recurrence and
projected Fisher matrices. It does not reconstruct or compare the reported
coefficient operands, terminal outcome, verdicts, table, Decimal confirmation,
forecast comparison, projected errors/totals/covariances/P1D/response, or the
saved-pair Fisher discrepancy. `validate_saved_evidence` authenticates only the
NPZ and instruction snapshot; it does not authenticate the recorded source
provenance, script hash, `checks.json`, or handoff.

Eight independent in-memory mutations of copied payloads were accepted:

- replacing `coefficient_sensitive` by `arithmetic_failure` and changing the
  verdict list;
- setting both reported coefficient changes to zero while retaining a true
  trigger;
- replacing the table's coefficients and errors;
- replacing the forecast comparison by zero change;
- changing projected errors, totals, covariance, P1D and response arrays;
- replacing the Decimal coefficients by zero;
- replacing the script hash; and
- replacing the source-report path and hash by unrelated values.

These are validator defects, not evidence that the unmodified W01 bundle is
false. The independent calculation above confirms its central result.

**Closure:** derive the expected source identity and provenance from
`load_source`, then reconstruct every reported attempt, coefficient comparison,
terminal/verdict, Decimal result, table row, projection product and forecast
comparison from authenticated arrays. Require the terminal reason to match the
first actual trigger or failure. Gate both saved single-spectrum Fisher
references. Authenticate an immutable diagnostic-script snapshot, final summary,
check inventory and handoff without a circular hash. Add a separate rejection
test for each mutation above, through the same full `--check-only` route used for
canonical evidence. Preserve truthful failed/capped and early-stop bundles.

### W01-R2 — medium: exact-zero scalar baselines can pass incorrectly

Location: `scripts/diagnose_legacy_forest_iterations.py:637-650`.

For scalar A and P_pixel, `_relative_change` returns `None` when a zero reference
is replaced by a nonzero value. The expression `abs(... or 0.0)` converts that
mismatch to zero discrepancy. A copied real source with `captured_a=0` therefore
passes `baseline_gate` against its nonzero recomputed A. Array weights already
use correct exact-zero handling.

The real W01 source has nonzero A and P_pixel, so this defect does not affect the
submitted numerical result.

**Closure:** use a zero-safe scalar comparison that accepts zero only against
exact zero and otherwise applies the per-quantity `5e-12` relative tolerance.
Add A-zero and P_pixel-zero mismatch controls plus matching-zero controls.

### W01-R3 — medium: baseline failures cannot produce the required bounded handoff

Locations: `scripts/diagnose_legacy_forest_iterations.py:910-915` and
`:1319-1348`; the incomplete test is at
`tests/test_legacy_forest_iterations.py:273-286`.

The baseline-failure test verifies only that evaluation raises before later
counts. `run_diagnostic` raises when count three is absent or the baseline gate
fails, so the CLI never reaches `write_evidence`. This loses the prescribed
attempt/failure record instead of producing the required source/arithmetic
finding and retained partial outcome. The same issue applies to an arithmetic
failure before the three-update gate.

The submitted real baseline succeeds exactly, so this does not alter its result.

**Closure:** make failed initialization, failed three-update arithmetic, and
baseline mismatch serializable terminal outcomes with the completed attempts
retained, no forecast projection, and no later evaluations. Add full CLI/evidence
tests for each route, including Decimal classification and successful
`--check-only` replay. A baseline mismatch must not receive
`reproduced_legacy`.

## Checks run

- Focused test: **23 passed in 0.43 s**.
- Ruff lint: **passed**.
- Ruff format check: **2 files already formatted**.
- Original canonical `--check-only`: **passed**.
- Independent source recurrence, Decimal and one-spectrum Fisher oracle:
  **passed**.
- Deliberate copied-payload semantic corruptions: **8/8 were incorrectly
  accepted**, establishing W01-R1.
- Exact-zero scalar mismatch probe: **incorrectly passed**, establishing
  W01-R2.

The smallest next action is W01 revision 2 only. The verified finite-update
sensitivity makes additional W01 counts unnecessary. After the repair is
independently reviewed, the user may stop this diagnostic sequence or choose the
next scientific question; this review does not prepare W02.

