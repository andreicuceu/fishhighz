# Step 11 revision 2: independent review

Date: 2026-09-13. Outcome: **passed independent review**. R1 and R2 are resolved;
no further required changes. User acceptance and progression remain pending.

Reviewed the current implementation, tests, example and handoff against the
revision-2 assignment. Only two production modules changed relative to the exact
revision-1 wheel: `fishhighz/adapters/legacy_inputs.py` and `fishhighz/survey.py`.
The existing runner, accepted numerical engines and historical report remain
intact. No production fixes were made during this review.

## Finding closure

- **R1, nonuniform redshift support:** the density reader now accepts explicit
  positive per-redshift-cell widths aligned with sorted redshift nodes, divides
  normalized raw counts by the correct row measure, and records the widths and
  policy. Uniform compatibility and explicitly selected legacy_first_spacing
  remain distinct; irregular nodes without a policy fail. No physical boundaries
  are inferred from centres. Tests and the example exercise unequal forecast
  bins, non-midpoint evaluation redshifts, gaps, area differences and caller
  order. Each volume uses its own bounds; source-table cells remain separate.
  All three formerly blocked LBG/LAE comparisons and the known negative LBG
  query now execute successfully under the documented legacy policy.
- **R2, dtype validation:** immutable input snapshots preserve non-object array
  dtypes until existing scientific validation. Complex/Boolean/string data are
  rejected without coercion; object arrays fail at snapshot setup. Valid
  scientific arrays become immutable float64 at the established boundary.
  Metadata values and dtypes are retained. All 20 original bypasses now raise
  ValueError, with warnings promoted to errors in the independent replay.

## New review evidence

All new artifacts are in `.validation/step11-review-r2/`. Old implementation and
review directories were not overwritten. Numerical commands used one thread
each for OMP, OpenBLAS and MKL on the login node.

| Check | New result |
| --- | --- |
| Live source identity | All 104 entries in implementation after.json match before review documentation edits. Independent comparison to the r1 wheel confirms only the two intended production changes; original r1 report hash matches its preserved copy. |
| `PATH="$PWD/.venv/bin:$PATH" scripts/check.sh` | **833 tests passed in 33.77s**, Ruff lint and format passed (99 files), exit 0. MUNGE diagnostics are retained in quick.log and did not prevent completion. |
| `.venv/bin/python .validation/step11-review-r2/independent.py` | Twelve newly generated irregular grids: 24 physical polynomial query checks and 24 selected raw-count measure checks, plus 12 independent legacy-divisor table checks. Maximum query relative error **4.663e-15**. All 20 original R2 rejection cases and three non-contiguous/endian metadata controls pass. |
| Copied reference.py rerun | **34 matching reference hashes, 33 passing newly evaluated reader comparisons, zero blocked**. Both unchanged negative QSO/LBG spline queries reject. Maximum comparison relative discrepancy **3.513e-15**; zero absolute allowance. |
| Existing installed-wheel probe rerun | All **32 source/wheel/installed modules**, METADATA and WHEEL match. All five examples, scalar/unequal-bin Fisher oracles, both density policies, 35 direct/orchestrated dtype rejection controls, optional-dependency blockers, Astropy geometry and retained Step 10 range checks pass. |

The independent density script uses randomized positive quadratic polynomials,
irregular centres and unrelated explicit widths. Expected off-grid values and
raw-count measures are computed independently, with rtol 5e-12 and 5e-13
respectively, atol 0. This supplements the implementation's fixed fixtures.

Reference output uses the exact saved Step 10 inputs/legacy values with matching
source/resource hashes, and newly evaluates revision-2 readers. It is a new
comparison against historical reference numbers, not a rerun of legacy forecasts.

The reviewed exact artifact is
`.validation/step11-r2-implementation/dist/fishhighz-0.1.0.dev0-py3-none-any.whl`,
SHA-256 `940eae9b77b9b4f9a1b07ce728a5d9addf25ee46a2e03228f4d60d004a36dce6`.
The probe ran from `/tmp` with `-I` in the implementation's existing wheel-env,
using copied scripts/examples and writing only new review artifacts. This is
an **existing-environment probe rerun**, not a new build or fresh installation.
Verified versions: Python 3.13.15, NumPy 2.5.3, SciPy 1.18.1, Astropy 8.0.1.

Scripts, exact commands represented by their logs, source/governance snapshots,
identity.json, production.diff, independent.json, reference.json and
wheel-probe.json retain the supporting evidence. No new real full forecast was
requested or run; it is not an outstanding acceptance gate.

## Handoff

The current assignment remains revision 2 with its acceptance contract preserved
and review status updated. Design, roadmap and package guidance record this pass,
distinct from user acceptance. No further repair instructions or Step 12 plan
were prepared. No production/test/example changes, dispatch, Slurm action,
commit or push occurred. Await the user's acceptance/progression request.
