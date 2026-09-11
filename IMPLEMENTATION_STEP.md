# Step 01: revise the package layout and validation workflow

Revision: 4 (2026-09-11).
State: changes requested after review of revision 3; ready for user dispatch.
Implementation root: `lib/fishhighz`.

Read [AGENTS.md](AGENTS.md), the
[roadmap](../../FISHHIGHZ_IMPLEMENTATION_PLAN.md), and the
[design](../../FISHHIGHZ_DESIGN.md). This is a revision of the existing scaffold
step, not a new implementation stage. Preserve working behavior and complete
only the changes and validation below.

## Reviewed starting point

The revision-3 scaffold and [implementation report](reviews/step-01.md) were
reviewed. All ten reported source/artifact hashes matched. The reviewer reran the
check script from outside the checkout (one pytest test, Ruff lint, and format
verification passed), installed both reported wheels in fresh offline environments,
and verified isolated imports and source-archive runner permissions.
See [review findings and evidence](reviews/step-01-review-r3.md).

Two items remain before this step can be accepted:

- **User request:** use `fishhighz/` directly under the package repository root,
  removing the extra `src/` directory. The previous layout followed revision 3
  correctly; this is a changed preference, not a defect in that implementation.
- **Reviewer finding:** repeated execution of the validation helper can reuse a
  previous wheel installation. `python -m venv` does not empty an existing
  environment, and pip can skip a changed wheel with the same version. Artifact
  selection by the first glob result can also pick an older build if several
  artifacts are present. The revision must eliminate this stale-evidence risk.

The current wheels themselves passed independent fresh-install checks. There
is no request to redesign the scaffold or introduce scientific functionality.

## Required revision work

### 1. Move to the direct package layout

The resulting maintained structure should include:

```text
lib/fishhighz/
    fishhighz/
        __init__.py
    tests/
        test_import.py
    scripts/
        check.sh
    pyproject.toml
    MANIFEST.in
    README.md
    .gitignore
    AGENTS.md
    IMPLEMENTATION_STEP.md
    reviews/
```

Move the existing minimal module into `fishhighz/__init__.py`, preserving its
behavior. Remove the obsolete `src/` tree after verifying it contains only the
moved source and generated package metadata/caches. If unexpected files are
present, preserve them and report the discrepancy. Do not broadly delete other
builds, environments, or review evidence to accomplish the move.

Update setuptools discovery in `pyproject.toml` to use the repository root and
explicitly include only `fishhighz` and its package descendants. Disable implicit
namespace discovery or otherwise ensure tests, scripts, reviews, and generated
validation directories cannot become import packages. Keep installed import name
`fishhighz`, version `0.1.0.dev0`, Python >=3.11, empty runtime dependencies, and
the existing small `dev` extra. No new scientific modules or dependencies.

Update the root-level egg-info ignore rule and any layout-dependent build,
manifest, README, or validation references. Keep the executable check runner and
single-source version definition. Review source-distribution and wheel contents;
metadata, tests, and the documented runner may be in the source archive, but only
the intended import package and distribution metadata belong in the wheel.

Reinstall the editable package after the move so its installed mapping no longer
points to `src`. Do not rely on running Python from the repository root: the new
layout makes an uninstalled checkout importable there, which can hide packaging
errors. Preserve isolated subprocess imports from outside the checkout without
PYTHONPATH or sys.path workarounds.

### 2. Make repeated validation test the current build

Update the retained validation helper under `.validation/`, or replace it with
an equally reproducible sequence recorded in the report. Keep this a small
validation aid rather than introducing a packaging test framework.

- Give each run a new validation-output directory and new wheel-install
  environments, including a separate environment for the wheel rebuilt from the
  source distribution. Environments must have no system site packages or dev
  extras. Do not treat forced reinstallation into an old environment as proof of
  a clean install.
- Produce/select artifacts in fresh run-specific output directories. Identify the
  wheel and source archive from that build explicitly and fail on unexpected
  multiplicity; do not select the first match from a directory containing older
  runs. Build from a clean source staging directory if needed to avoid stale
  build/egg-info contents. Validate the resulting contents against the new layout.
- Install the exact selected wheel with `--no-deps`, clearing PYTHONPATH and
  importing with isolated Python from an unrelated working directory. Confirm
  site-packages origin and matching installed metadata. Compare the installed
  package file bytes with the selected wheel payload so the evidence identifies
  the artifact actually exercised, not just its version number.
- Repeat the same checks for the wheel rebuilt from a separately extracted source
  archive. Avoid any dependence on adjacent checkout files.
- Retain logs, artifact paths/hashes, and helper/probe sources for each run.
  Preserve the previous revision's evidence; do not overwrite its final logs or
  artifacts. Run the updated workflow twice without manual cleanup, and demonstrate
  that both runs use distinct fresh environments and their own exact artifacts.

### 3. Refresh the handoff report

Before replacing `reviews/step-01.md`, preserve its existing revision-3 contents
as `reviews/step-01-r3.md` unless already archived. Keep the reviewer's
`reviews/step-01-review-r3.md` unchanged. Update the main handoff report to identify
revision 4, address both findings, and record the current evidence and remaining
limitations. Hash source/configuration/tests and the actual validation helper and
probe used, along with built artifacts. Do not hash a report into itself.

## Acceptance assessments

| Check | Required result |
| --- | --- |
| Layout and discovery | `fishhighz/__init__.py` is directly under the repository root; no obsolete `src/` remains. Wheel/source-archive contents reflect the new layout, and package discovery excludes unrelated directories. |
| Editable development workflow | Reinstall the current editable package in the development environment. Run the executable `scripts/check.sh` from the repository root and from an unrelated directory. Pytest, Ruff lint, and Ruff format verification all pass. |
| Fresh-process import | The existing quiet-import test passes outside the checkout with isolated Python, no PYTHONPATH injection, and neither Vega nor lyaforecast imported. Verify the editable install resolves to the new package directory. |
| Built-wheel independence | Install the exact newly built wheel with `--no-deps` into a fresh environment. Outside the checkout, verify quiet import, site-packages origin, installed metadata, and payload agreement with that wheel. |
| Source distribution completeness | Independently extract the new source archive, build its wheel with declared isolated build requirements, and repeat fresh installation/import/payload checks. Confirm the documented check runner is included and executable. |
| Repeat-run reliability | Run the updated validation workflow twice with unchanged source, without manual cleanup. Logs identify distinct new environments and output directories and exact artifacts for each run. Neither run skips installation because of an existing FishHighz installation. Byte-identical wheel hashes across builds are not required. |
| Scope and evidence | Updated report addresses both findings and includes commands, versions, exit statuses, artifacts, source/helper hashes, and limitations. Earlier evidence remains identifiable. No scientific implementation, sibling edit, invented metadata, commit, or publication. |

These retain the original installation/import checks while explicitly guarding
against source shadowing and stale builds during the layout migration. Use small
behavioral assessments; do not add tests that merely repeat packaging constants.
No multi-version Python matrix or scientific baseline is required for this step.

## Boundaries and completion

Only the scaffold, its development/validation aids, and its handoff report are in
scope. Do not implement forecasting APIs, add numerical runtime dependencies,
run DESI-2 forecasts, initialize Git, commit, push, publish, or advance the roadmap.
NumPy/Numba-friendly kernels and 2D P(k,mu) remain the agreed later numerical path.

The planning agent has updated the design/roadmap/AGENTS.md for the layout choice;
leave those documents and this plan unchanged during implementation unless the
user requests another planning revision. If an environment restriction prevents
a required check, retain the error and mark the assessment incomplete.

Complete the implementation and all assessments, write the handoff report, and
stop for the user's review. Passing checks is not user acceptance and does not
authorize another step.
