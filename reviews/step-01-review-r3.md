# Step 01 review of implementation revision 3

Reviewed 2026-09-11 against `IMPLEMENTATION_STEP.md` revision 3, before that plan
was replaced by revision 4. This is reviewer evidence, not user acceptance.

## Outcome

The maintained scaffold is small and satisfies the substantive revision-3 scope.
Its source and artifacts match the implementation report, and independent checks
passed. Step 01 remains open for the user's direct-layout request and one
validation-workflow correction. No production source was changed during review.

## Findings

### R1: repeated validation can exercise an old installed wheel

Priority: medium; fix before using the helper to validate the revised layout.

In `.validation/validate.py:51`, `assess` recreates a venv at a fixed existing
path and then calls pip without replacing an already installed same-version
package. Creating a venv at an existing path does not clear its site-packages.
A dry-run installation against the retained wheel environment reported:

```text
fishhighz is already installed with the same version as the provided wheel.
```

Thus a changed `0.1.0.dev0` wheel could be built and inspected while import checks
still run the old installation. The module-origin check verifies site-packages,
but does not identify which wheel supplied the installed bytes. The helper also
uses `next(glob(...))` for artifacts at lines 65, 67, and 74; that selection becomes
ambiguous if multiple builds are retained in the same directories.

Use fresh environments and output directories per run, select each exact new
artifact unambiguously, and compare installed package contents with its payload.
Demonstrate two consecutive executions without manual environment cleanup.
This finding does not invalidate the current artifacts: both were independently
installed successfully into fresh environments during this review.

### U1: user requests a direct package directory

The existing `src/fishhighz` structure correctly followed revision 3. The user now
requests `fishhighz/` directly under `lib/fishhighz`. Update package discovery,
egg-info ignore paths, documentation, editable installation, and validation.
Keep isolated installed-artifact checks so the new layout cannot conceal an
installation defect through imports from the working directory.

## Independent checks performed

- All seven source/configuration/test hashes and all three artifact hashes in
  `reviews/step-01.md` matched the files on disk.
- Ran `scripts/check.sh` using `.venv/bin` on PATH, cleared PYTHONPATH, and an
  unrelated temporary working directory: one pytest test passed; Ruff lint and
  formatting checks passed.
- Created separate new environments under `/tmp` and installed each of the
  reported direct/rebuilt wheels with `pip install --no-index --no-deps`.
  Both isolated import/metadata probes passed with site-packages origins. The
  environments were temporary and removed after these checks.
- Verified that both wheel payloads contain the current `fishhighz/__init__.py`
  bytes and that the source archive preserves the runner's executable mode.
- Read the retained build driver, source archive contents, build/installation
  evidence, package configuration, README, runner, smoke test, and handoff report.

Fresh review logs are in
[`../.validation/review-r3-1pgkatic/`](../.validation/review-r3-1pgkatic/):
`checks-outside.log`, `wheel.log`, and `rebuilt-wheel.log`.

The review reused the hashed distribution artifacts for independent installation;
it did not rerun the network-dependent distribution builds. Original build logs
and hashes remain in the implementer's report. Validation was on Python 3.13.15;
no claim is made that Python 3.11 or scientific behavior was tested.

The required changes and acceptance checks are consolidated in Step 01,
revision 4. No Step 02 assignment has been written.
