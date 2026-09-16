# Step 10 handoff — revision 2

Implemented the bounded R1 repair in `IMPLEMENTATION_STEP.md`, revision 2
(2026-09-13), following the user's dispatch. Ready for user and independent
review; acceptance and progression are not claimed. No Step 11 work occurred.
The original revision-1 handoff is preserved byte-for-byte in
[step-10-r1.md](step-10-r1.md). Its results remain explicitly historical.

## R1 repair

The old left-associated `r*w*w*variance` could lose an individual positive term
to underflow before multiplication by a large variance restored its mathematical
size. Another positive term then concealed that loss from the final P_pixel
positivity check.

`fishhighz/kernels/weights.py:_integrals` now evaluates the three contribution
products (`r*w`, then `r*w*w`, then `r*w*w*variance`) under a local NumPy
`under="raise"` policy before accumulating prefixes. Arithmetic order and the
weighting formulas are unchanged. The existing host exception boundary in
`fishhighz/weights.py` converts FloatingPointError to ValueError and now retains
the underlying arithmetic diagnosis alongside the field ID.

This uses the assignment's permitted **explicit rejection** option. It does not
try to expand float64's range. Inexact underflow in an individual contribution
is rejected even if the final total could have remained positive. Exact-zero
products from zero density, supplied weight or variance do not signal underflow
and remain supported. No exposed inputs are renormalized, no high-precision
runtime dependency is added, and no floor, clipping, jitter or algorithm change
is introduced. Prefix accumulation and coefficient evaluation retain their
previous arithmetic and validation policies.

Changed production files are only `fishhighz/kernels/weights.py` and
`fishhighz/weights.py`. README.md documents rejection of extreme normalizations.
Added `tests/test_weight_range.py` with ten regressions using independent
500-digit Decimal expressions for prefixes and coefficients. All 628 prior
tests and their tolerances are unchanged.

## Closure evidence

Evidence directory: `.validation/step10-r2-implementation/` (`$RUN` below).
The reviewer's original `range_regression.py` was copied without modification;
its four probes were replayed on source and the new installed wheel.

| Sample order | Weight scale | Repaired source and installed wheel |
| --- | --- | --- |
| Original | 1 | ValueError identifying field and underflow in multiplication |
| Reversed | 1 | Same explicit rejection |
| Original | 1e100 | Accepted; P_pixel = 1e-100 |
| Reversed | 1e100 | Accepted; P_pixel = 1e-100 |

Both accepted controls match the reviewer's Decimal coefficient with zero
relative discrepancy. Their I3 prefixes are `[1e50, 1e100]` and
`[1e100, 1e100]`, respectively. New ordinary tests independently check all
I1/I2/I3 prefixes and A/P_pixel with rtol=5e-13, atol=0, plus normalization
invariance between two accepted scales. They also test exact-zero density,
weight and variance in both row orders with another positive contribution, and
individual loss at the first and third product boundaries. Neither positive
totals nor conventional absolute tolerances can hide the R1 error.

- Focused Step 10 tests: **83 passed in 7.71 s** (`focused.log`).
- Final complete quick runner: **638 passed in 46.95 s**, Ruff lint passed,
  **85 files** formatted (`quick-final.log`). No prior tests were removed.
- Original independent R1 replay: all four probes passed on source and installed
  wheel (`range-source.json/log`, `range-wheel.json/log`). The unscaled cases
  are deliberate, clearly diagnosed rejections; rescaled cases must succeed.
- Installed scalar weight/noise and full-noise PSD/packing checks, amplitude
  Fisher, all four retained examples, Astropy geometry and blocked-Astropy/SciPy
  external weighting/noise subprocess passed (`wheel-probe.json/log`,
  `blocked-optional.log`). All **28 package Python modules** and METADATA/WHEEL
  match source/wheel/installed bytes; origins lie under the fresh environment.
- All **34 legacy source/resource hashes** still match
  (`reference-hash-check.json`). The **81 bounded comparisons** from revision 1
  and its independent review remain historical evidence, not a new execution.
  They were not rerun: this repair only changes operation-level validation;
  reference inputs, formulas and ordinary multiplication order are unchanged.
- Final whitespace, preservation and artifact-hash checks passed.

The unchanged weighted-noise example still gives amplitude information about
7324.26099945 and error 0.0116847141921, with one constrained amplitude and no
added prior. Full outputs are in the installed probe. This remains a synthetic
check, not a physical DESI-2 forecast. The archived handoff retains the bounded
real-input normalization, rejected spline samples, light-speed/response
conventions and limitations without overwriting those records.

## Preservation and exact artifact

Base commit remains `0d69786a06d5d676564a51fad14a7156951c2228`, with the previous
uncommitted/untracked work intact. `before.json` and `git-before.txt` snapshot
89 preexisting maintained files. Only README.md, the two weighting files and
this handoff changed; the other 85 preexisting files remain byte-identical.
The previous handoff is archived verbatim in both `reviews/step-10-r1.md` and
`$RUN/step-10-r1-handoff.md`. No earlier scientific modules, examples,
pyproject.toml, governance/planning documents, sibling resources, previous
validation bundles or review failure logs changed. Final source/preservation
and artifact manifests are `$RUN/after.json`, `preservation.json`,
`git-after.txt`, and `evidence-hashes.json`.

One exact new wheel was built into `$RUN/dist`, and installed into one fresh
`$RUN/wheel-env` with the existing template/cosmology extras. Its SHA-256 is:

`670a143f19fd9a507e6c0ecb07c43ba0c8a5bb988bdb2bf28393525fbc6b04ac`.

No earlier wheel or environment was overwritten. Initial package absence is
recorded in `fresh-environment.log`; installed module origins, metadata and byte
identity prevent a stale same-version install from passing. NumPy remains the
sole unconditional dependency, with no new extra. Python 3.13.15, NumPy 2.5.3,
Astropy 8.0.1 and SciPy 1.18.1 are recorded independently by the installed probe;
`versions.json` also records development tools.

## Commands and repairs during validation

Commands run from the package root unless stated otherwise. All numerical
scripts use `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`; the quick
runner exports these itself. Work stayed on the login node, without Slurm.
`$ROOT` is the absolute component root and `$RUN` the absolute revision-2 evidence
path.

```bash
.venv/bin/python -m pytest tests/test_weight_range.py tests/test_weights.py \
  tests/test_noise.py tests/test_weight_models.py -q
.venv/bin/python "$RUN/range_regression.py" source
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh
.venv/bin/python -m build --wheel --outdir "$RUN/dist"
.venv/bin/python -m venv "$RUN/wheel-env"
"$RUN/wheel-env/bin/python" -I -c \
  "import importlib.util; assert importlib.util.find_spec('fishhighz') is None"
"$RUN/wheel-env/bin/python" -m pip install \
  "$RUN/dist/fishhighz-0.1.0.dev0-py3-none-any.whl[templates,cosmology]"
# From /tmp with absolute paths:
"$RUN/wheel-env/bin/python" -I "$RUN/probe.py" \
  "$RUN/dist/fishhighz-0.1.0.dev0-py3-none-any.whl" "$ROOT" "$RUN"
"$RUN/wheel-env/bin/python" -I "$RUN/range_regression.py" wheel
git diff --check
```

The first quick run passed pytest but reported an import-order lint issue in the
new regression file. Ruff fixed that new import block; the final complete runner
passed. No existing tolerance changed. Sandbox DNS blocked isolated build and
optional dependency installation; approved retries succeeded. The first build
approval request was interrupted before execution, then retried successfully.
Logs are `build.log`, `build-approved.log`, `wheel-install.log` and
`wheel-install-approved.log`. Only one wheel artifact was produced. A stale
process-poll ID after that interruption failed harmlessly; subsequent logs and
completed checks establish the final outcomes. Intermittent NERSC MUNGE socket
messages accompanied successful checks and did not trigger scheduler actions.

All revision-2 acceptance checks are complete; no evidence is knowingly missing.
The real full forecast suite was not requested or run and is not an outstanding
gate. No commit, push, delegation, governance revision, broader scientific
change or Step 11 work occurred. Stop for user and independent review.
