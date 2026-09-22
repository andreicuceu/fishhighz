# Step 07 revision 1 review

Date: 2026-09-12.
Outcome: review passed; no required changes. Ready for progression when the user
requests Step 08. No implementation source or ordinary tests were changed by the
reviewer.

## Assessment

Reviewed the revision-1 assignment, handoff, template preparation/evaluation
modules, tests, README, dependency declarations, and saved acceptance evidence.
The implementation satisfies the bounded template-input assignment. It keeps
Astropy FITS input and SciPy coefficient preparation behind lazy imports, with
NumPy-only repeated evaluation and a NumPy-only base installation.

The loader selects the unique named PK binary table and hashes the same byte
snapshot it parses. It preserves signed PKSB and PK-PKSB samples, validates
metadata and the documented units, and converts h conventions explicitly.
Not-a-knot splines operate on linear amplitudes in log k; the derivative kernel
includes the 1/k factor. Domain checks precede logarithms and reject one-ULP
excursions. Explicit growth scales both components and derivatives once without
changing the prepared record. Numerical functions retain the intended array-only
boundary for later optimization.

The handoff accurately distinguishes interpolation tests from forecast accuracy.
Its initial oscillatory sampling fixture was refined to meet the derivative
target; the final tests retain a coarse-grid failure check and analytic oracles.
There is no outstanding repair request or reason to revise the implementation
requirements.

## Independent validation

From the package root:

```bash
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python .validation/step07-review-r1/verify.py \
  > .validation/step07-review-r1/verification.log 2>&1
git diff --check
```

- **408 quick tests passed in 31.97 seconds**, including 76 template tests, with
  no skips. Ruff lint and formatting passed (56 files already formatted).
- Twenty-four additional deterministic analytic cases used cubic polynomials
  in source log k, 4/5/9/17 irregular knots, three source/fiducial h combinations,
  and reference/nonreference redshift amplitudes. Powers and dP/dk agreed with
  independently evaluated polynomials; maximum absolute errors were
  `4.55e-13` and `4.37e-11`, respectively. The comparisons used
  `rtol=2e-11, atol=3e-10` to allow spline arithmetic and the small-k derivative
  factor. Unsorted/repeated queries, exact endpoints, interior knot neighbors,
  read-only inputs, owned outputs, exact sliced equality, and just-outside-domain
  rejection also passed.
- Both specified Vega FITS files were loaded read-only. Each retained 814 nodes,
  121 negative smooth samples, the specified redshift/H0, and the exact hash in
  the assignment. Full-power knot reconstruction had maximum absolute error
  `1.14e-13`; all 813 midpoint values and derivatives were finite. A second h_fid
  convention also reproduced the expected cubic power conversion. File hashes
  were unchanged after reading. No Vega/CAMB code or reference forecast ran.

The verifier and results are retained under `.validation/step07-review-r1/`.
Numerical commands used the three single-thread limits. NERSC MUNGE socket
messages appeared after successful subprocesses; the checks exited successfully.

## Source, artifact, and preservation evidence

All five inline handoff hashes match current files. The reported wheel hash is
`1af6778f6b68439001e3665f1b09eb69d0237b80db7fbbcb0000dc4232b72d57`, and all
19 package Python modules match the source. The 54-file pre-dispatch snapshot
independently confirms README.md and pyproject.toml were the only changed
preexisting maintained files. These checks preceded reviewer documentation/status
updates. Earlier scientific code, tests, reports, and governance were preserved.

The reviewer reran `templates-probe.py` from
`.validation/step07-r1-20260912T194343/` with its existing isolated
`wheel-env/bin/python -I`, from `/tmp`, with all three thread limits set to one.
Arguments were the exact reported wheel, absolute package source root, and the
new `.validation/step07-review-r1/` output directory. This preserves the original
generated fixture and evidence. `templates-probe.log` and `templates-probe.json`
record the successful rerun.

The installed probe verifies source/wheel/installed module and metadata identity,
optional dependency declarations, generated-FITS loading, unit/growth scaling,
domain errors, slicing, and evaluation with preparation-library calls disabled.
Imports resolve under the isolated environment's site-packages. This reruns the
reported fresh installation's probe; it does not claim a second fresh install.

The original base-only probe and fresh-environment evidence were inspected.
That probe verified actual absence of Astropy/SciPy before extras installation,
plus the external derivative/Fisher path. Its exact source/artifact identity
remains verified. The reviewer did not recreate that dependency-free environment;
the current installed environment contains extras. The rerun quick tests also
cover lazy imports and missing-extra diagnostics.

## Limits and progression

No real lyaforecast quick capture or full suite was required or run. Full mode
remains restricted to explicit user requests. These checks establish template
format, interpolation, units, amplitude, and packaging behavior; they do not
establish physical forecast accuracy, JIT performance, or a Python-version matrix.

Step 07 revision 1 and its handoff satisfy the assignment. Status documents now
record this review; no repair revision is requested. Await the user's instruction
before preparing Step 08. No commit or push was performed.
