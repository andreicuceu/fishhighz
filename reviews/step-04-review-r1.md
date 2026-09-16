# Step 04 revision 1 review

Outcome: review passed; no findings requiring changes. Ready to proceed when the
user requests the next step. Step 05 has not been drafted or implemented.

## Scope and provenance

Reviewed revision 1 of `IMPLEMENTATION_STEP.md`, the implementation handoff,
`fishhighz/covariance.py`, both kernel-package files, `tests/test_covariance.py`,
and the README covariance example against the accepted array/grid contracts.
No implementation code or tests were changed by the reviewer.

All five inline handoff hashes match the reviewed files. The exact wheel SHA-256
matches `0453ee2432d9fda20e99ca4a23a890847ede10738ed0696d765ba7848daf53e5`,
and all ten Python modules in that wheel match current source bytes. Independently
checked the 28 files in the implementer's pre-dispatch snapshot: README.md was
the only changed preexisting file, as reported. This check preceded the reviewer's
documentation/status updates. The accepted Step 03 source and tests are preserved.
The base commit remains `0d69786`; implementation is uncommitted.

## Independent validation

From the FishHighz root:

```bash
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python .validation/step04-review-r1/verify.py
git diff --check
```

- **186 tests passed in 28.13 seconds**, including all 40 covariance tests and
  the 146 retained tests. Ruff lint and formatting passed (34 files formatted).
- Sixty-three deterministic cases across one to six fields and full/deficient
  ranks matched an independent quadratic-estimator expression,
  `Cov_ab = 2 Tr(E_a T E_b T)/N`. This uses symmetric estimator matrices rather
  than the production pair lookup tables or its two-product implementation.
- Arbitrary selected subsets, field permutations, read-only inputs, whole versus
  sliced node batches, and composition ownership passed in those cases.
- All handoff hashes, source/wheel module bytes, and preservation checks passed.
  The independent verifier and JSON results are retained under
  `.validation/step04-review-r1/`.

Also reran the reported `.validation/step04/probe.py` using the reported
`.validation/step04/wheel-env/bin/python -I` from `/tmp`, with the exact wheel
and saved README example as arguments and all three thread limits set to one.
The probe confirmed quiet package import, all ten installed module identities,
the README example, analytic/kernel agreement, and singular/invalid-power
behavior. Module origins were under the isolated environment's site-packages.
Python was 3.13.15 and NumPy 2.5.3.

The implementation handoff and installation log document the fresh wheel install.
This review reran that existing installed environment's probe and independently
verified artifact/source identity; it does not claim a second fresh installation.

## Assessment

- The Wick covariance, selected-pair ordering, required dependencies, and
  positive-mu mode normalization follow the assignment. Signed cross powers and
  covariance entries are retained, and observed signal/noise are added once.
- Full active-field PSD validation catches multi-field indefiniteness and uses
  normalized correlations to account for disparate amplitudes. Its explicit
  tolerance and zero-auto policy match the handoff. Valid singular inputs pass
  without clipping, jitter, removal, or a false invertibility claim.
- Boundary validation and diagnostics are separate from the array-only kernel.
  The kernel fills and mirrors one triangle with one node-vector workspace,
  avoiding covariance-sized temporary expansions. Node slices work independently.
- Tests cover analytic limits, nontrivial cross noise, invalid inputs, overflow,
  array ownership, and singular cases. The README accurately limits its example
  to covariance from supplied arrays. No Fisher, survey, or model work was added.

## Limits and progression

No real lyaforecast quick capture or full suite was run; neither is required for
this synthetic covariance step. Full execution remains restricted to explicit
user requests. This review validates the NumPy reference on Python 3.13, not a
Python-version matrix, survey convergence, or compiled performance. Intermediate
arithmetic overflow remains an explicit documented error, and selected-covariance
solvability remains the responsibility of the later Fisher step.

The implementation and report satisfy Step 04 revision 1. No repair revision is
requested. Await the user's instruction to prepare Step 05.
