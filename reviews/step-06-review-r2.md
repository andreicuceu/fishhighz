# Step 06 revision 2 review

Date: 2026-09-12.
Outcome: review passed; finding R1 is resolved and no further changes are
required. Ready for progression when the user requests Step 07. No implementation
source or ordinary tests were changed by the reviewer.

## Repair assessment

Reviewed the revision-2 assignment, updated handoff, isolated revision diff,
convergence implementation, new tests, and README changes. The repair compares
exact, order-independent actual point sets during full-study preflight. It
checks every numerical provider/global column at every consecutive refinement,
before any model or Jacobian invocation. Only adjacent structural schedules are
retained; no model-value cache or new numerical allocation framework was added.

Errors identify provider, global parameter, refinement level, requested steps,
and repeated points, with an appropriate remedy. Moving only one point remains
a valid refinement. Single-step evaluation, mixed methods, tied bindings, and
existing numerical kernels are unchanged. The five added regression cases cover
the original exponential reproducer, later repeats, mixed providers in both
orders with zero calls, and a legitimate one-point refinement.

The [revision-1 finding](step-06-review-r1.md) is closed. Distinct stencils alone
still cannot prove derivative accuracy; the documentation retains that limit.

## Independent validation

From the package root:

```bash
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python .validation/step06-review-r2/verify.py
git diff --check
```

- **332 tests passed in 18.03 seconds**, including all retained Step 06 and
  earlier tests. Ruff lint and formatting passed (51 files already formatted).
- Sixteen independent repeated-stencil studies used positive and negative
  fiducials and refinement levels 1–4. Each failed with the expected diagnostic
  before any model call. Single-step derivatives still worked in those cases.
- Six positive/negative power-of-two boundary cases changed only one actual
  stencil point. All were accepted and returned successful convergence for a
  shifted-linear model. This checks that exact comparison avoids rejecting
  genuinely distinct but very close stencils.
- The verifier and results are retained under `.validation/step06-review-r2/`.

## Source, artifact, and preservation evidence

All three inline handoff hashes match the current files. The revised wheel hash
matches `54ea2f32b1613cf97ee1c1fbb9cb649a37245f9a0d643ac0c46d1684a3db86b8`;
all 17 package Python modules match current source. Comparing the revision-1
wheel with current source independently confirms `fishhighz/derivatives.py` is
the only changed package module.

The 53-file pre-dispatch snapshot shows exactly the four reported changes:
derivative implementation, derivative tests, README, and handoff. The archived
revision-1 handoff matches its original hash. These checks preceded the reviewer's
documentation/status updates. Earlier implementation, governance, review, and
reference evidence were preserved. The base commit remains `0d69786`; changes
are uncommitted.

The reviewer reran the extended `probe.py` from the reported evidence directory
`.validation/step06-r2-20260912T170452/`, using its existing isolated
`wheel-env/bin/python -I` from `/tmp`, with all three numerical thread limits
set to one. Arguments were the exact revised wheel, source root, and copied
standalone example. The output is saved as
`.validation/step06-review-r2/wheel-probe.json`.

The probe passed source/wheel/installed module and metadata identity, quiet
import, duplicate-stencil rejection with zero calls, the tied external derivative
and Fisher calculation, and the retained working example. Imports resolved into
the isolated environment's site-packages. Python was 3.13.15 and NumPy 2.5.3.
This reran the reported fresh installation's probe; it does not claim a second
fresh installation by the reviewer.

## Limits and progression

No real lyaforecast quick capture or full suite ran; neither is needed for this
repair. Full execution remains restricted to explicit user requests. Evidence
is synthetic and does not establish external cosmology accuracy, survey
convergence, JIT performance, or a Python-version matrix.

Step 06 revision 2 and its handoff satisfy the assignment. No additional repair
revision is requested. Await the user's instruction before preparing Step 07.
