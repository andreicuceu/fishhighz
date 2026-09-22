# Step 06 revision 1 review

Date: 2026-09-12.
Outcome: one required correction; do not progress to Step 07 yet. Revision 2 of
the current assignment describes the bounded repair. No implementation code or
ordinary tests were changed by the reviewer.

## Finding R1 — repeated effective stencils falsely pass convergence (P2)

In derivatives.py (local-only path: `../fishhighz/derivatives.py:296`), refinement preflight checks
each schedule independently but does not compare actual float64 points between
levels. The later comparison at line 320 (local-only path: `../fishhighz/derivatives.py:320`)
therefore treats repeated evaluations of the same stencil as convergence.

Reproducer using one numerical parameter and a plain provider returning exp(x):

```python
x = 1.5
h = 1.2 * np.spacing(x)
# Parameter("x", x, "target", step=h)
# check_convergence(..., atol=0, rtol=0)
```

The requested steps are `2.6645352591003756e-16` and
`1.3322676295501878e-16`, but both schedules use the same two points:
`(1.4999999999999998, 1.5000000000000002)`. Both finite differences return 4.0,
while the analytic derivative is `exp(1.5) = 4.4816890703380645`. The utility
returns `passed=True` with zero absolute/relative change. This is not an
independent refinement, and can falsely reassure a user about a rounding-dominated
derivative and its Fisher information. General warnings that convergence is not
proof of accuracy do not address this directly detectable failure.

During preflight, compare consecutive actual stencils for every numerical
provider/global column and every requested refinement. Reject unchanged effective
evaluation points with provider, parameter, refinement, requested-step, and actual
point context, before any external call. Use exact float64 point identity, not
`allclose`, which would also reject valid small refinements. Preserve ordinary
single-step differentiation, actual-offset coefficients, and family-change
reporting. Do not silently enlarge/shrink steps or return a successful assessment.

Add regressions for first-level and later-level duplicate stencils, including a
study with other valid columns/providers, and verify zero calls on preflight
failure. Retain legitimate zero-derivative convergence on distinct stencils.

## Other assessment and validation

The remaining reviewed behavior follows the assignment: complete exclusive pair
routing, original field indices, registry-aware bindings, analytic/mixed global
columns, tied global perturbations, bounds-aware second-order formulas, immediate
copies of provider output buffers, independent P1D invocation, and fixed-covariance
composition. Numerical kernels are separate from callable dispatch and metadata;
there is no persistent model-value cache or full stencil-history allocation.

Commands from the package root:

```bash
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python .validation/step06-review-r1/verify.py
git diff --check
```

- **327 quick tests passed in 17.28 seconds**, with Ruff lint/format passing
  (50 files already formatted). These tests do not cover finding R1.
- Nine independent cases compared polynomial derivatives across analytic,
  numerical, and mixed strategies at interior and boundary parameter states.
  Three reordered providers covered six required pairs, tied local slots,
  mutable output buffers, signed coefficients, and independent chain-rule
  expectations. Powers, Jacobians, and call counts passed.
- The independent verifier also records R1's actual points, derivatives, and
  false success in `.validation/step06-review-r1/results.json`. This reproduction
  records current broken behavior; it is not a passing regression for the repair.
- All seven handoff hashes match. The reported wheel hash matches
  `9a8b6ef69d330267cefc8b47a5725789fcfdb259adf8ebf4ffd32b85720b9c5e`,
  and all 17 wheel Python modules match source bytes. Of 45 files in the
  implementer's pre-dispatch snapshot, only README.md changed. This preservation
  check preceded the reviewer's governance/report edits.

The reviewer reran the reported `probe.py` with the existing isolated
`wheel-env/bin/python -I` from `/tmp`, all three thread limits set to one,
passing the exact wheel, source root, and copied example. Evidence is under
`.validation/step06-20260912T163621/`. Quiet import, source/wheel/installed
module identity, metadata identity, tied derivative/Fisher probe, and standalone
example passed. Python was 3.13.15 and NumPy 2.5.3. This verifies the reported
installation; it does not claim a second fresh install by the reviewer.

## Limits and next action

The handoff accurately describes its existing checks, but its completion claim
needs revision after R1 is fixed and tested. No real lyaforecast quick capture or
full suite ran; neither is required here. Full execution remains restricted to
explicit user requests. Synthetic checks do not establish external cosmology
accuracy, survey convergence, or JIT performance.

Implement the focused Step 06 revision-2 repair, rerun required quick checks and
one fresh installed-wheel probe, and provide an updated handoff for review.
No commit, push, or Step 07 work is authorized by this review.
