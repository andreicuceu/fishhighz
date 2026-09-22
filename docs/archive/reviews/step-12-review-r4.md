# Step 12 revision 4: independent performance review

Date: 2026-09-14. Outcome: **passes independent review for the optional compiled
performance path; no required code corrections identified**. User acceptance
is pending. The default NumPy-only matrix target remains unmet and is explicitly
disclosed, not waived or described as achieved. The compiled path meets both
revision-4 performance targets. This is not scientific acceptance of Step 12;
R2 metric/trial binding and R3 forest convergence remain deferred in the
[archived revision-3 assignment](step-12-instructions-r3.md). No Step 13 planning.

Reviewed the [implementation handoff](step-12-performance-r4.md), live production
changes, new tests, benchmark scripts/results, installed wheel evidence and both
figures against the exact revision-4 instructions reviewed (local-only path: `../.validation/step12-review-r4/IMPLEMENTATION_STEP.md`).
New review evidence is exclusively in `.validation/step12-review-r4/`.
No production code, tests, environments, historical evidence or reference inputs
were changed. No real forecast, CAMB calculation, installation, agent dispatch,
Slurm action, commit or push was performed during this review.

## Findings and performance limits

No new numerical or implementation defect was found in the reviewed scope.
The matrix fast paths preserve normalized acceptance thresholds and replay
exceptional/near-threshold batches through the scalar diagnostic. The optional
compiled solve preserves one-cell solve workspace and finite-intermediate checks,
with no fastmath. Missing optional imports and explicit underflow policies retain
NumPy execution. NumPy remains the default and sole unconditional dependency;
Numba/SciPy are an explicit optional extra.

The accuracy path reuses prepared factors while checking the resulting Fisher
against a direct NumPy solve using separately reconstructed Wick covariance.
The new tests perturb either side of that check. Public derivative/forecast
interfaces are unchanged. The three-payload cache is scoped to one study and
keys all refinement controls; independent coupled and failed-refinement studies
retain the baseline arrays, reports, controls and failures exactly.

**Performance limitation, disclosed:** at 2048 synthetic cells, factorization
plus contraction improves 2.77x with NumPy, below the 5x target. Warm optional
compilation achieves 17.56x. The composed evaluation improves 5.77x with NumPy
and 19.21x compiled, both above its 3x target. At 512 cells the corresponding
matrix gains are 2.82x/18.72x and composed gains 5.11x/15.08x. Meeting both goals
therefore requires explicitly selecting `FISHHIGHZ_FISHER_BACKEND=numba` with
`fishhighz[compiled]` installed. No further optimization was silently required
or performed to make the NumPy-only path meet the optional compiled result.
If the user requires both targets without compilation, that remains additional
performance work within Step 12 and needs a revised assignment.

Cold compilation, roughly 4–5 seconds for the observed signatures, is separate
from these warm gains and increases process memory. Real preparation and disk I/O
remain unmeasured. The handoff distinguishes these limitations appropriately.

## Source, archive and installed identity

Identity verification (local-only path: `../.validation/step12-review-r4/identity.json`) confirms:

- All 101 final source/test/example/script/configuration manifest entries match
  live files; all 139 baseline snapshot entries match their saved hashes.
- Exactly seven existing files changed from the preserved baseline: README,
  pyproject, `_information.py`, `covariance.py`, `fisher.py`,
  `validation/accuracy.py` and `validation/schema.py`. The added compiled module,
  performance tests and handoff are identified separately. Existing test files,
  scientific kernels outside the stated changes, and the r2 handoff are intact.
- All 49 packaged modules match live source and both installed environments
  byte-for-byte; installed metadata matches the wheel. Wheel SHA256:
  `95a74305d9a71e92bfd981c9ebd9ae9de6d87ac68bd48ce46861498f4fd32b2d`.
- Archived revision 3 remains byte-identical, SHA256
  `f74db52795b8b27cb3fd9839406e60236a4edde4ea6223df45df03d86f9bc752`.

The installed probes reused the implementer's existing environments; these are
fresh checks of existing installations, not fresh installation evidence. Separate
`-I` probes and examples execute from `/tmp/fishhighz-step12-review-r4`, verifying
site-packages imports, NumPy-only optional-dependency absence and lazy compilation.

## New tests and independent numerical checks

| New review check | Result | Evidence under the review directory |
| --- | --- | --- |
| Ordinary `scripts/check.sh` | 1163 passed, 25 optional-compiler skips in 55.93 s; Ruff lint/format passed | check.log |
| Installed compiled numerical/evidence regressions with `-I` | 358 passed in 30.63 s | installed-compiled.log |
| All new performance regressions repeated with installed compiled wheel from outside checkout | 129 passed | installed-outside-performance.log |
| Independent C/F oracles and preserved-baseline comparison | 32 SPD/scaling/backend combinations; maximum Fisher relative discrepancy 1.56e-15 | independent.py, independent.json, independent.log |
| Independent threshold/range/error comparisons | 28 cases with matching acceptance and first error context | independent.json |
| Compiled layout checks | C-order, Fortran-order, reversed and read-only arrays agree with reference | independent.json |
| Independent reduced convergence studies | Five scenarios retain arrays/reports exactly; one repeated evaluation removed in each | study_independent.py, study-independent.json |
| Installed standalone examples | Four NumPy-only examples and all six optional-dependency examples pass | installed-numpy.json, installed-numba.json, examples-*.log |
| Document/source preservation and whitespace | Passed | preservation.json |

The independent array fixtures cover one/two/three/five fields, two/twelve
parameter columns, exact inactive columns, signed cross terms, and observable
rescalings spanning 1e-120 to 1e120. Wick covariance is reconstructed directly;
Fisher information is checked against direct NumPy solve/einsum and the separately
loaded pre-optimization implementation. Rank/symmetry thresholds, subnormals,
overflow, invalid triangular factors, NaNs and batch-boundary diagnostics were
checked. Existing Step 10 underflow tests and original R1 evidence regressions
also pass. No tolerances were enlarged.

The five independent studies include constant information, unconverged weighting,
coupled refinement controls, magnitude-support failures, and a galaxy-only case.
Baseline/final evaluation counts are 15/14, 16/15, 37/36, 37/36 and 15/14. Their
exact matching outputs provide evidence that the cache removes repeated work
without reducing the scientific trial inventory or hiding failed refinements.

All numerical runs use one OMP/OpenBLAS/MKL thread. The ordinary suite uses the
package development environment. Independent compiled checks use the existing
compiled wheel environment; the fresh benchmark uses the original assessment
interpreter read-only. The ordinary log again includes environment-origin MUNGE
messages after pytest; pytest, Ruff and the runner return success. No scheduler
command was issued.

## Performance evidence and forecast estimate

The performance audit (local-only path: `../.validation/step12-review-r4/performance-audit.json`)
recomputes all reported speedups from saved sample medians and checks the linear
workload-model arithmetic. The benchmark uses matched artificial inputs and
step 0.001 for both implementations. Nested timing rows are correctly identified;
no component is added twice in the runtime model. Both
stage (local-only path: `../.validation/step12-r4-20260914T181308/stages.png`) and
scaling (local-only path: `../.validation/step12-r4-20260914T181308/scaling.png`) figures were inspected
and agree with their plotted data and labels.

One new bounded 512-cell compiled benchmark completed with reference imports
blocked. Its composed median is **26.327 ms**, versus the saved **26.213 ms**,
a 0.43% difference. This independently supports the reported warm timing at that
size. It is a synthetic benchmark, not a real 15x2pt run; no complete forecast
or real convergence study was recomputed.

The six historical report hashes and exact control inventories match: 106
successful distinct trials totaling 1318912 Fourier cells. The preserved
1167.369555-second primary sum and separate 1212.713899-second diagnostic sum
remain historical. The model's central warm primary estimates are **400.1 s
(6.67 min), 2.92x faster for NumPy**, and **290.8 s (4.85 min), 4.01x faster
compiled**. The stated scenario envelopes and the 12.612363-second saved-input
compatibility comparison are correctly qualified as unequal workloads.

These are conditional extrapolations, not independently measured real runtimes.
In particular, the 246.952-second residual is defined as historical time minus
a synthetic proxy; it is not an observed breakdown of real I/O, setup or failed
trials. Its constancy after optimization is an assumption. Two-point scaling
beyond 2048 cells, real forest preparation/conditioning, failed attempts, cold
setup and future convergence repairs can change the next runtime. The handoff
states these limits and makes no unsupported prediction for the 72 diagnostics.

## Disposition

No repair revision or next step was drafted. Governance status now records the
implemented and reviewed optional compiled path, the unresolved NumPy-only matrix
target and pending user acceptance. The exact reviewed assignment and original
handoff are preserved. R2/R3 scientific work remains archived. Stop for the user's
review; do not infer acceptance, restore revision 3 or run the real forecast.
