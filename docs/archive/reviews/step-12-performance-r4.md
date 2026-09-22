# Step 12 revision 4: performance implementation

2026-09-14. Performance changes are ready for independent review, with an
**incomplete NumPy-only matrix-throughput target**. Optional compilation exceeds
both performance targets. This is not scientific acceptance of Step 12. R2
metric/trial binding and R3 forest convergence remain deferred. No real 15x2pt,
other survey forecast, CAMB calculation, reference capture, convergence sweep,
Slurm action, sibling edit, commit or push was performed.

## Source identity and scope

All evidence is in the exclusive directory
step12-r4-20260914T181308 (local-only path: `../.validation/step12-r4-20260914T181308/`).
The baseline is the **dirty working tree**, not just Git HEAD
`0d69786a06d5d676564a51fad14a7156951c2228`. Source, tests, examples, scripts,
reports and package documents were copied before implementation into
baseline (local-only path: `../.validation/step12-r4-20260914T181308/baseline/`).
Baseline hashes (local-only path: `../.validation/step12-r4-20260914T181308/baseline-sha256.json`),
final hashes (local-only path: `../.validation/step12-r4-20260914T181308/final-source-sha256.json`),
Git state (local-only path: `../.validation/step12-r4-20260914T181308/final-git-status.txt`),
and preservation checks (local-only path: `../.validation/step12-r4-20260914T181308/preservation.json`)
identify the exact source. Existing reports, examples, tests, governance documents,
and the original assessment artifacts are retained. The separate r2 handoff
`reviews/step-12.md` is unchanged. The
[archived r3 assignment](step-12-instructions-r3.md) still has SHA256
`f74db52795b8b27cb3fd9839406e60236a4edde4ea6223df45df03d86f9bc752`.

Changes:

- `covariance.py` batches positive-auto field PSD checks in 256-cell blocks.
  Zero-auto, nonfinite, indefinite and near-threshold batches use the original
  scalar diagnostics. Covariance finiteness checking is vectorized.
- `_information.py` and `fisher.py` batch normalization, symmetry checks,
  eigenvalues and Cholesky. Eigenvectors remain in rank/null-space result
  diagnostics. Fast-path eligibility requires a minimum eigenvalue greater than
  twice the existing threshold; everything else replays the original scalar
  acceptance test. This changes eligibility, **not** the acceptance threshold.
  Failed batches replay in original cell order, preserving the first error.
- `kernels/_compiled_fisher.py` supplies an explicitly requested, lazy optional
  contraction. It keeps the reference row-wise forward substitution and matrix
  products, float64, `fastmath=False`, no parallel loops and one cell of solve
  workspace. Nonfinite intermediate/result status triggers scalar replay for
  exact diagnostic ordering. Returned arrays retain owned C-order storage.
- NumPy is still the default and sole required dependency. `fishhighz[compiled]`
  adds Numba and SciPy; SciPy is needed by Numba's matrix-product implementation.
  `FISHHIGHZ_FISHER_BACKEND=numba` selects compilation. Missing Numba, dtype
  conversion, or explicit NumPy underflow warning/error policy uses NumPy.
  Invalid backend names reject. There is no automatic compilation or disk JIT
  cache in the base installation. README documents selection and cold cost.
- `validation/accuracy.py` reuses prepared immutable factors and computes the
  Jacobian once. Private schema assembly consumes those factors. A separate
  direct NumPy solve on independently reconstructed Wick covariance checks both
  combined and individual-spectrum information. Public `schema.assemble`,
  `prepare_bin`, `run_bin` and external-callable signatures remain unchanged.
- Each `study` retains at most three payloads, protects its current selected
  control candidate from eviction, and preserves the existing trial-summary
  inventory. Keys include all six refinement controls, including the derivative
  step. Task, recipe, providers and options are fixed within this invocation;
  nothing is reused across public study/evaluate calls. Failures never enter
  the cache. Existing prepared/sample caches were not redesigned.

Fixed covariance/noise/weights, mean-only derivatives, Wick pair closure/order,
mode normalization, signed cross powers and wiggles, quadrature/cuts, h units,
response ownership, ap/at stencils, rank/convergence tolerances, bin independence,
ties and priors once are unchanged. No jitter, clipping, pseudoinverse, mode
removal, weighting-policy change or relaxed underflow guard was introduced.

## Measurements and profiling

The reference interpreter was used read-only to match the assessment:
`../lyaforecast/.validation/dev-env/bin/python`, Python 3.13.0, NumPy 2.3.5,
SciPy 1.15.3, Numba 0.66.0 and llvmlite 0.48.0. The node was login29; full BLAS
configuration are in each `assessment.json`; NumPy reports OpenBLAS 0.3.30.
Every numerical process sets OMP/OpenBLAS/MKL thread limits to 1. No reference
package is imported: benchmark import blockers reject CAMB, lyaforecast and Vega.

benchmark-matched.py (local-only path: `../.validation/step12-r4-20260914T181308/benchmark-matched.py`)
uses the same artificial five-field signal/noise, template and grids in both
implementations, with 15 spectra, two active parameters and 12 registry columns.
The final composed call invokes `AccuracyRecipe.evaluate` with synthetic
preparation, including geometry/grid-prepared model/noise covariance work.
It excludes real readers and forest sampling. The baseline executes the
preserved prepare/derivative/assemble/run_bin sequence. Both use step 0.001;
the earlier exploratory `final-*` runs used the accuracy default step and are
superseded by `matched-final-*` for equal-setting composed comparisons.
No smaller grids, fewer study trials or altered tolerances count as speedups.

Each benchmark has a 55-second alarm, 512 or 2048 cells, five timing samples
per main stage after warmup, three per batch-size comparison, and nine for one
cell. No 4096-cell or real-data benchmark was needed. Profiling and tracemalloc
are separate from unprofiled timing samples. Compilation precedes warm timings.
Each output directory must be new; do not overwrite the recorded artifacts.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# From the FishHighz root; choose a new output directory for every execution.
PYTHONPATH=. BENCH_CELLS=2048 BENCH_OUT=/new/bounded/output \
  FISHHIGHZ_FISHER_BACKEND=numba \
  ../lyaforecast/.validation/dev-env/bin/python \
  .validation/step12-r4-20260914T181308/benchmark-matched.py
# For the preserved baseline use its absolute directory as PYTHONPATH,
# BENCH_IMPLEMENTATION=baseline, and FISHHIGHZ_FISHER_BACKEND=numpy.
```

The initial matrix-only phase reduced the same 2048-cell composed fixture from
1.456 s to 0.470 s. Retaining the old duplicate production work but adding warm
compilation reduced it to 0.115 s. These staged runs are preserved in `baseline`,
`matrix` and `compiled`, with profiles. The final equal-setting stage table below
reports median [minimum, maximum] milliseconds. Nested rows must not be added:
`prepare` includes covariance/factorization; public `assemble` includes its own
factorization/contraction; `composed` includes preparation, J, private assembly
and its direct oracle. Validation and NPZ compression are outside composed time.

| Operation (2048 cells) | Preserved NumPy ms | Final NumPy ms | Final compiled warm ms |
|---|---:|---:|---:|
| field_psd | 115.787 [115.124, 142.384] | 4.671 [4.488, 4.877] | 4.378 [4.357, 4.580] |
| covariance | 123.528 [122.912, 134.514] | 7.624 [7.328, 7.773] | 6.894 [6.800, 6.917] |
| factor | 404.352 [400.167, 418.741] | 28.586 [28.438, 28.639] | 29.526 [27.173, 31.810] |
| fisher2 | 208.005 [199.445, 222.508] | 192.868 [183.851, 199.852] | 5.337 [5.220, 5.383] |
| fisher12 | 231.117 [207.331, 243.990] | 188.753 [186.028, 192.361] | 7.955 [7.660, 8.030] |
| prepare | 647.119 [639.213, 701.032] | 41.379 [40.520, 41.586] | 42.966 [41.500, 43.963] |
| derivatives | 8.810 [8.754, 8.840] | 9.849 [9.157, 10.302] | 9.996 [9.936, 10.115] |
| assemble | 603.776 [600.599, 607.595] | 235.383 [223.364, 236.941] | 46.035 [45.641, 50.054] |
| independent | 7.820 [7.804, 7.878] | 8.136 [7.899, 8.469] | 7.880 [7.870, 7.920] |
| validation | 26.730 [26.564, 26.809] | 28.561 [27.685, 28.698] | 27.241 [26.700, 27.291] |
| serialization_memory | 87.488 [85.941, 123.810] | 89.489 [88.176, 90.309] | 85.821 [85.598, 86.366] |
| composed | 1462.303 [1449.073, 1504.054] | 253.359 [252.949, 256.694] | 76.115 [76.034, 79.492] |
| three_spectra | 427.476 [426.289, 561.897] | 102.623 [99.335, 111.005] | 4.728 [4.644, 4.904] |
| one_cell | 0.356 [0.347, 0.380] | 0.334 [0.318, 0.423] | 0.226 [0.223, 0.252] |


Stage figure (local-only path: `../.validation/step12-r4-20260914T181308/stages.png`)
(PDF (local-only path: `../.validation/step12-r4-20260914T181308/stages.pdf`)) and
scaling figure (local-only path: `../.validation/step12-r4-20260914T181308/scaling.png`)
(PDF (local-only path: `../.validation/step12-r4-20260914T181308/scaling.pdf`)) use measured synthetic
data; both PNGs were visually checked. Batch sizes 1, 32, 128, 256 and 512 and a
three-spectrum subset are retained in the raw JSON. One-cell warmed usage does
not materially regress. Cold Numba cannot be treated as a small-call speedup.

| Cells | NumPy composed speedup | Compiled composed speedup | NumPy matrix speedup | Compiled matrix speedup |
|---|---:|---:|---:|---:|
| 512 | 5.11× | 15.08× | 2.82× | 18.72× |
| 2048 | 5.77× | 19.21× | 2.77× | 17.56× |

The composed target is ≥3×; factorization plus contraction is ≥5×.
NumPy misses the latter. Compilation meets both at both measured sizes.

| Implementation, 2048 cells | Traced peak bytes | Retained allocation count after call | Process peak RSS KiB | Profile calls |
|---|---:|---:|---:|---:|
| baseline | 23719178 | 132 | 108544 | 1401348 |
| final-numpy | 23966303 | 122 | 92160 | 124898 |
| final-numba | 23966592 | 127 | 206848 | 18364 |

The traced peak is approximately 24 MB, including the composed payload.
Retained allocation counts are snapshot counts, not total allocation events.
Process RSS includes imports/JIT and is not isolated kernel memory. Matrix
scratch is bounded to 256 cells; solve scratch is one selected-spectrum by
parameter array (15×12×8 = 1440 bytes for 12 columns), plus parameter contraction
buffers. The optional compiler increases process memory despite similar warm
array memory. Raw samples include shared-node noise; no profiler time enters
the throughput medians.

The second profile is
final compiled (local-only path: `../.validation/step12-r4-20260914T181308/matched-final-numba-2048/profile.txt`),
with final NumPy (local-only path: `../.validation/step12-r4-20260914T181308/matched-final-numpy-2048/profile.txt`)
and baseline (local-only path: `../.validation/step12-r4-20260914T181308/baseline-2048/profile.txt`)
for comparison. The baseline has roughly 1.4 million profiled calls; the optimized
compiled evaluation has about 18 thousand. Two derivative evaluations become one;
11 model calls become 6, two factorizations become one, and two production
contractions become one production contraction plus a direct solve oracle.

The NumPy contraction still spends about 0.18 s in per-cell checks and Python
row solves. The all-cell NumPy prototype was not adopted because it violates
the public one-cell solve-workspace contract.

After compilation, normalized eigenvalue/factor work remains the largest
numerical stage. Model/derivatives and independent evidence checks now contribute
material fractions. Their extra complexity is not justified by the remaining
short synthetic costs. Template/forest-weight JIT, inactive-column pruning,
geometry reconstruction, and broad immutable-copy changes were left alone.
The full writer still validates twice: removing that repetition would require
preserving its distinct malformed/scientifically-failed error ordering. The
measured checker is about 27 ms, while NPZ compression alone is about 86 ms;
no serialization/error-predicate redesign was retained. Disk I/O is unmeasured.

First compiled assembly is about 3.8 s versus approximately 0.046 s warm. The
separate cold probe (local-only path: `../.validation/step12-r4-20260914T181308/cold/cold.json`)
records an additional readonly/non-C-contiguous signature: first composed call
0.985 s, second 0.103 s in that process. Compilation is signature dependent;
allow roughly 4–5 s for these two observed signatures in a fresh process, plus
unmeasured real setup. At 2048 cells, the contraction-only cold crossover is
roughly 20–25 calls; at 512 cells it is roughly 80–90. These are estimates from
cold excess divided by warm savings, not universal thresholds. No JIT cost is
charged to default NumPy.

## Conditional model for the next historical workload

The saved reference is **1167.369555 s** for six primary accuracy records,
**12.612363 s** for six compatibility records and, separately,
**1212.713899 s** for 72 requested accuracy diagnostics. These are unequal
workloads. Compatibility begins from saved powers/Jacobians. Shared background
setup is outside the record sums, and no single fixed-setting real accuracy
forecast was separately timed.

Historical controls (local-only path: `../.validation/step12-r4-20260914T181308/historical-controls.json`)
were read from six small JSON reports, retaining their hashes. No saved-case
recomputation or real recipe construction was performed.

| Bin | Successful distinct trials | Sum of trial cells | Recorded unresolved entries |
|---|---:|---:|---:|
| 0 | 17 | 208896 | 4 |
| 1 | 17 | 208896 | 4 |
| 2 | 18 | 225280 | 2 |
| 3 | 18 | 225280 | 2 |
| 4 | 18 | 225280 | 2 |
| 5 | 18 | 225280 | 2 |
| Total | 106 | 1318912 | 16 |

The failure entries include repeated magnitude-level failures; they are not a
complete timed execution trace. Final grids remain 16384 cells, with original
k/mu/z/magnitude/weight/step controls and all stopping/refinement rules. The model
does not subtract trials. It also gives no extra credit for avoiding a repeated
final payload, because historical repeated-evaluation timing is not isolated.

The central model fits `t(N)=a+b*N` through the two measured composed medians,
then sums it over the **106 recorded successful controls**. It leaves the
historical time not represented by that synthetic extrapolation unchanged.
This is a workload proxy, not a measured decomposition of the real run. It
extrapolates beyond 2048 cells and assumes the artificial model's scaling is
informative; real preparation caches, forest inputs and conditioning can differ.

| Model term | Preserved NumPy | Final NumPy | Final compiled |
|---|---:|---:|---:|
| Per-evaluation intercept [ms] | 39.601 | 18.685 | 9.579 |
| Slope [µs/cell] | 694.679 | 114.587 | 32.488 |
| 106-trial synthetic proxy [s] | 920.418 | 153.111 | 43.865 |

The unchanged residual is 246.952 s = 1167.369555 s minus
the 920.418 s baseline proxy. It contains all unassigned real costs and
synthetic-to-real mismatch, not independently measured setup or I/O. Keeping
it unchanged gives the central estimates below. No nested stage timing is
added to the proxy, and no reduction in the historical inventory is assumed.

The following decomposition states what is and is not timed; it prevents adding
nested matrix timings to the composed proxy a second time.

| Contribution | Synthetic evidence | Treatment in real-runtime model |
|---|---|---|
| Shared background/template/input setup before records | Excluded | Unknown additive cold setup; outside 1167.37 s |
| Per-trial preparation, model and J | Included in composed; separate prepare/J rows above | Included once in successful-trial proxy; real reader/weight difference unmeasured |
| Field PSD, selected C, normalization and Cholesky | Included within preparation | Included once; not separately added to proxy |
| Production F and summaries | Included within composed assembly | Included once; factor reuse removes duplication |
| Separate direct C/F oracle | Included in final composed | Included once; not inferred from cached production F |
| Writer mathematical validation | About 27 ms per 2048-cell synthetic payload | Historical real total unknown; retained within unchanged residual |
| Compression and filesystem output | About 86 ms in-memory NPZ; disk unmeasured | Retained within unchanged residual |
| Failed trials and repeated successful evaluations | Original controls/errors retained; no real timings per attempt | Retained within unchanged residual; no invented subtraction |
| Optional JIT | Measured separately by signature | Roughly 4–5 s additive cold cost for observed signatures |

Conservative and optimistic scenarios retain the assessment's hypothetical
accelerated fractions 0.70 and 0.90, respectively, but use the newly measured
2048-cell composed acceleration in
`T=1167.369555*((1-f)+f/s)`. The central estimate instead uses the control-inventory
proxy above. Fractions are assumptions, not fitted real fractions or confidence
limits. Full arithmetic is in `summary.json` and `workload-model.json`.

| Backend | Scenario | Primary records [s] | Minutes | Speedup |
|---|---|---:|---:|---:|
| numpy | conservative | 491.8 | 8.20 | 2.37× |
| numpy | central | 400.1 | 6.67 | 2.92× |
| numpy | optimistic | 298.8 | 4.98 | 3.91× |
| numba | conservative | 392.7 | 6.55 | 2.97× |
| numba | central | 290.8 | 4.85 | 4.01× |
| numba | optimistic | 171.4 | 2.86 | 6.81× |

These are warm-record predictions. Add the separate cold setup/JIT terms.
The compiled scenario envelope is approximately 3.0–6.8×, with the
inventory-based central estimate approximately 4×; the NumPy estimate is
approximately 2.4–3.9×. These are conditional scenarios, not error bars.

This replaces the preliminary 2.7–7.7× scenario range with evidence from guarded
production code. Unknown real input preparation, I/O and failed-trial costs,
shared-node noise, extrapolation beyond measured cell counts, and future changed
convergence choices can put the next runtime outside this envelope. A future
convergence repair may change the workload and invalidate this prediction.
**The next real runtime has not been measured.**

Even the compiled scenarios remain well above the historical 12.61 s
compatibility sum. Compatibility also benefits from common covariance/Fisher
optimizations: the equal-array and public-assembly rows measure that benefit
synthetically. Its saved-input advantage and absence of 106 refinement trials
remain. No credible fresh compatibility total or diagnostic total follows from
these data. The separate 72 diagnostics retain their historical 60 completed and
12 failed requests; their real preparation/failure mixture was not timed here,
so no separate prediction is claimed for their 1212.71 s sum.

## Numerical checks, installation and remaining limits

The new self-contained regressions cover seeded SPD five-field/all-pair and
three-pair selections, signed cross terms, 2/12 columns, exact inactive zeros,
positive observable-unit rescaling and direct NumPy C/F/constrained-error
oracles. Existing analytic one/two-field, bin-combination, ties and priors-once
tests are retained. Relative norm limits remain 5e-12 or stricter; the measured
2048-cell Fisher residual against the separate direct solve is about 2.1e-15.
No tolerance was enlarged.

Batch-boundary tests place errors at cells 0, 255, 256, 257 and 511, followed by
another error at 512. They compare exact scalar error text for negative/zero
rows, asymmetry, indefinite/near-rank-threshold matrices and subnormals. Both
backends retain structural lower-triangle, nonfinite, solve-overflow and
accumulation-overflow rejection. Explicit NumPy underflow warning/error policies
are compared directly. `test_weight_range.py` and its mixed-product/exact-zero/
representable-subnormal controls are unchanged and pass.

The reduced study probe (local-only path: `../.validation/step12-r4-20260914T181308/study_probe.py`)
compares the preserved implementation with the final source. Its baseline/final
JSON records contain identical controls, matrices, failure entries, metrics and
pass state; actual evaluations fall from 15 to 14.

Cache tests count evaluations over all supported controls, check independent
trial matrices, verify hits/evictions and the three-payload bound, preserve failed
trials and ordering, and prove no reuse across separate study calls. A new test
perturbs either production information or the direct oracle and requires
rejection. Existing R1 evidence mutation tests pass; the known R2 binding defect
is **deferred**, not newly passed. Public NumPy-only execution is tested with
optional imports blocked, including an explicit unavailable-Numba request.

Actual final checks:

| Check | Result | Evidence |
|---|---|---|
| `PATH="$PWD/.venv/bin:$PATH" scripts/check.sh` | 1163 passed, 25 optional-compiler skips; Ruff lint/format passed | `check-complete.log` |
| Installed compiled wheel, explicit NumPy/Numba regressions plus covariance/Fisher, weights, evidence, results, parameters and forecast tests | 358 passed | `wheel-compiled-complete.log` |
| Installed NumPy-only wheel, all new performance regressions | 104 passed, 25 optional-compiler skips | `wheel-numpy-complete.log` |
| Earlier broader NumPy-only installed numerical/evidence checks | 228 passed, 26 skips, 1 SciPy-specific test deselected | `wheel-numpy-final.log`; final changed paths rerun above |
| Installed standalone examples | Four in NumPy-only environment; all six in compiled environment | `wheel-*-examples-complete.log` and verification JSON |
| Final source-to-wheel and wheel-to-installed byte comparisons, metadata, `-I` import location | Passed in both environments | verification JSON and final SHA256 manifest |
| `git diff --check` | Passed | final local check |

All 1059 original tests are retained; 129 new parametrized cases give 1188 total
cases in the ordinary suite. Skips concern the optional compiler in the base
development environment; that path is explicitly tested in the compiled wheel.
The compiled installed command, from the external temporary test directory, was:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
FISHHIGHZ_FISHER_BACKEND=numba /absolute/path/to/compiled-env/bin/python -I \
  -m pytest tests/test_fisher.py tests/test_covariance.py \
  tests/test_weight_range.py tests/test_step12_revision2.py \
  tests/test_step12_performance.py tests/test_results.py \
  tests/test_parameters.py tests/test_forecast.py -q
```

The exact environment path is the artifact directory's `compiled-env`; all test
fixtures/examples/scripts were copied to the `/tmp` directory stated below.
Build command: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
.venv/bin/python -m build --wheel --outdir
.validation/step12-r4-20260914T181308/dist`. Local wheel installations use
`pip install --no-index --no-deps --force-reinstall` against that exact wheel.


The exact final wheel is
fishhighz-0.1.0.dev0-py3-none-any.whl (local-only path: `../.validation/step12-r4-20260914T181308/dist/fishhighz-0.1.0.dev0-py3-none-any.whl`).
Its digest is in
final-wheel-sha256.txt (local-only path: `../.validation/step12-r4-20260914T181308/final-wheel-sha256.txt`).
Every packaged module was compared byte-for-byte with current source and with
both installed copies. Metadata/import paths and example inventories are in
`wheel-numpy-verification.json` and `wheel-numba-verification.json`. Tests and
examples execute under `/tmp/fishhighz-step12-r4-20260914T181308` with `python -I`.
The NumPy environment contains NumPy and test dependencies, without Numba/SciPy;
the separate compiled environment pins NumPy 2.3.5, SciPy 1.15.3, Numba 0.66.0
and llvmlite 0.48.0. No sibling/shared environment was modified.

Initial validation issues are preserved in their logs: a no-isolation build
lacked setuptools; sandbox DNS blocked isolated build/install dependencies,
then approved isolated downloads succeeded. The first compiled test exposed
Numba-owned buffer metadata, corrected by an owned host result copy. An initial
NumPy-only installed test selection included a SciPy-specific template test and
omitted a standalone script fixture; the fixture was copied, and the
SciPy-specific test passed in the compiled environment. An initial payload-cache
regression found eviction of the selected candidate; the protected three-entry
cache now passes. Earlier exploratory wheels remain in `initial-dist` and
`cache-repair-dist`; only `dist` is the final wheel.

The ordinary check log contains environment-origin MUNGE socket diagnostics
following pytest completion; pytest and Ruff exit successfully. No scheduler
command was issued. Performance completion remains conditional on accepting the
optional compiled path: the NumPy-only 5× matrix target is explicitly incomplete.
The 3× composed target passes for both backends. R2/R3 and scientific acceptance
remain with the user; the archived scientific assignment was not resumed.
