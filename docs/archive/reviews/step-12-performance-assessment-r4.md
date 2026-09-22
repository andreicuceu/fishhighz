# Step 12 revision 4: preliminary performance assessment

Date: 2026-09-14. Planning assessment only; no production optimization or real
forecast was executed. The user has suspended the convergence/evidence repair
assignment and requested performance optimization first. The exact previous
instructions are [archived revision 3](step-12-instructions-r3.md), SHA256
`f74db52795b8b27cb3fd9839406e60236a4edde4ea6223df45df03d86f9bc752`.
[Revision 4](../notes/IMPLEMENTATION_STEP.md) is proposed for user approval/dispatch.
R2 evidence binding and R3 forest convergence remain unresolved; see the
[unchanged review](step-12-review-r2.md). Step 12 is not scientifically accepted.

## Historical timing and unequal work

Read-only extraction from
the r2 manifest (local-only path: `../.validation/step12-r2-implementation/profiles-streamed/manifest.json`)
is saved with its hash in
historical-timings.json (local-only path: `../.validation/step12-performance-plan-r4/historical-timings.json`).
These are old record timings, not measurements made during this assessment.

| 15x2pt workload | Historical time | Included work |
| --- | ---: | --- |
| Six compatibility primary records | 12.612363 s | Saved upstream powers/Jacobians; FishHighz covariance, Fisher, checking/output |
| Six accuracy primary records | 1167.369555 s (19 min 27.37 s) | 17,17,18,18,18,18 distinct successful refinement trials, repeated evaluations, checking/output |
| 72 accuracy diagnostic requests | 1212.713899 s (20 min 12.71 s) | Separate diagnostic evaluations, including failed requests |

Accuracy primary bins took 177.36, 175.31, 201.99, 204.74, 203.84 and 204.13 s.
Its final grid has 16384 cells (128 k intervals, order 4, 32 mu nodes), versus
5000 for compatibility. Each redshift bin has two active ap/at parameters in a
12-column global registry. Shared background setup is outside these record
sums; a fresh process may incur additional setup. The ratio 92.56 compares
unequal workloads. The time for one fixed-setting real accuracy evaluation is
unknown. Reaching compatibility's saved-input time for the entire convergence
study is an aspiration, not a justified prediction from this evidence.

## Short targeted measurements

Used the existing read-only reference interpreter to match the historical
NumPy/SciPy environment: Python 3.13.0, NumPy 2.3.5, SciPy 1.15.3; available Numba
0.66.0 and llvmlite 0.48.0. Host login29, OMP/OpenBLAS/MKL thread limits all 1.
No environment was installed or changed. The script had a 55-second wall-time
alarm and completed successfully. Peak process RSS was 137216 KiB, including
imports, profiling, arrays and JIT; this is not isolated kernel memory.

The fixture contains five fields, 15 spectra, 2048 Fourier cells, two active
parameters in a 12-column registry, an artificial smooth-plus-sinusoidal wiggle
template, constant background functions and synthetic full noise. Field names
and bin identities use the 15x2pt selection only. No real survey/model inputs
are evaluated; an import blocker rejects CAMB, lyaforecast and Vega. A separate
256-sample/12-update forest-weight fixture uses artificial positive inputs.

Reproduce from the package root, writing only to a new artifact directory if
repeating later (the recorded script currently writes beside itself):

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=. \
  ../lyaforecast/.validation/dev-env/bin/python \
  .validation/step12-performance-plan-r4/assess.py
```

Script (local-only path: `../.validation/step12-performance-plan-r4/assess.py`),
raw timings and numerical checks (local-only path: `../.validation/step12-performance-plan-r4/assessment.json`),
log (local-only path: `../.validation/step12-performance-plan-r4/assessment.log`), and
cProfile table (local-only path: `../.validation/step12-performance-plan-r4/profile.txt`)
are new synthetic evidence. Times below are unprofiled medians after warmup,
usually three samples; composed evaluation uses two and weighting/JIT use five.
Nested rows overlap and must not be summed as independent costs.

| Operation, 2048 cells unless stated | Median time |
| --- | ---: |
| Field-power PSD validation | 121.97 ms |
| Gaussian covariance including field PSD | 130.98 ms |
| Selected covariance validation and factorization | 434.90 ms |
| Fisher from factors, 2 parameters | 182.27 ms |
| Fisher from factors, 12 parameters | 185.54 ms |
| Kaiser model | 1.56 ms |
| Full derivative evaluation | 8.86 ms |
| prepare_bin | 656.53 ms |
| run_bin from prepared state | 199.19 ms |
| Evidence assemble | 608.68 ms |
| Payload mathematical checker | 27.12 ms |
| Direct independent covariance/Fisher contraction | 7.91 ms |
| In-memory compressed NPZ | 86.30 ms |
| Grid construction | 0.72 ms |
| Forest weights: 256 samples, 12 updates | 0.46 ms |
| Composed prepare/derivatives/assemble/run_bin check | 1478.02 ms |

At 512 cells, factorization took 102.09 ms and Fisher contraction 45.66 ms:
approximately linear scaling over this small interval. Disk I/O, real reader
interpolation, CAMB, a complete refinement study and actual forecast runtime
were not measured. NPZ timing is compression into BytesIO, not filesystem I/O.
The composed fixture mirrors the duplicate production path in accuracy.evaluate;
it is not an invocation of the real AccuracyRecipe.

The profiled composed call took 2.055 s because profiling adds overhead; use
1.478 s for throughput. Within that profile, two factor_covariance calls took
1.145 s cumulatively, two fisher_from_factors calls 0.458 s, and two field PSD
validations 0.394 s. These three disjoint call families account for about 97%
of this deliberately simple fixture. The real workflow has additional costs,
so this fraction must not be transferred directly to its 1167-second runtime.

## Bottlenecks verified in live code

1. factor_covariance (local-only path: `../fishhighz/fisher.py:31`) loops over cells, invoking
   inspect_information (local-only path: `../fishhighz/_information.py:11`) separately on every
   15-by-15 block. It repeatedly copies/normalizes arrays, builds triangle indices,
   computes eigenvectors unused by factorization, then calls Cholesky. Batched
   validation/eigenvalues/factorization is the first priority. Retain eigenvectors
   for the result/null-space path and all normalized rank/symmetry conventions.
2. fisher_from_factors (local-only path: `../fishhighz/fisher.py:82`) validates/copies per-cell
   factors and Jacobians, then uses a Python forward-substitution loop and a
   small contraction for each cell. Similar costs for 2 and 12 columns show
   Python/validation overhead dominates here. Fused guarded loops are a strong
   targeted-JIT candidate; inactive-column optimization is secondary initially.
3. Field PSD validation (local-only path: `../fishhighz/covariance.py:60`) loops over small
   field matrices. In the synthetic fixture it dominates Wick-product work:
   bare Wick assembly takes only 2.75 ms. The preparation path validates noise
   and total field power separately; these are different arrays and cannot
   simply be treated as one check.
4. AccuracyRecipe.evaluate (local-only path: `../fishhighz/validation/accuracy.py:507`) prepares
   covariance/factors, evaluates J, calls evidence assembly (another factorization
   and Fisher contraction), then run_bin (another derivative evaluation and
   contraction). Preserve a separate numerical oracle while avoiding repeated
   identical production work. study (local-only path: `../fishhighz/validation/accuracy.py:544`)
   retains summaries but recomputes arrays when evaluate is requested again
   for an existing control key. Prepared-state cache size is only three;
   repeated control dependencies warrant measured, bounded reuse.
5. Evidence writing (local-only path: `../fishhighz/validation/evidence.py:189`) calls the full
   payload checker twice, first without and then with a scientific-pass
   requirement. Separate those predicates without dropping mathematical checks.
   Large transient arrays, metadata and compression may become significant after
   the matrix loops are optimized; profile them again at that point.

Forest convergence failure is not evidence that forest iteration is the runtime
bottleneck. Its synthetic 12-update calculation is cheap compared with matrix
handling. Real input preparation remains unmeasured. Likewise, Kaiser/template
JIT is not justified as the first optimization by these results.

## Scratch throughput experiments and their limits

No production functions were replaced. Prototypes in the assessment script
operate only on already valid SPD arrays and omit much of the public error
contract; the batched solve allocates all-cell workspaces. They establish
potential throughput, not acceptable replacements or full-API speedups.

| Prototype | Median | Comparison with current public operation |
| --- | ---: | ---: |
| Batched normalized eigenvalues plus Cholesky | 22.23 ms | 19.6x factorization |
| Batched NumPy forward solve and contraction | 2.47 ms | 73.9x Fisher contraction |
| Batched generic NumPy solve and contraction | 4.79 ms | 38.0x Fisher contraction |
| Fused Numba forward solve/contraction, warm | 0.390 ms | 467.9x Fisher contraction |

First Numba compilation plus execution took 1.306 s, so it does not benefit a
single tiny cold call. No fastmath or parallelism was enabled. Relative norm
residuals versus the current reference were 2.28e-15 for factors, 5.05e-15 for
batched Fisher and 4.93e-15 for compiled Fisher, all below 5e-12 on this fixture.
Malformed, near-singular, subnormal and overflow cases were not tested in these
prototypes. The implementation plan requires those tests and preserves the
public one-cell solve-workspace contract before claiming a production speedup.

## Preliminary expected improvement

A guarded matrix-stage gain of 10–30x is plausible from the throughput
experiments but is not yet measured in production. To avoid applying the 97%
synthetic profile fraction to real data, use explicit hypothetical fractions
f of the historical primary runtime that benefit. With stage acceleration s,
`T_new = 1167.369555 * ((1-f) + f/s)`; this assumes unchanged workload and leaves
all other historical costs unchanged. Additional cold setup is separate.

| Scenario (assumptions, not fitted fractions) | f | s | Overall speedup | Estimated six-bin primary time |
| --- | ---: | ---: | ---: | ---: |
| Conservative | 0.70 | 10 | 2.70x | 432 s (7.20 min) |
| Central | 0.85 | 20 | 5.19x | 225 s (3.75 min) |
| Optimistic | 0.90 | 30 | 7.69x | 152 s (2.53 min) |

This is a provisional scenario envelope, not a statistical interval or measured
end-to-end forecast gain. Real model/reader/setup costs, compression/filesystem
costs, the control-dependent number of cells, failed trial paths, shared-node
noise and strict numerical checks limit extrapolation. Removing duplicate work
may improve on the simple model, but must not be double-counted in stage gains.
A future scientific repair that changes trial counts invalidates equal-work
comparison. The 72 diagnostics require their own model; no separate prediction
is supported here. Compatibility may itself become faster after shared kernels
are optimized. Even these scenarios leave the historical convergence workload
well above its 12.61-second saved-input comparison.

Revision 4 therefore starts with matrix validation/factorization, contraction
and repeated calculation, then mandates another profile and an updated detailed
performance report. Synthetic goals are 5x guarded factorization/contraction
and 3x composed evaluation, with scientific/error parity. They are targets for
implementation, not results of this planning turn.

## Preservation and checks performed

The bounded assessment completed, its synthetic production-path cross-check and
three prototype residual checks passed, and blocked reference imports remained
absent. Historical timing metadata was read without re-running any forecast.
The 1059-test result belongs to the earlier r2 review; no new full pytest suite
or wheel build was needed for these planning-only edits. Production modules,
tests, examples and the r2 reports remain unchanged. Governance documents and
one active assignment were synchronized, with new assessment/archive files.
See preservation-check.json (local-only path: `../.validation/step12-performance-plan-r4/preservation-check.json`)
for the final source/archive identity and document checks. No implementation,
real forecast, dispatch, installation, commit or push occurred.
