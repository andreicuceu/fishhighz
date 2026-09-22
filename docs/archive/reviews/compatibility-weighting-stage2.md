# Compatibility weighting stage 2

Implemented 2026-09-17 against stage 2, including the coordinator's clarification
that a doubled-count candidate must survive the remaining saved trajectory.
This is a bounded, signed-input iteration and magnitude-spacing study, not a
forecast, a continuum result, or adoption of a weighting prescription.

## Scientific result

The full-sample historical recurrence converges numerically in all 12 noise-entering
forest-auto contexts on all three grids. The intrinsic and moment-aliasing full
sums converge in all six QSO autos and LBG autos in bins 2–5; LBG bins 0–1 do not
qualify. Prefix aliasing exhibits transient plateaus, later drift and appreciable
magnitude-spacing sensitivity. Literal intrinsic-prefix arithmetic fails at
updates 15–17 throughout the study; this is a floating-point failure, not proof
of mathematical nonconvergence. All 54 distinct pair preparations were retained,
including the 42 additional contexts that do not enter additive saved noise.

## Calculation and convergence definition

For each of 54 contexts, all five variants were attempted through update 96 on
107, 213 and 425 endpoint-inclusive magnitude nodes: 810 attempted trajectories.
The 107-node run uses the saved literal arrays. The refined grids retain the same
endpoints, full rectangular endpoint weights dm, redshifts and auxiliary scalars.
No accuracy-profile quadrature or interpolated saved weights were substituted.
All finite iterates and I1/I2/I3/A/P_pixel are saved, with checkpoints
0, 3, 6, 12, 24, 48 and 96 where reached. Failure records retain the last valid
state, failed update and reason. The forward fixed-point residual at t is
||F(w_t)-w_t||_inf/||w_t||_inf, available whenever update t+1 was executed;
at the cap or failed transition it is explicitly unavailable (NaN in NPZ only).

Amplitude is ||w||_inf. Signed normalized shape is w/||w||_inf. Amplitude,
shape, A and P_pixel changes use their own relative scales, without absolute
floors. Three consecutive steps, including relative vector residuals, must pass
rtol=1e-3 or 1e-4. A candidate t must pass an actual t-to-2t comparison of
all four measures and the weight vector, with 2t<=96. Every subsequent saved
state must also remain within tolerance of that candidate. This last condition
rejects early plateaus: for example bin-0 mixed-forest prefix aliasing at 107
passes an early doubled comparison but later differs by ~82%; only a later
candidate is supported by the complete bounded trajectory. Some refined QSO
prefix cases have no surviving candidate. No state beyond 96 was evaluated.
These are finite-count numerical criteria, not an existence/convergence proof.

Signed densities, signals and moments are retained. On 107 nodes, 106 of 270
variant/context trajectories contain at least one negative full moment; none
contains an exactly zero full moment in its retained valid states. Per-moment
counts are saved. The classifier separately represents exact zero amplitude,
invalid arithmetic, decay with stable normalized quantities, continued decay,
drift/oscillation, unconfirmed candidates and transient plateaus. No all-zero
amplitude trajectory occurred. Stable A/P_pixel alone never qualifies.

## Counts and iteration ranges

Confirmed finite nonzero trajectories below survived the complete remaining
saved tail. Autos are the 12 contexts supplying additive covariance noise.
Ranges are actual confirmation updates at rtol=1e-4 (twice the selected candidate).
All 54 includes the autos; the per-context table below separates preparations.

| Variant | Nodes | Autos 1e-3 / 1e-4 | All 54 1e-3 / 1e-4 | Auto confirmation range, 1e-4 |
| --- | ---: | ---: | ---: | --- |
| prefix_intrinsic | 107 | 0 / 0 | 0 / 0 | None |
| prefix_intrinsic | 213 | 0 / 0 | 0 / 0 | None |
| prefix_intrinsic | 425 | 0 / 0 | 0 / 0 | None |
| sum_intrinsic | 107 | 10 / 10 | 40 / 40 | [14, 40] |
| sum_intrinsic | 213 | 10 / 10 | 41 / 41 | [14, 40] |
| sum_intrinsic | 425 | 10 / 10 | 40 / 39 | [14, 40] |
| prefix_aliasing | 107 | 11 / 11 | 22 / 22 | [12, 94] |
| prefix_aliasing | 213 | 11 / 11 | 20 / 20 | [12, 72] |
| prefix_aliasing | 425 | 10 / 9 | 18 / 17 | [14, 72] |
| sum_aliasing | 107 | 10 / 10 | 39 / 39 | [12, 34] |
| sum_aliasing | 213 | 10 / 10 | 37 / 36 | [12, 34] |
| sum_aliasing | 425 | 10 / 10 | 34 / 34 | [12, 34] |
| sum_historical | 107 | 12 / 12 | 46 / 46 | [10, 30] |
| sum_historical | 213 | 12 / 12 | 46 / 46 | [10, 30] |
| sum_historical | 425 | 12 / 12 | 43 / 43 | [10, 30] |

Original-grid confirmed bins per pair context (out of 6), identical counts at
both tolerances; selected counts themselves are recorded for each tolerance.

| Context | prefix_intrinsic | sum_intrinsic | prefix_aliasing | sum_aliasing | sum_historical |
| --- | ---: | ---: | ---: | ---: | ---: |
| lya(qso)_lya(qso) | 0 | 6 | 6 | 6 | 6 |
| lya(qso)_qso | 0 | 4 | 2 | 5 | 5 |
| lya(qso)_lbg | 0 | 4 | 1 | 4 | 4 |
| lya(qso)_lae | 0 | 3 | 0 | 2 | 1 |
| lya(qso)_lya(lbg) | 0 | 4 | 5 | 4 | 6 |
| qso_lya(lbg) | 0 | 5 | 1 | 5 | 6 |
| lbg_lya(lbg) | 0 | 5 | 1 | 5 | 6 |
| lae_lya(lbg) | 0 | 5 | 1 | 4 | 6 |
| lya(lbg)_lya(lbg) | 0 | 4 | 5 | 4 | 6 |

At 107 nodes and 1e-4: intrinsic sums have 40 converged, 3 decaying with stable
normalized quantities, 4 still decaying in shape/coefficients and 7 drifting;
moment-aliasing sums have 39 converged, 2 stable-normalized decays, 4 other decays,
8 drifting and 1 unconfirmed; historical sums have 46 converged and 8 drifting;
prefix aliasing has 22 converged and 32 drifting; intrinsic prefixes have 54 invalid.
Full status counts at both tolerances and all grids are in summary.json.

The intrinsic-prefix failures are overflow or division by a zero prefix moment.
For example bin-0 QSO auto at t=14 has min|J1_prefix|=2.61e-316 while its
full-sample I1=4.58e-4 is finite; update 15 overflows the noise division.
Stopping that branch preserves literal arithmetic and its last valid state;
it does not establish that all full-sample normalized quantities diverge.

## Original interpolation and magnitude sensitivity

The signed density sampler reproduces Tracer's raw magnitude mask, target-density
normalization above z=2.15, first-dz/first-dm conversion, quadratic
RectBivariateSpline, original magnitude sentinel and angular/velocity conversion.
Raw SNR tensors are read with SNRReader(smoothing='none'); the original Gaussian
sigma=10 smoothing and linear 3D interpolation then precede the original bright
clamp, out-of-range noise sentinel and post-scaling SNR floor. The floor is the
inherited input policy, not a new density/weight regularization. Auto/galaxy-cross
contexts use the corresponding forest SNR; the mixed forest context uses LBG
density and the raw geometric-mean QSO/LBG SNR before smoothing.

The mixed raw wavelength axes are not identical. The original reader truncates
to the shorter row count, pairs equal row indices and keeps the primary QSO
wavelengths; this study preserves that behavior rather than realigning tables.
This is a compatibility limitation, not a new physical interpolation choice.
All 54 original-node density and variance vectors reproduce bitwise exactly.
Each bin JSON hashes all 26 actual raw files (two density and 24 SNR files) and
its source report. No upstream package import, CAMB evaluation, P3D query or
installation was needed. Provenance is anchored to the live raw files and their
exact match at saved nodes; no fresh upstream forecast was captured.

Maximum relative shifts among autos with confirmed 1e-4 endpoints on both grids
are below, separately for A and P_pixel. They describe spacing sensitivity of
this interpolation and recurrence, not a continuum proof.

| Variant | 107→213: A / P_pixel | 213→425: A / P_pixel |
| --- | --- | --- |
| sum_intrinsic | 5.68e-4 / 4.36e-4 | 3.61e-5 / 8.61e-5 |
| sum_aliasing | 3.82e-4 / 4.42e-4 | 2.89e-5 / 9.66e-5 |
| sum_historical | 3.56e-3 / 2.73e-3 | 1.94e-3 / 1.17e-3 |
| prefix_aliasing | 2.19e-2 / 1.81e-2 | 1.08e-2 / 1.43e-2 |

Individual values and three-update comparisons are in summary.json. The largest
three-update auto spacing shift is prefix_aliasing A, 0.233 on 213→425.
For the diagnostic signed cross contexts, grid sensitivity can be much larger:
among full-sample intrinsic cases confirmed on 107/213, maximum P_pixel change
is 9.46 (946%); aliasing sums 2.85 and historical sums 5.31. These contexts do not
supply additive noise in the saved covariance and must not be conflated with the
auto statistics. Tightening 1e-3 to1e-4 changes confirmed original-grid auto
A/P_pixel by at most 1.01e-5/4.13e-6 (intrinsic),3.26e-7/1.90e-7
(aliasing),1.02e-7/2.52e-8 (historical),3.79e-6/3.39e-7 (prefix aliasing).
No BAO-error tolerance claim follows without the later stage-4 comparison.

## Proposal for stage 3

Support opt-in stopping initially only for the three full-sample variants in
positive-signal forest-auto contexts, with explicit nonconvergence outside
demonstrated cases; do not infer eligibility of negative-signal cross preparations. Historical sums have all 12 autos
eligible; intrinsic and moment-aliasing sums have 10/12 (all QSO bins, LBG bins 2–5).
These same contexts qualify at both tolerances on every tested grid. Do not
infer eligibility of the extra 42 signed contexts from the autos: retain their
individual statuses. For initial controls use minimum 3 updates, three stable
successive steps, rtol=1e-4, cap 96, and actual doubled-count confirmation.
The observed auto candidate/confirmation ranges are7–20/14–40 (intrinsic),
6–17/12–34 (aliasing),5–15/10–30 (historical) on all three grids.

This proposal concerns finite-grid iteration only. Historical autos still have
0.194% A and 0.117% P_pixel spacing sensitivity in the final refinement, both
exceeding 0.1%, despite excellent iteration stability. Prefix aliasing should remain fixed-count-only initially: its apparent
plateaus can fail later, two refined QSO autos lack1e-4 confirmation within 96,
and LBG bin 2 remains drifting. Nine of 12 autos pass both tolerances on all grids,
but that restricted empirical observation does not validate a general early-stop
rule. Intrinsic prefixes remain fixed-count-only under this literal arithmetic.
The subsequent implementation/review must verify proposed stopping against the
saved longer trajectories; no adaptive API was added in stage2.

## Evidence and focused checks

- `fishhighz/validation/weight_convergence.py`: relative metrics, trajectories and
  bounded classification; reuses stage-1 recurrence source unchanged.
- `scripts/study_compatibility_weight_convergence.py`: raw sampling and per-bin
  trajectories; `scripts/summarize_compatibility_weight_convergence.py`: derived
  tables/classification and plots from saved arrays, without new updates.
- `.validation/compatibility-weighting-stage2/bin-{0..5}.npz`: every retained
  weight vector, magnitude/density/variance input, moment/coefficient trajectory,
  amplitude, step changes and available residuals for all three grids.
- Corresponding JSON: statuses/counts at both tolerances, checkpoint coefficients,
  failure reasons and signed-moment counts, raw/source hashes and reader checks.
  `summary.json`: all status/range counts, tolerance and grid comparisons.
- `weights-bin-{0..5}.png` and `convergence-bin-{0..5}.png`: 9×5 panels per bin,
  covering every context and variant on 107 nodes. Signed weights use symlog;
  missing later checkpoints on failed branches are not substituted.

Commands from lib/fishhighz, using OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1,
MKL_NUM_THREADS=1 (each bin run separately):

```text
.venv/bin/python scripts/study_compatibility_weight_convergence.py   .validation/step12-r5-20260914T191855Z/profiles-checked   ../lyaforecast .validation/compatibility-weighting-stage2 --bin N
.venv/bin/python -m pytest tests/test_compatibility_weights.py tests/test_weight_convergence.py -q
.venv/bin/ruff check fishhighz/validation/weight_convergence.py   scripts/study_compatibility_weight_convergence.py   scripts/summarize_compatibility_weight_convergence.py tests/test_weight_convergence.py
git diff --check
```

21 tests passed in 0.15s, including algebraic subthreshold decay with exactly stable
coefficients, tiny-amplitude relative metrics and a false early plateau. Ruff
and whitespace checks passed. The raw study ran in short single-thread per-bin
chunks. The local venv lacks Matplotlib; plotting and derived summaries used the
existing read-only NERSC interpreter
`/global/common/software/nersc/pe/conda-envs/26.8.0/python-3.14/nersc-python/bin/python`,
with `PYTHONPATH=.`, `MPLCONFIGDIR=/tmp/fishhighz-stage2-mpl`, and
`XDG_CACHE_HOME=/tmp/fishhighz-stage2-cache`, invoking the summary script with
`--bin N` and once without that argument. It does not recompute weight updates.
Known post-exit MUNGE warnings and initial font-cache warnings did not alter
successful checks. No full suite, fresh environment, forecast, Slurm operation,
commit, push, production change or later stage was executed.

Ready for the authorized independent scientific review.
