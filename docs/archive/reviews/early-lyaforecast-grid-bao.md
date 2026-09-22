# Early lyaforecast magnitude-grid BAO diagnostic

The historical full-sum prescription (`sum_historical`) passes the diagnostic
0.1% BAO-error benchmark for the **213→425-node refinement** in every individual
spectrum and the full joint forecast in user bins 2–6. The largest absolute
fractional changes are 0.0843233% (individual) and 0.00562013% (joint), both
transverse in user bin 2. However, the original **107→425** comparison reaches
0.286030% for the LBG-forest auto: the original grid is not uniformly within
0.1% of the finest tested grid. The joint difference is at most 0.0214464%.
These finite-grid comparisons do not establish a continuum error bound.

Only user bins 2–6 (internal indices 1–5, mean redshifts 2.3525, 2.5875, 2.8225,
3.0575, 3.2925) enter any calculation, table or figure. The first bin is excluded.
No other weighting variant is evaluated.

## Retained calculation

The diagnostic script (local-only path: `../scripts/check_early_lyaforecast_grid_bao.py`) loads
immutable Step-12 revision-5 compatibility records 002, 004, 006, 008 and 010,
and the corresponding saved stage-2 trajectories and stage-3 stopping metadata.
It selects finite-nonzero confirmed historical-sum states at rtol=1e-4 on each
of 107, 213 and 425 magnitude nodes. It checks source hashes, bin/context/grid,
classification, candidate and confirmed counts, exact coefficient equality to
stage 3, and its reviewed stopped-to-later trajectory discrepancy below 1e-4.
No raw interpolation or new weight iteration is performed. Confirmation counts
on all three grids are QSO 16,16,14,12,10 and LBG 28,30,26,24,22 in bin order.

Only the two forest-auto A and P_pixel coefficients change. The existing
`validation.compatibility_bao.forecast` subtracts the original forest noise and
adds the changed noise, retaining the intrinsic covariance-redshift signal.
Original mean Jacobians, k/mu nodes, mode counts, volumes, response convention,
galaxy noise, density/input policies, signed measures and all cross powers are
fixed. The same covariance contraction gives 15 individual two-parameter BAO
forecasts and the full joint forecast, including inter-spectrum covariance.
No plot cut enters either Fisher contraction. Original arrays remain immutable.

## Results

Entries below are signed percentages at the maximum absolute change, pooling
both error components. Individual maxima are transverse LBG-forest autos except
for 213→425 in user bins 3–6, where they are transverse QSO-forest autos. Joint
maxima are transverse in bins 2–4 and parallel in bins 5–6.

| User bin | Individual 107→213 | Individual 213→425 | Individual 107→425 | Joint 107→213 | Joint 213→425 | Joint 107→425 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | −0.201877% | −0.0843233% | −0.286030% | −0.0158272% | −0.00562013% | −0.0214464% |
| 3 | −0.0241031% | +0.00390517% | −0.0210464% | −0.00373492% | +0.000791930% | −0.00294302% |
| 4 | −0.0310716% | +0.00331090% | −0.0304381% | −0.00293040% | +0.000167288% | −0.00276311% |
| 5 | −0.0297454% | +0.00254420% | −0.0276404% | −0.00168030% | +0.000198320% | −0.00148199% |
| 6 | −0.0287418% | +0.00186488% | −0.0273496% | −0.00621090% | +0.000384098% | −0.00582682% |

For the most sensitive LBG-forest auto (user bin 2), transverse errors are
0.06744376343338815, 0.06730761019907051 and 0.06725085418211424 on the three
grids. Joint transverse errors are 0.014290396915373416, 0.01428813514680011 and
0.014287332134498343. Coefficient sensitivity alone does not determine the
BAO-error sensitivity because A and P_pixel can change in opposite directions.
In the 107→425 comparison, both components of all five spectra containing
the LBG forest in user bin 2 exceed 0.1% (10 components): the four cross-spectra
change by approximately −0.112% to −0.117% parallel and −0.140% to −0.143%
transverse; the auto changes by −0.241641% and −0.286030%, respectively.
All other individual and joint components remain below 0.1% for this comparison,
and every component passes the 213→425 benchmark. The 0.1% threshold here is a
diagnostic benchmark, not a stipulated absolute grid-accuracy guarantee. No grid was extended after examining these results.

Fractional changes CSV (local-only path: `../.validation/early-lyaforecast-grid-bao/fractional_changes.csv`)
and JSON summary (local-only path: `../.validation/early-lyaforecast-grid-bao/summary.json`) retain
both sigma components and their ratio operands for all 15 spectra plus joint,
five bins and three comparisons (480 component comparisons). Per-bin JSON files
`bin-1.json` through `bin-5.json` retain all 240 forecasts, Fisher matrices,
selected coefficients/counts, source digests and checks.

Parallel figure (local-only path: `../.validation/early-lyaforecast-grid-bao/parallel-refinement.png`)
and transverse figure (local-only path: `../.validation/early-lyaforecast-grid-bao/transverse-refinement.png`)
each show 15 individual panels and one joint panel. Both error components are
masked together if either exceeds 0.2; ratios require both operands to pass;
NaNs preserve gaps. All points in these five bins survive this plotting cut.
The saved extrema include every valid point irrespective of the cut.

## Focused verification and reproduction

All 107-node errors and Fisher matrices reproduce the stage-4
`sum_historical/converged_1e-4` results **bitwise**. Galaxy-only individual errors
and Fisher matrices remain **bitwise unchanged** on all three grids. An independent
scalar Wick construction at three Fourier cells in user bin 2 on the 425-node
grid agrees exactly; per-cell linear solves reproduce the vectorized Fisher
contraction to relative 2.12e-15. All selected states match the independently
reviewed stage-3 finite-tail evidence. The focused existing BAO tests pass
(3 tests, 0.27 s), including analytic covariance/Fisher and paired masking checks.
Ruff lint and formatting pass. No new production API or defaults were changed.

From the package root, the exact numerical commands are:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for bin in 1 2 3 4 5; do
    PYTHONPATH=. .venv/bin/python scripts/check_early_lyaforecast_grid_bao.py --bin "$bin"
done
PYTHONPATH=. MPLCONFIGDIR=/tmp/fishhighz-early-grid-mpl XDG_CACHE_HOME=/tmp/fishhighz-early-grid-cache /global/common/software/nersc/pe/conda-envs/26.8.0/python-3.14/nersc-python/bin/python scripts/check_early_lyaforecast_grid_bao.py --summarize
.venv/bin/python -m pytest tests/test_compatibility_bao.py -q
.venv/bin/ruff check scripts/check_early_lyaforecast_grid_bao.py
.venv/bin/ruff format --check scripts/check_early_lyaforecast_grid_bao.py
```

First-bin execution took 2.06 s; the remaining per-bin executions took
1.93–2.28 s, all single-threaded on the login node. Known post-exit MUNGE warnings
did not alter successful test status. Existing evidence and unrelated dirty
files were preserved. No full suite, Slurm operation, commit, push, governance
edit, accuracy reassessment or stage-5 decision was made.

Independent scientific review and user acceptance remain separate. Stop for review.
