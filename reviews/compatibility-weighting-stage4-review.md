# Scientific review: compatibility weighting stage 4

## Scientific summary and verdict

**PASS.** No scientifically consequential correction is needed before the user
reviews Stage 4. The saved-input calculation compares all five prescribed
weighting variants at exactly three updates and uses converged coefficients only
where finite nonzero convergence was established on the 107-node compatibility
grid. Unavailable forest-auto states are propagated to the dependent individual
spectra and full joint result without substituting a capped, failed or partial
state.

The covariance calculation changes only forest-auto noise. It retains the saved
intrinsic signal, galaxy noise, cross powers, BAO Jacobian, Fourier modes and
parameter conventions, applies the response once, and rebuilds every affected
Wick covariance entry. Individual forecasts use their own diagonal variance;
the joint forecast contracts the full inter-spectrum covariance. The numerical
result therefore supports the assigned finite-iteration comparison. It does not
establish magnitude-grid convergence, physical optimality or adoption of a new
weighting profile.

## Convergence availability and state use

All five variants retain 90 individual and six joint results at three updates.
For the `rtol=1e-4` converged comparison, the independently recovered
individual/joint availability counts are:

| Variant | Individual / 90 | Joint / 6 |
| --- | ---: | ---: |
| `prefix_intrinsic` | 36 | 0 |
| `sum_intrinsic` | 80 | 4 |
| `prefix_aliasing` | 85 | 5 |
| `sum_aliasing` | 80 | 4 |
| `sum_historical` | 90 | 6 |

I compared every available forest-auto preparation with the saved Stage-2
107-node coefficient trajectory at its reported update. The three-update states
match trajectory index 3 exactly. Converged states match the saved confirmed
count exactly and retain the reviewed Stage-2 classification. The prefix states
are correctly labelled as Stage-2 tail-tested evidence; their use here does not
make them eligible for the Stage-3 adaptive API. Every unavailable auto has null
coefficients and no reported update, and every dependent forecast has null
Fisher/errors with an explicit missing-auto reason. The joint result is present
if and only if both forest autos are available.

The 36 available `prefix_intrinsic` individual results are precisely the
galaxy-only spectra and therefore contain no forest-weight dependence. Across
all variants and modes, 540 galaxy-only Fisher matrices and errors are bitwise
equal to the immutable saved individual outputs.

## Independent numerical and table checks

The six bin files contain 1,440 forecast rows. All 1,260 available Fisher
matrices are finite and rank two; the remaining 180 rows are unavailable only
because a required converged forest auto is absent. I recomputed all source
SHA256 bindings and checked the saved coefficient operands rather than relying
on the output metadata alone.

The three-update compatibility control reproduces the saved baseline in every
bin. The maximum relative discrepancy is `4.67e-17` for the joint Fisher matrix
and `9.11e-17` for joint BAO errors; individual errors agree exactly. The direct
three-cell Wick construction agrees exactly with the rebuilt covariance, and
the independent per-cell Fisher solves agree with the vectorized contraction to
`1.89e-15` or better. Together with inspection of the reconstruction path, this
checks that the response and forest noise are each applied once and that the
saved signal and Jacobian remain fixed.

Tightening the stopping tolerance from `1e-3` to `1e-4` was evaluated for all
390 results available at both tolerances, including the five results hidden by
the plotting cut. The largest fractional uncertainty change is
`8.669219007106932e-7`, well below the required `0.001`. The 1,440 saved
fractional-change rows reproduce their table operands exactly. At convergence,
available joint errors differ from the three-update compatibility baseline by
`-7.1497%` to `+0.8956%`; the within-variant convergence changes are
`-0.0371%` to `+0.9450%`. These signed changes are scientific results, not
evidence that one variant should become a default.

The paired plotting cut is applied when either uncertainty exceeds 0.2 or is
unavailable/nonfinite. Both components are then omitted, and a fractional ratio
is present only when both operands survive the cut. All rows, including hidden
valid values and reasons, remain in the tables. The eleven figures preserve
gaps and the joint panel shows only the bins with a complete full-covariance
forecast.

## Focused validation and scope

I ran the dedicated Stage-4 unit tests independently:

```text
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    .venv/bin/python -m pytest tests/test_compatibility_bao.py -q
3 passed in 1.22s
```

The tests cover the paired 0.2 cut and ratio operands, unaffected galaxy
information when a forest auto is unavailable, and a correlated joint Fisher
calculation that differs from the sum of individual Fisher matrices. The
coordinator separately reports all 33 focused tests and Ruff lint/format checks
passing for the four Stage-4 files.

This verdict is restricted to the saved six-bin compatibility calculation and
its finite trajectories through update 96. It does not resolve the previously
measured magnitude-spacing sensitivity, justify physical accuracy, alter
production defaults or constitute user acceptance. I ran no new forecast, full
suite, upstream capture, Slurm action, commit, push or Stage-5 analysis.

Ready for user review. Stop for review.
