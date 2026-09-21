# Scientific review: compatibility weighting stage 2

## Scientific summary and verdict

**PASS.** No scientifically consequential correction is needed before stage 3.
The saved calculation supports a deliberately narrow adaptive-stopping study for
the three full-sample variants in forest-auto contexts. At `rtol=1e-4`, on each
of the 107-, 213- and 425-node grids, `sum_intrinsic` and `sum_aliasing` converge
in 10/12 autos, while `sum_historical` converges in 12/12. The two excluded cases
for each intrinsic/moment-aliasing sum are the LBG autos in bins 0 and 1.

This is finite-grid evidence through 96 updates. It is neither an asymptotic
convergence proof nor evidence that iteration convergence removes magnitude-grid
sensitivity. In particular, the historical sum remains sensitive at the final
213-to-425 refinement by as much as 0.194% in `A` and 0.117% in `P_pixel`, despite
convergence of its iteration on every auto grid.

## Recurrence and convergence classification

I inspected the five stage-1 recurrences, the trajectory/classification source,
the raw-input study, the summarizer, focused tests, six JSON/NPZ pairs and the
stage-2 summary. The implementation retains weight amplitude and signed inputs;
no floor, clipping, normalization or substitute state is introduced. Every
finite state stores full-sample `I1`, `I2`, `I3`, `A` and `P_pixel`. Arithmetic
failure terminates the branch with its last finite state and reason.

The convergence measures use relative infinity-norm changes with no absolute
denominator floor. They separately test amplitude, signed normalized shape,
`A`, `P_pixel` and the complete weight vector. A candidate requires three stable
updates, an actual comparison to update `2t`, and agreement with every later
saved state through the cap. Stable coefficients alone cannot qualify.

An independent check of the decisive false plateau gives the following result
for bin 0, mixed forest, 107 nodes, `prefix_aliasing`, at `rtol=1e-3`: candidate
`t=12` agrees with `t=24` in amplitude, shape, `A`, `P_pixel` and the full vector,
all within tolerance, but `A` later differs from its candidate value by 0.8183 at
`t=37`. The tail condition therefore correctly rejects that early candidate.
The first candidate supported by the complete trajectory is `t=38`, confirmed
at update 76.

The amplitude-decay classification is also effective on the saved data. For the
bin-0 LBG auto with `sum_intrinsic`, the amplitude falls from 0.9657 to
`1.50e-132`; its final normalized shape and coefficients are stable at roundoff,
but its final relative amplitude change and residual are both 0.9591. It has no
candidate at either tolerance and is classified as decay, not finite nonzero
convergence.

I reconstructed the reported inventory and counts directly from the evidence:
54 distinct bin/context preparations, 12 forest autos, five variants and three
grids give 810 records. At `rtol=1e-4`, the reported auto/all-context counts and
auto confirmation ranges agree with the JSON/NPZ arrays. The 107-node signed
sample has 106/270 trajectories with at least one negative full moment and none
with an exactly zero full moment. All 54 `prefix_intrinsic` trajectories fail
arithmetically at updates 15--17; this is correctly reported as finite-precision
failure rather than mathematical nonconvergence.

The independently recomputed auto-grid maxima agree with the handoff. For
107-to-213 and 213-to-425 nodes respectively, the maximum relative shifts in
`(A, P_pixel)` are `(5.68e-4, 4.36e-4)` and `(3.61e-5, 8.61e-5)` for
`sum_intrinsic`; `(3.82e-4, 4.42e-4)` and `(2.89e-5, 9.66e-5)` for
`sum_aliasing`; `(3.56e-3, 2.73e-3)` and `(1.94e-3, 1.17e-3)` for
`sum_historical`; and `(2.19e-2, 1.81e-2)` and `(1.08e-2, 1.43e-2)` for
`prefix_aliasing`, using only cases confirmed on both grids.

## Raw-input reconstruction and scope

The raw sampler reproduces the compatibility preparation rather than introducing
a new input policy. It applies the saved magnitude mask and target-density
normalization, first-spacing `dz` and `dm`, quadratic density spline, magnitude
sentinel, and division by `2.998e5/(1+z_q)`. It forms raw SNR tensors without
prior smoothing, applies the legacy Gaussian smoothing before linear 3D
interpolation, then applies the bright clamp, out-of-range sentinel and SNR floor.

For the mixed forest, it uses the LBG density and forms the QSO/LBG geometric mean
row by row before smoothing. The QSO and LBG wavelength arrays have equal lengths
here but differ by 0.05 Angstrom at corresponding rows; retaining the primary QSO
axis exactly matches the legacy implementation. This is a preserved compatibility
limitation, not a physical wavelength realignment. Direct `np.array_equal` checks
reproduce all 54 original-node density vectors and all 54 variance vectors
bitwise; the stored relative reader discrepancies are exactly zero. The evidence
also binds the 26 raw files and six source reports by digest.

The additional 42 cross-preparation contexts are useful signed algebraic
diagnostics, but they do not enter additive forest noise in the saved covariance.
Their sometimes very large grid sensitivity cannot be promoted to an auto-noise
or BAO conclusion. Conversely, no BAO-error tolerance follows from coefficient
stability before the stage-4 Fisher comparison.

## Stage-3 restriction

Stage 3 may proceed with opt-in stopping for `sum_intrinsic`, `sum_aliasing` and
`sum_historical` in forest-auto contexts only, with minimum 3 updates, three
stable successive steps, `rtol=1e-4`, cap 96 and actual doubled-count
confirmation. Positive signal is a necessary domain restriction here, not a
guarantee of convergence: the four low-redshift LBG variant/bin cases identified
above must report nonconvergence at the cap. Tests should compare every returned
stopped state to the later saved trajectory through update 96. Prefix variants
and all non-auto contexts remain ineligible, and failure or cap exhaustion must
not silently return a substitute iterate.

## Checks run

With one BLAS/OpenMP thread:

```text
.venv/bin/python -m pytest tests/test_weight_convergence.py -q
3 passed in 0.18 s

.venv/bin/python - <<'PY'
# Independent JSON/NPZ inventory, bitwise input replay, status/count/range,
# false-plateau, amplitude-decay and grid-sensitivity calculations above.
PY
810 records; 54 contexts; 12 autos; 108/108 input arrays bitwise equal
```

The known post-exit MUNGE socket message followed successful Python execution.
No forecast, full suite, source correction, Slurm action or later-stage work was
performed in this review.
