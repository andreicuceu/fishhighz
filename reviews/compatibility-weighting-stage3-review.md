# Scientific review: compatibility weighting stage 3

## Scientific summary and verdict

**PASS.** No scientifically consequential correction is needed before the user
reviews Stage 3. The implementation is restricted to the three full-sample
recurrences in the two positive-signal forest-auto contexts, and the returned
states satisfy the stated finite-trajectory stopping prescription. Prefix
variants, non-auto contexts and nonpositive signals are explicitly ineligible.
The existing exactly-three-update path and production defaults are unchanged.

This result supports conditional stopping only on the tested compatibility
magnitude grids through update 96. It does not establish asymptotic convergence,
physical optimality or magnitude-grid convergence. In particular, the
historical full-sum recurrence retains the Stage-2 213-to-425-node shifts of
0.194% in `A` and 0.117% in `P_pixel`, despite converging in iteration.

## Stopping rule and status semantics

The defaults are `min_updates=3`, three stable transitions, `rtol=1e-4` and a
cap of 96. Each transition tests relative weight amplitude, signed normalized
shape, `A`, `P_pixel` and the complete weight vector, without an absolute scale
floor. A candidate is retained only while every intervening transition and
every candidate-to-current comparison remain within tolerance. It is returned
only after the recurrence has actually reached twice the candidate count. This
implements the required protection against the transient plateaus found in
Stage 2.

The API distinguishes all four required outcomes:

- `converged` returns the state at the actual doubled count and records the
  candidate and candidate-to-stop changes;
- `capped` returns the finite state at update 96 and does not label it
  converged;
- `arithmetic_failure` returns the last state for which both weights and
  coefficients were finite;
- `ineligible` performs no recurrence and returns no substitute weights.

`updates` counts completed weight transitions, while `state_updates` identifies
the returned finite weight/coefficient state. This distinction is correct when
an update succeeds but coefficient evaluation fails. `last_step` contains the
five preceding-transition changes. `forward_residual` is explicitly `None`, so
no unreported update is evaluated after the returned state. There is no silent
fallback to update 3 or to the last attempted state.

## Independent numerical checks

I reconstructed the stopping decision from selected Stage-2 NPZ trajectories
using an independent implementation of the five relative measures and candidate
logic. The decisive `rtol=1e-4` cases agree with the Stage-3 evidence:

- Bin 4, 107-node LBG auto, `sum_intrinsic`: candidate 20, returned update 40.
  Returned weights and coefficients are bitwise equal to saved state 40. The
  maximum stopped-to-later changes through update 96 are
  `[8.71e-10, 8.50e-8, 8.30e-8, 3.40e-8, 8.54e-8]` for amplitude, shape, `A`,
  `P_pixel` and vector respectively, all below `1e-4`.
- Bin 4 LBG auto reaches the reported longest `sum_aliasing` result at candidate
  17/update 34 on all three grids. Bin 2 LBG auto reaches the longest checked
  `sum_historical` result at candidate 15/update 30 on all three grids.
- The twelve expected `rtol=1e-4` caps are exactly the bin 0 and bin 1 LBG autos
  for `sum_intrinsic` and `sum_aliasing` on all three grids. Each retains saved
  state 96; the last-transition changes remain too large to claim convergence.
- A forced coefficient failure after one successful weight transition returns
  `arithmetic_failure`, `updates=1`, `state_updates=0`, the seed weights and seed
  coefficients, with no fabricated `last_step` or forward residual. An update
  failure before completion correspondingly leaves both counts at zero.

The saved inventory contains 1,620 results: 192 `converged`, 24 `capped` and
1,404 `ineligible`. All converged results have `updates = 2*candidate`. At
`rtol=1e-4`, actual returned counts span 14--40 for `sum_intrinsic`, 12--34 for
`sum_aliasing` and 10--30 for `sum_historical`; at `rtol=1e-3` they span 10--26,
10--26 and 8--24 respectively. The validation script binds each bin to its
Stage-2 JSON/NPZ evidence, checks bitwise equality at the returned count, and
tests every later saved state through update 96. It also checks all 810
exactly-three-update trajectories.

## Scope of the verdict

The implementation and saved evidence are sufficient for the bounded Stage-3
claim. No BAO forecast or stopping-tolerance effect on BAO uncertainties has
been tested here; that scientific criterion belongs to the separately
authorized Stage 4. This review performed no forecast, production change, full
suite, Slurm action, commit or push.
