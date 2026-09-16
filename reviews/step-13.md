# Step 13 revision 2: normalized cumulative-weight limit

Revision 2 is implemented and ready for independent/user review. The exact
finite recurrence and all fixed-grid results remain verified. The continuum
question is now classified honestly as **unresolved within the stated bounds**
for `lya(qso)` and `lya(lbg)` in all six saved 15x2pt bins. No production
weighting, estimator, input policy, arithmetic guard or forecast result changed.

The revision-1 handoff is preserved byte-for-byte as
[step-13-r1.md](step-13-r1.md), and the original derivation remains
[step-13-weight-limit-r1.md](step-13-weight-limit-r1.md). The corrected technical
assessment is [step-13-weight-limit-r2.md](step-13-weight-limit-r2.md). Numerical
evidence is under
[step13-r2-20260915T010410Z](../.validation/step13-r2-20260915T010410Z/).

## S13-R1: asymptotic conclusion

The diagnostic continues to evolve the exact nonlinear recurrence with retained
log amplitude and relative shape. At fixed finite dimension,
`F(w) <= M w` and `rho(M) < 1` establish amplitude decay; the triangular
eigenvector gives the unique-dominant fixed-grid normalized limit. Each saved
batch now records the exact one-step nonlinear correction
`M w - F(w) = (M w)^2/(1 + M w)`. The 1024-update state is explicitly not used
as an error-bounded asymptotic approximation, so no unsupported accumulated
remainder claim remains.

The continuous source measure is
`dR = rho_z(m) (1+z_source)/c dm` on bounded ordered magnitude support, with
positive saved variance. Revision 2 distinguishes the fixed-grid iteration-first
limit from refinement-first and joint limits. The five available quadratures are
finite and non-nested; they do not verify a non-atomic refinement family,
uniform nonlinear remainder control or limit interchange. The conditional
backward-product/tail-measure concentration argument is derived in the technical
report, but its asymptotic assumptions cannot be established from these saved
orders.

All 24 signed order-doubling ratios for the fixed-grid asymptotic
`P_pixel` are positive, 0.934--0.982. The order-16/64 exponent remains
0.9615--0.9798 and effective-measure exponent -0.9079 to -0.2671. This is strong
empirical grid dependence and fails the screening target, but it is not promoted
to a continuum theorem. Actual distribution means, widths, faint-tail
probabilities and source measures replace the revision-1 dominant-diagonal
coordinate heuristic. All twelve population decisions are unresolved.

Analytic/synthetic controls retain the scalar finite-noise case and the genuine
constant-coefficient linearized Volterra divergence. Added controls include the
bounded sequence that the old rising-trend rule misclassified, a false late
plateau, violation of the constant-variance proportionality, repeated and
near-dominant spectra, and high-precision nonlinear trajectories.

## S13-R2: source and attempt binding

Initialization and finalization independently reconstruct each expected
bin/field/order from the immutable r5 diagnosis and accuracy records. Magnitudes,
masses, variance, `L`, `Delta_v`, `S` and `B` are compared separately at
`5e-12` relative tolerance with zero absolute tolerance. Source/report hashes,
profile matches, attempted checkpoints, caps, all historical comparisons,
finite-trajectory comparisons and final summary operands are recomputed.

All 60 batches and 840 checkpoints completed. Four profile/source matches use
order 32 (bins 0--1, both fields), eight use order 64 (bins 2--5), and every
batch has its separate diagnosis-source comparison. Historical discrepancies
are at most `9.49e-14` in `A` and `1.22e-13` in `P_pixel`. Maximum trajectory
time was 0.801 s; all batches satisfy the 30 s work bound.

The offline checker reconstructs the complete final bundle and compares the
summary, decision table and three SVGs byte-for-byte. The two preserved r1 review
probes now reject. Additional full-bundle probes reject swapped orders,
self-consistent replacement arrays with old provenance, falsified attempt
status, altered historical operands and an altered final conclusion. Truthful
complete, failed, capped and identical-input controls pass their declared
contracts. Results are in `mutation-results.json`.

## S13-R3: scalar numerical range

All recurrence, spectrum and eigenvector uses of `L/Delta_v` now evaluate
`log(L) - log(Delta_v)` before any ratio is formed. Related analytic scalar
operations use log arithmetic or explicit ordinary-value availability. Both
reciprocal extreme one-cell controls give `w0=1/2`, `w1=1/3`, `d=1` and
`A=P_pixel=1`. A direct diagnostic-versus-Decimal comparison passes after every
ordinary float64 weight has underflowed, at two Decimal precisions.

The unchanged production path still rejects the Step 10 mixed-product case and
the saved approximately `8.12e-566` product. Exact zeros and genuinely
out-of-range coefficients retain explicit branches/availability flags.

## Validation and identity

- Focused mathematical/diagnostic suite: **75 passed**.
- Required quick command: **1318 passed, 25 optional-compiler skips** in
  177.78 s; Ruff lint and format passed. Known post-pytest MUNGE warnings were
  nonfatal.
- Isolated exact-wheel affected tests outside the checkout: **136 passed**.
- All six installed examples passed with isolated imports and NumPy backend.
- All **53 source/wheel/installed Python modules** and installed METADATA/WHEEL
  bytes match. Relative to revision 1, 52 modules are unchanged and only
  `fishhighz/validation/weight_limit.py` differs.
- Wheel SHA256:
  `1d616dcb203f3f7fc1d2efba78533df74532c3b401cf214efb2754e58d0b31d4`.
- Interpreter: Python 3.13.15; validation backend: NumPy 2.5.3, one thread.

The first isolated wheel build failed under sandboxed DNS while obtaining the
declared setuptools build dependency; the approved retry succeeded. Both logs
are preserved. The inherited explicit-Numba/SciPy-blocked subprocess limitation
and unmet NumPy matrix-speed target were not assigned and remain unchanged.

## Changed-file boundary and scientific status

The implementation changes only:

- `fishhighz/validation/weight_limit.py`;
- `scripts/diagnose_weight_limit.py`;
- `tests/test_weight_limit.py`;
- the stale Step 13 paragraph and new diagnostic subsection in `README.md`;
- the revision-2 handoff/derivation and preserved revision-1 handoff;
- the unique revision-2 evidence directory.

Production kernels, `prepare_forest_weights`, survey readers, noise,
covariance/Fisher/model code, strict defaults and Step 12 validators are
unchanged. No real model, reader interpolation, CAMB/P3D/P1D call, forecast,
sensitivity suite, Slurm job, commit or push occurred.

The relation
`N_F = [A P1D W^2 + P_pixel] d_deg^2/a_v` shows why the strong positive saved
grid trend matters, but coefficient diagnosis alone does not establish Fisher
convergence. Step 12 retains six unconverged accuracy bins and 12 unavailable
diagnostics. No production prescription is adopted and no forecast rerun is
requested automatically. Stop here for independent/user review; do not advance
to Step 14.
