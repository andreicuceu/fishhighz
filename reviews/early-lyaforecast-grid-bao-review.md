# Scientific review: early lyaforecast magnitude-grid BAO diagnostic

## Scientific summary and verdict

**PASS.** No scientifically consequential correction is needed before the user
reviews this extra diagnostic. For `sum_historical` at `rtol=1e-4`, the final
213-to-425-node refinement changes every available individual BAO uncertainty by
less than 0.1%; the largest absolute change is 0.0843233%. The largest joint
change is 0.00562013%. Both occur in the transverse component of user bin 2.

The result does not show that the original 107-node grid is accurate to 0.1%.
The 107-to-425 change reaches 0.286030% for the transverse LBG-forest auto in
user bin 2, while the corresponding maximum joint change is 0.0214464%. Thus
the last tested refinement passes the diagnostic benchmark, but the original
grid fails it for all five individual spectra involving the LBG forest in user
bin 2. These three finite grids do not supply a continuum error bound or justify
extending the grid automatically.

## Scope and stopped-state binding

The calculation uses only internal bins 1--5, corresponding to user bins 2--6
and saved Step-12 records 002, 004, 006, 008 and 010. Internal bin 0 does not
appear in the primary records, derived table or figures. The inventory contains
240 forecasts: 15 individual spectra plus the full joint contraction on three
grids in each of five bins.

I rebound all 30 forest-auto preparations to the Stage-2 trajectory arrays and
Stage-3 stopping metadata. Every preparation is a finite-nonzero
`sum_historical` state at `rtol=1e-4`; its reported coefficient vector is
bitwise equal to the trajectory row at `confirmed_at`, and its update count is
twice the saved candidate. The counts are QSO
`16, 16, 14, 12, 10` and LBG `28, 30, 26, 24, 22` in bin order, identically on
107, 213 and 425 nodes. Every reported stopped-to-later discrepancy is below
`1e-4`. No update-96, last-valid or otherwise hidden iterate enters the BAO
calculation.

Only the two forest-auto aliasing and pixel-noise coefficients change. The
immutable signal, cross powers, galaxy noise, response, observed Jacobian,
Fourier modes, volume and two-parameter convention are retained. The existing
reassembly rebuilds all affected Wick entries, and the joint forecast contracts
the full inter-spectrum covariance.

## Independent numerical checks

I reconstructed all 480 component ratios from the 240 per-bin primary records.
They agree exactly with the JSON summary and CSV, including the comparison
operands and signs. There are no values hidden by the 0.2 plotting cut in this
selected five-bin sample. The global extrema are:

| Comparison | Largest individual change | Largest joint change |
| --- | ---: | ---: |
| 107 to 213 | -0.201877% | -0.0158272% |
| 213 to 425 | -0.0843233% | -0.00562013% |
| 107 to 425 | -0.286030% | -0.0214464% |

For 107 to 425, all ten components of those five LBG-forest spectra exceed the
0.1% benchmark: the four cross-spectrum shifts span -0.112293% to -0.143092%,
and the auto shifts are -0.241641% and -0.286030%. For 107 to 213, five
components across four LBG-forest spectra exceed 0.1%; for 213 to 425, none do.

As an independent covariance/Fisher check, I reconstructed the 425-node
LBG-forest auto in user bin 2 directly from the retained total power, subtracting
the saved baseline noise and adding the selected historical-sum noise. Using
the analytic auto variance `2 T^2 / N_modes` over the saved Fourier cells gives

```text
sigma_parallel   = 0.058228755260148854
sigma_transverse = 0.06725085418211425
```

These agree with the saved forecast to machine precision; the maximum relative
Fisher discrepancy is `1.22e-15`. The corresponding 107- and 213-node errors
give transverse shifts of -0.201877% and -0.0843233% for the two adjacent
refinements, confirming that the final increment is smaller while the total
107-to-425 displacement remains above the benchmark.

All 107-node Fisher matrices and errors are bitwise equal to the reviewed
Stage-4 `sum_historical/converged_1e-4` results. The 90 galaxy-only grid records
(six spectra, five bins and three grids) remain bitwise equal to their immutable
saved values. All recorded source digests match the current bound inputs and
calculation source.

## Tables, figures and limits

The two 16-panel figures contain all 15 individual spectra and the full joint
result. Their ratio construction applies the paired 0.2 uncertainty cut to both
components and requires both grid operands; NaNs would preserve gaps. No point
in this selected scope is cut. The saved tables retain every valid value
independently of plotting.

The implementation handoff reports the focused BAO tests and Ruff checks
passing; the coordinator separately reports 33 focused tests plus lint, format
and whitespace checks passing. I did not repeat the five forecast executions.
This review is restricted to the authorized saved-array diagnostic. It does not
select a production default, change governance or the accuracy profile, include
the first bin or another weighting prescription, establish physical optimality,
or constitute user acceptance. No full suite, Slurm action, commit or push was
performed.

Ready for user review. Stop for review.
