# Independent scientific review: S2 numerical refinement

## Scientific verdict

**PASS.** I find no scientifically consequential correction required before the
coordinator proceeds to the authorized S3 comparison using the existing 18-record
bundle. All six revised-accuracy bins satisfy every prescribed isolated and
combined finite-refinement criterion. All 11 required forest preparations converge
at both tested tolerances, and bin 1 contains only the three authorized QSO-forest
auto, QSO auto and cross spectra.

This verdict qualifies the reported calculation at the finite controls that were
tested. It does not establish a continuum error bound, select a preferred weighting
profile, perform S4 attribution, or constitute scientific acceptance of the
forecast.

## Independent reconstruction of the accuracy qualification

I independently read the frozen manifest, report and NPZ operands with
`.validation/s2-review/check_accuracy_qualification.py`. The check recomputes BAO
errors from each saved Fisher matrix, reconstructs every joint and per-spectrum
refinement metric, verifies the selected trial identities and adjacent operands,
and checks the adaptive stopping records without calling the forecast runner.

The evidence contains six accuracy records and 78 successful trials, 13 per bin,
with no failed or unresolved trials. Every bin stops at the same controls:
`k_intervals=128`, `mu_order=32`, `magnitude_order=16`, `z_order=32`,
`step=2.5e-4`, and `weight_rtol=1e-5`. The recorded levels are exactly
32/64/128 in k, 8/16/32 in mu, 4/8/16 in magnitude, 8/16/32 in volume,
1e-3/5e-4/2.5e-4 in derivative step, and 1e-4/1e-5 in weight tolerance.
For each coordinate, the saved metric operands are the final state and the stated
immediately preceding level; the combined lower state changes all six controls
together.

Recomputing the stored contract metrics gives the following extrema over all bins
and all seven refinement families:

- Fisher relative change: `2.4529926003e-5`, combined refinement in bin 4;
- joint BAO-error change: `1.2685029479e-5`, combined refinement in bin 4;
- individual BAO-error change: `1.1410945542e-5`, LBG auto in the same bin and
  refinement;
- volume relative change: `2.0991652860e-11`, volume refinement in bin 6.

All are below the declared limits of `1e-3` for Fisher and BAO errors and `1e-6`
for volume. Direct high/low signed component ratios reproduce the handoff values:
the largest individual change is `-1.1410815334e-5` and the largest joint change is
`-1.2684868571e-5`, both for the transverse component in bin 4. The small difference
between these values and the contract metrics is solely the direction of the ratio;
it has no bearing on qualification.

The stopping inventory independently reproduces the reported update counts. QSO
forests use 18/22 updates in bin 1 and then 16/20, 16/18, 14/16, 12/14 and 10/12
for `rtol=1e-4/1e-5`; LBG forests use 28/36, 30/36, 26/30, 24/28 and 22/26 in
bins 2--6. Counts are invariant among the retained trials at a given bin and
tolerance. At final tolerance, the largest preceding-step residual is
`1.0675332e-11` and the largest doubled-candidate confirmation residual is
`1.5855689e-6`, below `1e-5`. Every returned state has finite positive weights,
positive A, non-negative pixel noise, and no candidate approaches the 96-update
cap.

## Fixed-compatibility grid diagnostic and evidence scope

The second review script, `.validation/s2-review/check_bin1_grid.py`, independently
reconstructs the selected bin-1 covariance cell by cell from the saved response,
P1D, coefficients and three Wick elements. Direct NumPy contractions reproduce all
three individual Fishers and the selected joint Fisher. The 213-to-425-node change
is at most `4.1283341874e-5` individually and `1.5121609606e-5` jointly; the QSO
galaxy auto is unchanged. The final increment reverses sign and exceeds the
preceding increment, so this diagnostic supports the stated finite comparison but
does not imply monotonic or continuum convergence.

I also verified that the reused bins 2--6 evidence and its numerical operands retain
the recorded hashes. Those historical forecasts were not rerun or relabelled. The
accuracy bundle retains three selected spectra in bin 1 and all 15 spectra in bins
2--6, with every selected individual spectrum included in the refinement tests.

Both independent saved-operand checks pass with one numerical thread. I ran no real
forecast, broader suite, Slurm action, commit or push. S2 is numerically qualified
within its stated scope; scientific acceptance remains with the user.
