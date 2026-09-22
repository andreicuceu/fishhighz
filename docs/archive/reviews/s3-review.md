# Independent scientific review: S3 three-profile comparison

## Scientific verdict

**PASS.** I find no scientifically consequential correction required in the
three-profile BAO comparison. The tabulated endpoint differences, coverage,
qualification labels and 1% trigger inventory agree with an independent
reconstruction from the frozen S2 Fisher operands. The S3 products support the
reported comparison and are ready for user review.

This verdict confirms the endpoint calculation and its numerical interpretation.
It does not attribute the differences to particular physical or numerical choices,
authorize S4, select a preferred profile, or constitute scientific acceptance of
the revised forecast.

## Independent Fisher and table reconstruction

I used `.validation/s3-review/check_comparison.py` to invert all 252 available
individual and joint 2x2 Fisher matrices directly. The derived radial and transverse
uncertainties, correlation coefficients, ranks and ellipse areas reproduce
`three-profile-tables.json`. The script also reconstructs all three profile ratios
from the raw NPZ operands rather than from the S3 summary.

Each profile has 78 available individual results and six joint results, together
with 12 explicit bin-1 exclusions. All 84 selected matrices per profile have rank
two. There are no failed, unconstrained or plot-hidden valid results, and the 0.2
paired-component plotting cut leaves the numerical tables and ratio checks intact.
The tables preserve `qualified_finite_refinement` for all six accuracy records and
`reproducibility_control` for the full- and fixed-compatibility records. They also
retain the distinct legacy and early-lyaforecast recipe identities, saved weighting
records and accuracy trial contract; no forest convergence state is relabelled as
accuracy qualification.

For each profile, I recomputed the bin-1 joint Fisher by solving the saved 3x3
selected covariance before inversion. It agrees with the saved joint Fisher and is
not the sum of the three individual Fishers: the relative joint-versus-sum
differences are 0.5460, 0.5380 and 0.5453 for full compatibility, fixed
compatibility and accuracy, respectively. Bins 2--6 retain all 15 spectra. The
comparison therefore uses the required observable set at every endpoint.

## Scientific comparison

The accuracy/fixed endpoint is reproduced exactly. Accuracy lowers all six joint
radial errors by 2.6693--3.8610% and all six joint transverse errors by
0.6665--1.6089%. The largest individual endpoint change is a reduction of
8.856208% in the bin-1 Lyα(QSO) × QSO radial uncertainty. The full-compatibility
comparisons also reproduce the reported extrema, including the much larger
LBG-forest-auto changes; these are endpoint differences and do not isolate their
cause.

The independent 1% inventory finds 111 of 168 accuracy/fixed components crossing
the declared trigger, spanning 62 of 84 results: 56 individual results and all six
joint results. Every crossing is a reduction. All joint radial components cross;
joint transverse components cross in bins 1, 2, 3 and 6. The per-bin component and
result counts agree exactly with the S3 handoff, and none is hidden by the plotting
cut. The largest reviewed S2 finite-refinement effect, 0.00127%, is far smaller than
these endpoint differences, so the trigger is not numerically ambiguous at the
tested controls.

There is no rank or availability change between profiles. For accuracy/fixed, the
largest absolute correlation shift is 0.03624 and the largest major-axis rotation
is 0.04519 rad (2.59 degrees); neither reveals an unreported orientation
discontinuity. Individual ellipse-area ratios span 0.8814--1.0166, while the six
joint ratios are 0.9545--0.9745, consistent with the reported 2.55--4.55% joint
area reductions.

The weights, A and pixel-noise panels are correctly described as saved endpoint
diagnostics. Their differences do not identify which change produced a BAO shift.
The handoff makes no causal claim and records that no attribution switches were
run. This is the required interpretation because the 1% threshold diagnoses where
an S4 study would be informative; it does not itself supply an attribution.

The independent check output is
`.validation/s3-review/comparison-check.json`. All 17 frozen output-file hashes,
including 12 PNGs, the complete table and summary, reproduce the artifact manifest;
the separate S3 handoff hash also agrees. I ran no forecast, source modification,
full suite, Slurm action, commit or push. Scientific acceptance and any decision to
assign S4 remain with the user.
