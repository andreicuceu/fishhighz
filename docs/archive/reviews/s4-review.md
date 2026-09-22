# Independent scientific review: S4 direct attribution

## Scientific verdict

**PASS.** I find no scientifically consequential correction required. The direct
controls support the stated result: the revised mean-squared mixed
forest--galaxy damping convention is the largest isolated contribution to the
smaller accuracy-profile joint BAO errors in all six bins, while the grouped
full-AP derivative prescription consistently opposes part of that reduction.
The calculation is numerically resolved within its finite controls and the
reported tables reproduce the saved Fisher operands.

This result identifies the dominant measured profile difference for the stated
DESI-2 selection. It does not show that the revised mixed-pair reconstruction
physics is preferable, provide an additive causal partition, establish a
continuum limit, select a forecasting profile, or constitute scientific
acceptance.

## Endpoint closure and controlled interpretation

I checked the source inventory against the live Kaiser, response, weighting and
legacy-derivative implementations. Fixed compatibility retains the captured
legacy Fourier measure, mean/derivative prescription and covariance construction,
while revised accuracy changes the spectrum/decomposition, mean redshift,
response convention, damping, AP derivative and numerical measures. Both
endpoints use the same early-lyaforecast recurrence and the same declared
representative angular/velocity mode. A weighting-method change is therefore not
an explanation of the endpoint difference.

The reconstructed endpoints close in all six bins. The largest componentwise
endpoint-error residual is 1.2866421e-6 in fractional units, or 0.000128664%,
well below the 0.1% numerical target. Accuracy Fisher closure ranges from machine
precision to 1.96e-8 relative; fixed closure remains below 1.00e-6 in the saved
matrix norms. The fixed residual correctly includes the small difference between
captured covariance-side interpolation of CAMB growth rate at geometric redshift
and the live exact-redshift reconstruction. The controlled redshift switch
separately changes the mean/J evaluation from arithmetic to geometric centre,
including bias, growth rate, damping and mean response.

The mixed-pair switch changes only the forest--galaxy squared width from the
legacy unreconstructed forest value to the mean of the forest and reconstructed
galaxy auto-width squares. Its forward joint effects span -4.10543% to -2.37603%
radially and -2.38308% to -0.954074% transversely. Reverse effects span +2.46025%
to +4.30386% and +0.981388% to +2.47954%, respectively. The grouped AP-operator
switch gives +1.01977% to +1.17409% and +1.12319% to +1.30083% forward, with the
opposite-sign reverse effects. Thus the dominant damping change is partially
cancelled rather than additively reproduced by the operator and smaller controls.

The operator result is correctly limited to a grouped estimator change: AP
prefactor, remapped coordinates/RSD/damping, fixed observed response and central
differentiation replace the projected backward-k derivative. No unsupported
subterm shares are claimed. The Fourier-only forward switch is also correctly
described as a continuous extension of the grid-dependent legacy estimator,
because it retains the native backward interval and lower-bound treatment on new
nodes.

## Independent numerical checks

I independently inverted every joint and individual Fisher matrix in all 364
trials. The largest error reconstruction discrepancy is 6.66e-16, every selected
matrix has rank two, all saved forest-weight solutions report convergence, and
the forecast selection is exactly three spectra in bin 1 and all 15 spectra in
bins 2--6. The largest tested joint or individual error refinement is 0.0107825%;
the largest joint Fisher Frobenius change is 0.000214949. These are finite GL
endpoint/single-switch refinements, not a convergence statement for the native
legacy estimator or every interaction.

For bin 1, I reconstructed the complete selected Wick covariance and Fisher
contraction from saved `total`, `modes` and `observed_j` arrays for both endpoints
and for forward/reverse mixed damping and AP-operator witnesses. Joint and
individual Fisher/error agreement is 1.47e-14 or better. All reconstructed
selected covariances are positive definite. The joint Fisher differs materially
from the sum of individual Fishers, confirming that the three-spectrum covariance,
rather than an individual-information sum, is used. The four witness summaries
are bitwise identical to their main-run summaries.

The isolated volume control follows sigma proportional to V^(-1/2) to
2.22e-15 relative and changes errors by only +0.00640% to +0.01035% forward, so it
cannot explain the observed reductions. No rank change occurs in any retained
single-switch control. The largest tested two-switch interaction is 0.00413745 in
natural-log units, or 0.413745 log-percent, for an individual result; the largest
joint interaction is 0.281285 log-percent. The separately reconstructed residual
after summing all forward log effects remains below the 1% attribution trigger.
These nonzero values justify the report's nonadditivity caveat but do not overturn
the dominant-driver conclusion.

## Tables, evidence and scope

I reconstructed every endpoint, forward/reverse switch and interaction value in
`.validation/s4/attribution-summary.json` directly from
`.validation/s4/direct-r1/summary.json`; the maximum numerical discrepancy is
zero at stored precision. The interaction statistic is consistently reported as
the natural logarithm; 0.00413745 is explicitly distinguished from 0.413745
log-percent. The published joint ranges, endpoint geometry, rectangular-425
checks, refinement maxima and no-rank-change statement therefore agree with the
raw trial records.

Independent review artifacts are in `.validation/s4-review/`: the six-bin
endpoint rerun, endpoint and witness contraction checks, all-trial inversion and
selection check, and final-table comparison. I ran no Slurm action, full package
suite, commit, push, S5 calculation or production-source modification. S4 is
complete for user review; scientific acceptance and any subsequent profile or
package decision remain with the user.
