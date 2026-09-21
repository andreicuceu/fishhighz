# Independent scientific, software and documentation review: S5

## Verdict

**PASS.** I find no consequential scientific, software or documentation
correction outstanding. The new guide and example consistently expose the
user-selected S2--S4 accuracy prescription without changing its numerical
implementation or relabeling the retained evidence. The example is a useful
standalone API calculation and is clearly separated from the numerically
qualified DESI-2 result.

This review confirms internal consistency, practical installation and the stated
finite evidence scope. It does not independently validate the mixed-pair
reconstruction physics, establish a continuum limit or constitute scientific
acceptance of the forecast.

## Scientific prescription

I checked [RESEARCH_BASELINE.md](../RESEARCH_BASELINE.md) against the live
weighting, response, Kaiser, covariance, survey and accuracy-profile source and
against the S1--S4 handoffs. It records the full-sample early-lyaforecast
recurrence without weight renormalization, the representative mode
`(2.4 deg^-1, 0.00035 s/km)`, and the final coefficients with the correct units.
The final forest noise uses mode-dependent intrinsic P1D and the field response
only in the aliasing term; pixel noise remains unsmoothed.

The guide also records the physical FWHM conversion, exact fiducial template
growth, CAMB growth rate, damping normalization, and mean-squared auto widths for
cross spectra. Its AP prescription leaves the smooth component at identity and
maps the complete wiggle contribution, including the AP prefactor, Kaiser
factors and damping. Weights, response, noise, integration grid, volume and
covariance remain fixed while differentiating the selected mean; no covariance
derivative is added.

The fit is explicitly two AP parameters per bin, with no broadband parameters,
nuisance parameters or priors. Bin 1 retains only the Ly-alpha(QSO) auto, QSO
auto and their cross in the full three-spectrum covariance; bins 2--6 retain all
15 spectra. The `full-compatibility`, `fixed-compatibility` and `accuracy`
identities have distinct reproduction, weighting-isolation and research roles.
Strict public nonnegative inputs are separated from the validation-only density
floor and literal signed recurrence paths.

The numerical record correctly distinguishes the public/recipe default
`rtol=1e-4` from the final S2/S3 evidence state at `rtol=1e-5`. It gives the
retained k, mu, redshift, magnitude and finite-difference controls and preserves
the finite-test interpretation of the S2/S4 refinement maxima. D1 reference-mode
sensitivity and D2 early-versus-W12 comparison remain deferred.

## Independent example and covariance checks

I reconstructed the example's selected Wick covariance cell by cell from its
saved total powers and mode counts, without using production covariance
assembly. Contracting that covariance with the observed mean Jacobian reproduces
the public joint Fisher matrix to `4.50e-17` relative. Its off-diagonal
inter-spectrum covariance is nonzero. The joint Fisher differs from the sum of
the three single-spectrum Fishers by `1.3562509043` in the reported maximum
relative metric, directly demonstrating that the joint result is not an
independent-information sum.

The synthetic example varies exactly `ap` and `at`, adds no prior, and returns a
rank-two result with errors `(0.0468765060, 0.0225581179)`. Early weights converge
after 10 updates from candidate 5 under the stated `(rtol, min_updates,
stable_steps, max_updates)=(1e-4,3,3,96)` controls. An independent scalar damping
check reproduces the mean-squared cross-width rule, and the recorded response
matches the physical resolving-power FWHM conversion exactly. Repeating the
derivative calculation with a different batch size and step scale leaves the
prepared grid, modes, response, noise, total power, covariance factors and forest
weights bitwise unchanged.

The example constructs its template, geometry and survey arrays in memory and
imports no Vega, lyaforecast or CAMB package. Its Einstein--de Sitter background,
biases and damping widths are labeled illustrative; its numerical output is not
presented as a DESI-2 forecast or as inheriting the S2 qualification.

## Software, documentation and distribution

The README now lists all five public weighting choices and explains adaptive
versus fixed-count behavior, explicit failure, auxiliary routing and caller
ownership of nondefault stopping-control provenance. Public docstrings describe
the same method set and available preparation/result metadata. Three validation
error messages now enumerate the already-supported adaptive methods. The package
module description no longer claims that forecasting is unavailable.

I checked 10 local documentation links and the source-distribution manifest.
The research guide and examples are included in the source archive. Coordinator
checks build the wheel from the source distribution, compare all 61 installed
modules byte for byte with the wheel and current source, and run the extracted
example with the installed package under isolated Python from `/tmp`. Its JSON
matches the source execution exactly and no neighboring scientific package is
imported. The installation inherits the recorded read-only NumPy/SciPy/Astropy
environment, so this is an installation and checkout-independence check rather
than a newly solved minimal environment.

The coordinator's preservation check finds every package module's executable AST
unchanged after removing docstrings and normalizing the three enumerated error
messages. All 18 retained decision, S1--S4 report and numerical-evidence files
also remain byte-identical. Thus S5 changes documentation, messages, examples,
tests and distribution contents, not the accepted numerical recipe or historical
evidence.

Independent review evidence is in `.validation/s5-review/review-check.json` and
its replayable `check_review.py`. My focused command

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pytest -q \
  tests/test_research_bao_example.py tests/test_full_sum_weights.py \
  tests/test_survey_composition.py
```

passes all 14 tests. The coordinator additionally records 101 focused source
tests, two installed-package tests, Ruff, formatting, Markdown, whitespace,
distribution and preservation checks. No real-survey forecast, broad reference
suite, Slurm action, commit or push was performed. S5 is complete for user
review; scientific acceptance remains with the user.
