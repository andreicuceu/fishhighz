# Step 12 revision 2 handoff

Implementation and the prescribed bounded full validation are complete and ready
for review. **The scientific acceptance gate fails:** 44/78 primary records pass
(39 compatibility plus 5 accuracy); 34 accuracy bins remain unconverged. All
468 diagnostic requests were attempted, with 422 completed and 46 unavailable
because of guarded numerical underflow. `execution_finished=true` and
`complete=false` are intentional. This is not an acceptance claim.

## Preservation and scope

The revision-1 handoff is preserved byte-for-byte as [step-12-r1.md](step-12-r1.md).
The original dirty/untracked checkout at `0d69786a06d5d676564a51fad14a7156951c2228`
was inventoried before edits. Governance, `pyproject.toml`, all six examples,
production scientific modules, adapters, sibling source/data and historical
baseline environments were preserved. No Step 13 work, delegation, commit,
push, scheduler action or change to a scientific tolerance was made.

Evidence lives under `step12-r2-implementation` (local-only path: `../.validation/step12-r2-implementation/`).
`before.json`, `status-before.txt`, the assignment snapshot, `changed-files.json`
and the final preservation audit bind this work to the reviewed tree. Changes
are confined to the validation layer, its scripts/tests, README and handoff.

## Evidence contract and R1 repair

Schema 2 binds exact recipe, case, bin edges, profile, observed-field identities,
selected/required pair closure and ordering, parameter IDs, grid nodes/measures,
response ownership and numerical thresholds. The writer and offline checker
reconstruct Wick covariance and Fisher information from saved total powers,
observed Jacobians and modes. The checker uses direct solves independently of
the writer's factor path, inspects information sign/rank and recomputes named
covariance, errors and correlations. Unavailable errors have an explicit mask;
tables/plots encode them as null rather than finite error bars.

The five original R1 probes reject. Consistently rehashed corruptions exercise
wrong dimensions, negative information, inconsistent errors/covariance/correlation,
pairs, parameter order, bins, nodes/J alignment, pass status and false convergence.
The standalone attribution checker additionally reconstructs stage differences,
branch reference values and branch differences from its saved matrices. Loading
that script and checking saved attribution matrices also pass an import blocker
for FishHighz model modules, CAMB and lyaforecast; model setup imports are lazy.
Valid finite, null-direction and roundoff controls pass. Scientific success
requires boolean `passed=True`; execution completion alone does not pass.
Historical schema 1 has an explicitly limited inspector and cannot satisfy this
gate. `scripts/check_desi2_full.py` additionally requires the exact 78 primary and
468 diagnostic requests, so a narrowed full-tier selection cannot pass the full
assignment gate.

Real reports are individual hash-bound JSON files, loaded one at a time. This
repairs an I/O scaling problem found during execution: inline sampled-input
reports had grown the partial manifest beyond 76 MB. The interrupted bundles and
wheels remain intact. The corrected wheel differs only in `validation/evidence.py`,
`profiles.py` and `plots.py`; all 45 remaining package modules, including schema,
scientific code and numerical helpers, are byte-identical. Twenty-nine completed
records were revalidated and their C/F reassembled with the corrected wheel.
Their original powers/Jacobians/convergence-input producer and verified module
identities remain explicit in `cached_numerical_inputs`. No reference rerun or
silent numerical restudy was substituted for this reuse.

## Reproducible commands and environments

Commands run from the component root with `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1` and `MKL_NUM_THREADS=1`. Reference/profile work uses the
existing read-only `../lyaforecast/.validation/dev-env/bin/python`; the actual
imported lyaforecast source directory is checked against the requested checkout.
It uses NumPy 2.3.5, SciPy 1.15.3, Astropy 8.0.1 and CAMB 2.0.1. The final fresh
`wheel-env-streamed` uses Python 3.13.15, NumPy 2.5.3, SciPy 1.18.1 and Astropy
8.0.1. Every resolved source/resource and actual import origin is inventoried.

The final wheel is
`dist-streamed/fishhighz-0.1.0.dev0-py3-none-any.whl`, SHA-256
`be9015e0a0f945357018f92079a76d57581ea48d5a8f94862abec7ea63d421dc`.
All 48 source/wheel/installed modules and installed METADATA/WHEEL bytes match.
Six examples, optional-dependency blockers and schema/profile regressions run
outside the checkout with `-I`. The real interpreter explicitly imports this
wheel through the saved `runtime/streamed_wheel_entry.py`; no historical
reference environment was modified.

The exact commands, logs and copied drivers are retained in the evidence directory.
The main invocations are:

```bash
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh
git diff --check

# Fresh actual NewForecast calls, seven serial cases:
PYTHONPATH=. ../lyaforecast/.validation/dev-env/bin/python \
  scripts/capture_desi2_reference.py --suite full --reference ../lyaforecast \
  --output .validation/step12-r2-implementation/reference-01

# Real profiles; the saved entry driver inserts the exact wheel on sys.path:
../lyaforecast/.validation/dev-env/bin/python -I \
  .validation/step12-r2-implementation/runtime/streamed_wheel_entry.py \
  scripts/desi2_profiles.py run --suite full --reference ../lyaforecast \
  --template ../vega/vega/models/Planck18/Planck18_z_2.406.fits \
  --reference-bundle .validation/step12-r2-implementation/reference-01 \
  --wheel .validation/step12-r2-implementation/dist-streamed/fishhighz-0.1.0.dev0-py3-none-any.whl \
  --reuse-completed .validation/step12-r2-implementation/profiles-final \
  --output .validation/step12-r2-implementation/profiles-streamed
```

For offline regeneration, use `scripts/plot_desi2.py`, `audit_desi2.py`,
`attribute_desi2.py` and `check_desi2_plots.py` as documented in README.
`attribute_desi2.py --replot SAVED_ATTRIBUTION --output NEW_DIRECTORY` verifies
stored matrices and renders saved values without loading a template or evaluating
a model. Every output directory is exclusive. Plot caches are confined to the
new validation directory.

Build/install logs retain the initial sandbox DNS failures and successful approved
retries. Development logs also retain the corrected test-fixture/import-format
failures, interrupted preliminary profile runs, and the fixed arithmetic-mean
diagnostic setup. MUNGE/font-cache messages are environment output, not submitted
scheduler work or scientific failures. Current successful checks supersede these
attempts without deleting their history.

## Completed implementation checks

- `check-handoff-final.log`: **1059 passed in 51.39 s**, Ruff checks passed,
  134 files already formatted. The original 1017 tests are retained.
- `pytest-offline-installed-final.log`: **226 Step 12 tests passed in 18.15 s** in the fresh
  installed environment outside the source tree with isolated imports. The final
  attribution-difference corruption regression also passed in that environment
  (`attribution-differences-installed.log`, 11.64 s).
- `probe-streamed-installed.log`, `step12_probe-streamed-installed.log`,
  `desi2_synthetic-streamed-installed.log` and `revision2_probe-streamed-installed.log`:
  exact wheel payload, six examples, optional-dependency blockers, independent
  numerical oracles and all 78 synthetic profile records passed.
- `policy-final.log` / `policy-final.json`: **104 actual Tracer/Spectrograph checks**,
  maximum relative discrepancy 3.1516201241062103e-16. Adapter/scientific bytes are
  unchanged between that wheel and the final artifact.
- `external-streamed.log` and `external-streamed/manifest.json`: actual caller-built
  PowerSpectrum, explicit amplitude wrapper, analytic trace oracle and schema-2
  reconstruction passed using the final exact wheel. The earlier finite-difference
  amplitude calculation was 248739.3526956913 versus oracle 248739.352695692.
- `reference-01.log`: all seven fresh actual NewForecast runs and all 39 bins
  completed. `reference-historical-comparison.json` records exactly zero change
  in sigma_ap, sigma_at and correlation relative to the accepted historical full
  reference. Historical full/quick offline checks passed with their schema-1
  limitations explicitly retained.
- `preservation-final.json`: 43 protected files, the exact r1 handoff, all 48
  source/wheel modules and 75 input hashes checked. `git diff --check` passed.

Recorded task times are in `timings-final.json`: the seven fresh reference
calls total 814.43 s. The final bundle records 4230.69 s for primary records
(including 29 reassemblies) and 6148.12 s for diagnostics. Earlier interrupted
attempts retained 5 records / 325.74 s and 29 records / 1961.60 s, respectively;
the latter supplied the 29 verified cached inputs. These are recorded task-time
sums, not a claim of total end-to-end wall time.

The final strict offline check and full gate both exit 1 with “partial or
scientifically failed validation”, after matching the exact 78/468 inventory.
The plot-data scan checks all 500 available numerical payloads successfully;
`plot-check-verified.log` passes all 720 plotted data series across 78 profiles
and 37 files. The initial plot checker failed on Matplotlib object-array data;
using converted numeric artist data fixed the checker, with the failed log retained.
`audit-final.json` passes. All 39 detailed attribution chains and both branches
per bin pass reconstruction; saved attribution replotting also completed offline.

The independent-bin comparison checks all 14 case/profile combinations, preserving
separate per-bin parameters and their marginal errors. Compatibility agrees with
reference blocks; only the galaxy-only accuracy case has all bins converged.
See `combined-bins/manifest.json` and its saved matrices. No shared-parameter
sum was substituted for the assigned independent-bin forecast.


## Scientific definitions

**Compatibility** uses actual captured upstream observed means, total powers and
Jacobians, including endpoint-inclusive k sums, ten mu midpoints, polynomial-8
peak extraction, backward dlogk derivative with the first sample zero, legacy
pair-specific inputs and original reconstruction rules. Mean redshift is the
arithmetic centre, while covariance uses the geometric centre. FishHighz forms
C/F independently; upstream final F is only a comparator. Literal field PSD
and selected-covariance conditioning are audited separately. Production field-PSD
validation was not relaxed.

**Accuracy** uses the accepted `prepare_bin`/`run_bin` and `KaiserModel` engines,
fixed observed cuts [0.01,0.5] h/Mpc, independent per-bin wiggle-only ap/at and
identity smooth scaling. Mean/noise preparation uses the geometric evaluation
redshift; plots retain arithmetic centres and original widths. CAMB is prepared
once at every actual model redshift, including the template and damping-reference
redshifts (template 2.406; damping reference 2.3). Power growth is the squared
sigma8 ratio. The configured galaxy reconstruction factor is 2.0; forests use
1.0. Biases, beta, f, auto widths,
noise, response, weights, volume and cuts are held fixed during differentiation.
The derivative includes inverse coordinates, Q, remapped RSD and damping.

The template uses its existing signed `PK-PKSB` wiggle and `PKSB` smooth columns.
The available template/background cosmological metadata comparison is saved in
`template-background.json`: H0, baryon/CDM/neutrino densities, curvature, w, TCMB
and ns match exactly. Current CAMB sigma8 at the template redshift is
0.3003665595810009 versus FITS 0.300381905975677, a relative difference
-5.1089610828114473e-5. The header does not encode complete solver/primordial
provenance, so this is a limited metadata comparison.
No template decomposition or amplitude was retuned. Resolving power is converted
to Gaussian sigma with `c/(R*2*sqrt(2*ln(2)))`, c=299792.458 km/s. Supplied pixel
and Poisson noise are not response multiplied; aliasing carries W squared.
Damping crosses use the mean squared auto widths. Sampling noise is independently
prepared per observed field; no sample-overlap noise was invented.

The magnitude domain is fixed and partitioned at spline/support/SNR boundaries
and density zeros. All raw/effective density, variance, query masks and quadrature
samples are retained. Geometry, k, mu, magnitude, finite-difference and weighting
controls and their actual trial F/pair-F arrays are saved. Initial levels and the
at-most-two extra quadrature refinements follow the assignment. Weights use
3/6/12 updates, with 24 tested when needed; no replacement optimizer or relaxed
arithmetic guard was introduced.

The raw first-spacing convention, negative-density extension and SNR floors/clamps
remain explicit assumptions. Floor variants, support-removal controls and artificial
±10% widths are recomputed consistently; inactive policies have zero-mask proof.
Width diagnostics scale evaluated densities at fixed normalized counts, including
the negligible floor contribution. They are not measured cell widths or systematic
error bars. Per-bin/per-field contributions and signed diagnostic changes are in
the audit and CSV tables.

## Full primary comparison

All 78 primary records contain finite, independently reconstructed information.
Compatibility passes in all 39 bins: maximum relative F discrepancy
1.9564536459043957e-15, combined-error discrepancy 1.4432899320127035e-15,
and selected-pair error discrepancy 1.887379141862766e-15. Accuracy passes the
full convergence criterion in 5/39 bins. The other 34 are **unresolved**, even
where a finite curve differs little from the reference.

The table gives signed ranges of accuracy/reference uncertainty changes over
every original bin and the largest absolute correlation difference. Accuracy/
compatibility differences are equal at the displayed precision; complete tables
retain all three comparisons without rounding. Positive percentages mean larger
forecast errors. Every case has a material flag, including the galaxy-only case
through correlation (and its highest-bin parallel error).

| Case | Compatibility passed | Accuracy converged | Δσ ap (%) | Δσ at (%) | max absolute Δcorr |
|---|---:|---:|---:|---:|---:|
| `lbg_lae_3x2pt` | 5/5 | 5/5 | -1.114 to -0.186 | -0.612 to +0.413 | 0.01584 |
| `lya_lbg_lae_3x2pt` | 5/5 | 0/5 | +7.908 to +15.752 | +18.490 to +23.945 | 0.01896 |
| `lya_lbg_lae_6x2pt` | 5/5 | 0/5 | +0.216 to +4.855 | +1.209 to +6.470 | 0.01925 |
| `lya_qso_2x2pt` | 6/6 | 0/6 | -4.474 to -1.205 | -3.056 to +0.158 | 0.02280 |
| `lya_qso_lbg_lae_15x2pt` | 6/6 | 0/6 | -3.307 to +2.365 | -0.610 to +3.318 | 0.01814 |
| `lya_qso_lbg_lae_4x2pt` | 6/6 | 0/6 | -7.248 to -3.901 | -3.503 to -0.769 | 0.01983 |
| `lya_qso_lbg_lae_8x2pt` | 6/6 | 0/6 | -6.897 to +6.285 | -2.420 to +12.563 | 0.01783 |

**Convergence:** k, mu, volume and derivative-step checks pass in all 39 bins.
Weight stabilization fails in all 34 forest-containing bins; magnitude checks
also fail in 23 bins. Across all bins, the final-two maxima are:

| Control | Failed bins | relative F | combined error | selected-pair error | relative volume |
|---|---:|---:|---:|---:|---:|
| k | 0 | 1.07566e-05 | 3.53635e-06 | 5.01253e-06 | 0 |
| mu | 0 | 1.01913e-06 | 7.76037e-07 | 1.1481e-06 | 0 |
| magnitude | 23 | 0.00153344 | 0.000809159 | 0.00336268 | 0 |
| volume | 0 | 2.09918e-11 | 1.04965e-11 | 1.04974e-11 | 2.09917e-11 |
| step | 0 | 1.4145e-05 | 9.06593e-06 | 8.37277e-06 | 0 |
| weights | 34 | 0.523135 | 0.200678 | 0.385003 | 0 |
| combined | 34 | 0.524084 | 0.200947 | 0.385924 | 2.09917e-11 |

All primary final controls use 128 k intervals (order 4), mu order 32, volume
order 32 and step 2.5e-4. Magnitude orders are 16/32/64 and retained finite weight
counts are 6 or 12; the galaxy-only bins use order 16 and 12 updates. Complete
per-bin choices, attempted levels, failures and trial matrices are retained in
the primary reports and `audit-final.json`. Several higher-count trials fail
underflow or overflow. The largest final-two weight-study combined-error change
is 20.07%, and the selected-pair change reaches 38.50%. No convergence claim can
be inferred from a finite six- or twelve-update output.

### Figures and numerical tables

| Case | Main PNG / PDF | Selected spectra PDF | Sensitivity / detailed attribution |
|---|---|---|---|
| `lbg_lae_3x2pt` | PNG (local-only path: `../.validation/step12-r2-implementation/plots-final/lbg_lae_3x2pt.png`), PDF (local-only path: `../.validation/step12-r2-implementation/plots-final/lbg_lae_3x2pt.pdf`) | Pairs (local-only path: `../.validation/step12-r2-implementation/plots-final/lbg_lae_3x2pt-pairs.pdf`) | Sensitivity (local-only path: `../.validation/step12-r2-implementation/sensitivities-final/lbg_lae_3x2pt.png`), Attribution (local-only path: `../.validation/step12-r2-implementation/attribution-final/lbg_lae_3x2pt.png`) |
| `lya_lbg_lae_3x2pt` | PNG (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_lbg_lae_3x2pt.png`), PDF (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_lbg_lae_3x2pt.pdf`) | Pairs (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_lbg_lae_3x2pt-pairs.pdf`) | Sensitivity (local-only path: `../.validation/step12-r2-implementation/sensitivities-final/lya_lbg_lae_3x2pt.png`), Attribution (local-only path: `../.validation/step12-r2-implementation/attribution-final/lya_lbg_lae_3x2pt.png`) |
| `lya_lbg_lae_6x2pt` | PNG (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_lbg_lae_6x2pt.png`), PDF (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_lbg_lae_6x2pt.pdf`) | Pairs (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_lbg_lae_6x2pt-pairs.pdf`) | Sensitivity (local-only path: `../.validation/step12-r2-implementation/sensitivities-final/lya_lbg_lae_6x2pt.png`), Attribution (local-only path: `../.validation/step12-r2-implementation/attribution-final/lya_lbg_lae_6x2pt.png`) |
| `lya_qso_2x2pt` | PNG (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_2x2pt.png`), PDF (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_2x2pt.pdf`) | Pairs (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_2x2pt-pairs.pdf`) | Sensitivity (local-only path: `../.validation/step12-r2-implementation/sensitivities-final/lya_qso_2x2pt.png`), Attribution (local-only path: `../.validation/step12-r2-implementation/attribution-final/lya_qso_2x2pt.png`) |
| `lya_qso_lbg_lae_15x2pt` | PNG (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_lbg_lae_15x2pt.png`), PDF (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_lbg_lae_15x2pt.pdf`) | Pairs (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_lbg_lae_15x2pt-pairs.pdf`) | Sensitivity (local-only path: `../.validation/step12-r2-implementation/sensitivities-final/lya_qso_lbg_lae_15x2pt.png`), Attribution (local-only path: `../.validation/step12-r2-implementation/attribution-final/lya_qso_lbg_lae_15x2pt.png`) |
| `lya_qso_lbg_lae_4x2pt` | PNG (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_lbg_lae_4x2pt.png`), PDF (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_lbg_lae_4x2pt.pdf`) | Pairs (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_lbg_lae_4x2pt-pairs.pdf`) | Sensitivity (local-only path: `../.validation/step12-r2-implementation/sensitivities-final/lya_qso_lbg_lae_4x2pt.png`), Attribution (local-only path: `../.validation/step12-r2-implementation/attribution-final/lya_qso_lbg_lae_4x2pt.png`) |
| `lya_qso_lbg_lae_8x2pt` | PNG (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_lbg_lae_8x2pt.png`), PDF (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_lbg_lae_8x2pt.pdf`) | Pairs (local-only path: `../.validation/step12-r2-implementation/plots-final/lya_qso_lbg_lae_8x2pt-pairs.pdf`) | Sensitivity (local-only path: `../.validation/step12-r2-implementation/sensitivities-final/lya_qso_lbg_lae_8x2pt.png`), Attribution (local-only path: `../.validation/step12-r2-implementation/attribution-final/lya_qso_lbg_lae_8x2pt.png`) |

The all-case summary (local-only path: `../.validation/step12-r2-implementation/plots-final/all-cases.png`) marks unresolved
accuracy bins. `plots-final/tables.json` and its CSV contain all 39 combined-bin
and 234 selected-pair-bin results, including null/availability and pass labels.
`all-comparisons.csv` retains both profiles/reference and accuracy/compatibility
differences; `pair-convergence.csv` records every selected spectrum/control;
`sensitivity-comparisons.csv` retains signed diagnostic changes. All values are
linked to source hashes, original bin bounds and arithmetic plot centres.
`audit-final.json` includes normalized field-PSD checks, selected covariance/F
conditioning, raw-policy counts and effective density contributions.

The field audit finds **zero material PSD violations** in either profile.
Minimum normalized field eigenvalues are 0.13415 (compatibility) and 0.13854
(accuracy). Maximum selected-covariance condition numbers are 1.39450e7 and
1.16413e7; maximum Fisher condition numbers are 3.58388 and 3.38806.
The actual legacy auto-power minimum is positive (32.5343). Thus an actual
legacy PSD failure is not an explanation for the observed differences.

All seven primary PNGs were visually inspected and are byte-identical to their
final versions; the selected-spectrum PDFs contain all 41 required pages.
The all-case overview and representative detailed attribution/sensitivity figures
were inspected. The corrected `sensitivities-final` figures explicitly break lines
at failed variants, with 422 finite effects and 46 unavailable gaps independently
checked. They supersede the original supplementary diagnostic renderings, which
could bridge failed bins. `visual-inspection-final.json` records the inspected files.

## Attribution and limitations

Every available bin receives an endpoint-reproducing common-node array chain.
The detailed chain separates volume, diagonal sampling/aliasing/pixel noise,
observed signal, template-model/polynomial peak, signed PK-PKSB peak, full wiggle
mapping derivative and final quadrature. The interpolation residual is measured;
these sequential effects depend on order and are not independent additive errors.
Separate branches isolate the exact legacy forest/galaxy cross-damping rule and
negative-density flooring within the literal three-update reference recipe.
The latter reproduces captured legacy weights/coefficient arrays first; it is a
validation-only signed-array diagnostic, not a change to strict readers.

One-setting growth, evaluation-redshift, resolution, reconstruction,
rectangular-magnitude and three-update diagnostics supplement the input-policy
studies in every bin. Derived weights/noise are recomputed consistently for these
settings. The detailed audit also reconstructs primary intrinsic power, response,
observed J, total power and final per-field weights/noise from the recorded fixed
model and samples. All stage and branch matrices can be checked offline from
saved T/J/modes; all plotting values are bound to numerical tables and files.

## Diagnostic completion and retained-input limits

All 468 diagnostic requests were attempted: **422 completed; 46 failed**. Every
failure is `underflow encountered in multiply` in the existing weighting
integrals/coefficients. These are unavailable comparisons, not zero effects.
Successful diagnostics establish finite, internally consistent calculations at
the recorded controls; they do not repair the primary convergence failures.

| Case | Completed / requested | Unavailable variants and zero-based bins |
|---|---:|---|
| `lbg_lae_3x2pt` | 60/60 | None |
| `lya_lbg_lae_3x2pt` | 55/60 | floor_low: 0,1,2,3,4 |
| `lya_lbg_lae_6x2pt` | 55/60 | floor_low: 0,1,2,3,4 |
| `lya_qso_2x2pt` | 66/72 | rectangular_magnitude: 0,2,3,5; floor_low: 1,4 |
| `lya_qso_lbg_lae_15x2pt` | 60/72 | floor_low: 0,1,2,3,4,5; remove_bright: 0,1,2,3,4; width_plus: 2 |
| `lya_qso_lbg_lae_4x2pt` | 66/72 | rectangular_magnitude: 0,2,3,5; floor_low: 1,4 |
| `lya_qso_lbg_lae_8x2pt` | 60/72 | floor_low: 0,1,2,3,4,5; remove_bright: 0,1,2,3,4; width_plus: 2 |

The failures comprise 26 `floor_low` (1e-22), 8 `rectangular_magnitude`,
10 `remove_bright`, and 2 `width_plus` variants. Central 1e-20 policies stay
unchanged. `coverage-final.json` and `diagnostic-coverage-final.json` retain the
exact errors and inventory. All zero-mask proofs, successful variants and their
signed combined/selected-pair differences remain in the numerical tables.

Successful diagnostic ranges below are signed combined uncertainty changes
relative to the retained primary accuracy result, over both parameters/all bins.
They are conditional on the recorded, often unconverged controls; failed variants
are excluded and remain explicitly unavailable.

| Variant | Minimum change (%) | Maximum change (%) |
|---|---:|---:|
| `floor_low` | +0 | +0 |
| `floor_high` | +0 | +0 |
| `remove_bright` | +0 | +7.34914 |
| `remove_sentinel` | +0 | +0 |
| `width_minus` | -10.1076 | -3.9009 |
| `width_plus` | +3.84041 | +10.3384 |
| `growth_eds` | -0.390701 | +1.13835 |
| `legacy_resolution` | +0 | +0.700252 |
| `rectangular_magnitude` | -7.08046 | +0.345898 |
| `no_reconstruction` | +2.02288 | +18.6235 |
| `arithmetic_mean` | -0.162555 | +0.0195095 |
| `three_weights` | -20.0678 | +0 |

Zero floor/sentinel effects in completed variants do not bound the missing
floor-low comparisons or establish physical support outside the tables.
The artificial width changes and three-update weighting shifts are substantially
larger than the formal quadrature tolerance. Reconstruction is also a material
model choice; these effects are not independent systematic error bars.

### Measured attribution

`attribution-summary.json` and `attribution-stages.csv` retain every bin/stage,
including signed errors, correlation and Fisher changes. Ranking by the largest
absolute combined-error step gives:

| Case | Largest sequential effect (number of bins) | Range across both parameters/bins for that stage (%) |
|---|---|---|
| galaxy 3x2pt | full wiggle mapping derivative (5) | +1.199 to +1.530 |
| forest/galaxy 3x2pt | per-field noise (5) | +17.461 to +29.548 |
| forest/galaxy 6x2pt | noise (3), template model with polynomial peak (2) | noise +1.058 to +7.711; model -3.606 to -1.094 |
| forest/QSO 2x2pt | model (5), noise (1) | model -6.339 to -3.189; noise +0.431 to +4.973 |
| 15x2pt | model (5), noise (1) | model -4.640 to -1.763; noise +0.256 to +5.744 |
| 4x2pt | model (6) | -9.066 to -3.990 |
| 8x2pt | model (3), noise (3) | model -9.331 to -3.976; noise +0.819 to +16.881 |

These ordered replacements explain the finite endpoints, with maximum relative
endpoint residual 4.46102e-15. They do not establish a unique physical attribution:
the model/polynomial stage groups template and model conventions, and the final
quadrature stage also retains a measured Jacobian interpolation discrepancy of
0.004269–0.004498. That interpolation residual is not silently assigned to physics.
The exact legacy cross-width branch changes combined errors by up to 8.54952%;
negative-density flooring within the separately reproduced legacy recipe changes
them by up to 0.185758%. These branches have different baselines and must not be
added to the main chain.

Primary intrinsic power, response, observed Jacobian and total-power reconstruction
have zero relative discrepancy in all 39 bins. Reconstructed forest noise agrees
to at most 1.55654e-16. Consequently the saved accuracy outputs agree with the
stated implementation; the remaining scientific obstacle is weight stabilization
(and magnitude convergence in 23 bins), alongside the documented input/model
assumptions. A small reference/accuracy difference in some bins does not resolve
those failures. Further changes to weighting or input policy require review of
this evidence, not an automatic extra refinement or changed tolerance.

The results quantify numerical reproducibility within the specified model. They
do not establish calibrated coverage outside source tables, physical raw cell
widths, sample overlap, mode mixing, redshift errors, new P1D physics or within-bin
signal/density evolution. No equivalence to a complete realistic survey model or
joint optimality of the cumulative weighting rule is claimed. Remaining findings
are for independent/user review; Step 12 is not accepted and Step 13 was not started.
