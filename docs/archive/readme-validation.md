# Historical README validation record

This archive preserves the earlier README account. Historical instructions do
not authorize new calculations or change the adopted research prescription.

## Historical DESI-2 validation and explicit legacy compatibility (Step 12)

This section records the earlier Step-12 comparison and remains available for
reproduction. Its fixed-reference `inverse_variance` accuracy profile is not the
current research baseline. The current three-profile identity and selected
S2--S4 prescription are given in [RESEARCH_BASELINE.md](../research/RESEARCH_BASELINE.md).

`python examples/desi2_synthetic.py` runs all seven explicit selections on tiny
synthetic arrays, checking their correlated Fisher matrices against an independent
five-field covariance. It requires only NumPy. `fishhighz.validation.cases`
provides `CASE_IDS`, `recipe(case)`, `selection(case)`, `bins(case)` and
`verify_inventory(directory)`. Recipes preserve original tracer indices and
reference upper-triangle pair order, including covariance-required unselected
spectra. Every original option is checked; this is not a general INI translator.

Strict `DensityReader` and `SNRReader` remain the default. Separately construct
`LegacyDensity(reader, negative_policy="reject" | "floor_negative")` and
`LegacySNR(reader)` from validated tables. Their `sample` methods return immutable
`values`, `raw` and `provenance` records. `sample_legacy_forest(...)` converts
these into normalized arrays and plain provenance for `ForestInput`:

- Density strictly outside the raw magnitude domain becomes **1e-20**; closed
  endpoints remain inside. Redshift uses the reference spline domain extension,
  with original and effective coordinates recorded. This is not validated
  physical extrapolation. Exact zeros and positive densities below the floor
  remain unchanged. Negative in-domain interpolants fail with `reject`;
  `floor_negative` replaces only negatives by **1e-20**, an explicit extension
  beyond legacy behavior. Compatibility preserves the reference normalization
  reduction order, including zero-masked forest rows.
- SNR clamps bright magnitudes to the first node. Faint magnitudes or source
  redshift/wavelength outside closed domains return variance **1e20**, independently
  of pixel/exposure scaling. Inside, the **1e-10 SNR floor** applies after
  `sqrt(pixel_width_angstrom)*sqrt(exposure_count/file_exposure_count)` scaling.
  Exposure time must match the file. Arithmetic failures and invalid inputs fail.
- Diagnostics include source hashes, normalization/units, original/effective
  queries, raw values, reason masks/counts and constants. Out-of-range SNR raw
  placeholders are zero and explicitly marked unevaluated. Diagnostics are
  snapshots, with no counters updated during derivatives. Each forest owns its
  own source, SNR, response, weights and independent P1D.

Revision 2 historically provided two explicit validation profiles for the seven
original cases (39 bins, 78 primary records). They share observed cuts and
selected spectra, but represent different estimators:

- **Maximum compatibility** captures actual upstream mean/total powers and BAO
  derivatives, then independently constructs Wick covariance and Fisher
  information in FishHighz. It preserves legacy endpoint sums, polynomial peak
  extraction, backward derivative, pair-specific noise and redshift conventions.
  Literal field-matrix PSD failures are reported separately; selected covariance
  must still be positive definite. This path is validation-only.
- **Historical maximum accuracy** uses the accepted `prepare_bin`/`run_bin` and `KaiserModel`
  engines: geometric evaluation redshift, integrated volume, Gauss–Legendre k/mu
  and composite magnitude quadrature, CAMB growth at actual redshifts, physical
  FWHM resolution, accepted auto/cross damping and full wiggle-mapping derivatives.
  Each forest now uses fixed inverse-variance weights at the declared
  `q_star=0.00035 s/km`: `B_star=P1D(z_eval,q_star)*W_field(q_star)^2` and
  `nu=B_star/(B_star+l_p*v)`. The intrinsic P1D, field response and `B_star` are
  recorded per forest; no P3D signal or auxiliary coordinates enter this choice.
  Noise and weights are fixed during differentiation. An explicitly selected
  `accuracy_method="legacy"` retains cumulative compatibility/replay behavior.
  “Maximum” is bounded by this existing model and the prescribed convergence
  tests; the fixed mode is a convention, not a fitted or optimized scale.

The default is the 15x2pt bin [2.47,2.705]. Full execution must be explicitly
selected and authorized. The interpreter must already provide the requested
reference checkout, CAMB, and FishHighz optional extras. The scripts never modify
reference code or install dependencies automatically. Use fresh output paths:

```bash
python scripts/desi2_validation.py preflight --reference /path/to/lyaforecast
# Separately authorized fresh reference capture:
python scripts/capture_desi2_reference.py --suite full \
  --reference /path/to/lyaforecast --output /path/to/new-reference
# Both profiles; omit --suite full for the one-bin quick default:
python scripts/desi2_profiles.py run --suite full \
  --reference /path/to/lyaforecast --template /path/to/template.fits \
  --reference-bundle /path/to/new-reference --wheel /path/to/exact.whl \
  --output /path/to/new-profiles
python scripts/desi2_profiles.py check --output /path/to/new-profiles
# Strict full-assignment inventory and scientific gate:
python scripts/check_desi2_full.py /path/to/new-profiles
# Offline tables/plots, including failed convergence (Matplotlib optional):
python scripts/plot_desi2.py --bundle /path/to/new-profiles \
  --output /path/to/new-plots
# Sensitivity panels with explicit gaps and labels for failed variants:
python scripts/plot_desi2_sensitivities.py --plots /path/to/new-plots \
  --output /path/to/new-sensitivity-plots
python scripts/audit_desi2.py --bundle /path/to/new-profiles \
  --output /path/to/new-audit.json
python scripts/combine_desi2_bins.py --bundle /path/to/new-profiles \
  --output /path/to/new-independent-bin-check
python scripts/attribute_desi2.py --bundle /path/to/new-profiles \
  --template /path/to/template.fits --output /path/to/new-attribution
# Verify saved attribution matrices and redraw without model evaluation:
python scripts/attribute_desi2.py --replot /path/to/new-attribution \
  --output /path/to/redrawn-attribution
python scripts/check_desi2_plots.py --plots /path/to/new-plots \
  --reference /path/to/new-reference
```

`--profiles`, `--cases` and `--bins` select explicit subsets; subsets never count
as the full 78-record comparison. The Python entry point is
`fishhighz.validation.profiles.run(...)`. `wheel_identity` verifies every imported
module against the exact supplied artifact; use an installed copy of that wheel
or explicitly place that wheel on the interpreter's import path.

Schema 2 binds case/profile/bin/field/pair/parameter identities and reconstructs
nodes, Wick covariance, Fisher information and rank-aware errors offline. It
recomputes recorded convergence metrics and requires boolean scientific success.
Execution finished and scientific acceptance are distinct. Real reports live in
individually hashed JSON files, keeping the manifest bounded as diagnostics grow.
After an interrupted validation run, `--reuse-completed /path/to/bundle` can reuse
immutable numerical inputs only when every scientific/schema module and the
actual reference interpreter, source and resources match. Report-I/O/controller
changes are the only permitted code differences. C/F are reassembled with the
current exact wheel; both the original producer and reuse provenance are retained. Unconstrained errors
use explicit availability masks; plots/tables use nulls. Historical schema-1
bundles can only be inspected with `evidence.inspect_legacy`, which returns a
limited result and cannot pass this gate. No pickle or external model imports
are needed to read evidence.

Fixed-reference accuracy studies use trial-contract version 2 and vary k, mu,
volume, magnitude and derivative step independently, then test their combined
lower controls. Iterations and a `weights` convergence metric are explicitly
inapplicable. Historical version-1 evidence retains the six-control cumulative
schedule; missing method metadata means legacy, never fixed-reference. Cache
reuse checks the requested method, and `three_weights` remains a legacy-only
diagnostic. The existing cumulative rule may underflow at higher counts; such
historical failures and unresolved single-pair convergence remain findings,
without altered tolerances or a new optimizer. Floors/clamps and
`legacy_first_spacing` remain labeled input policies.
Floor changes, support removal and artificial ±10% cell widths are diagnostics,
not calibrated systematic uncertainties. Additional one-setting physical swaps
and common-node array chains record order dependence and unresolved grouped
causes. Significant unexplained differences prevent an unqualified pass.

The historical r1 one-bin study used 107 rectangular magnitude nodes, three weight
updates, EdS power growth and legacy c/R resolution. It is preserved separately
and is not the revision-2 accuracy result. All six examples and the explicit
legacy adapters remain available.

The callable `IntrinsicP3D(provider, *, routes, n_fields, h_source, h_fid,
k_domain, z_domain)` wraps a caller-prepared object's **intrinsic**
`compute_p3d_hmpc(z,k,mu,corr)`. It validates all paired queries and explicit routes
before calling a provider that might clamp. It preserves signed crosses and
arbitrary pair order, maps `k_source=k_fid*h_fid/h_source` and
`P_fid=P_source*(h_fid/h_source)**3`, and adds no response/noise. It has zero
local parameters; the real validation uses a separate `A*P_external` wrapper
and independent matrix amplitude-Fisher oracle. That amplitude example is not
the primary BAO forecast; P1D remains independently chosen.

Run the bounded actual external-model example in the identified reference
environment with the exact FishHighz wheel imported:

```bash
python scripts/external_desi2.py --reference /path/to/lyaforecast \
  --template /path/to/template.fits --wheel /path/to/exact.whl \
  --output /path/to/new-external-amplitude
```

This caller constructs the actual reference `PowerSpectrum` once and saves a
separately typed, independently reconstructed amplitude record. It does not run
`NewForecast` or replace the full BAO comparison.

The native `Forecast` interface provides configuration, CLI and serialization
independently of these historical validation adapters. General lyaforecast INI
translation is not implemented. The closed Step 13 study provides a
validation-only amplitude/shape diagnosis of the cumulative weighting
recurrence; it is not a production weighting mode. See
`reviews/step-13.md` for its bounded scientific conclusion and
`reviews/step-12.md` for the preceding forecast evidence.

### Step 12 revision 5: six-bin 15x2pt diagnosis

Revision 5 retains the revision-4 NumPy default and optional compiled contraction.
New evidence uses schema 3: every accuracy metric references distinct actual
refinement trials, including their controls, combined and individual-spectrum
Fisher matrices, volumes and failed outcomes. The offline checker replays the
bounded control sequence from saved summaries without evaluating a model.
Schema-2 bundles remain historical, with explicitly limited convergence checks;
they cannot satisfy the new 15x2pt acceptance gate.

Run only the explicitly selected case for this assignment, using the exact newly
built wheel in an isolated environment and the preserved upstream reference:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export FISHHIGHZ_FISHER_BACKEND=numba # optional fishhighz[compiled]
python scripts/desi2_profiles.py run --suite full \
  --cases lya_qso_lbg_lae_15x2pt --bins 0 1 2 3 4 5 \
  --reference ../lyaforecast \
  --template ../vega/vega/models/Planck18/Planck18_z_2.406.fits \
  --reference-bundle .validation/step12-r2-implementation/reference-01 \
  --wheel /absolute/path/to/exact.whl --output /new/exclusive/bundle
python scripts/check_desi2_15x2pt.py /new/exclusive/bundle
python scripts/plot_desi2.py --bundle /new/exclusive/bundle --output /new/plots
python scripts/check_desi2_plots.py --plots /new/plots \
  --reference .validation/step12-r2-implementation/reference-01 \
  --cases lya_qso_lbg_lae_15x2pt
```

The scoped gate requires 12 primary records, 180 individually selected-spectrum
results and 72 diagnostic requests. The historical seven-case gate is separate.
Plots retain failed convergence and unavailable diagnostics. Compatibility
reconstructs the legacy estimator from saved upstream powers and Jacobians;
accuracy uses the existing physical-response, growth and full-wiggle-derivative
prescription. Neither profile introduces sample-overlap noise or calibrated
coverage beyond the retained input policies. Cumulative weighting stabilization
remains a scientific limitation; a changed prescription requires user review.

The cumulative-update diagnosis can be reproduced from the saved six-bin
samples without a CAMB solve:

```bash
python scripts/recheck_desi2_15x2pt.py --source /historical/schema2/bundle \
  --output /new/derived/evidence
python scripts/diagnose_desi2_weights.py --bundle /six-bin/bundle \
  --reference ../lyaforecast \
  --template ../vega/vega/models/Planck18/Planck18_z_2.406.fits \
  --output /new/weight-diagnosis
python scripts/plot_desi2_weight_diagnosis.py --source /new/weight-diagnosis \
  --output /new/weight-figures
```

This diagnostic holds the initial weight function fixed while refining magnitude
integration, then compares the retained finite cumulative updates. It does not
replace the primary weights. Decimal products in the diagnosis are evaluated to
80 significant digits to identify unrepresentable operations; they are not used
in forecast calculations. Historical failure controls reconstructed by exact
controller replay are explicitly distinguished from newly recorded attempts.

For the scoped independent-bin check and controlled comparisons:

```bash
python scripts/combine_desi2_bins.py --bundle /six-bin/bundle \
  --output /new/independent-bins --cases lya_qso_lbg_lae_15x2pt
python scripts/attribute_desi2.py --bundle /six-bin/bundle \
  --template ../vega/vega/models/Planck18/Planck18_z_2.406.fits \
  --output /new/attribution
python scripts/plot_desi2_sensitivities.py --plots /new/plots \
  --output /new/sensitivity-figures
```

After resuming a run that reassembles primary Fisher matrices, bind saved
sensitivity comparisons to those final matrices before plotting:

```bash
python scripts/bind_desi2_sensitivities.py --source /resumed/bundle \
  --output /new/checked-comparisons
```

This derives a new bundle, preserving numerical arrays and prior reports. It
recomputes comparison metrics rather than relaxing checks for nearly zero effects.

### Step 12 revision 6: weak-spectrum evidence checks

Revision 6 compares each individual-spectrum Fisher matrix and derived summary
on its own scale, including primary/final-trial and replayed operands. Scaling
before subtraction and norm evaluation avoids range loss for very weak or strong
information. The checker reconstructs convergence metrics and the verdict from
saved study arrays; detached operands cannot certify convergence. Evidence
consistency remains 5e-12, separately from the unchanged scientific thresholds.

This repair uses synthetic regressions and offline inspection only. The preceding
real-run commands document revision 5 and require explicit authorization for any
new execution. Saved revision-5 compatibility passes in six bins; all six accuracy
records remain unconverged, with 12 of 72 diagnostics unavailable. Scientific
weighting decisions and Step 12 acceptance remain open. See the
[current handoff](reviews/step-12.md) and [preserved revision-5 handoff](reviews/step-12-r5.md).

### Step 13 revision 2: cumulative-weight limit diagnosis

`fishhighz.validation.weight_limit` evolves the exact nonlinear recurrence in a
range-safe log-amplitude/log-shape representation. It records ordinary-value
availability, the fixed-grid linearized spectrum, the normalized sightline
distribution, tail measure and coefficient concentration without changing
`prepare_forest_weights` or its underflow guards. The explicit controller
`scripts/diagnose_weight_limit.py` reads only preserved Step 12 arrays and binds
every population, quadrature order and attempted checkpoint to those immutable
inputs. `--check-finalized` independently reconstructs all final tables and
figures before accepting a saved bundle.

The five available magnitude orders show strong, consistently positive growth of
the fixed-grid asymptotic pixel-noise coefficient, but they are finite and
non-nested. They do not establish the required non-atomic continuum measure,
uniform nonlinear remainder control or an exchange of iteration and refinement
limits. Revision 2 therefore classifies all twelve saved forest populations as
unresolved in the continuum limit. No forecast was reassembled and no production
prescription was adopted; Step 12 scientific acceptance remains open.
