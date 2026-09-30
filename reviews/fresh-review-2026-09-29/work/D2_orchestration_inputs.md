# D2: orchestration and survey-input layer (implementation audit)

Scope: survey_config.py, public.py, forecast.py, results.py, accuracy.py, cli.py, parameters.py,
survey.py, magnitude.py, resources.py, adapters/legacy_{inputs,compat}.py, bundled INI and data
README, docs/user/*.md, docs/research/RESEARCH_BASELINE.md. Reference: lyaforecast (tracer.py,
spectrograph.py, weights.py, covariance.py, fisher.py, forecast_new.py), read-only.

Method: code reading plus preparation-level calls on the bundled inputs (readers, prepare_survey,
prepare_bin for all six bins, _pair_spec preparation for selected pairs in bins 1, 3, 6).
No Forecast.run(), run_bin, or pytest was executed. Scratch scripts: /tmp/claude-81394/.../scratchpad/D2/insp*.py.
The seed map (fishhighz/reviews/.../map_fishhighz.md) does not exist at the given path; the
reviews/ tree did not exist, so no seed line refs were used. All line refs below are from this
audit.

## Summary

No implementation error that changes the forecast numbers was found in the orchestration or
survey-input layer. Inputs (z bins, z_eval, z_source, n_bar, rho, SNR variance, L, pixel width,
volume, a_v, d_deg) reproduce independent recomputations to <=1e-9 (except the 0.05% volume
sample noted below, which is a check-tolerance effect, not a code difference). One
scientifically material input-coverage issue (M1) and one confirmed CLI defect (m1) are
reported. All recipe values in accuracy.py/INI_DEFAULTS agree with RESEARCH_BASELINE.md except the
rtol statement (m4).

## Findings

### M1 (major, input coverage; behaviour inherited from lyaforecast, silent in outputs): lya(lbg) source density is extrapolated at constant value in bins 5 and 6

- Location: survey_config.py:1409-1425 (sample at z_source), adapters/legacy_compat.py:110-131
  (`redshift="legacy_spline_extension"`), data/lbg_matched_dndzdr.txt.
- Defect: the LBG/LAE table spans z = 2.38, 2.60, 2.83, 3.07, 3.29 only. The Lya(LBG) background-source
  density is queried at z_source = (1+z_eval)*1215.67/sqrt(1040*1205) - 1, which is 3.404 (bin 5) and 3.660
  (bin 6). Both exceed 3.29. RectBivariateSpline with bbox clamps, so dN/dz is the 3.29-row value in both
  bins (total 113.5 deg^-2 per unit z, versus 524.9 at z_src=3.149 in bin 4). `redshift_extension` is
  True in the per-node provenance masks (verified for bins 5, 6) but nothing is warned, and the saved
  settings.json does not carry it (only ini identity/hash; ForestInput.provenance stays in the prepared state).
- Evidence (prepare_bin output, A = aliasing coefficient in deg^2, smaller = more forest information):

  | bin | z_src | lya(qso) dN/dz, A | lya(lbg) dN/dz, A, extended |
  |---|---|---|---|
  | 4 | 3.149 | 30.8, 0.0615 | 524.9, 0.00328, no |
  | 5 | 3.404 | 15.0, 0.114 | 113.5, 0.0164, yes |
  | 6 | 3.660 | 5.7, 0.274 | 113.5, 0.0142, yes |

  In bins 5-6 the LBG forest carries ~20x (bin 6) the source density of the QSO forest and lower A, so
  the joint forest information (and lya(lbg)xlbg/lae/qso crosses) is set by the extrapolated value.
  A ~ 1/rho, so a factor-2 error in the assumed LBG source density is a factor-2 change in the LBG-forest
  aliasing noise power in those bins. The table row at 3.29 already drops by 4.5x from 3.07, so a plateau
  beyond the last row is not a neutral choice. Bin-2 galaxy LBG/LAE n_bar (z_eval=2.350<2.38) are likewise
  clamped to the 2.38 row (n_bar_LBG = 8.3e-6 h^3/Mpc^3 versus 4.9e-5 in bin 3); this one is a lower-side clamp
  on a steeply rising dN/dz.
- Expected impact: bins 5-6 joint errors depend on an unsupported LBG source density. Not an
  arithmetic error; matches lyaforecast, and INI comments document "legacy spline boundary extension".
  Recommend reporting sensitivity of bin 5-6 joint sigma to rho_LBG(z_src>3.29) (e.g. 0.5x, 0x) and surfacing
  the flag in the saved output. I did not run the Fisher to quantify sigma.

### m1 (minor, confirmed): CLI `--joint-only` always fails on the default (bao) recipe

- Location: cli.py:16-21 calls `forecast.run(individuals=False)`; public.py:489-490 raises
  `ValueError("individuals=False requires full_shape mode")` for any mode not handled by the earlier
  full_shape / bao_marginalized dispatches (public.py:469-488). The seed's "cli.py:491" is not a cli.py
  line (file has 27 lines); the raising line is public.py:490.
- Evidence: `Forecast().run(individuals=False)` on the bundled INI raises exactly that message before
  prepare() (executed, instantaneous). So `fishhighz-forecast --output X --joint-only` on any BAO-mode INI
  aborts after INI parse, with no result. ini.md advertises `--joint-only` under bao_marginalized/full_shape
  only; cli.md omits the flag from its argument table and states there are no such flags.
- Impact: no wrong numbers; the flag is unusable in the default mode, and there is no scientific reason
  the BAO path cannot skip the per-spectrum loop (public.py:525-551 is separable).

### m2 (minor, documented compatibility convention): `legacy_first_spacing` inflates LBG/LAE dN/dz by up to 7%

- Location: legacy_inputs.py:153-161 (widths = dz[0]), legacy_compat.py:103 (division by widths[:,None]*dm).
- The LBG/LAE z spacings are 0.22, 0.23, 0.24, 0.22 but every row is divided by 0.22. Relative to
  central-difference cell widths (0.220, 0.225, 0.235, 0.230, 0.220) the density in rows 2.60, 2.83, 3.07
  is high by 2.3%, 6.8%, 4.5%. Integrated over cells: total dN/dz*dz is 4.8% above the 380 deg^-2 target
  (LBG) and 2.6% above 430 (LAE); both n_bar and the LBG-forest source density inherit this. QSO table
  is uniform (0.1), no effect. Matches lyaforecast; INI comments state the convention.
- Impact: few-percent overestimate of LBG/LAE n_bar and lya(lbg) rho in bins 3-4; sign is toward smaller
  errors.

### m3 (minor, documented extension): negative-density floor raises n_bar and I1 slightly

- The quadratic tensor spline goes negative in the tails (min ~ -34 per cell for QSO). floor_negative
  replaces those by 1e-20 (lyaforecast keeps them). Measured on bins 1, 3, 6: n_bar(floor)/n_bar(raw) =
  1.0044-1.0060 for qso/lbg/lae; I1(floor)/I1(raw) = 1.0004-1.0083 (largest for lya(lbg)). Negative
  roots are inserted into the magnitude partition, so the floor is integrated at the correct location.

### m4 (minor, doc/code): recipe identity says rtol=1e-4, native default is 1e-5

- RESEARCH_BASELINE.md:8-9,144-148,242-244 state the recipe identity and "public adaptive solver default"
  remain rtol=1e-4 and that callers "must record an override" to reproduce the final S2/S3 state.
  The native Forecast path already uses 1e-5 (accuracy.py:25 CONTROLS weight_rtol; INI_DEFAULTS
  weighting_rtol=1e-5, numerical.weight_rtol=1e-5), which ini.md:38-40 says explicitly. accuracy.py:10
  STOPPING rtol=1e-4 is the value carried in `accuracy_settings()` provenance and is not what the native
  run uses. Only prepare_forest_weights (weights.py:215) defaults to 1e-4. Converged weights are recorded
  in either case; effect on results is at the 1e-3 % level quoted in the baseline.

### m5 (minor, documented deviation from lyaforecast): resolving power interpreted as FWHM

- survey_config.py:1379-1381 uses resolving_power_fwhm_to_sigma(R) = c/(2.3548 R) = 50.92 km/s for R=2500.
  lyaforecast (covariance.py:177-183) uses sigma = c/R = 119.9 km/s. RESEARCH_BASELINE and the annotated INI
  state the FWHM interpretation. Effect on W^2 (k*mu/a_v, bin 2): k=0.2 h/Mpc: 0.991 (fishhighz) vs 0.950;
  k=0.5: 0.944 vs 0.725. Negligible at the BAO-scale, a few percent at the upper k limit. This is a recipe
  choice, not a bug; flagged so that comparisons to lyaforecast are not read as pure implementation checks.

### m6 (minor, fragility): run() converts ap_step with hard-coded 0.001

- public.py:493 `step_scale = ap_step/0.001` relies on the registry step at survey_config.py:1208
  (`step=0.001`). Consistent today (effective central-difference half-width 2.5e-4, verified via
  derivatives._stencil), but the two constants are unlinked.

### m7 (minor, output completeness)

- SurveyResult.save (public.py:140-194) writes per-record sigma only for individual and joint records;
  combined info is saved only as a 12x12 Fisher matrix. Since the AP parameters are separate per bin, the
  combined Fisher is block-diagonal and adds no information beyond the per-bin joint entries; results.md
  states this. Nothing wrong, but no cross-bin combination (shared alpha) is offered by the facade.

## Doc/code mismatches (all minor)

1. cli.md:15-19 argument table omits `--joint-only`; text says no such flags exist; cli.py:16-18 defines it (and it fails in bao mode, m1).
2. python.md `run()` table omits the `individuals` argument (public.py:465); only ini.md mentions it, and only for full_shape/bao_marginalized.
3. RESEARCH_BASELINE.md:256-258 says "General INI/CLI translation and result serialization also remain outside the current research-ready
   Python milestone"; INI parsing, CLI, and JSON/NPZ saving are implemented (survey_config.py, cli.py, public.py:140).
4. RESEARCH_BASELINE.md rtol statements versus native default (m4).
5. RESEARCH_BASELINE.md:8 gives the recipe id "early-lyaforecast-2026-09-18" with rtol 1e-4 while accuracy.py REVISION is identical but CONTROLS carry 1e-5 (same as m4).
6. data/README.md: inventory sizes and all 30 SHA-256 values (5 files, 24 SNR tables incl. per-magnitude hashes) verified to match. No mismatch.

## Checks passed (one line each)

- INI_DEFAULTS/CONTROLS vs baseline: k=[0.01,0.5] h_fid/Mpc, 128 intervals x GL4, mu order 32, z order 32, magnitude order 16, ap_step 2.5e-4, weight method early_lyaforecast, representative mode (2.4 deg^-1, 0.00035 s/km), min/stable/max updates 3/3/96, beta_F 1.45, Sigma_perp amplitude 3.26 with sigma8(z_eval)/sigma8(2.3)/sqrt(r), Sigma_par=(1+f)Sigma_perp, r=2 (galaxy/QSO) and 1 (forest), template z 2.406, damping ref z 2.3: all match (except m4).
- Node count: 512 k x 32 mu = 16384 nodes per bin (verified in prepared specs).
- z edges 2.0-3.41 in six 0.235 bins equal lyaforecast's linspace; z_eval = sqrt((1+zmin)(1+zmax))-1 = 2.1153, 2.3504, 2.5856, 2.8207, 3.0558, 3.2909; used consistently for geometry, CAMB sigma8/f, bias, response wavelength, galaxy density.
- Bin 1 selection = lya(qso)xlya(qso), lya(qso)xqso, qsoxqso; prepared bin 1 has fields active {lya(qso), qso}, 3 selected = 3 required pairs, only qso in `galaxies`, only lya(qso) in `forests`; bins 2-6: 15 selected/15 required, three galaxies and two forests.
- Excluded list = 12 pairs in bin 1, none in bins 2-6 (public.py:552-580 logic read; run not executed).
- Individual spectra (public.py:278-339): each gets its own PairSelection with all needed auto terms (required pairs [[0,0],[0,4],[4,4]] for lya(qso)xlya(lbg); [[0,0],[0,1],[1,1]] for lya(qso)xqso), providers filtered to required pairs, same ForestInput objects: forest weights, A, P_pixel are bit-identical to the joint bin (checked bins 1, 3, 6).
- Joint uses the full selected-spectrum covariance (forecast.py:133-140, gaussian_covariance over selection); reduction via `fix_except` is a principal submatrix (results.py:236-258), correct because other-bin rows are exactly zero; marginalization of (ap,at) inside each record (public.py:47-54, 2x2 inverse); ordering sigma_ap,sigma_at,corr correct.
- Multi-bin combination: combine_results sums data Fisher (results.py:261-297), prior once; run_forecast checks non-overlap, common registry and h_fid.
- Saved vs returned: settings.json floats are the same Python floats (JSON round trip exact); individual_fisher/joint_fisher/combined_fisher are data_fisher of the same objects in record order (public.py:184-193); combined is the full 12x12.
- Density normalisation: strict z>z_norm_min, target sum over masked table, magnitude bounds; matches lyaforecast for lya(qso) (16.1-26.75 mask, z>2.15), qso (z>2.15, no mask), lbg/lae (all z, no mask); table magnitude ranges (QSO 16.1-24.9, LBG/LAE 21.75-26.75) lie inside survey bounds so masks are inert.
- dN/dz dm -> per deg^2 per km/s: rho = dN/dz dm * (1+z)/c with z=z_source for forests (density_per_velocity called with z_source, survey_config.py:1422) and z_eval for galaxies (noise.py:27); matches lyaforecast weights.py:175.
- Galaxy n_bar in (h/Mpc)^3 = sum(q rho) a_v/d_deg^2 (noise.py:32): recomputed independently as (integral dN/dz dm)/(c D_M^2/H (pi/180)^2 h^3) from the raw tables: ratio 1 +/- 4e-10 for qso, lbg, lae in all six bins (qso 4.49e-5 ... 5.98e-6; lbg 8.3e-6 ... 3.0e-5; lae 1.16e-4 ... 9.0e-5). lyaforecast's Riemann sum on 107 points is 0.4-0.6% lower (discretisation only).
- Geometry: a_v=H/(h(1+z)), d_deg = h D_M pi/180 recomputed, exact; volume matches (z_max-z_min)*Omega*dV/dz at z_eval to 0.02-0.04% (curvature across bin, expected); units (Mpc/h)^3.
- SNR: 1/(SNR^2 * pix * N_exp/N_file), pix in Angstrom (0.8), SNR per Angstrom tables (NEXP=4 in header, num_exposures=4 -> factor 1); smoothing sigma=10 samples, reflect, truncate 4 along wavelength, same as lyaforecast; clamp 1e-10 (variance 1e20), sentinel 1e20 for m>24.75 or z/lambda outside table, bright clamp at 19.25; QSO and LBG variances at m = 17, 19.25, 20.1, 22.3, 24.7, 24.75, 24.76, 25.5 equal an independent scipy recomputation (ratio 1 exactly).
- Effective SNR magnitude limit: tables end at 24.75; forest weights fall to ~1e-20 above; QSO density ends at 24.9 (no loss beyond), LBG density extends to 26.75 but carries no forest information above 24.75; galaxy n_bar integrates the full 21.75-26.75. Consistent with lyaforecast (large noise for rmag > table maximum).
- Forest geometry: z_source = lambda_obs/sqrt(1040*1205)-1 with lambda_obs=1215.67(1+z_eval) (2.383 ... 3.660), L = c ln(1205/1040) = 44147 km/s, pixel width c*0.8/lambda_obs (63.3 ... 46.0 km/s); same as lyaforecast _get_zq_bin/_get_forest_length/_get_pix_kms.
- Magnitude quadrature (magnitude.py): partition contains 16.1, 26.75, density magnitude nodes and spline knots (within bounds), quadratic zero-crossing roots (algebra verified: p(0)=y0, p(1/2)=y1, p(1)=y2), and all SNR magnitude nodes (19.25 clamp kink, 24.75 sentinel step); 162-168 intervals x GL16 = 2592-2688 nodes; no interval straddles a density/SNR breakpoint.
- Forest x galaxy and forest x forest crosses: response products W_i W_j (forest x galaxy = W_forest, lya(qso) x lya(lbg) = W^2, both forests share pixel/resolution), noise zero off diagonal (independent sampling, as lyaforecast cross branch); every forest's weights use its own auto P3D and P1D at the representative mode.
- Bias inputs: analytic lya/qso match lyaforecast constants; LBG/LAE tabulated bias with linear extrapolation as interp1d; beta_F constant 1.45 equals lyaforecast (alpha=0).
- Forest-weight solver in the native path: converged in 12-36 updates per forest per bin (status "converged"), final steps ~1e-14.
- resources.py, parameters.py, survey.py: no defects affecting the default path.

## Not checked / limits

- No Fisher matrices were computed (no run_bin/Forecast.run); Fisher information ordering (joint >= individual) and derivative accuracy were not tested.
- Model, response kernels, weights recurrence internals, Fisher/covariance kernels were read only where needed to confirm conventions; they belong to other audit scopes.
- bao_marginalized.py and full_shape paths were not audited beyond confirming they do not alter the default (mode=bao) dispatch or parsing.
- The archive and other review files were not read, per the barrier.
