# Independent review of fishhighz against lyaforecast: DESI-2 15×2pt BAO forecast

**Scope.**
- fishhighz (fh): HEAD 456129b, `data/desi2_accuracy.ini`.
- lyaforecast (lf): NewForecast, `examples/desi2/lya_qso_lbg_lae_15x2pt.ini`.
- Survey: 5000 deg², six bins in 2.0 < z < 3.41, fields lya(qso), qso, lbg, lae and lya(lbg).
- Papers: ME07 (astro-ph/0607122), MW11 (1102.1752), FR14 (1308.4164).

**Notation.**
- σ∥ ≡ σ(α∥) and σ⊥ ≡ σ(α⊥) are joint 15×2pt errors.
- "Combined" means the sum of the six per-bin 2×2 Fisher matrices.
- Unless stated, Δ = σ/σ_fh,baseline − 1, in percent, quoted as ∥/⊥.

Every number traces to `work/E_results.md`, `work/E3_corrections.md`, `work/F1_verification.md` or `appendix/runs_summary.csv`. The detailed review is in `report.md`.

## 1. High-level summary of findings and main conclusions

fishhighz is scientifically sound for DESI-2 BAO forecasts within the linear, Gaussian, fixed-damping Fisher framework it implements.

- **Numerical convergence.** Every k, μ, magnitude, AP-step and weight-tolerance refinement changes σ by < 10⁻³ % (E §2.4).
- **Forest noise model.** It is ME07 Eqs. 13, 16 and 17 (= FR14 Eqs. 23–25 = MW11 Eq. 5), implemented correctly. An independent implementation reproduces the production A and P_pix to all printed digits (F1).
- **Wick/Fisher layer.** With lf's arrays as input it reproduces lf to 1.8 × 10⁻¹⁵ (`fh_full_compat`). This validates the covariance and Fisher arithmetic, not fh's physics (F1 F5).

**Headline comparison.**
- σ_fh/σ_lf = 0.942/0.965 (combined joint: 0.007191/0.005020 against 0.007625/0.005205).
- Against a grid-converged lf reference the ratio is 0.941/0.960 (E3 §6).
- Per bin, fh is lower by 2.2–10.3 % in σ∥ and by −0.35 to 7.2 % in σ⊥.

**Attribution.** The difference is fully attributed to identified conventions:
- The cumulative chain closes onto lf to +0.35/+0.55 % (+0.17/+0.02 % against the converged lf), and to ≤ 0.7 % in every bin.
- Two changes dominate:
  - lf's prefix-I1 forest weights, a genuine lf bug: +3.1/+3.0 %.
  - fh's mixed forest–galaxy damping rule, an analysis-plan choice: +3.1/+1.5 %.

**Headline caveats.**
1. **Mixed damping.** The rule Σ²_ij = (Σ²_i + Σ²_j)/2 is realistic only if forest × galaxy crosses are measured against a reconstructed galaxy catalogue. DESI Lyα×QSO BAO does not reconstruct the QSOs. With lf's rule, σ∥ is 3.3 % larger.
2. **Weights.** The baseline "early lyaforecast" weights are neither ME07's nor minimum variance.
   - MW11 Eq. 11 inverse-variance weights lower σ by 1.8/1.9 %, and the forest autos by 7–10 %.
   - The baseline depends on the arbitrary weighting reference mode at ±1.3 %.
3. **Larger uncertainties.** Each of the following effects, in at least one direction, is larger than the entire lf→fh difference:
   - reconstruction efficiency;
   - galaxy redshift errors;
   - survey density normalisation;
   - forest bias values;
   - continuum-fitting loss of low-k∥ modes;
   - the k⊥ Nyquist or beat-coupling limit.

   The fh–lf difference of ≈ 6/4 % is therefore not the dominant uncertainty of the DESI-2 forecast (§4).
4. **Bin-6 LBG forest.** Its noise rests on constant extrapolation of the LBG density table beyond z = 3.29. A range-averaged density truncated at the table edge raises bin-6 σ by +21.7/+16.8 % (combined +2.1/+1.5 %).

**Bugs.**
- No fh bug changes a forecast number. The fh defects are a CLI crash, a stale validation-test fixture and stale provenance text.
- In lf, the prefix-I1 weight recurrence is a genuine bug of scientific size for the LBG forest (lya(lbg) auto σ +42/+53 %).
- The lf resolution convention (σ = c/R) and backward-difference normalisation are minor genuine bugs.

## 2. Scientifically meaningful changes lyaforecast → fishhighz

"Effect" is the isolated toggle: Δ of fh with the lf convention substituted (joint, combined). Where applicable, the step of the cumulative chain, measured against lf, is given in brackets.

| change | lyaforecast | fishhighz | paper basis | effect on joint σ∥/σ⊥ | verdict | reason |
|---|---|---|---|---|---|---|
| Forest weights: I1 normalisation | prefix I1(<m) in each source's noise, `weights.py:199-205, 339-358` | full-sample scalar I1, `kernels/full_sum_weights.py:31-52` | ME07 Eqs. 14, 19; FR14 Eqs. 26–29; MW11 Eq. 2 | +3.12/+3.02 (+3.12/+3.02), combined weight toggle `legacy3`; lya(lbg) auto +42/+53 | justified | ME07 and FR14 integrate I1 to m_max. The prefix under-weights bright sources and inflates LBG-forest N_F by ×1.6–2.1. |
| Forest weights: aliasing in the weighting signal | S = P3D only (FR14-literal), despite the code comment, `weights.py:150-154` | S′ = S + B/(L I1), `full_sum_weights.py:45-46` | ME07 §II B (P_S "including the aliasing term"); FR14 Eq. 29 text; MW11 Eq. 11 | inside `legacy3`; ME07-literal McDonald −0.19/−0.17; MW11 inverse variance −1.83/−1.85 | questionable | The early form is in none of the papers. N_eff exceeds the MW11 minimum by 0.2–22 %. |
| Forest weights: convergence | 3 fixed updates of a map with no non-trivial fixed point, `weights.py:114-116` | converged to rtol 10⁻⁵ in 12–36 updates | ME07 §II B ("a few iterations") | inside `legacy3`; rtol 10⁻⁴ or 10⁻⁶ gives \|Δ\| < 10⁻⁶ | justified | In lf, N_eff/N_min = 1.28, 1.98, 16.2 after 1, 3, 40 updates (LBG bin 4), i.e. a transient. |
| Resolution | σ = c/R = 120 km/s, `covariance.py:182-183`, `spectrograph.py:301` | σ = c/(2√(2 ln 2) R) = 50.9 km/s, `response.py:51-54` | ME07 §II A: FWHM = 2.355 σ_R = λ/R | +0.23/+0.05 (+0.19/+0.04) | justified | lf over-smooths by 2.355 in σ; negligible for BAO. |
| Cross-pair damping | any pair containing Lyα unreconstructed, `fisher.py:156-161` | Σ²_ij = (Σ²_i + Σ²_j)/2, `kernels/kaiser.py:48-55` | FR14 Eq. 16 and "no reconstruction" for Lyα; SE07 §2; neither specifies crosses | +3.31/+1.65 (+3.11/+1.47); lya(qso)×qso +8.5/+4.4 | justified with caveat | Product of propagators √(D_iD_j). Valid only for crosses with a reconstructed catalogue. |
| BAO derivative | −[μ², 1−μ²] ∂(D P_w)/∂ln k at fixed μ, no Q, backward FD, `fisher.py:135-167` | full AP remap of the wiggle term (Q, μ′, D(k′)), central FD, `kernels/kaiser.py:13-78` | FR14 §4.1.1 (wiggle-only derivative, damping outside); SE07 Eq. 5 | lf-like −0.16/−0.52 (−0.17/−0.54); SE07 form +0.43/−0.27; Q = 1: +0.05/−1.17 | justified with caveat | At fixed Σ, the Q, μ′ and D(k′) terms carry ≤ 0.7 % per bin of non-BAO information (forest autos ≤ 2.7 %). |
| Wiggle/no-wiggle split | 8th-order polynomial in ln k fitted to the full model, `fisher.py:169-198` | Vega template PK − PKSB, `models/templates.py` | FR14 §4.2 (subtract a smoothed spectrum) | information ratio 0.999–1.000 at k > 0.05 (B row 25); not isolated in a forecast, bounded by the closure residual | justified | Equivalent for the BAO information. |
| Growth | EdS ((1+2.3)/(1+z))², `power_spectrum.py:68-77` | CAMB [σ8(z)/σ8(2.406)]², `survey_config.py:1285-1286` | linear theory | +0.19/+0.25 (+0.21/+0.26); bin 1 −0.26/−0.37, bin 6 +0.55/+0.70 | justified | Removes an EdS error of up to 1.4 % in P at z = 3.3. |
| Evaluation redshift | model at arithmetic z_c, covariance at λ-geometric z | single z_eval = √((1+z_min)(1+z_max)) − 1, `survey_config.py:928` | ME07 §II A (central z) | −0.08/−0.08 (−0.06/−0.06) | justified | Removes an internal inconsistency. |
| Volume | point evaluation at z_mean, `covariance.py:237-242` | 32-node GL in z, `geometry.py:162-164` | FR14 Eq. 20 | V_lf/V_fh = 1.0002 | justified | Negligible. |
| Mode quadrature | 500 linear k (backward FD scaled by dk/k) × 10 μ midpoints | 128 × GL4 in k × GL32 in μ | FR14 Eq. 20 (small-bin limit) | lf converged grid +0.19/+0.53 on lf (E3); fh refinements < 10⁻³ % | justified | lf carries a ~0.5 % numerical floor, with partial cancellation between k and μ. |
| Density normalisation | first-spacing dz, z > 2.15, signed spline | identical, plus a floor on the negative interpolant, `legacy_compat.py:122` | – (input policy) | floor: n̄ +0.44–0.60 %, σ < 0.3 %; per-row cell widths (shared) +0.91/+1.35 | floor justified; first spacing questionable (shared) | The floor removes unphysical negative densities; the first spacing over-counts LBG/LAE by 4.8/2.6 %. |
| Floors/clamps | z-clamp through `bbox` | same clamp; same SNR sentinels | – | shared LBG-forest clamp; range-averaged truncation +2.12/+1.51 (bin 6 +21.7/+16.8) | questionable (shared) | Constant extrapolation of a table that falls 5.9× between z = 3.07 and 3.29. |
| Magnitude quadrature | rectangle rule, 107 nodes, `survey.py:26-29` | GL16 on density, spline-root and SNR breakpoints (2592 nodes), `magnitude.py:11-58` | ME07 Eqs. 14–18 | lf with 54/214 nodes: +0.75/+1.15, +0.17/+0.16; lf is non-monotone at 10⁻⁴ | justified | Removes lf grid noise (lya(lbg) auto up to 5.5 %). |
| Bin-1 selection | all 15 spectra with clamped LBG/LAE | lya(qso)², lya(qso)×qso, qso² (user directive) | – (input coverage) | lf restricted identically (`joint_fhsel`); all 15 spectra with clamped densities −1.55/−1.81 | justified | No LBG/LAE table support below z = 2.27; lya(lbg) alone adds only −0.05. |
| Multi-bin combination | per-bin output only | independent bins, summed Fisher, `results.py:261-297` | FR14 Eqs. 11, 20 | none per bin | justified | Cross-bin correlations at k∥ < 0.035 h Mpc⁻¹ are neglected, not double counted. |

**Headline ratio and closure.**
- fh/lf = 0.942/0.965 against `lf_baseline`, and 0.941/0.960 against the grid-converged lf (`lf_conv_k2000`: num_k 2000, n_μ 40, 214 magnitudes).
- The chain from fh (legacy-3 weights, c/R, lf damping, EdS, arithmetic z, lf-like AP) reaches +0.35/+0.55 % of `lf_baseline` and +0.17/+0.02 % of the converged lf. The per-bin residual is ≤ +0.7 %.
- The largest residual is in the lya(lbg) auto, +3.4/+5.1 %. It comes from fh's S and B auxiliary spectra, not from the update rule.
- Weights and damping account for ≈ 6.2/4.5 % of the chain; the other steps sum to < 0.5 %.
- Over 20 matched variations the fh/lf ratio is stable at 0.943 ± 0.010 / 0.965 ± 0.008. It moves by up to 2.9 % only when density, reconstruction or the magnitude limit changes (E §3, fig. 5).

## 3. Bugs found

### 3.1 fishhighz

| ID | defect | file:line | impact on σ | F1 verdict | class |
|---|---|---|---|---|---|
| fh-1 | `--joint-only` always raises in BAO mode | `cli.py:21` → `public.py:489-490` | none (software) | CONFIRMED (B1a) | bug |
| fh-2 | stale fixture: `test_independent_accuracy_check[production, oracle]` raises AttributeError, so the independent accuracy check is void | `tests/test_step12_performance.py:161,174` vs `validation/accuracy.py:644` (commit 815f661) | none on forecasts | CONFIRMED genuine (B1c); the other 5 of the 7 pytest failures are environmental (`sys.executable -I`; fishhighz is not installed in `vega_test`) | bug (validation) |
| fh-3 | STOPPING rtol = 10⁻⁴ while CONTROLS use 10⁻⁵; RESEARCH_BASELINE says "remain 1e-4" | `accuracy.py:10` vs `:25`; `docs/research/RESEARCH_BASELINE.md:239-242` | \|Δσ\| < 10⁻⁶ | CONFIRMED (B1b) | provenance |
| fh-4 | volume element lacks the curvature factor for Ω_k ≠ 0 | `geometry.py:163-164` | none (flat baseline) | not tested by F1 (D1 F3) | latent |

The claim that per-category `k_max` masking is missing in BAO mode (D1 F2) is refuted: the parser rejects those keys in BAO mode (`survey_config.py:527-541`; F1 B1d).

**Design choices (not bugs) with measured consequence:**
- early weights instead of MW11 inverse variance: 1.8 %;
- arbitrary weighting reference mode: ±1.3 %;
- mixed damping: 3.3/1.7 %;
- full-AP wiggle derivative: ≤ 0.7 % per bin.

**Input-data limitations shared with lf:**
- LBG/LAE table z coverage is 2.27–3.40 and is clamped outside it (+2.1/+1.5 %).
- First-spacing normalisation (+0.9/+1.35 %).
- The z_norm_min = 2.15 threshold sits on a cell centre (a value of 2.1 gives +2.0/+2.6 %).
- The quadratic tensor spline raises the LBG source density by 8.5 % at z_src = 2.89.
- The SNR-table EXPTIME semantics are undocumented.

### 3.2 lyaforecast

| ID | defect | file:line | impact on σ | F1 verdict | class |
|---|---|---|---|---|---|
| lf-1 | cumulative prefix I1(<m) in the per-source noise, iterated 3 times on a map with no non-trivial fixed point; bright LBGs are suppressed (the weight peaks at r = 24.07, where N/B = 17) | `weights.py:199-205, 339-358, 114-116` | joint +3.1/+3.0 %; lya(lbg) auto +42/+53 %; LBG N_F ×1.63–2.09, dependent on the iteration count | CONFIRMED (W1) | bug, major for the LBG forest |
| lf-2 | resolution σ = c/R instead of c/(2.355 R) | `covariance.py:182-183`, `spectrograph.py:301` | +0.23/+0.05 % | CONFIRMED (F3) | bug, minor for BAO |
| lf-3 | backward difference normalised by dk/k: offset by half a bin and inflated by ≈ dk/(2k) | `fisher.py:163-165`, `power_spectrum.py:46` | lf ≈ 0.9 % optimistic at num_k = 500 | CONFIRMED (F4) | bug, minor |
| lf-4 | the comment "weights include aliasing as signal" is contradicted by `signal_power = self._p3d_w` | `weights.py:150-154` | inside lf-1; the sign varies with S L I1/B | CONFIRMED (W2) | code/comment inconsistency; the implementation is FR14-literal |
| lf-5 | signal model at the arithmetic z_c, covariance at the λ-geometric z_mean | `survey.py:78-82`, `covariance.py:134-138` | −0.08 % | E toggle | inconsistency, negligible |

**lf approximations that fh corrects (not bugs):**
- EdS growth (+0.19/+0.25 %);
- 10 μ midpoints (σ∥ ≈ 0.7 % pessimistic);
- a 107-node magnitude grid (non-monotone at 10⁻⁴);
- the negative spline overshoot is kept (−0.5 % in n̄);
- the damping-envelope derivative is kept inside ∂P/∂α, contrary to SE07 and FR14.

The "any Lyα pair unreconstructed" damping rule is an analysis choice. It is the realistic one for unreconstructed QSOs.

## 4. Scientifically relevant approximations and assumptions in fishhighz

In the last column, "larger" means the effect exceeds the net lf→fh difference (5.8/3.5 %) in at least one direction. "Comparable" means it exceeds the largest single convention step, which is ≈ 3 %.

| approximation (fh) | paper statement | measured Δσ∥/Δσ⊥ (joint, combined) | vs lf→fh |
|---|---|---|---|
| No galaxy redshift-error or FoG damping | SE07 Eq. 27; absent in FR14 Eq. 15 | σ_z = 2 h⁻¹Mpc: +2.55/+0.73; 4 h⁻¹Mpc: +9.14/+2.68 (bin 5: +12.9) | **larger** (∥) |
| Reconstruction as a scalar R = 2 on galaxy Σ | FR14 r(n̄P) interpolation, §4.1.1 | off: +10.42/+6.38; FR14 r per tracer: −2.82/−1.61; R = 3.43: −5.07/−2.85 | **larger** |
| No continuum-fitting loss of low-k∥ forest modes | ME07 §II B (first 2N_q,los k∥ modes dropped); MW11 §4.4, §5 | cut at k∥ > 0.015 ≈ 2π/L: +0.54/+2.84; at > 0.03: +1.09/+5.86; forest autos σ⊥ +9.8/+23.5 | **larger** (⊥ at 0.03) |
| No k⊥ Nyquist cut; information not capped | ME07 §II B; FR14 after Eq. 29 ("definitely must be correct"); MW11 Eq. 10 (beat coupling caps it instead) | with n_2D (k_Nyq = 0.389→0.078 h Mpc⁻¹): +0.40/+1.23 (bin 6: +1.7/+3.0); with MW11 n̄_eff: +1.29/+5.52 (bin 6 σ⊥: +10.3) | **larger** (⊥, n̄_eff); smaller (n_2D) |
| Forest bias b_F(2.33) = −0.1352, β_F = 1.45 | FR14 quotes none (MEO03 Table I) | eBOSS DR16 (−0.117, 1.669): +1.80/+3.07; DESI DR1 (−0.108, 1.743): +3.20/+4.71 (forest autos +15/+30). Caveat: these are pure-forest values after HCD separation. | **comparable/larger** (⊥) |
| Galaxy biases | FR14 b(z)D(z) = const | QSO bD = 1.2: +1.19/+2.14; LAE bD = 0.89: −3.04/−5.12 | comparable (⊥) |
| Survey density taken as tabulated; no targeting efficiency | FR14 ×0.8 (DESI), ×0.73 (BOSS) | density ×0.5: +36.3/+49.1; z_norm_min = 2.1: +2.03/+2.59; per-row cell widths: +0.91/+1.35 | **larger** (scale) |
| Early FKP-like weights at one reference mode, instead of MW11 optimal weights | ME07 Eq. 19; FR14 Eq. 29 (k = 0.07, μ = 0.5); MW11 Eq. 11 | inverse variance: −1.83/−1.85 (forest autos −7 to −10); reference mode ×0.5 / ×2 / ME07 mode: +1.34 / −1.28 / −1.24 | smaller |
| Fixed Σ; no amplitude or broadband marginalisation | SE07 §2 and FR14 §4.1.1 (damping outside ∂); FR14 Eqs. 19–20 | +1.3–1.5 % once A_w, Σ, b, β and the broadband are marginalised (C2; analytic, not a forecast) | smaller |
| Gaussian covariance, diagonal in k; MW11 (4/N) w̄² P_F P_F′ and beat-coupling terms dropped | MW11 Eqs. 9–10; ME07 §II B (O(N_q⁻¹ᐟ²) corrections) | not computed; inverse-variance weights raise A by up to ×1.45, so part of their gain may lie in these terms | unknown (bracketed by the Nyquist row) |
| Source density taken at a single z_src (midpoint rule in z_s) | ME07 §II A (central z) | both forests range-averaged: −1.12/−0.90 (lya(qso) auto bin 6: −23.8/−25.4) | smaller |
| LBG density extrapolated at constant value beyond z = 3.29 | – | range-averaged and truncated at 3.40: +2.12/+1.51 (bin 6: +21.7/+16.8) | smaller (bin 6: **larger**) |
| Forest range 1040–1205 Å | ME07: 1041–1185 Å; FR14: 985–1200 Å (BAO) | 1040–1185 Å: +1.94/+1.94 | smaller |
| P1D (PD2013 with a low-k plateau) | ME07 Eq. 13; MW11 Eq. 6 | ×0.8 / ×1.2: −1.10 / +0.93; no plateau: −0.15/−0.70 | smaller |
| Independent redshift bins; single z_eval per bin | FR14 Eqs. 11, 20; ME07 §III (a two-bin split gives slightly smaller errors) | not computed; second order, O((Δz/(1+z))²) ≈ 0.5 % | smaller |
| Linear theory: Kaiser, no forest D_NL (FR14 Eq. 22) | FR14 Eqs. 15, 21 | not computed; D_NL acts at k∥ ≳ 0.3, where the wiggles are already damped | smaller (expected) |
| Shared-sample noise neglected for lya(qso)–qso and lya(lbg)–lbg | MW11 Eq. 34 | correlation coefficient ≲ 1 % if same-object pairs are excluded (C2 §5) | smaller |
| SNR tables: √N_exp scaling, ambiguous EXPTIME semantics | MW11 Eqs. 3–4 | N_exp = 2 / 8: +5.07/+4.84 and −5.53/−5.86; the forest noise is 34–86 % aliasing | unresolved input |

**Conclusion.** The fh–lf convention differences are resolved to ≲ 0.5 % and are no longer the limiting uncertainty. The inputs that now control the DESI-2 Lyα BAO precision at the 5 % level are analysis-plan choices, and must be fixed before that precision is quoted:
- the reconstruction model, including whether forest × galaxy crosses use reconstructed catalogues;
- the galaxy redshift errors;
- the continuum k∥ loss;
- the Nyquist or beat-coupling treatment;
- the forest bias values;
- the targeting efficiency and density normalisation.

