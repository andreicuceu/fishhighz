# Detailed review: fishhighz vs lyaforecast for the DESI-2 15×2pt Lyα + galaxy BAO forecast

**Scope.**
- Codes: fishhighz (fh) at HEAD 456129b, with an unmodified working tree; lyaforecast (lf) NewForecast, read-only.
- Case: DESI-2 15×2pt, 5000 deg², six bins of width 0.235 in 2.0 < z < 3.41.
- Fields: lya(qso), qso, lbg, lae, lya(lbg).
- Paper abbreviations: ME07 (astro-ph/0607122), MW11 (1102.1752), FR14 (1308.4164), SE07 (astro-ph/0701079).

**Notation.**
- σ∥ ≡ σ(α∥) and σ⊥ ≡ σ(α⊥).
- "Joint" is fh's selection: in bin 1, lya(qso)², lya(qso)×qso and qso²; in bins 2–6, all 15 spectra. The lf joint is restricted to the same selection (`joint_fhsel`).
- "Combined" is the sum of the per-bin 2×2 Fisher matrices, with α∥ and α⊥ common to all bins. This is a review construct applied identically to both codes.
- Δ = σ/σ_baseline − 1, in percent, quoted as ∥/⊥.
- Code paths are given relative to `lib/fishhighz/fishhighz/` or `lib/lyaforecast/lyaforecast/`.

## 1. Scope and method

### 1.1 Question

This review asks whether fh is a scientifically sound replacement for lf for DESI-2 BAO forecasts. To answer it, the review:
- identifies every scientifically meaningful change from lf;
- attributes the change in σ to each difference;
- measures the size of the approximations shared by both codes against the three papers.

### 1.2 Agents and work products (`work/`)

| stage | file | content |
|---|---|---|
| paper digests | A1_ME07, A2_MW11, A3_FR14 | Equation digests from the arXiv sources, including stated approximations and ambiguities |
| inventory | B_inventory, map_* | 29-row stage-by-stage lf/fh inventory with verified file:line references; small synthetic calculations |
| scientific assessment | C1_weights, C2_signal_fisher, C3_survey_noise | Derivations and preparation-level calculations on the production inputs, plus requested runs |
| implementation audit | D1_numerical_core, D2_orchestration_inputs | Independent re-implementations of kernels, geometry, quadrature and Wick/Fisher; input-layer audit |
| numerical runs | E_results (+ `appendix/runs_summary.csv`, figures) | 107 completed full forecasts with attribution toggles, variations and follow-ups |
| verification | F1_verification | Adversarial re-test of every claim against code, papers and saved runs. **Authoritative where it conflicts with earlier files.** |
| corrections | E3_corrections | 12 corrective runs requested by F1 |

### 1.3 Information barrier

No file under `lib/fishhighz/reviews/` outside this directory, and nothing in `lib/fishhighz/docs/archive/`, was read. F1 checked the work files and confirms that the barrier holds. The only archived number used, the McDonald bin-2 amplitude of 0.8668 at 96 updates, is quoted in `docs/research/FOREST_WEIGHTING_DECISION.md:318`.

### 1.4 Runs and resources

**Environment.**
- Harness and raw outputs: `/pscratch/sd/a/acuceu/fishhighz-review/` (`harness/`, `runs/<case>/result.json`).
- Conda environment `vega_test`: Python 3.13.0, CAMB 2.0.1.
- One numerical thread per process.
- fh was placed on `PYTHONPATH`, because it is not installed in `vega_test`.
- Toggles were injected by monkeypatch (`toggles.py`, `toggles2.py`, `toggles3.py`; described in `TOGGLES.md`). No package file was modified.

**E set (82 production cases and 26 follow-ups).**
- 108 cases were launched and 107 completed: 79 fh, 26 lf and 2 fh compatibility profiles.
- `fh_t_w_mcdonald` failed by design at the 96-update cap; `fh_r3_mcd2000` is its retry.
- The runs used a Perlmutter interactive CPU-node allocation, with up to 48 concurrent single-thread cases.
- Total cost was 3.44 core-h. Wall time was 88–161 s per fh case and 126–142 s per lf case; peak RSS was ≤ 578 MB.

**E3 set (12 corrective cases).**
- Run on allocation 59090748 (node nid200232); all 12 completed with status ok.

**Tests.**
- `pytest`: 1576 passed and 7 failed. Five failures are environmental; two are genuine and form one defect (fh-2, §3).
- `ruff`: clean (`appendix/pytest_ruff.log`).

**Line numbers.** File:line references cited in the deliverables were re-checked against the code while writing them.

### 1.5 Not run or not resolved

1. **Σ-marginalised BAO forecast.** It needs width derivatives in `KaiserModel`, i.e. a package change. The +1.3–1.5 % estimate is C2's analytic calculation on prepared bins.
2. **MW11 covariance terms.** No forecast includes the off-diagonal or beat-coupling terms of MW11 Eq. 10.
3. **SNR table semantics.** Whether EXPTIME is per exposure or total in the SNR tables was not confirmed; the `desi_quicklya.py` defaults are not in the repository.
4. **Other fh modes.** `full_shape` and `bao_marginalized` were not audited; only `mode = bao` was.
5. **Missing physics.** No continuum-distortion model beyond k∥ masks, and no HCD or metal contributions to P1D or to the signal.

## 2. Forest noise and weights in common notation

### 2.1 Quantities (observed coordinates: deg and km/s)

- Source density: ρ(m) = dN/(dz dm deg²) · (1+z_s)/c, per deg² per km/s per magnitude.
- Forest length: L = c ln(1205/1040) = 44147 km/s. The pixel width is l_p = Δv = c · 0.8 Å/λ_obs, i.e. 63 → 46 km/s.
- Pixel noise: σ_N²(m) = [SNR_Å² Δλ (N_exp/N_exp,file)]⁻¹. The one-sightline white-noise power is N(m) = l_p σ_N²(m), MW11's P_N,n (Eq. 3).
- Aliasing numerator: B(k∥) = P1D(k∥) W²(k∥).
- Weighting signal: S = P3D(k, μ) W² a_v/d_deg², where a_v = H/(h(1+z)) and d_deg = h D_M π/180.
- Unit conversion: 1 deg² km/s = d_deg²/a_v (h⁻¹Mpc)³ (D1 check 2).
- Weight moments: I1 = ∫dm ρw, I2 = ∫dm ρw², I3 = ∫dm ρσ_N²w² (ME07 Eqs. 14, 15, 18; FR14 Eqs. 26–28).

### 2.2 Observed power

Averaging the weighted field δ_w = (w/w̄)(δ_F + δ_N) (ME07 Eq. 10) over Poisson sightline positions gives

 P_obs(k) = W² P3D(k) + A B(k∥) + P_pix, with A = I2/(L I1²) and P_pix = l_p I3/(L I1²).

- These are ME07 Eqs. 13, 16 and 17, and FR14 Eqs. 23–25.
- The mapping to MW11 Eq. 5 is n̄⁻¹ \overline{w̃²} = A and n̄⁻¹ P̄_N = P_pix.
- The three papers agree exactly.
- fh implements this in `noise.py:78-90`; A and P_pix come from `kernels/weights.py:23-37` and `kernels/full_sum_weights.py:92-102`. lf uses the same final expressions (`covariance.py:292-294, 352-363`).

### 2.3 Minimum-variance weights

The weights enter only through

 N_eff(w) ≡ A B + P_pix = ∫ρw²(B + N) / [L (∫ρw)²].

**Cauchy–Schwarz bound.** N_eff ≥ [L ∫ρ/(B+N)]⁻¹, with equality iff w ∝ 1/(B + N). This is MW11 Eq. 11. The minimum is N_eff,min = B/n̄_eff, with n̄_eff = L ∫ρν and ν = B/(B+N) (MW11 Eqs. 12–13).

**Optimality within fh's model** (F1 W3):
- ∂P3D/∂α is independent of w.
- T → T + D, with D ⪰ 0 independent noise, cannot increase the Fisher information of the multi-tracer Wick likelihood.
- Therefore the information at each (k, μ) is monotone in N_eff.

**k dependence.** The optimum depends on k∥ only, through B(k∥):
- B(k∥ = 0.1, 0.2, 0.3)/B(ref) = 0.97–0.99, 0.86–0.92 and 0.76–0.83.
- The per-node optimum lies within 0.4 % of the scalar weight.
- MW11 finds a 2 % loss from using w(k∥ = 0).

### 2.4 FKP-type prescriptions in one form

ME07 Eq. 19 and FR14 Eq. 29 give w = (P_S/P_N)/(1 + P_S/P_N), with P_N(m) = σ_N² l_p/(I1 L). Equivalently, w_i = c/(c + N_i) with c = P_S L I1. The prescriptions differ only in c:

| prescription | c | source |
|---|---|---|
| minimum variance (`inverse_variance`; also the seed) | B | MW11 Eq. 11 |
| McDonald (`sum_aliasing`) | S L I1 + B I2/I1 | ME07 §II B: P_S is the total flux power "including the aliasing term" |
| FR14-literal (`sum_intrinsic`) | S L I1 | FR14 text after Eq. 29 |
| early lyaforecast (fh baseline, `sum_historical`) | S L I1 + B | none of the three papers |
| lf NewForecast | S L I1(<m_i) | none (prefix bug) |

**Consequences for the production inputs** (C1, F1 W1–W4):

1. **Early weights.**
   - At the reference mode, S L I1/B is 3.7, 2.9, 1.9, 1.25, 0.68 and 0.30 for the QSO forest in bins 1–6, and 1.8, 21, 19, 2.5 and 3.6 for the LBG forest in bins 2–6.
   - The resulting excess N_eff,early/N_min is 1.216, 1.144, 1.089, 1.042, 1.013 and 1.002 (QSO) and 1.018, 1.196, 1.202, 1.053 and 1.076 (LBG).
   - The early P_S = S + B/(L I1) equals the total observed power at the optimal weights, since B/(L I1) = N_eff,min at w = ν. It is not a moment-consistent aliasing term. This is harmless for the covariance, which is exact for any w.
2. **lf prefix.**
   - For the LBG forest, I1(<m)/I1 is 0.02–0.25 over r = 22.5–23.5, where N ≳ B already. The weight then peaks at r = 24.07, where N/B = 17.
   - In the MW11 weighting, 44 % of the information lies at r < 23.5.
   - After 3 updates, N_F(lf3)/N_F(early) is 2.09, 1.76, 1.65, 1.76 and 1.63 (LBG bins 2–6) and 0.94–1.07 (QSO).
   - The map has no non-trivial fixed point. In LBG bin 4, N_eff/N_min is 1.28, 1.64, 1.98, 4.73 and 16.2 after 1, 2, 3, 10 and 40 updates. lf's LBG-forest noise is therefore an iteration-count transient.
3. **McDonald and FR14-literal recurrences.**
   - Both maps are homogeneous of degree 1 in the weight amplitude, and they collapse iff Λ < 1.
   - LBG bin 2 has Λ_McD = 1.04, which is critical slowing: the solver converges at 624 updates (it returns state 2t, with t = 312), with amplitude 0.906. It also has Λ_FR14 = 0.87, i.e. collapse: amplitude 2.7 × 10⁻⁵ at 96 updates and 2 × 10⁻⁶² at 1000.
   - A and P_pix are invariant under w → aw, and the collapse limit is noise-only inverse variance. The "unavailable" McDonald cases therefore reflect the amplitude criterion and the cap, not a missing covariance.
4. **Reference mode.**
   - (2.4 deg⁻¹, 3.5 × 10⁻⁴ s/km) maps to k ≈ 0.051–0.052 h Mpc⁻¹ with μ = 0.69–0.81.
   - S there is 1.9–2.6× its value at ME07/FR14's (k = 0.07, μ = 0.5).
   - lf's comment "~0.035 h/Mpc" corresponds to k⊥ = k∥, not to μ = 0.5.
5. **Noise regime.**
   - Aliasing is 34, 46, 51, 63, 74 and 86 % of the lya(qso) noise at (0.14, 0.6) in bins 1–6, and 5–21 % of the lya(lbg) noise.
   - The forest autos are noise dominated: P/N is 0.43 → 0.05 for lya(qso) and 0.17–0.47 for lya(lbg).
   - The weak response to N_exp therefore reflects aliasing, not cosmic-variance limitation (F1 E-x).
6. **Sightline densities** (E3 §4).
   - QSO forest: n_2D = 63.5, 43.2, 31.3, 18.8, 9.7 and 3.9 deg⁻²; MW11 n̄_eff = 20.6, 19.7, 14.4, 10.6, 6.6 and 3.1 deg⁻².
   - LBG forest: n_2D = 74–436 deg⁻².

## 3. Per-stage comparison

### 3.1 Geometry and cosmology

| item | lf | fh | assessment |
|---|---|---|---|
| evaluation redshift | model at z_c = (z_min+z_max)/2 (`survey.py:78-82`); covariance, noise and V at z_mean = √(λ_min λ_max)/1215.67 − 1 (`covariance.py:134-138`) | z_eval = √((1+z_min)(1+z_max)) − 1 for every quantity (`survey_config.py:928`) | fh is consistent; the mismatch is ≤ 0.0022 in z; toggle −0.08/−0.08 |
| volume | Ω c ln(λ_max/λ_min)(dr/dθ)²/(dv/dr) (`covariance.py:202-213, 237-242`) | Ω h³ Σ₃₂ w c D_M²/H (`geometry.py:162-164`) | V_lf/V_fh = 1.0002; fh is flat-only (fh-4) |
| a_v, d_deg | `cosmoCAMB.py:95, 140-141` | `geometry.py:165-166` | identical (D1 check 1) |
| linear P | CAMB at z_ref = 2.3, np.interp (`cosmoCAMB.py:78-80`) | Vega template at z = 2.406, PK/PKSB, cubic in ln k, no extrapolation (`models/templates.py:140-229`) | agree to ±0.1 % (B row 5) |
| growth | EdS (`power_spectrum.py:68-77`) | [σ8(z_eval)/σ8(2.406)]² (`survey_config.py:1285-1286`; `kernels/kaiser.py:71`) | P_lf/P_fh = 1.004 → 0.987 over bins 1–6; toggle +0.19/+0.25 |

The volume error from single-z evaluation is O((Δz/(1+z))²) ≈ 0.5 %. ME07 §III found that a two-bin split gives slightly smaller errors, so single-z is adequate at ≲ 1 % (C2 §4).

### 3.2 Survey inputs

| item | lf | fh | assessment |
|---|---|---|---|
| density normalisation | Lyα: magnitude cut, then z > 2.15 (`tracer.py:148-151`); QSO z > 2.15, LBG/LAE all rows (`tracer.py:201-204`); dz from the first spacing; signed spline | first spacing (`adapters/legacy_inputs.py:153-161`); normalisation and spline (`adapters/legacy_compat.py:100-108`); negative interpolant floored at 10⁻²⁰ (`:122`) | identical except the floor, which adds 0.44–0.60 % to n̄ and < 0.3 % to σ |
| first spacing | LBG/LAE rows all divided by 0.22 | same | LBG/LAE spacings are 0.22–0.24; totals +4.8/+2.6 %; rows 2.60, 2.83 and 3.07 high by 2.3, 6.8 and 4.5 % |
| z_norm_min = 2.15 | strict `>` (`tracer.py:150, 203`) | same | the stored row value 2.1499999999999995 is excluded, also in exact arithmetic. The threshold sits on a cell centre, a 16 % ambiguity: 2.1 gives +2.03/+2.59, 2.2 is bit-identical (F1 S2) |
| z coverage | RectBivariateSpline with bbox, clamps (`tracer.py:173-177, 224-228`) | same; the `redshift_extension` flag is not saved (D2 M1) | LBG/LAE rows span only 2.38–3.29. lya(lbg) z_src = 3.404 and 3.660 in bins 5–6, so the clamp gives 113.5 deg⁻² per unit z (F1 S1) |
| 2D spline | kx = ky = 2 | same | ringing at sharp magnitude cut-offs (−0.44 to −0.60 %); along z it peaks at 769 deg⁻² per unit z at z = 2.93, against nodal values 700/665 (C3 §1) |
| magnitude quadrature | rectangle rule with 107 nodes (`survey.py:26-29`) | GL16 on breakpoints (`magnitude.py:11-58`; `survey_config.py:1398-1403`) | integrals agree to 1.0000 (B row 13); lf grid sensitivity in §4.4 |
| SNR | `spectrograph.py:102-132, 254-277` | `adapters/legacy_inputs.py:299-382`; `legacy_compat.py:151-256` | identical: σ = 10 samples along λ, bright clamp at r = 19.25, sentinel 10²⁰ beyond 24.75 |
| forest z_src, L | `covariance.py:197-200, 230-233` | `survey_config.py:1388-1397, 1409-1427` | identical; single-z_src midpoint rule |
| galaxy n̄ | `weights.py:275-293`; `covariance.py:384-390` | `noise.py:19-33` | same local approximation; fh exact to 4 × 10⁻¹⁰ against an independent calculation; lf Riemann sum −0.4 to −0.6 % |

**Data facts** (C3 §1, §3):
- QSO counts are non-zero only for r ∈ [16.5, 23.5]; LBG/LAE counts only for r ∈ [21.75, 24.25].
- Variations of r_max ≥ 24.75 are therefore null tests: fh changes by ≤ 8 × 10⁻¹⁵. lf changes by 2.5 × 10⁻⁴ (combined) and by up to 7 × 10⁻³ in the lya(lbg) auto.
- The pixel variations are null by construction, since σ² Δv is independent of Δλ.

### 3.3 Forest weights

| item | lf | fh | assessment |
|---|---|---|---|
| reference mode | (2.4 deg⁻¹, 3.5 × 10⁻⁴ s/km), z-independent (`weights.py:80-82`); S from smooth P3D | same mode (`accuracy.py` REFERENCE); S from the full auto including damped wiggles (`weights.py:88-159`, S and B at 128-140) | S differs by < 1 %; the mode itself is not ME07/FR14's (§2.4) |
| seed | (B/Δv)/(B/Δv + σ²) (`weights.py:129-132`) | same (`kernels/full_sum_weights.py:14-18`) | the seed is the MW11 optimum at the reference k∥ |
| update | S/(S + σ²Δv/(I1(<m)L)) (`weights.py:150-154, 339-358`); cumsum I1 (`:199-205`) | S′/(S′ + σ²Δv/(I1 L)), S′ = S + B/(I1 L) (`full_sum_weights.py:31-52`) | §2.4 |
| iterations | 3 fixed (`weights.py:114-116`) | rtol 10⁻⁵, minimum 3, 3 stable, cap 96, doubled-count confirmation (`full_sum_weights.py:111-213`); raises if not converged (`weights.py:320-324`) | fh needs 12–36 updates |
| alternatives | – | `mcdonald`, `inverse_variance` (`weights.py:336-341`), `legacy` | the INI parser accepts only `early_lyaforecast` |

Doc cross-check: see §6, items 1–8.

### 3.4 Signal

| item | lf | fh | assessment |
|---|---|---|---|
| Kaiser and biases | ∏ b(1+βμ²) (`analytic_biases.py:92-121`; `tracer.py:49-61`) | forest b(1+βμ′²), galaxy b + fμ′² (`kernels/kaiser.py:40-45`; `models/biases.py:85`; `survey_config.py:1288-1301`) | identical values: b_F(2.33) = −0.1352 with (1+z)^2.9; β_F = 1.45; b_QSO = 3.54 ((1+z)/3.33)^1.44; b_LBG = 3.3; LAE linear through 1.76 at z = 2.5 and 2.42 at z = 3.0, extrapolated |
| AP derivative | −[μ², 1−μ²] ∂(D P_w)/∂ln k at fixed μ (`fisher.py:135-167`) | P = G[K_iK_j(μ)P_sm + Q K_iK_j(μ′) D_ij(k′) P_w(k′)], Q = 1/(α∥α⊥²) (`kernels/kaiser.py:13-37, 58-78`; `models/kaiser.py:249-297`) | see below |
| wiggle split | polynomial (`fisher.py:169-198`) | PK − PKSB | information ratio 0.999–1.000 |
| damping | Σ⊥ = 3.26 σ8(z)/σ8(2.3), Σ∥ = (1+f)Σ⊥, /√R for pairs without Lyα (`fisher.py:156-161, 217-235`) | same widths per field; forest R = 1, galaxy R = 2; cross Σ² is the mean of the two squared widths (`survey_config.py:1303-1314`; `kernels/kaiser.py:48-55`) | Σ⊥ is 0.8 % above FR14 (9.4 σ8(z)/0.9 gives 3.234 at z = 2.3); R = 2 corresponds to FR14 r = 0.707 |
| response | σ = c/R (`covariance.py:177-183`; `spectrograph.py:298-303`) | σ = c/(2.355 R) (`response.py:51-54`; `survey_config.py:1376-1385`); the same pixel sinc | ME07 §II A FWHM definition |
| P1D | PD2013 (`analytic_p1d_PD2013.py:21-35`) | identical to 7 × 10⁻¹⁶ (`models/p1d.py:15-56`) | plateau below 5.7 × 10⁻⁴ s/km |
| fiducial total power in covariance | undamped smooth P3D + noise (`covariance.py:352-363`) | W_iW_j P_ij with damped wiggles + N (`forecast.py:162-165`) | < 0.3 % |

**Decomposition of the fh derivative** (C2 §1). ∂/∂ln α∥ of the wiggle term is the sum of four pieces:
- (a) −μ² KKD ∂P_w/∂ln k: the BAO phase shift.
- (b) −KKD P_w: the Q term, degenerate with the wiggle amplitude.
- (c) ∂(KK)/∂μ² · (−2μ²(1−μ²)) D P_w: the Kaiser remap.
- (d) +k∥²Σ∥² KKD P_w: the damping remap, exactly degenerate with ∂/∂ln Σ∥.

SE07 §2 and FR14 §4.1.1 keep only (a). Terms (b)–(d) partially cancel. Measured effect of the alternatives:
- SE07 form: +0.43/−0.27 combined, with per-bin σ∥ +0.23 to +0.71 and forest-auto σ⊥ up to +2.7 % (bin 6).
- lf-like form: −0.16/−0.52.
- After marginalising A_w, Σ, b, β and a per-pair broadband (285 parameters), fh, lf and SE07 agree to 0.1 %. This rests on C2's analytic calculation and is rated PLAUSIBLE by F1.

**Cross damping** (C2 §2; F1 F2):
- In the propagator picture, the cross wiggle damping is G_iG_j = √(D_iD_j). Arithmetic mean widths reproduce this exactly.
- For a reconstructed galaxy leg, the relative displacement argument gives Σ²_cross ≈ (Σ²_u + Σ²_rec)/2.
- This holds only if the cross is measured against the reconstructed catalogue.

### 3.5 Noise and covariance

| item | lf | fh | assessment |
|---|---|---|---|
| forest noise | `covariance.py:352-363` | `noise.py:78-90` | identical expression |
| galaxy shot noise | `covariance.py:384-390` | `noise.py:36-40` | identical |
| cross noise | none (`covariance.py:394-418`) | none; independent sampling (`noise.py:98-141`) | shared. The lya(qso)×qso same-object term has a correlation coefficient ≲ 1 %, and standard analyses exclude own-forest pairs (C2 §5) |
| Wick covariance and modes | (P_imP_jn + P_inP_jm)/N_modes, N_modes = V k² dk dμ/(2π²) (`fisher.py:121-131`; `covariance.py:244-249`) | same (`kernels/covariance.py:6-24`; q_mode at `grids.py:84-87`) | equals FR14 Eq. 20 and SE07 Eq. 1; D1 confirms with an FFT box count and a Monte Carlo test of the Wick formula |

FR14's high-noise Lyα+QSO combination σ = (σ_Lyα⁻¹ + σ_QSO⁻¹)⁻¹ is 10–30 % optimistic relative to the exact 3-spectrum Fisher at DESI-2 densities (C2 §5).

### 3.6 Derivatives and Fisher

| item | lf | fh | assessment |
|---|---|---|---|
| parameters | (α∥, α⊥) per bin (`fisher.py:30, 54-73`) | ap_i, at_i per bin (`survey_config.py:1206-1212`); `fix_except` principal block (`public.py:506`) | equivalent |
| finite difference | backward, scaled by dk/k (`fisher.py:163-165`; `power_spectrum.py:46`) | central, step 2.5 × 10⁻⁴ (`derivatives.py:59-104`; `public.py:493`) | fh converged to 10⁻⁶ (D1 check 16); lf ≈ 0.9 % optimistic (F4) |
| contraction | Σ dPᵀC⁻¹dP; `np.linalg.inv` (`fisher.py:58-71, 257`) | Cholesky Σ(L⁻¹J)ᵀ(L⁻¹J) (`fisher.py:127-203`; `forecast.py:269-272`) | agree to 1.6 × 10⁻¹⁵ (D1 check 18) |
| weights and ∂C | fixed; no ∂C/∂θ | same | FR14 Eq. 11 convention |

### 3.7 Selection and combination

| item | lf | fh | assessment |
|---|---|---|---|
| pairs | all 15 pairs in every bin (`forecast_new.py:60-85`) | bin 1: 3 spectra (INI `[pairs bin 1]`; `survey_config.py:900-929`) | fh's selection follows a standing user directive. The LBG/LAE tables have no support below z = 2.27; lf fills bin 1 with the clamped z = 2.38 row, i.e. a fictitious 107 LAE deg⁻² |
| individual spectra | own covariance including autos (`forecast_new.py:266-282`) | `public.py:278-339` | identical definitions; weights bit-identical to the joint (D2) |
| combination | per bin only | block-diagonal all-bin Fisher (`results.py:261-297`) | independent bins (FR14 Eqs. 11, 20). Radial modes with k∥ < 2π/Δχ ≈ 0.035 are correlated across bins; this is neglected, second order |

**Bin-1 cost** (E; F1):
- Adding the table-supported lya(lbg) and its crosses changes bin 1 by −0.54/−0.59 % (combined −0.05/−0.05).
- Adding all 15 spectra with clamped LBG/LAE changes bin 1 by −14.3/−17.9 % (combined −1.55/−1.81). This gain comes from the clamp.

## 4. Numerical results

### 4.1 Baseline (E Table 1)

| bin | z | fh σ∥ | fh σ⊥ | lf σ∥ (fhsel) | lf σ⊥ (fhsel) | fh/lf − 1 ∥ % | fh/lf − 1 ⊥ % |
|---|---|---|---|---|---|---|---|
| 1 | 2.00–2.23 | 0.02417 | 0.01811 | 0.02473 | 0.01804 | −2.24 | +0.35 |
| 2 | 2.23–2.47 | 0.01930 | 0.01410 | 0.02047 | 0.01462 | −5.71 | −3.60 |
| 3 | 2.47–2.71 | 0.01671 | 0.01172 | 0.01800 | 0.01237 | −7.13 | −5.25 |
| 4 | 2.71–2.94 | 0.01425 | 0.00936 | 0.01501 | 0.00966 | −5.10 | −3.19 |
| 5 | 2.94–3.17 | 0.01562 | 0.01054 | 0.01636 | 0.01080 | −4.48 | −2.35 |
| 6 | 3.17–3.41 | 0.02076 | 0.01616 | 0.02313 | 0.01741 | −10.25 | −7.18 |
| combined | | 0.007191 | 0.005020 | 0.007625 | 0.005205 | −5.77 | −3.53 |

- lf's native bin-1 joint, with all 15 spectra and clamped LBG/LAE, is 0.02134/0.01488.
- Against the grid-converged lf (`lf_conv_k2000`), fh/lf − 1 = −5.95/−4.04 % (E3 §6).

Key spectra, combined (fh/lf − 1, ∥/⊥ %):

| spectrum | fh/lf − 1 |
|---|---|
| lya(qso) auto | +1.90/+3.53 |
| lya(lbg) auto | −27.94/−32.17 |
| qso auto | −0.39/+0.15 |
| lya(qso)×qso | −6.89/−2.14 |
| lbg auto | −0.14/+0.38 |

Per bin, the lya(qso) auto goes from +3.1/+5.1 (bin 1) to −8.8/−10.3 (bin 6). The lya(lbg) auto is −24 % to −42 % in every bin.

Forest weighting integrals, lf/fh:

| forest | quantity | lf/fh, bins 1–6 |
|---|---|---|
| lya(qso) | A | 1.07–1.17 |
| lya(qso) | P_pix | 0.86 → 0.51 |
| lya(lbg) | A | 0.62–1.05 |
| lya(lbg) | P_pix | 1.59–2.08 |

### 4.2 Compatibility profiles (F1 F5)

`fh_full_compat` reproduces `lf_baseline` to 1.8 × 10⁻¹⁵ in every joint and key spectrum. This is trivial with respect to the physics:
- The adapter (`adapters/desi2_compatibility.py:1-7`) captures lf's total power, mode counts, signals, damping widths and `legacy_jacobian`.
- fh supplies only the Wick assembly and the Fisher contraction.

`fh_fixed_compat` uses early weights on lf's arrays and numerics. It sits at −2.52/−2.48 % from lf and +3.45/+1.09 % from `fh_baseline`, so lf's non-weight numerics and inputs account for +3.5/+1.1 %.

### 4.3 Attribution (E §2; fig. 2)

**Isolated toggles** (vs `fh_baseline`, combined):

| toggle | Δσ∥/Δσ⊥ % | largest per-bin effect |
|---|---|---|
| legacy-3 weights | +3.12/+3.02 | bin 1 −1.5; bin 6 +8.3/+6.9; lya(lbg) auto +41.8/+53.0; lya(qso) auto −2.7/−3.6 |
| McDonald, cap 2000 | −0.19/−0.17 | ≤ 0.55 per bin |
| c/R resolution | +0.23/+0.05 | lya autos +0.5/+0.2 |
| lf damping rule | +3.31/+1.65 | per-bin σ∥ +2.46 to +4.30; lya(qso)×qso +8.46/+4.43; autos 0 |
| EdS growth | +0.19/+0.25 | bin 1 −0.26/−0.37; bin 6 +0.55/+0.70 |
| arithmetic z | −0.08/−0.08 | ≤ 0.26 |
| AP Q = 1 | +0.05/−1.17 | lya(qso) auto σ⊥ −2.15 |
| AP SE07 | +0.43/−0.27 | per-bin σ∥ ≤ +0.71; lya(qso) auto σ⊥ up to +2.7 (bin 6) |
| AP lf-like | −0.16/−0.52 | – |

**Cumulative chain** (vs `lf_baseline`):

| step | combined Δ vs lf % | step % |
|---|---|---|
| fh baseline | −5.77/−3.53 | – |
| + legacy-3 weights | −2.83/−0.62 | +3.12/+3.02 |
| + c/R | −2.65/−0.58 | +0.19/+0.04 |
| + lf damping | +0.38/+0.89 | +3.11/+1.47 |
| + EdS | +0.59/+1.15 | +0.21/+0.26 |
| + arithmetic z | +0.53/+1.10 | −0.06/−0.06 |
| + lf-like AP (`fh_cum6_lf`) | +0.35/+0.55 | −0.17/−0.54 |
| (alternative endpoints) Q = 1 / SE07 | +0.56/−0.04 ; +0.97/+0.81 | |
| `fh_cum6_lf` vs converged lf | +0.17/+0.02 | |

Per-bin closure is +0.14 to +0.53 (∥) and +0.28 to +0.69 (⊥).

Along the chain, the lya(lbg) auto ends at +3.44/+5.05 %. There, fh's legacy-3 P_pix lies 6–7.5 % above lf's in bins 3–6, with A within 3 %. The residual therefore lies in S and B, the auxiliary spectra, and not in the update rule.

### 4.4 Variations (Δ vs own-code baseline, combined; E tables)

| group | case | fh | lf | comment |
|---|---|---|---|---|
| survey | area 3000 / 8000 deg² | +29.10 / −20.94 | same | exact A⁻¹ᐟ² scaling to 1.7 × 10⁻¹⁵ |
| survey | r_max 23.75 | +10.00/+18.55 | +11.73/+20.15 | |
| survey | r_max 24.25 | +1.08/+1.54 | +0.32/+0.39 | |
| survey | r_max 25.75 / 26.25 | 0 | ≤ 2.5 × 10⁻⁴ | null tests |
| survey | densities ×0.5 | +36.25/+49.14 | +40.25/+52.66 | |
| survey | densities ×2 | −20.92/−27.45 | −22.76/−28.56 | |
| survey | N_exp 2 | +5.07/+4.84 | +4.48/+4.16 | lya(qso) auto +14.6/+18.9; weak response because of aliasing |
| survey | N_exp 8 | −5.53/−5.86 | −5.57/−5.56 | |
| instrument | R 2000 / 4000 | +0.03 / −0.03 | +0.10 / −0.11 | |
| instrument | pixel 0.4 / 1.6 Å | ≤ 10⁻⁴ | ≤ 10⁻⁴ | null |
| instrument | 1040–1185 Å | +1.94/+1.94 | +1.92/+1.86 | |
| instrument | 1050–1205 Å | +0.64/+0.62 | +0.55/+0.51 | |
| analysis | k_max 0.3 / 0.4 | +1.10/+1.62 ; +0.04/+0.09 | +1.06/+1.76 ; +0.08/+0.06 | |
| analysis | reconstruction off | +10.42/+6.38 | +7.72/+5.02 | |
| analysis | β_F 1.2 / 1.7 | +2.27/+0.49 ; −2.11/−0.48 | – | |
| analysis | b_Lyα ×0.85 / ×1.15 | +4.06/+3.90 ; −3.58/−3.69 | – | |
| analysis | weighting mode ×0.5 / ×2 / ME07 | +1.34/+1.30 ; −1.28/−1.30 ; −1.24/−1.26 | – | forest autos 4–9 % |
| numerics | k/μ/magnitude order, AP step, weight rtol | < 10⁻³ % | – | fh converged |
| numerics | lf num_k 250 / 1000 | – | −0.87/−0.78 ; +0.45/+0.40 | lf FD dk/k |
| numerics | lf n_μ 5 / 20 | – | +2.19/+0.79 ; −0.53/−0.19 | |
| numerics | lf n_mag 54 / 214 | – | +0.75/+1.15 ; +0.17/+0.16 | non-monotone |

**Expectation checks** (E §2):
- A⁻¹ᐟ² area scaling: pass.
- Monotonicity in density, N_exp, R, k_max, reconstruction, area, r_max, forest range and σ_z: all pass for fh (14/14 sub-checks each). lf fails the r_max sequence in bin 1 at the 10⁻⁴ level (107-node grid).
- fh pixel null tests: ≤ 7.6 × 10⁻⁴.

**fh/lf ratio over 20 matched variations** (E §3; fig. 5):
- R∥ = 0.9432 ± 0.0101 and R⊥ = 0.9651 ± 0.0078.
- The largest movers are density ×0.5 (−2.85/−2.31 %), reconstruction off (+2.50/+1.29 %), density ×2 (+2.39/+1.56 %) and r_max 23.75 (−1.55/−1.33 %).
- The driver is the lya(lbg) auto, whose ratio ranges over 0.62–0.83. This is the prefix-I1 effect in the noise-dominated LBG forest.

### 4.5 Follow-ups (E §4; Δ vs `fh_baseline`, combined)

**Weights.**
- Inverse variance: −1.83/−1.85 %.
  - Per-bin joint: −5.1/−5.4 in bin 1, −0.4/−0.3 in bin 5.
  - Forest autos: lya(qso) −6.8/−9.2 and lya(lbg) −7.4/−10.0. The lya(qso) auto in bin 1 is −10.0/−13.1.
  - A rises by up to ×1.45 (QSO bin 1) while P_pix falls by up to ×0.50.
- Inverse variance with B at k∥ = 0.1 is identical to within 3 × 10⁻⁵. This is a weak test, since B changes by only 0.6–2.6 % (F1 W5).

**Reconstruction.**
- FR14 r(n̄P) per tracer: −2.82/−1.61 %. The QSO auto moves from −5.2 % (bin 1) to +5.8 % (bin 6).
- Global R = 3.43: −5.07/−2.85 %.

**Galaxy σ_z.**
- 2 h⁻¹Mpc: +2.55/+0.73 %.
- 4 h⁻¹Mpc: +9.14/+2.68 %; bin 5 +12.9/+3.3; QSO auto +19.7/+4.4.

**Forest k∥ cut** (applied to spectra that contain a forest).

| k_c (h Mpc⁻¹) | joint Δ | forest autos σ⊥ | forest autos σ∥ |
|---|---|---|---|
| 0.01 | +0.37/+1.95 | +6.5 | +1.4 |
| 0.015 | +0.54/+2.84 | +9.8 | +2.1 |
| 0.02 | +0.73/+3.91 | +14.2 | +3.0 |
| 0.03 | +1.09/+5.86 | +23.5 | +5.2 |

**Source redshift averaging** (both forests, clamped tables): −1.12/−0.90 %. The lya(qso) auto changes by +0.4 (bin 1) to −23.8/−25.4 (bin 6). The lya(lbg) auto changes by +17.6/+23.5 (bin 3) and −30.9/−36.7 (bin 5).

**Density normalisation.**
- z_norm_min = 2.1: +2.03/+2.59 %; bin-1 joint +7.7/+10.7; QSO auto +14.0/+15.0.
- z_norm_min = 2.2: bit-identical to the baseline.

**P1D.**
- ×0.8: −1.10/−1.10 %; lya(qso) auto −7.0/−8.4.
- ×1.2: +0.93/+0.91 %.
- No plateau: −0.15/−0.70 %.

### 4.6 Corrective runs (E3) and the earlier estimates they replace

| question | earlier result | corrected (E3) |
|---|---|---|
| k⊥ Nyquist | `fh_r6_nyq`: null. It was mis-specified with n2d = 90 deg⁻², giving k_Nyq = 0.37–0.40 (F1 A1c REFUTED) | Per-bin n_2D (k_Nyq = 0.389, 0.303, 0.246, 0.183, 0.127, 0.078 for lya(qso)): +0.40/+1.23; bin 6 +1.66/+3.00; lya(qso) auto σ⊥ +2.08. Both forests: +0.40/+1.24. MW11 n̄_eff: QSO forest only +0.76/+2.68; both forests +1.29/+5.52 (bin 6 +3.07/+10.27) |
| LBG density beyond the table | `fh_dzero_340` +3.76/+2.69; `fh_dzero_329` +279/+365 in bin 6. The latter is a node-edge artefact that also removed LBG/LAE galaxies (F1 S1) | range-averaged over z_s and truncated at 3.40 (lya(lbg) only): +2.12/+1.51; bin 5 −1.90/−1.68; bin 6 +21.67/+16.75; lya(lbg) auto +11.4/+12.2. Range-averaged with the clamp: −0.55/−0.39. The clamp lowers the combined σ by ≈ 2.7/1.9 % relative to the truncated treatment |
| QSO bias in b·D terms | literal b = 1.2: +5.40/+9.40 (mis-specified, F1 F6) | b_QSO D = 1.2: +1.19/+2.14; QSO auto +16.2/+19.9 |
| LAE bias in b·D terms | literal b = 0.89: +6.09/+9.45 (mis-specified) | b_LAE D = 0.89: −3.04/−5.12; bin 2 −7.45/−12.79; LAE auto −22/−28 |
| forest bias | eBOSS b_F = −0.117, β_F = 1.67: +1.80/+3.07 (valid) | DESI DR1 b_F = −0.108, β_F = 1.743: +3.20/+4.71; forest autos +15.2/+29.6 |
| LBG/LAE cell widths | global rescale: +1.10/+1.49 | per-row widths: +0.91/+1.35; bins 3–5 +0.8 to +1.7 in ∥, bins 1 and 6 unchanged |
| lf grid convergence | lf numerical floor ±1 % (F4) | num_k 1000: −0.04/+0.33. num_k 2000, n_μ 40, 214 magnitudes: +0.19/+0.53. fh/lf is then 0.941/0.960 and the closure is +0.17/+0.02 % |

**Other corrected earlier estimates.** The values in the second column are superseded; the verdicts are F1's.

| earlier estimate | superseded by |
|---|---|
| resolution change ~1 % (B row 9) | +0.23 % |
| lf μ midpoints make lf ~2 % optimistic in σ⊥ (B row 23) | lf is pessimistic in σ∥ |
| AP toy fh/lf = 1.009/1.015 (B row 26) | not reproduced; forecast gives 0.999/0.997 |
| mixed damping 1–2 % (B row 8) | +3.3 % |
| early-weight excess 0.1–1 % (D1 F1) | 0.2–22 % in N_F, 1.8 % in σ |
| forest autos "cosmic-variance limited" (E) | refuted; they are noise dominated |
| weighting-mode insensitivity of 3 × 10⁻⁵ (W5) | a conflation; the FKP baseline carries a ±1.3 % mode dependence |
| continuum cut effect σ∥ < 1 % (C3 toy) | forest autos +2–5 % |
| "no Nyquist cut needed" (C3) | the papers disagree (FR14 vs MW11); E3 gives the size |
| lf k grid irrelevant (C2) | lf's dk/k derivative scaling makes it ≈ 0.9 % |
| McDonald ≈ 300 updates (C1) | 624, since the returned state is 2t |

### 4.7 Figures (`appendix/figures/`)

- fig1: per-bin joint σ for fh baseline, lf, fixed-compat and the full chain, with residuals.
- fig2: cumulative chain waterfall (combined, relative to lf), with the Q = 1 and SE07 endpoints.
- fig3: Δσ/σ of all fh variations (symlog), coloured by group.
- fig4: lya(qso) and lya(lbg) auto σ per bin for the early, inverse-variance, McDonald, legacy-3, ME07-mode and lf weightings.
- fig5: σ_fh/σ_lf over the matched variations.

The figures predate E3 and do not include the E3 cases.

## 5. Verification (F1, condensed)

| ID | claim | F1 verdict | resolution |
|---|---|---|---|
| W1 | lf prefix I1 is a bug; LBG N_F 1.6–2.1×; the recurrence diverges | CONFIRMED | dominant lf→fh change |
| W2 | lf omits aliasing from S despite the comment; FR14-literal | CONFIRMED | – |
| W3 | MW11 Eq. 11 is minimum variance; early is 0.2–22 % above it; −1.8 % joint | CONFIRMED (range corrected) | caveat: A ×1.45 moves weight into the unmodelled MW11 Eq. 10 terms |
| W4 | McDonald is slow, not collapsing (624 updates) | mechanism CONFIRMED; "decision record incorrect" DOWNGRADED to incomplete | §6 item 6 |
| W5 | weighting mode changes σ by 3 × 10⁻⁵ | REFUTED as stated | early weights: ±1.3 % |
| S1 | LBG clamp in bins 5–6; truncation +279 % | clamp CONFIRMED (shared); 279 % REFUTED as an attribution | E3 +2.12/+1.51 |
| S2 | first spacing +4.8/+2.6 %; z_norm_min excluded "only via float" | first CONFIRMED; second REFUTED (threshold on a cell centre) | E3 per-row +0.91/+1.35 |
| S3 | single-z_src biases QSO-forest noise by up to 24 % | CONFIRMED | −1.12/−0.90 joint |
| F1 | full-AP derivative within ≤ 0.5 % of SE07/FR14; 0.1 % after marginalisation | CONFIRMED (≤ 0.7 % per bin); marginalised PLAUSIBLE | no forecast |
| F2 | mixed damping; lf rule +2.5–4 % | numbers CONFIRMED; physics PLAUSIBLE, conditional | analysis-plan choice |
| F3 | FWHM convention correct; lf c/R a bug; +0.23 % | CONFIRMED | – |
| F4 | lf μ grid pessimistic, k grid optimistic | CONFIRMED | E3: converged lf +0.19/+0.53 |
| F5 | chain closes to +0.35/+0.55 %; 10⁻¹⁵ compat match | CONFIRMED; the match is trivial for the physics | E3: +0.17/+0.02 |
| F6 | b_F, β_F ~15 % off; literal-b runs mis-specified | CONFIRMED (HCD caveat) | E3 b·D and DESI DR1 runs |
| A1a | no σ_z damping: +2.6/+9.1 % | CONFIRMED | – |
| A1b | no continuum k∥ cut: σ⊥ +2–6 % | CONFIRMED | – |
| A1c | no Nyquist cut: null | REFUTED (mis-specified) | E3 +0.40/+1.23 to +1.29/+5.52 |
| B1a | `--joint-only` crash | CONFIRMED | fh-1 |
| B1b | stale rtol | CONFIRMED, no numerical effect | fh-3 |
| B1c | 7 pytest failures | 5 environmental, 2 genuine | fh-2 |
| B1d | latent per-category k_max in BAO | REFUTED as reachable | – |
| E-x | forest autos signal limited | REFUTED | – |

## 6. Documentation mismatches

**fishhighz.**

1. `docs/methods/weights.md:46` calls early "the recommended research method". This is unsupported: N_F lies 0.2–22 % above the minimum, and inverse variance gives −1.8 % in σ. The statement that inverse variance "does not establish a multi-mode BAO optimum" is over-cautious, because within fh's covariance the information is monotone in N_eff.
2. `docs/methods/weights.md` calls I1/I2/I3 "prefix integrals". In the full-sum methods only the totals enter.
3. The `weights.py` and `kernels/weights.py` docstrings say the cumulative formulas "implement ME07". The prefix (legacy) recurrence is not ME07's, since ME07 Eq. 14 integrates to m_max.
4. `docs/research/FOREST_WEIGHTING_DECISION.md` §3 attributes the aliasing-in-P_S prescription to FR14 Eqs. 23–29. FR14 does not state it; only ME07 does.
5. Decision record §2 and §6 do not note two facts. First, the S → 0 limit of the early form is the MW11 minimum-variance weight. Second, at the chosen mode S L I1/B = 0.3–21, far from that limit.
6. Decision record §5 Step 2 labels McDonald "unavailable" in LBG bin 2. On the production inputs it converges at 624 updates. The §3 and §5 caveats ("does not prove the asymptotic behavior"; "unresolved within the tested cap") are correct as written, so F1 rates the record incomplete rather than incorrect.
7. `docs/research/RESEARCH_BASELINE.md:239-242` states the default rtol is 10⁻⁴. The native path uses 10⁻⁵ (fh-3).
8. RESEARCH_BASELINE presents (2.4 deg⁻¹, 3.5 × 10⁻⁴ s/km) as "the representative mode". It is not ME07/FR14's (0.07, 0.5), and it controls a ±1.3 % effect on the joint σ.
9. RESEARCH_BASELINE lists the AP prefactor, the remapped Kaiser factor and the remapped damping in the derivative. It does not state that SE07 §2 and FR14 §4.1.1 exclude these terms as non-BAO information (≤ 0.7 % per bin).
10. RESEARCH_BASELINE writes the reconstruction as Σ/√r_i with r_i = 2. FR14's "50 % reconstruction" (Σ × 0.5) corresponds to r = 4 in fh units; r = 2 corresponds to FR14 r = 0.71.
11. RESEARCH_BASELINE and `docs/methods/kaiser.md` state the cross-damping rule without its physical basis, and without the condition that crosses use reconstructed catalogues.
12. `docs/methods/survey.md` says the strict readers reject extrapolation and have no floors. The native INI path uses `adapters/legacy_compat` (z clamp, floor, bright clamp, sentinel, first spacing), documented only in `INI_DEFAULTS` (`accuracy.py:36-87`).
13. `docs/methods/geometry.md` claims support for curved FLRW, but the volume element lacks the curvature factor (fh-4).
14. `docs/user/cli.md:15-19` omits `--joint-only` (fh-1). `docs/user/python.md` omits the `individuals` argument of `run()`. RESEARCH_BASELINE.md:256-258 calls INI/CLI/serialisation out of scope, although all three are implemented.
15. `data/README.md` and the INI comments omit several facts:
    - the LBG/LAE tables cover only z = 2.27–3.40 (5 rows);
    - the counts vanish for r > 23.5 (QSO) and r > 24.25 (LBG/LAE), so r_max ≥ 24.75 and `density_magnitude_bounds` are inert;
    - the SNR-header `--ab-magnitude` lines are stale;
    - EXPTIME is ambiguous between per exposure and total: the LBG tables were generated with `--total-exptime 4000 --nexp 4`, the QSO tables with `--nexp 4` only.
16. `grids.py:6-7` has a stale docstring saying no Fisher is implemented.

Documentation that agrees with the code: `covariance.md`, `fisher.md` (half-sphere 1/(2π²) normalisation), `kaiser.md` (smooth part never damped; D at wiggle coordinates), and `survey.md` `sample_forest_readers`.

**lyaforecast.** In `weights.py`:
- the comment "weights include aliasing as signal" contradicts the code (lf-4);
- the "(McDonald & Eisenstein 2007)" docstrings on the prefix recurrence are incorrect;
- the comment `kt_w_deg = 2.4  # ~0.035 h/Mpc` corresponds to k⊥ = k∥, not to μ = 0.5.

## 7. Limitations of this review

1. **Fisher level only.** Every result is a Fisher forecast in a fixed Planck18 fiducial cosmology. No mock, likelihood or correlation-function validation was attempted.
2. **Harness toggles.** Toggles and corrections were implemented as harness monkeypatches, not package options. F1 verified the σ_z and mask implementations; the others are checked through the chain closure and the null and monotonicity tests.
3. **Unforecast estimates.** The Σ/broadband-marginalised offset (+1.3–1.5 %) is analytic. The MW11 Eq. 10 covariance terms were not computed. The inverse-variance gain may therefore be partly offset.
4. **External inputs.** The density and SNR tables were taken as given. Targeting efficiency, purity, redshift success and the EXPTIME semantics cannot be established from the repository.
5. **Bias references.** Bias comparisons use eBOSS DR16 and DESI DR1 pure-forest values (F1, from the arXiv sources). DESI DR2 was not consulted. The HCD-inclusive large-scale bias (≈ −0.164) exceeds fh's, so "fh optimistic" is not established.
6. **Barrier cost.** The barrier excluded the archived S1–S5 evidence, including the author's S4 attribution. Cross-checks against it are limited to what `docs/research/` quotes.
7. **Common-α combination.** The "combined" error with a common α is a review construct. fh's native combined Fisher keeps per-bin α and is block-diagonal.
8. **Figures.** They do not include the E3 cases.

## 8. Recommendations (scientifically consequential items only)

| # | recommendation | expected impact on joint σ | smallest check |
|---|---|---|---|
| 1 | Adopt MW11 Eq. 11 inverse-variance weights as the baseline, or report them alongside early | −1.83/−1.85 %; forest autos −7 to −10 %; removes the ±1.3 % reference-mode dependence and the iteration and cap issues | already run (`fh_r1_invvar`). Add an evaluation of the MW11 Eq. 10 (4/N) w̄² P_F P_F′ term with the inverse-variance A (×1.45 in QSO bin 1), to show the gain is not an artefact of the diagonal covariance |
| 2 | Make the forest × galaxy damping rule an explicit analysis-plan input; use lf's rule for Lyα×QSO without QSO reconstruction | +3.31/+1.65 % if switched | decision only (`fh_t_damp_lf` exists) |
| 3 | Include a transverse information limit: Nyquist with n_2D or n̄_eff, or the MW11 beat-coupling covariance | +0.40/+1.23 % (n_2D) to +1.29/+5.52 % (n̄_eff); bin 6 σ⊥ up to +10 % | compute the MW11 Eq. 10 off-diagonal terms for lya(qso) in bins 5–6 to choose between the two cuts |
| 4 | Replace constant extrapolation of the LBG density beyond z = 3.29 by truncation or a documented luminosity-function extrapolation; average over z_s; record the extension flag in the saved output | +2.12/+1.51 %; bin 6 +21.7/+16.8 % | done (E3); needs an input decision or an extended table |
| 5 | Include the continuum-fitting k∥ loss, cutting at k∥ ≳ 2π/L or using MW11's continuum noise | +0.54/+2.84 % at k∥ > 0.015 | done (`fh_kpar0p015`); the mask becomes a package option |
| 6 | Add per-tracer galaxy redshift-error damping (SE07 Eq. 27) | +2.55 to +9.14 % in σ∥ for σ_z = 2–4 h⁻¹Mpc | done; needs measured σ_z per tracer (DESI QSO, LBG, LAE) |
| 7 | Update b_F and β_F to DESI DR1/DR2 values, with a consistent HCD treatment | +3.20/+4.71 % (DESI DR1 pure forest) | done; add an HCD-bias variant |
| 8 | Use an n̄P-dependent reconstruction r (FR14), per tracer or for the joint galaxy field | −2.82/−1.61 % (per tracer) to −5.07/−2.85 % (R = 3.43); +10.4/+6.4 % if off | done; a decision on the reconstruction plan |
| 9 | Implement Σ-marginalised BAO (width derivatives) and report it as a companion | +1.3–1.5 % (analytic) | a one-bin forecast with Σ∥ and Σ⊥ free per field |
| 10 | Define the density targets physically (z range, efficiency); use per-row cell widths; move z_norm_min off a cell centre | +0.91/+1.35 % (widths); +2.03/+2.59 % (z_norm_min = 2.1) | done; documentation and input decision |
| 11 | Repair the fh-2 fixture so that the independent accuracy check runs | none on σ; restores validation coverage | update `tests/test_step12_performance.py:161,174` to set `p3d.selection` |

