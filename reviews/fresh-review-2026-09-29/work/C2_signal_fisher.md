# C2: signal model, BAO Fisher formalism, geometry and covariance (fishhighz vs lyaforecast)

Scope: DESI-2 15×2pt BAO forecast (`fishhighz/data/desi2_accuracy.ini`, 5 fields, 6 bins 2.0–3.41). "fh" = fishhighz, "lf" = lyaforecast NewForecast. Paper references:
- ME07 = McDonald & Eisenstein 2007. Equation text verified against the arXiv LaTeX source astro-ph/0607122.
- MW11 = McQuinn & White 2011. From the A2 digest only, not re-verified.
- FR14 = Font-Ribera et al. 2014. Eq. numbers 15, 16, 19, 20, 21–29 verified from the arXiv LaTeX source 1308.4164.
- SE07 = Seo & Eisenstein 2007. From the arXiv LaTeX source astro-ph/0701079. Equation numbers are my count of numbered rows in that source:
  - (3)–(4): damped BAO.
  - (5): full Fisher with damping outside the derivative.
  - (26): 2-D fitting formula.
  - (27): R(k,μ) with redshift errors.

Quantitative checks are in `/tmp/claude-81394/.../scratchpad/C2/`, in the scripts `prep.py`, `c2lib.py`, `check_fd.py`, `deriv_variants.py`, `marg.py`, `margjoint.py`, `grid_check*.py`, `np_recon.py`, `damp_rules.py`, `fog.py`, `combo.py`, `cuts*.py`, `crossnoise.py`, `growth.py` and `bias_cmp.py`. They use the real prepared fh bins: `Forecast(desi2_accuracy.ini).prepare()`, including the converged forest weights, noise, responses, GL nodes and mode counts. On top of those bins they use an independent analytic re-implementation of the wiggle model. I validated this analytic Jacobian against a central finite difference of the fh kernels (`kernels/kaiser.py`), with a maximum relative difference of 6e-7. The fiducial power is reproduced exactly. No forecast was run, i.e. no `Forecast.run`.

## Summary

The fh signal/Fisher core is a defensible upgrade of lf. The key worry is that the full AP remap of the wiggle component retains non-BAO information. It is correct in principle but numerically negligible.
- At fixed nuisances, fh differs from a Seo–Eisenstein/FR14 derivative by ≤0.5% for the 15×2pt joint fit, and by ≤1% for galaxy autos.
- Forest autos differ by up to −2.6% (σ⊥, Lyα(QSO) bin 6).
- After marginalising the wiggle amplitude, Σ∥/Σ⊥, biases, RSD and a per-pair broadband (285 parameters), fh, lf and the SE07 form agree to 0.1% for the joint fit. All three give errors that are 1.3–1.5% larger than the unmarginalised fh baseline.

Scale of the individual changes:
- The largest single change is the mixed forest–galaxy damping rule, at 2.5–4% on σ∥ in the joint fit. It is physically better motivated than lf's rule.
- The reconstruction model is the most uncertain input. A fixed R=2 corresponds to FR14's r=0.71 and is neither conservative nor optimistic uniformly. FR14's r(n̄P) gives −1 to −7%.
- Resolution (FWHM convention), growth (CAMB σ8 instead of EdS) and the geometric-mean z_eval are correct choices with ≤1% effect.
- One inventory claim is wrong in sign. lf's 10 μ-midpoints make lf pessimistic in σ∥, by 0.3–2.3%, not optimistic in σ⊥.
- One inventory toy ratio is not reproduced: σ_fh/σ_lf = 1.009/1.015.

Two approximations shared by both codes each have a larger effect than any lf→fh change:
- The absence of galaxy redshift-error/FoG damping: +2 to +13% σ∥ for σ_z = 2–4 h⁻¹Mpc.
- The absence of low-k∥ forest continuum loss: +6% σ⊥ for a forest k∥ > 0.01 cut.

## 1. BAO derivative

**fh.** The observable is `P = G[K_iK_j(μ) P_sm(k) + Q K_iK_j(μ') D_ij(k∥',k⊥') P_w(k')]` (`kernels/kaiser.py:58-78`, `models/kaiser.py:249-297`).
- Q = 1/(α∥α⊥²).
- The mapped coordinates are k∥' = kμ/α∥ and k⊥' = k√(1−μ²)/α⊥.
- The smooth component is not dilated.
- There are no nuisances.
- The derivative is a central finite difference with step 2.5e-4.

Analytically, ∂P_w-term/∂lnα∥ is the sum of four pieces:
- (a) −μ² KKD dP_w/dlnk. This is the phase shift, i.e. BAO.
- (b) −1·KKD P_w. This is Q.
- (c) ∂(KK)/∂μ² · (−2μ²(1−μ²)) D P_w. This is the Kaiser remap.
- (d) +k∥²Σ∥² KKD P_w. This is the damping remap.

For α⊥ the corresponding factors are −(1−μ²), −2, +2μ²(1−μ²) and k⊥²Σ⊥².

Terms (b)–(d) are all proportional to P_w, i.e. in phase with the wiggles and nearly orthogonal to (a).
- Term (d) is exactly degenerate with ∂/∂lnΣ, because ∂lnD/∂lnα∥ = −∂lnD/∂lnΣ∥ (derived).
- Term (b) is exactly degenerate with a free wiggle amplitude.

**lf.** The derivative is −[μ², 1−μ²] d(D·P_w)/dlnk at fixed μ (`fisher.py:135-167`).
- The Kaiser factor is outside the derivative and dμ/dα is neglected, consistent with SE07.
- It includes the damping-envelope derivative P_w dD/dlnk. SE07 §2, text before its eq. (5), explicitly excludes this: "we take the exponential factor … outside of the derivatives … equivalent to marginalizing over a large uncertainty in Σ⊥ and Σ∥". FR14 (text after its eq. 16) does the same: "damping factors, along with the RSD factor, are taken outside the Fisher matrix derivatives".
- lf de-wiggles with an 8th-order polynomial fitted to the full model.
- lf uses a backward difference normalised by dk/k_i instead of ln(k_i/k_{i−1}).

**Quantitative** (`deriv_variants.py`, `marg.py`, `margjoint.py`). The table gives σ/σ(SE-like), where "SE-like" is term (a) only, on the fh grid with prepared covariance.

| case | fh (a–d) | fh, Q=1 | lf-like d(DP_w)/dlnk |
|---|---|---|---|
| joint 15×2pt, bin 3, ∥/⊥ | 0.995/1.002 | 0.995/0.989 | 0.993/0.997 |
| joint, bin 6 | 0.997/1.002 | 0.998/0.991 | 0.996/0.998 |
| QSO auto, bins 2–6 | 0.999–1.001 | 0.996–1.000 (⊥) | 0.996–0.998 |
| Lyα(QSO) auto, bin 3 | 0.991/0.991 | 0.991/0.970 | 0.992/0.996 |
| Lyα(QSO) auto, bin 6 | 0.991/0.974 | 0.991/0.962 | 0.995/0.997 |

The terms partially cancel. Q (b) raises σ⊥, while the Kaiser (c) and damping (d) remaps lower it. Q=1 alone lowers the joint σ⊥ by about 1.1%, and forest autos by 2–3%.

Marginalisation test, 15×2pt, bins 3 and 6. The table gives σ/σ(SE-like, no nuisances).

| nuisances | SE-like | fh | lf-like |
|---|---|---|---|
| none | 1.000/1.000 | 0.995/1.002 | 0.993/0.997 |
| A_w | 1.003/1.005 | 0.996/1.002 | 0.997/1.003 |
| Σ∥,Σ⊥ per field | 1.013/1.013 | 1.009/1.012 | 1.013/1.013 |
| A_w+Σ+b_i+(f,β_i) | 1.014/1.014 | 1.012/1.013 | 1.014/1.014 |
| + per-pair broadband (19 terms/pair; mult. k^n μ^{2l} on P_sm + additive k^{−1..1}μ^{0,2,4}) | 1.014/1.015 | 1.013/1.014 | 1.014/1.015 |

For single tracers the pattern is the same:
- Galaxies: fh equals SE-like to ≤0.1% once A_w and Σ are free.
- Forests: fh stays 0.3–0.5% below SE-like, from the Kaiser remap term (c). This term is wiggle-shaped and is not absorbed by smooth broadband terms.

**Verdicts.**
- Full AP remap of the wiggle component with Q, μ' and remapped damping: **questionable in principle, justified in practice**.
  - Term (d) uses the damping anisotropy as a ruler at fixed Σ, which SE07 and FR14 explicitly avoid.
  - Term (c) is the "artificial (non-BAO-distance) breaking of the degeneracy" that FR14 (Lyα paragraph, after its eq. 29 block) avoids by dividing the noise by (1+βμ²)² rather than putting the RSD factor in the derivative.
  - Term (b) is BAO-amplitude information.
  - The net effect is ≤0.5% on the joint, ≤1% on galaxies, and ≤2.6% on σ⊥ of individual forest autos, where the sign favours fh.
  - After marginalisation the three prescriptions coincide.
  - Stating that fh's 2-parameter result is "BAO-only" is acceptable only with this ~0.5% caveat.
  - D1's statement that the wiggle-only remap is "equivalent to marginalising a free broadband shape" is not exact. Fixing the smooth component is not marginalisation. It is numerically equivalent to 0.1% only once Σ and A_w are also free.
- Smooth component unscaled: **justified**. It is the standard BAO-only construction: the broadband carries no dilation information by construction. It is not the FR14/SE07 construction, which removes broadband distance information by Ω_b subtraction. The effect is equivalent within 0.1% (table above).
- No broadband/bias/Σ marginalisation: **questionable**, but uniform and small. All prescriptions are 1.3–1.5% optimistic relative to the marginalised result, dominated by Σ. So SE07's "outside the derivative ≡ marginalising Σ" statement is itself only accurate to about 1.3%. The `bao_marginalized` mode frees only b_i and β_Lyα. My test shows that this changes σ by <0.1%, so it cannot decide this question.
- Central FD: **justified**. D1 found it converged to 1e-6. It also removes lf's backward-difference half-bin offset and dk/k mis-scaling. lf-exact (polynomial de-wiggling + backward FD on the lf grid) gives σ 0.1–0.9% smaller than lf-analytic, i.e. lf was slightly optimistic.
- Inventory toy check (row 26: σ_fh/σ_lf = 1.009 ∥, 1.015 ⊥ for a Lyα auto): **not reproduced**.
  - Like-for-like (same grid and covariance), I find fh/lf = 0.999/0.995 for Lyα(QSO) in bin 3, and 1.003/1.003 for QSO.
  - Including lf's grid and FD (fh on the fh grid over lf-exact on the lf grid), I find 0.991/0.993 for Lyα(QSO) bin 3 and 1.003/1.011 for the joint 15×2pt in bin 3.
  - The inventory number probably mixes the derivative change with its (sign-incorrect) μ-grid effect; see §8.

## 2. Damping and reconstruction

- **Σ⊥ = 3.26 σ8(z)/σ8(2.3), Σ∥ = (1+f)Σ⊥** (`survey_config.py:1305-1314`). FR14 eq. 16 and the following text give Σ⊥ = 9.4 σ8(z)/0.9 and Σ∥ = (1+f)Σ⊥. These follow from SE07's Σ0 = 12.4 h⁻¹Mpc for σ8 = 0.9 with G(0) = 0.758 (SE07 §2). With CAMB Planck18, σ8(2.3) = 0.3097, so FR14 gives 3.234 at z = 2.3. fh's normalisation is therefore 0.8% larger than FR14's, i.e. 0.8% more damping, and the z-scaling is identical. Across the bins, Σ⊥,unrec(z_eval) = 3.45, 3.21, 3.01, 2.83, 2.67, 2.52 h⁻¹Mpc. Verdict: **justified** (equivalent to FR14).
- **Reconstruction Σ/√R, R=2, for galaxies only.** FR14 ("50% reconstruction") multiplies Σ by r(n̄P), with r = 0.5 at high density and r interpolated from x = n̄P(0.14, 0.6)/0.1734. fh's R=2 is r = 1/√2 = 0.707, which FR14 assigns at x ≈ 1, i.e. n̄P ≈ 0.17. The table gives n̄P(0.14, 0.6) computed with fh's own signal and noise (`np_recon.py`) and FR14's r in brackets.

| bin (z_eval) | QSO | LBG | LAE | Σ n̄P (gal) | Lyα(QSO) P/N | Lyα(LBG) P/N |
|---|---|---|---|---|---|---|
| 1 (2.115) | 0.33 (0.61) | – | – | – | 0.43 | – |
| 2 (2.350) | 0.27 (0.65) | 0.06 (0.89) | 0.22 (0.68) | 0.54 (0.55) | 0.39 | 0.17 |
| 3 (2.586) | 0.19 (0.69) | 0.29 (0.63) | 0.27 (0.65) | 0.75 (0.54) | 0.27 | 0.41 |
| 4 (2.821) | 0.14 (0.74) | 0.94 (0.53) | 0.27 (0.65) | 1.34 (0.52) | 0.20 | 0.47 |
| 5 (3.056) | 0.09 (0.80) | 0.84 (0.53) | 0.22 (0.67) | 1.15 (0.52) | 0.12 | 0.20 |
| 6 (3.291) | 0.06 (0.89) | 0.13 (0.75) | 0.28 (0.64) | 0.47 (0.57) | 0.05 | 0.25 |

  - These values bracket FR14's DESI-QSO value of 0.09 at z ≈ 2.5 (FR14 §4.3).
  - Per tracer, FR14's rule is more aggressive than R=2 for LBG in bins 4–5, and less aggressive for QSO at z > 3 and for LBG in bin 2.
  - A joint reconstruction using all galaxy tracers in the volume (Σ n̄P ≈ 0.5–1.3) gives r ≈ 0.52–0.57 in every bin.
  - Joint 15×2pt σ relative to fh (`damp_rules.py`):
    - FR14 r per tracer: 0.955–0.998 (∥), 0.977–1.000 (⊥).
    - FR14 r for the combined field: 0.935–0.970 (∥).
    - No reconstruction: 1.074–1.127 (∥), 1.045–1.070 (⊥).
  - Verdict: **questionable** as a physical model, but it is a stated input. R=2 is an FR14-compatible intermediate value, not an n̄P-dependent prescription. I recommend adopting FR14's r(n̄P) for a physically motivated baseline, applied to the combined galaxy field if the reconstruction is joint. Also note that all FR14 reconstruction numbers derive from low-z calibrations (Padmanabhan 2012; White 2010), and none of them is validated at z > 2 with Σ ≈ 3 h⁻¹Mpc.
- **No reconstruction for forests**: **justified**. FR14 (Lyα paragraph) states "damping factors, with no reconstruction". No forest reconstruction method exists.
- **Cross-pair rule Σ²_ij = (Σ²_i+Σ²_j)/2.**
  - In the propagator (IR-resummation / Zel'dovich) picture, each field's wiggles are multiplied by G_i = exp(−k²σ²_i/2), where σ²_i is the single-point displacement variance and Σ²_i = 2σ²_i. The cross wiggles are then damped by G_iG_j = √(D_iD_j), which is exactly the fh rule.
  - For a reconstructed galaxy field, the residual displacement is Ψ − S. The forest–galaxy pair displacement variance is ⟨(Ψ₁−Ψ₂+S₂)²⟩ ≈ Σ²_u − ⟨S²⟩, while the reconstructed auto has ≈ Σ²_u − 2⟨S²⟩ (dropping ⟨Ψ₁S₂⟩ at 100 h⁻¹Mpc separation). Hence Σ²_cross ≈ (Σ²_u + Σ²_rec)/2, again the fh rule. This is my derivation, not in FR14 or SE07.
  - lf's "any Lyα pair unreconstructed" ignores the fact that the galaxy leg has been displaced back, so it is too pessimistic.
  - Effect on the joint fit, lf rule vs fh: σ∥ +2.5 to +4.1%, σ⊥ +1.0 to +2.3%. For the Lyα(QSO)+QSO block: +2.4 to +3.8% / +1.2 to +2.2%. This is consistent with the RESEARCH_BASELINE/S4 attribution.
  - Caveat: the rule assumes the reconstruction displacement field is available over the forest volume and that the cross-correlation is measured post-reconstruction. DESI Lyα×QSO BAO analyses do not reconstruct the QSOs.
  - Verdict: **justified physically**. It is the realistic baseline only if cross-correlations with reconstructed galaxy catalogues are planned.
- **Damped wiggles in the fiducial covariance**: **justified**. It is the consistent choice for the nonlinear observable. FR14 uses the undamped linear P in the covariance (text after eq. 16). The effect is <0.3% (inventory row 22).

## 3. Resolution

ME07 §II A (verified in source) states "resolution σ_R (FWHM = 2.355 σ_R = λ/R)". Under the standard spectroscopic definition R ≡ λ/Δλ_FWHM, which is the DESI instrument convention, the Gaussian LSF has σ_v = c/(2√(2ln2)R) = 50.9 km/s at R = 2500 (`response.py:51-54`, `survey_config.py:1379-1381`). lf's σ = c/R = 119.9 km/s is 2.355× too broad. Verdict: fh is **justified**, and lf's convention is a bug.
- FR14 does not specify a resolution convention (A3 digest).
- The residual approximations are a Gaussian, wavelength-independent LSF, and a single R = 2500 for both QSO and LBG spectra.
- The effect on BAO is ≈1% on σ∥ (inventory row 9). This is a planned toggle.

## 4. Growth, evaluation redshift, volume, template range

- **Growth.** fh uses G = [σ8(z_eval)/σ8(2.406)]² applied to both components. This is exact in linear theory for scale-independent growth, with a ≲0.1% neutrino scale dependence for Σmν = 0.06 eV. lf uses EdS ((1+2.3)/(1+z))² on P(2.3).
  - My ratio P_lf/P_fh = 1.003, 0.998, 0.994, 0.991, 0.988, 0.986 for bins 1–6. This includes lf's arithmetic z_c and agrees with inventory row 6 to 1e-3.
  - Verdict: fh **justified** (lf bug fix). The effect is ≤0.7% in σ.
- **z_eval.**
  - fh evaluates signal, noise, weights and damping at the geometric mean in (1+z).
  - lf evaluates the model at the arithmetic z_c and the covariance at the λ-geometric z. The mismatch is 0.0016–0.0022, giving σ8² differences of 0.07–0.14%.
  - Verdict: **justified**, removing an internal inconsistency. The effect is negligible.
- **Single-z vs integrated.**
  - fh integrates the volume (GL in z) but evaluates the power at a single z_eval. The volume ratio to lf is 1 + 2e-4.
  - The bin is 0.235 wide, so the first-order evolution cancels at the geometric centre and the second-order error is O((Δz/(1+z))²) ≈ 0.5%.
  - ME07 §III found that splitting into two z bins gives "slightly smaller" errors. A Fisher density evaluated at z_eval is therefore adequate at the ≲1% level.
  - Verdict: **justified approximation**.
- **Template.** The template is a Vega FITS file at z = 2.406, with domain 1e-4–1147 h⁻¹Mpc. The AP-displaced k' = k(1 ± 2.5e-4) is far inside the domain, and out-of-domain queries raise rather than extrapolate (D1 item 14). P_template matches CAMB at 2.406 to 0.1% (inventory row 5). Verdict: **justified**.

## 5. Covariance

- **Normalisation spot-check** (analytic).
  - The number of independent modes of a real field with μ ∈ [0,1] is N_ind = Vk²ΔkΔμ/(4π²).
  - fh uses N = Vk²ΔkΔμ/(2π²) = 2N_ind, and Cov = (T_imT_jn + T_inT_jm)/N. For an auto this reduces to P²/N_ind, the correct exponential-variance result.
  - It is identical to FR14 eq. 20, ⟨ΔP_ijΔP_mn⟩ = 2π²/(Vk²ΔkΔμ)(P_imP_jn + P_inP_jm) for 0 < μ < 1. It is also consistent with SE07 eq. (1): d³k/(2(2π)³) over μ ∈ [−1,1].
  - D1 confirmed this with an FFT box count and a Monte Carlo test.
  - Verdict: **correct**.
- **Multi-tracer Wick covariance, joint 15×2pt**: **justified**. FR14 eq. 20 uses this form for galaxies, and FR14 §4.3 states that a full Lyα×QSO multi-tracer Fisher "is possible" but uses a high-noise approximation instead.
  - In the high-noise limit FR14 combines errors as an inverse sum, σ = (σ_Lyα⁻¹ + σ_QSO⁻¹)⁻¹. For the DESI-2 densities this is 10–30% optimistic relative to the exact joint 3-spectrum Fisher: joint/FR14 = 1.30/1.21 (bin 2), 1.24/1.17 (bin 3), 1.14/1.10 (bin 5). FR14 itself calls the approximation "borderline".
  - The joint fit is 6–21% better than a naive independent sum of the two autos (`combo.py`).
- **Galaxy shot noise only on autos; no forest–galaxy cross noise.**
  - The QSO's own forest pixels lie 24–430 h⁻¹Mpc in front of it, all at r⊥ = 0. The same-object term in ⟨δ_q W δ_F⟩ is therefore not a Poisson delta function. It is the analogue of the MW11 eq. 5 aliasing term: N_qF(k∥) = n̄_2D⁻¹ Π_qF(k∥), where Π_qF = ∫dr∥ ξ_qF(r⊥=0, r∥) cos(k∥r∥) over the forest range.
  - With linear Kaiser in bin 3, |Π_qF| = 0.01–0.13 h⁻¹Mpc for k∥ = 0–0.2 (`crossnoise.py`). This is ≲10% of the forest's own aliasing numerator P1D, which is about 0.5 h⁻¹Mpc in comoving units.
  - After normalising by the QSO shot noise and the forest aliasing, the correlation coefficient of the cross noise is ≲1%.
  - Standard analyses exclude QSO–own-forest pairs because continuum-fitting distortion dominates them. This removes the mean of this term at negligible cost in signal.
  - Verdict: **justified**, at the ≲1% level, for lya(qso)×qso and, by the same argument, lya(lbg)×lbg. It is a shared lf/fh approximation (inventory row 21), and neither code models the continuum distortion (§9).
- **Shared sample and forest cross noise** between lya(qso) and lya(lbg): both codes set it to zero. This is justified because the two forest samples are distinct sightlines.

## 6. Combination across bins and bin-1 selection

- **Independent bins, summed Fisher.** This follows FR14 (text after eq. 11: independent experiments add) and FR14 eq. 20 (bins independent). It is **justified** for the BAO dilations, which are separate per-bin parameters (`fix_except(ap_i, at_i)`).
  - The bin depth at z = 2.6 is about 180 h⁻¹Mpc, so radial modes with k∥ < 2π/Δχ ≈ 0.035 are correlated across bin boundaries.
  - In the continuum-integral Fisher this correlation is second order. Any cross-bin information is neglected, not double counted.
- **Bin-1 selection** (only lya(qso)², lya(qso)×qso, qso²).
  - The LBG and LAE density tables start at z = 2.38. Bin 1 is [2.0, 2.235], and lya(lbg) sources for bin-1 forests sit at z_src ≈ 2.38, at the table edge.
  - lf fills bin 1 with the clamped z = 2.38 densities, together with an LAE bias extrapolated to 1.25.
  - Physically, an LAE at z < 2.2 has observed Lyα at < 3890 Å, at the edge of the DESI blue-arm throughput.
  - Verdict: **justified** as an input-coverage decision, not a density argument. At the bin-2 density, LAE with n̄P = 0.22 would not be negligible.

## 7. Biases and β

The fh and lf values are identical (inventory row 7). The table compares them with FR14 bD = const using CAMB D(z) = σ8(z)/σ8(0) (`bias_cmp.py`).

| z_eval | b_F (fh) | b_QSO fh / FR14 1.2/D | b_LAE fh / FR14-HETDEX 0.89/D | b_LBG fh |
|---|---|---|---|---|
| 2.35 | −0.138 | 3.57 / 3.19 | 1.56 / 2.37 | 3.3 |
| 2.59 | −0.168 | 3.94 / 3.41 | 1.87 / 2.53 | 3.3 |
| 3.06 | −0.240 | 4.70 / 3.85 | 2.49 / 2.85 | 3.3 |
| 3.29 | −0.282 | 5.10 / 4.06 | 2.80 / 3.01 | 3.3 |

- **Forest b_F ∝ (1+z)^2.9 and β_F = 1.45 constant.**
  - FR14 takes b_F and β_F from MEO03 Table I without quoting them (A3 digest). ME07 marginalises β.
  - From memory, and not verified here: eBOSS DR16 and DESI DR1 measure b_F ≈ −0.11 to −0.12 and β_F ≈ 1.6–1.7 at z_eff ≈ 2.33, including HCD modelling. fh's b_F(2.33) = −0.135 is about 15% larger in amplitude, and β_F is about 15% smaller. The line-of-sight amplitude b_F(1+β_F) is similar (0.33 vs 0.31); the transverse forest power is about 30% higher.
  - Verdict: **questionable input** to be checked against the DESI DR2 fits.
- **QSO b(z)** is 10–25% above FR14's DESI bD = 1.2. It is close to the Laurent et al. 2017 fit I recall (3.3–5.7 across the bins, not verified). Verdict: **justified**, and consistent with current measurements.
- **LBG b = 3.3 constant.** This is the value ME07 uses for WFMOS LBGs (§III validation, verified in source). A magnitude-limited sample should show some evolution of b with z. Verdict: acceptable.
- **LAE b(z)**, linear through 1.76 at z = 2.5 and 2.42 at z = 3.0 and extrapolated, is 7–35% below FR14's HETDEX bD = 0.89. The steep slope (1.3 per unit z) makes the extrapolated bins 2 and 6 input-driven. Verdict: **questionable at the extrapolated ends**.
- **Galaxy b + fμ² with CAMB f(z_eval)** is FR14 eq. 15. Verdict: **justified**.

## 8. Mode grids

fh uses 128 k-intervals with 4 GL nodes each and 32 GL nodes in μ. lf uses 500 linear k points with full dk weights and 10 μ midpoints. I evaluated the same derivative and covariance on both grids (`grid_check*.py`). The k grid is irrelevant: 500 vs 2000 points, rectangle vs trapezoid, all agree to 1e-4. The 10-midpoint μ rule is not.

| | σ(lf grid)/σ(fh grid), ∥/⊥ |
|---|---|
| QSO bins 3, 6 | 1.008/1.003 |
| LAE | 1.009–1.010/1.003–1.004 |
| Lyα(QSO) bin 3, 6 | 1.014/1.006, 1.023/1.010 |
| joint 15×2pt | 1.007–1.009/1.003 |
| 40 midpoints | ≤1.0014 |
| 200 midpoints | 1.0000 |

- The midpoint rule underestimates ∫μ⁴(...)dμ, e.g. Σμ⁴Δμ = 0.1983 vs 0.2, and more so for the steep (1+β_Fμ²)⁴ forest weight. lf is therefore **pessimistic** in σ∥ by 0.3–2.3%.
- The inventory's claim (row 23) that "lf σ⊥ is about 2% optimistic from μ midpoints" is **wrong in sign and axis** for the real prepared noise.
- Verdict: fh GL is **justified**.

## 9. Approximations in this part of fishhighz

Each item is listed with its source statement.

1. **Linear theory plus Kaiser.** There is no forest D_NL(k,μ), such as FR14 eq. 22 (McDonald 2003) or Arinyo-i-Prats 2015, and no galaxy FoG or redshift-error damping.
   - For the forest, D_NL affects k ≳ 0.3 along the line of sight, where the wiggles are already damped.
   - For galaxies, SE07 eq. (27) prescribes R(k,μ) = (1+βμ²)² exp(−k²μ²Σ_z²). Adding a Gaussian σ_z (`fog.py`, bins 3 and 5) changes σ as follows:
     - QSO only, 2 h⁻¹Mpc: joint σ∥ +0.2–0.5%.
     - QSO only, 4 h⁻¹Mpc: +0.8–1.6%, and +20% on the QSO auto.
     - All galaxies, 2 h⁻¹Mpc: +2–3.5%.
     - All galaxies, 4 h⁻¹Mpc: +7–13% σ∥ and +2–3% σ⊥.
   - DESI QSO redshift errors plus nonlinear velocities (σ_v ≈ 3–4 h⁻¹Mpc in DESI Lyα×QSO fits, from memory) and LAE/LBG line offsets make this the largest unmodelled effect in this part of fh. It is shared with lf.
2. **Gaussian covariance, diagonal in k.** This follows FR14 eq. 20 and MW11 eq. 9/10. The forest beat-coupling and (4/N)P_FP_F' terms of MW11 eq. 10 are dropped. MW11 states these "cap" the information at high n̄. There is no super-sample covariance.
3. **Fixed damping widths.** Σ is outside the parameter vector and inside the derivative through term (d). This departs from the SE07 §2 prescription and is worth +0.4–1.3% (§1).
4. **Reconstruction as a scalar R=2 on galaxy Σ.** This is not n̄P-dependent and differs from FR14's r(n̄P) (§2).
5. **No broadband, bias, A_w or Σ marginalisation.** The result is 1.3–1.5% optimistic (§1).
6. **Independent z bins; single-z signal evaluation per bin.** Evolution across a bin is ignored, as ME07 §II A also does.
7. **Flat sky and plane-parallel.** There is a single line of sight per mode and a single (a_v, d_deg) per bin. This is exact in the Fourier Fisher limit used by FR14 eq. 19. Wide-angle effects are negligible for 5000 deg² at χ ≈ 4 Gpc/h at BAO scales.
8. **Fiducial covariance fixed.** There is no ∂C/∂θ term, a Gaussian-likelihood approximation following FR14 eq. 11.
9. **Continuum fitting not modelled.** ME07 §II B drops the first 2N_q,los k∥ modes; MW11 omits the first three line-of-sight modes. fh keeps all k∥ ≥ 0.
   - A forest-only k∥ > 0.01 h/Mpc cut raises the Lyα auto σ⊥ by 6.2–6.5%. At k∥ > 0.02 the increase is 14%. The joint fit rises by 1.5%/7.9% at k∥ > 0.01 when the cut is applied to all pairs, which is an upper bound.
   - ME07 found "no noticeable difference", but its survey was Nyquist-limited.
   - This is a shared lf/fh approximation that is relevant for transverse BAO.
10. **k-range [0.01, 0.5] with no Nyquist k⊥ cut** for the forest. ME07 §II B cuts k⊥ at the effective sightline Nyquist frequency as a conservative choice. At k_max = 0.3 (0.4) the joint σ rises by 0.9–1.0 (0.0)% in ∥ and 1.3–1.4 (0.1)% in ⊥.
11. **No cross noise between forest and galaxies; independent sampling.** This is justified at the ≲1% level (§5).
12. **Linear growth G on the template.** The template is linear total matter (δ_tot), with no neutrino-specific δ_cb (inventory row 5).
13. **Galaxy n̄ at z_eval** (local approximation) with the full volume. The forest weights use a single reference mode (D1).
14. **Gaussian LSF and top-hat pixel.** R is identical for QSO and LBG spectra.

## 10. Documentation mismatches

- `RESEARCH_BASELINE.md` says the derivative "includes the AP prefactor, remapped template, remapped Kaiser factors, and remapped fixed-width damping". That matches the code, but the document does not state that terms (b)–(d) are non-BAO information that SE07 §2 and FR14 §4.1.1 exclude. It should quote the ≤0.5% (joint) and ≤2.6% (forest) size, or offer the SE form as an option.
- D1 ("Modelling choices") calls the wiggle-only remap "equivalent to marginalising a free broadband shape". That is not exact: fixing the smooth component is not marginalisation. It is equivalent to 0.1% only when A_w and Σ are also marginalised (§1).
- `RESEARCH_BASELINE.md` gives the reconstruction as `r_i = 2` in Σ/√r_i. Readers familiar with FR14 will read "50% reconstruction" as Σ×0.5, i.e. R = 4 in fh units. The document should state that R=2 corresponds to FR14 r = 0.71.
- `RESEARCH_BASELINE.md` and `kaiser.md` describe the cross-damping rule correctly but give no physical justification. The propagator argument in §2 could be added.
- `kaiser.md` states "Smooth power is never damped" and D is evaluated at the wiggle coordinates. Both match the code, `kernels/kaiser.py:48-78`.
- `geometry.md` says `prepare_astropy_geometry` supports curved FLRW through the transverse distance, but the volume element lacks the 1/√(1+Ω_kH0²D_M²/c²) factor (D1 F3). The baseline is flat, so there is no impact.
- `covariance.md` and `fisher.md` agree with the code and with FR14 eq. 20. The normalisation statement "half the sphere with 1/(2π²)" is correct.
- Inventory B rows 23 and 26 are superseded by §8 and §1. The μ-grid effect has the opposite sign, and the toy fh/lf derivative ratio is not reproduced.

## Requested runs

Each run is a full 15×2pt forecast with per-bin and combined σ(α∥), σ(α⊥) and ρ against the baseline. The expected sizes come from the single-bin estimates above.

1. **Planned toggles, with expectations to check against:**
   - lf cross-damping rule: joint σ∥ +2.5–4%, σ⊥ +1–2%.
   - Reconstruction factor 1: σ∥ +7–13%, σ⊥ +4–7%.
   - AP wiggle-only Q=1: joint σ⊥ −1%; Lyα autos σ⊥ −2 to −3%.
   - Seo–Eisenstein-like derivative: joint within ±0.5%; forest autos σ⊥ up to +2.7%.
   - Resolution c/R: σ∥ +~1%.
   - Arithmetic z and EdS growth: ≤0.7%, opposite signs in bin 1 and bins 2–6.
   - k_max 0.3: σ +1%/+1.4%. k_max 0.4: ≤0.1%.
2. **New, reconstruction physics.**
   - (a) Per-tracer, per-bin FR14 r(n̄P) using the table in §2. Set reconstruction_factor_i = 1/r_i², for example QSO bins 1–6: 2.69, 2.37, 2.10, 1.83, 1.56, 1.26.
   - (b) If only one global factor is supported, run R = 1/0.54² ≈ 3.4, the joint galaxy-field reconstruction, as the optimistic bracket.
   - Expected σ∥ change: −1 to −7%.
3. **New, galaxy redshift errors/FoG.** Apply a Gaussian exp(−k²μ²σ_z²/2) per galaxy leg, in both the signal and the covariance, with σ_z = 2 and 4 h⁻¹Mpc for QSO, LBG and LAE. This needs a galaxy InstrumentResponse with a Gaussian σ in comoving units; a response in h⁻¹Mpc, not km/s, is sufficient. Expected σ∥: +2–13%.
4. **New, forest continuum proxy.** Drop the forest legs at k∥ < 0.01 and < 0.02 h/Mpc. This needs a per-pair node mask, the same machinery missing for `k_max_<category>` (D1 F2). Expected forest σ⊥: +6–14%.
5. **New, Σ-marginalised BAO.** Free Σ∥ and Σ⊥ per field (plus A_w), with the existing `bao_marginalized` biases. Expected σ: +1.3–1.5% uniformly. This run decides whether the fixed-Σ baseline should carry a marginalisation correction. The current `bao_marginalized` mode (b_i, β_Lyα only) will not change σ by more than 0.1% and is therefore not informative for this question.
6. **New, bias inputs.** Run b_F(2.33) = −0.117 with β_F = 1.67 at the same (1+z)^2.9 evolution. Separately, run LAE bD = 0.89 (FR14 HETDEX) and QSO bD = 1.2 (FR14). These bracket the forest and galaxy signal amplitudes.
