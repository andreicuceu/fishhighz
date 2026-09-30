# C1: forest weighting, aliasing and pixel noise in fishhighz vs lyaforecast

Scope: forest weights, aliasing and pixel-noise coefficients only. No BAO forecast was run. The numbers come from the real DESI-2 inputs of the bundled `desi2_accuracy.ini`. These are the production `prepare_survey`/`prepare_bin` states: 16-point Gauss–Legendre magnitude quadrature, `floor_negative` densities, 2592 magnitude nodes. Scripts and raw output are in `/tmp/claude-81394/.../scratchpad/C1/` (`extract.py`, `analyse.py`, `modes.py`, `mcd_test.py`, `extra.py`, `nyq.py`, `*.out`). My `sum_historical` recurrence reproduces the production A and P_pix to all printed digits in all 11 forest autos. Paper equations were checked against the arXiv LaTeX sources: ME07 astro-ph/0607122 (numbering of the arXiv source; Eqs. 13–19 are §II.B), FR14 1308.4164 (Eqs. 23–29), and MW11 1102.1752 (Eqs. 5, 11–13).

## Summary

The final noise model is the ME07/FR14/MW11 expression, and it is implemented correctly: A = I2/(L I1²), P_pix = l_p I3/(L I1²), N_F = A·P1D·W² + P_pix.

The lyaforecast prefix I1(<m) is a bug relative to all three papers. It is confirmed. The prefix recurrence has no non-trivial fixed point. After 3 updates it gives an arbitrary transient in which bright LBGs are suppressed. As a result, N_F(LBG forest) is 1.6–2.1× the fishhighz value. The QSO forest changes by only −6% to +7%. The omission of aliasing from the lyaforecast weighting signal is also confirmed. It contradicts ME07 and the code comment, but it matches a literal reading of FR14 Eq. 29.

The main new finding concerns all FKP-type prescriptions: early (fishhighz baseline), McDonald (ME07) and FR14-literal. They are all w ∝ 1/(c + l_pσ²) with c = P_S L I1. The minimum-variance per-sightline weight has c = P1D. This follows from ME07's own Eq. 13 and is MW11 Eq. 11. At the chosen reference mode, S·L·I1/B = 0.3–21, so the early weights raise N_F above the minimum by 2–22%. The largest excesses are QSO forest bins 1–3 and LBG forest bins 3–4. An SE07-type proxy for the forest-auto σ improves by up to 10%. The inverse-variance option, which already exists, reaches the per-mode optimum to ≤0.4%. D1 F1's estimate of 0.1–1% is refuted. The McDonald "non-convergence" is a genuine property of the ME07 amplitude recurrence (collapse iff Λ<1, defined below), plus a 96-update cap. It does not affect the scale-invariant coefficients.

## Derivations

### Common notation (observed coordinates: angle in deg, line-of-sight velocity in km/s)

- ρ(m) is the source density per deg² per km/s of source velocity per magnitude: ρ = dN/(dz dm deg²)·(1+z_s)/c. On the quadrature, r_i = ρ_i q_i.
- L = c ln(λ_max^rf/λ_min^rf) is the forest length (km/s). n_2D = L∫ρ dm is the number of sightlines per deg² that pierce a given absorber position. ME07's "L_q l² dn_q/dm" is the probability that a cell is probed (text after ME07 Eq. 14).
- N(m) = l_p σ_N²(m) (km/s) is the 1D white pixel-noise power of one sightline. It follows the "P = σ²v" rule of ME07 §II.B and is MW11's P_N,n (Eq. 3). σ_N is the δ_F noise per pixel. The SNR tables are "S/N per Å for mean quasar with mean forest", so no extra ⟨F⟩⁻² is needed (MW11 Eq. 3 uses continuum S/N).
- B(k∥) = P1D(k∥) W²(k∥) (km/s) and S = P3D(k,μ) W² a_v/d_deg² (deg² km/s). W is the pixel-top-hat × Gaussian response.
- I_n = ∫dm ρ wⁿ for n = 1, 2, and I3 = ∫dm ρ σ_N² w² (ME07 Eqs. 14, 15, 18; FR14 Eqs. 26–28).
- Conversion to comoving units: 1 deg² km/s → d_deg²/a_v (Mpc/h)³, with d_deg = h D_M π/180 and a_v = H/(h(1+z)).

### Observed power

The weighted field is δ_w = (w/w̄)(δ_F + δ_N) (ME07 Eq. 10). Averaging over Poisson sightline positions (ME07 Eqs. 11–13; FR14 Eq. 23; MW11 Eq. 5) gives

P_obs(k) = P3D(k) W² + A·P1D(k∥) W² + P_pix,

with A ≡ P_w^2D = I2/(L I1²) in deg² (ME07 Eq. 16, FR14 Eq. 24) and P_pix ≡ P_N^eff = l_p I3/(L I1²) in deg² km/s (ME07 Eq. 17, FR14 Eq. 25).

The mapping to MW11 Eq. 5 uses n̄ = n_2D, w̃ = w/⟨w⟩, n̄⁻¹ \overline{w̃²} = A and n̄⁻¹ P̄_N = P_pix. The three papers agree exactly.

In comoving units, A d_deg² = 1/n_eff,2D, with n_eff = L I1²/I2 the density of "aliasing-effective" sightlines. The 3D signal term carries no weight factor: the mean weight is spatially uniform, because sources of all magnitudes are Poisson-mixed on the sky. Only the self-pair (aliasing) and noise terms depend on w.

### Moments, effective densities and the minimum-variance weight

Weights enter the covariance only through N_eff(w) = A B + P_pix = ∫ρw²(B+N) / [L(∫ρw)²]. By Cauchy–Schwarz,

N_eff ≥ 1/[L ∫ρ/(B+N)], with equality iff w ∝ 1/(B+N).

This is MW11 Eq. 11 (w̃_n ∝ [P_los + P_N,n]⁻¹). It gives N_eff,min = B/n̄_eff with n̄_eff = L∫ρν, ν = B/(B+N) (MW11 Eqs. 12–13).

In fishhighz the covariance is Wick in T = W W P + N. Adding a positive diagonal noise to T is a PSD increase of the Wick covariance, so the Fisher information of every auto and cross spectrum is monotone in N_F at each node. The minimum-N_F weights are therefore the minimum-variance weights for the BAO parameters within this model. This holds up to the weak k∥ dependence of B(k∥), and the latter is tested below ("opt" column). MW11 §4.2 (Eq. 34) makes the same point for crosses.

The effective densities are ME07's n_p^eff = I1 L/l_p (text after Eq. 18), the aliasing density n_eff = L I1²/I2, and MW11's n̄_eff.

### Per-source weights (all prescriptions in one form)

ME07 Eq. 19 and FR14 Eq. 29 give w = (P_S/P_N)/(1 + P_S/P_N), with P_N(m) = σ_N² l_p/(I1 L). Hence

w_i = c/(c + N_i), c ≡ P_S L I1.

In shape this is an inverse-variance weight with "aliasing floor" c. The optimum has c = B. The prescriptions differ only in c:

| prescription | c | source |
|---|---|---|
| minimum variance (`inverse_variance`, seed) | B | MW11 Eq. 11; S→0 limit of ME07 Eq. 19 with the early P_S |
| McDonald (`sum_aliasing`) | S L I1 + B I2/I1 | ME07 §II.B: P_S = "total flux power ... at k=0.07, μ=0.5, including the aliasing term", i.e. S + B·P_w^2D via Eq. 16 |
| FR14-literal (`sum_intrinsic`) | S L I1 | FR14 text after Eq. 29: "signal power at some typical wavenumber" (no aliasing stated) |
| early lyaforecast (fishhighz baseline) | S L I1 + B | in none of the three papers; early lyaforecast py/forecast.py per the decision record |
| lyaforecast NewForecast | S L I1(<m_i) (source-dependent) | none; see below |

Contents and location of P_S:
- ME07 includes aliasing through the moment expression B I2/(L I1²), at comoving k = 0.07 h/Mpc, μ = 0.5 (k∥ = 0.035, k⊥ = 0.061). Whether P1D is resolution-smoothed there is not stated.
- FR14 uses the same mode. Aliasing is not mentioned.
- MW11 has no P_S. Aliasing enters as a per-sightline noise, on the same footing as P_N,n, at each k∥ (Eq. 11). MW11 recommends a scalar version, σ_los² on ≳10 Mpc (§2.4, §6).

Identity for the early form: at w = ν (the MW11 weights), B/(L I1) = N_eff,min exactly. The early P_S = S + B/(L I1) is therefore the total observed power including pixel noise, evaluated at optimal weights, rather than ME07's "flux power including aliasing". It is not a moment-consistent aliasing term. That does not matter for the covariance, which is exact for any w.

## lyaforecast assessment

(a) Prefix I1(m): confirmed as a bug relative to all three papers.
- In `weights.py`, `compute_int_1` computes `np.cumsum` and `_compute_noise_power_m` sets `noise_power = σ²/(I1(<m) L/Δv)`.
- ME07 Eq. 14 integrates to m_max, and "P_N(m) = σ_N²(m) l_p/I1 L_q" uses that single scalar. FR14 Eqs. 26 and 29 do the same. MW11 normalises over all N sightlines. The comment "cumsum so we can plot as a function of magnitude" describes the diagnostic A(m_lim), P_pix(m_lim). That use is valid only with weights recomputed per limit. It is not valid inside the per-source noise.
- Physical consequence: a source of magnitude m sees the density of brighter sources only. Its noise level is overestimated by I1/I1(<m), so bright, high-S/N sightlines are under-weighted, and the weight profile becomes non-monotonic.
- QSO forest, bin 4: I1(<m)/I1 = 0.02, 0.10, 0.30, 0.62 at r = 19, 20, 21, 22. But N ≪ S L I1(<m) there, so only r ≲ 19.5 is suppressed (w_lf/w_max = 0.73 at r = 19).
- LBG forest: the counts are steep (4%, 18%, 65% of sources brighter than 23.0, 23.5, 24.0). I1(<m)/I1 = 0.02–0.25 over r = 22.5–23.5, where N is already comparable to or larger than B (44–270 km/s vs B = 39 km/s). The lyaforecast weight peaks at r ≈ 24, where N/B ≈ 17. This is anti-optimal: in the MW11 weighting, 84% of the information sits at r < 24 and 18% at r < 23.
- This is why the effect is large for the LBG forest and small for the QSO forest. Result: N_F,lf3/N_F,early = 1.6–2.1 (LBG) and 0.94–1.07 (QSO) (table below). This agrees with B_inventory row 17.

(b) Aliasing omitted from the weighting signal: confirmed. `signal_power = self._p3d_w`, despite the comment "weights include aliasing as signal".
- This contradicts ME07 but is the literal FR14 Eq. 29 reading.
- Its effect is not uniformly adverse. Removing the +B from c moves the weights towards the optimum when S L I1 > B (QSO forest, z < 2.9: N_F,fr14/N_F,early = 0.93–0.96). It moves them away when S L I1 < B (QSO forest, bins 5–6: 1.006–1.063).
- It also opens the collapse channel (see the McDonald discussion below; Λ_FR14 = 0.87 in LBG bin 2).

(c) Fixed 3 iterations. The prefix map diverges: the bright nodes collapse and the collapse cascades to fainter magnitudes. In bin 4, N_eff/N_eff,min is:

| updates | 1 | 2 | 3 | 10 | 40 |
|---|---|---|---|---|---|
| LBG forest | 1.28 | 1.64 | 1.98 | 4.7 | 16 |
| QSO forest | 1.005 | 1.008 | 1.012 | 1.03 | 1.23 |

lyaforecast's forest noise is therefore a transient that depends on the iteration count (≈ +25% per update for the LBG forest). It is also mildly grid dependent: 2.05–2.11 × optimum over 162–2592 nodes in bin 3. It is not a converged estimator.

## fishhighz assessment (verdicts)

| item | verdict | basis |
|---|---|---|
| Full-sample scalar I1 in the update (prefix fix) | justified | ME07 Eqs. 14, 19; FR14 Eqs. 26, 29; MW11 Eq. 2 normalisation |
| Final A = I2/(L I1²), P_pix = l_p I3/(L I1²); aliasing × W², P_pix unsmoothed; d_deg²/a_v conversion | justified | ME07 Eqs. 13, 16, 17; FR14 Eqs. 23–25; MW11 Eqs. 5–7; ME07 §III (aliasing suppressed by resolution) |
| Early update S + B/(L I1) | questionable | Not in any paper. Its "inconsistency" with the final A is harmless, because the covariance is exact for any w. But it is an FKP heuristic that is not minimum-variance: c = S L I1 + B vs the optimum B. |
| Early as baseline, on robustness grounds | unjustified as an optimality choice | N_F is 2–22% above N_eff,min (below). The robustness argument concerns the amplitude degree of freedom, not the covariance. `inverse_variance` is paper-backed (MW11 Eq. 11), iteration-free and within ≤0.4% of the per-mode optimum. |
| McDonald option | justified as the ME07 reading; suboptimal | ME07 §II.B text; within 0–1.5% of early |
| Per-source noise l_p v_i/(L I1) | justified (ME07's stated choice) | ME07 text after Eq. 19 ("not unambiguously defined") |
| Seed B/(B + l_p v) | justified; it is the MW11 optimum at the reference k∥ | MW11 Eq. 11 |
| Reference mode (2.4 deg⁻¹, 3.5×10⁻⁴ s/km), fixed in observed coordinates | questionable | Not ME07/FR14's (k = 0.07, μ = 0.5). It maps to k = 0.051–0.052 h/Mpc with μ = 0.69→0.81 (bins 1→6), i.e. k∥ = 0.036→0.041 and k⊥ = 0.037→0.030. S is 1.9–2.6× the ME07-mode value. This choice maximises the FKP regulariser. For the early form, N/N_min at [fh mode, ME07 mode, k = 0.15 & μ = 0.5] is QSO bin 1 [1.22, 1.08, 1.01] and LBG bin 4 [1.20, 1.06, 1.004]. The comoving drift itself is minor. |
| S including damped wiggles | immaterial | <1% in S at k = 0.05 |
| Convergence criterion (amplitude, shape, A, P_pix; rtol 1e-5; doubled count; cap 96) | numerically sound; physically over-strict | A and P_pix are invariant under w → a w. Requiring amplitude convergence is what makes the McDonald/FR14 collapse "unavailable". |
| Weights held fixed in the mean derivatives | justified | The weights define the estimator, fixed at the fiducial model. The additive noise has no α dependence. The Tr[C⁻¹∂C C⁻¹∂C] term is omitted, as in FR14 Eq. 11. |

### Decision record §3 (homogeneous fixed points; McDonald collapse)

The homogeneous algebra is correct: McDonald has a positive fixed point iff s + B > N, while early tends to B/(B+N). For magnitude-dependent populations I derived and verified the following criterion.

Linearise at w → 0. The map is homogeneous of degree 1 in the weight amplitude a, because P_N ∝ 1/I1 has degree −1 while the aliasing term B I2/(L I1²) has degree 0. The eigenvector is u_i = N_min/N_i. Collapse occurs iff

Λ ≡ [S L Σ r u + B Σ r u²/Σ r u]/N_min < 1

(for FR14-literal, Λ has the first term only).

Production inputs (non-negative, floored densities):

| sample | Λ_McD | Λ_FR14 | outcome |
|---|---|---|---|
| LBG forest bin 2 | 1.04 | 0.87 | FR14-literal collapses (amplitude 3×10⁻⁵ at update 96, 10⁻⁶² at update 1000; A and P_pix stable). McDonald is critically slow: amplitude 0.910 at update 96, then converges to 0.906 by ≈300 updates. |
| QSO forests | 34–169 | – | no collapse |
| LBG forests bins 3–6 | 2.3–3.9 | – | no collapse |

With S at the ME07 mode or at k = 0.15, McDonald collapses in LBG bins 2–5.

Conclusions on the collapse:
- It is a real property of the literal ME07 recurrence in sparse, noisy samples. It is not caused by signed densities (these inputs are non-negative) or by the wiggles in S. It is caused by the representative mode only insofar as S enters Λ.
- In the collapse limit the normalised weights tend to 1/N_i, i.e. noise-only inverse variance, and A and P_pix have well-defined limits. With the weights renormalised each update, N_eff/N_eff,min = 1.00–1.20 (McDonald) and 1.00–1.19 (FR14-literal). The "unavailable" LBG bins are therefore an artefact of the amplitude criterion and the cap, not of a missing covariance.
- I could not test the 107-node signed grid of the archived study, which is behind the barrier. Its bin-2 amplitude of 0.8668 at update 96 is consistent with the same slow approach.

## Quantitative checks

N_F in (Mpc/h)³ at nodes nearest to the listed (k, μ); ratios to early in parentheses. P_F is the forest-auto signal including W² at that node. Variant labels: early = production; mcd = McDonald converged (LBG bin 2 at 400 updates, beyond the production cap); ivar = B/(B+N) with B at the reference mode; opt = w ∝ 1/(B(q) + N) at each node's own q (per-mode lower bound); lf3 = lyaforecast prefix, S only, 3 updates; fr14 = sum_intrinsic converged.

The "σ proxy" is √(I_early/I) with I = Σ modes·e^{−(6k)²}·[P/(P+N)]² over k > 0.02. It is a weights-only effective-volume proxy for the forest auto, not a forecast.

| bin | forest | quantity | early | mcd | ivar | opt | lf3 | fr14 |
|---|---|---|---|---|---|---|---|---|
| 2 | lya(qso) | A [deg²] | 0.02697 | 0.02724 | 0.03449 | – | 0.02893 | 0.02853 |
| 2 | lya(qso) | P_pix [deg² km/s] | 0.6896 | 0.6667 | 0.3601 | – | 0.5752 | 0.5766 |
| 2 | lya(qso) | N_F(0.10,0.48), P_F=26 | 56.5 | 55.8 (0.987) | 49.4 (0.874) | 49.4 (0.874) | 53.4 (0.945) | 53.1 (0.939) |
| 2 | lya(qso) | N_F(0.20,0.90), P_F=24 | 54.6 | 53.8 (0.986) | 46.9 (0.860) | 46.9 (0.859) | 51.3 (0.941) | 51.0 (0.935) |
| 2 | lya(qso) | N_F(0.30,0.90), P_F=11 | 52.4 | 51.6 (0.985) | 44.1 (0.842) | 44.0 (0.840) | 49.0 (0.935) | 48.7 (0.930) |
| 2 | lya(qso) | σ proxy / early | 1 | 0.992 | 0.924 | 0.924 | 0.967 | 0.964 |
| 2 | lya(lbg) | A / P_pix | 0.01114 / 2.811 | 0.01876 / 2.637 | 0.01435 / 2.686 | – | 0.00702 / 6.233 | 0.0202 / 2.635 |
| 2 | lya(lbg) | N_F(0.10,0.48), P_F=26 | 134.0 | 133.9 (0.999) | 131.7 (0.982) | 131.7 (0.982) | 280.0 (2.089) | 135.2 (1.008) |
| 2 | lya(lbg) | N_F(0.30,0.90), P_F=11 | 132.3 | 131.0 (0.990) | 129.5 (0.979) | 129.5 (0.978) | 278.9 (2.108) | 132.1 (0.998) |
| 2 | lya(lbg) | σ proxy / early | 1 | 0.998 | 0.987 | 0.987 | 1.745 | 1.005 |
| 4 | lya(qso) | A / P_pix | 0.06148 / 1.421 | 0.0624 / 1.344 | 0.06942 / 0.9575 | – | 0.06976 / 0.9872 | 0.06895 / 0.9765 |
| 4 | lya(qso) | N_F(0.10,0.48), P_F=42 | 189.0 | 187.0 (0.989) | 181.4 (0.960) | 181.4 (0.960) | 183.5 (0.971) | 181.4 (0.960) |
| 4 | lya(qso) | N_F(0.30,0.90), P_F=18 | 167.2 | 164.9 (0.986) | 156.8 (0.938) | 156.4 (0.935) | 158.8 (0.950) | 157.0 (0.939) |
| 4 | lya(qso) | σ proxy / early | 1 | 0.992 | 0.971 | 0.971 | 0.978 | 0.971 |
| 4 | lya(lbg) | A / P_pix | 0.00328 / 1.451 | 0.003286 / 1.446 | 0.00549 / 1.099 | – | 0.003518 / 2.465 | 0.003298 / 1.437 |
| 4 | lya(lbg) | N_F(0.10,0.48), P_F=42 | 78.1 | 77.9 (0.997) | 65.0 (0.832) | 65.0 (0.832) | 128.7 (1.649) | 77.5 (0.992) |
| 4 | lya(lbg) | N_F(0.30,0.90), P_F=18 | 76.9 | 76.7 (0.997) | 63.0 (0.819) | 63.0 (0.819) | 127.5 (1.657) | 76.3 (0.992) |
| 4 | lya(lbg) | σ proxy / early | 1 | 0.998 | 0.905 | 0.905 | 1.344 | 0.996 |
| 6 | lya(qso) | A / P_pix | 0.2737 / 2.896 | 0.2754 / 2.765 | 0.2783 / 2.561 | – | 0.3179 / 1.513 | 0.3181 / 1.336 |
| 6 | lya(qso) | N_F(0.10,0.48), P_F=66 | 1106 | 1105 (0.999) | 1104 (0.998) | 1104 (0.998) | 1185 (1.071) | 1176 (1.063) |
| 6 | lya(qso) | N_F(0.30,0.90), P_F=29 | 910 | 907 (0.997) | 904 (0.994) | 903 (0.992) | 957 (1.052) | 948 (1.042) |
| 6 | lya(qso) | σ proxy / early | 1 | 0.999 | 0.998 | 0.998 | 1.059 | 1.052 |
| 6 | lya(lbg) | A / P_pix | 0.01421 / 3.349 | 0.01452 / 3.265 | 0.01862 / 2.763 | – | 0.01418 / 6.006 | 0.01514 / 3.134 |
| 6 | lya(lbg) | N_F(0.10,0.48), P_F=66 | 229.6 | 226.2 (0.985) | 213.4 (0.929) | 213.4 (0.929) | 372.6 (1.623) | 221.3 (0.964) |
| 6 | lya(lbg) | N_F(0.30,0.90), P_F=29 | 219.4 | 215.7 (0.983) | 200.0 (0.912) | 199.6 (0.910) | 362.4 (1.652) | 210.4 (0.959) |
| 6 | lya(lbg) | σ proxy / early | 1 | 0.990 | 0.952 | 0.952 | 1.393 | 0.976 |

Other bins, as ivar/early N_F at (0.1, 0.5) and σ proxy:

| bin | forest | N_F ivar/early | σ proxy |
|---|---|---|---|
| 1 | QSO | 0.823 | 0.895 |
| 3 | QSO | 0.919 | 0.946 |
| 3 | LBG | 0.836 | 0.904 |
| 5 | QSO | 0.987 | 0.989 |
| 5 | LBG | 0.950 | 0.965 |

The size of the effect is set by S L I1/B at the reference mode:

| forest | bins | S L I1/B |
|---|---|---|
| QSO | 1–6 | 3.7, 2.9, 1.9, 1.25, 0.68, 0.30 |
| LBG | 2–6 | 1.8, 21, 19, 2.5, 3.6 |

Further diagnostics:
- Sightline densities: n_2D = 64, 43, 31, 19, 10, 3.9 deg⁻² for the QSO forest in bins 1–6, and 74–436 deg⁻² for the LBG forest. MW11 n̄_eff is 21, 20, 14, 11, 6.6, 3.1 (QSO) and 7–30 (LBG) deg⁻².
- Nyquist k⊥ = π√n_2D/d_deg for the QSO forest is 0.39, 0.30, 0.25, 0.18, 0.13, 0.078 h/Mpc in bins 1–6. The same forest-auto proxy puts 2%, 7% and 21% of its information at k⊥ > k_Nyq in bins 4–6.
- 8–10% of the proxy information comes from k∥ < 0.02 h/Mpc, which continuum fitting would remove.

## Approximations in the fishhighz forest-noise model (keyed to paper statements)

1. One representative mode for the weights. ME07: "shortcut of using one constant for P_S"; FR14 Eq. 29. In the FKP forms the mode choice changes N_F by up to 20%. For inverse-variance weights the mode is immaterial (≤0.4%; "opt" column), because P1D is on its plateau (floor at 5.7×10⁻⁴ s/km). MW11 Eq. 11 and §2.2 note the weak k∥ dependence.
2. Scalar weight per source, independent of k∥. MW11 Eq. 11 has w̃_n(k∥), but the paper argues constancy (Table 1; §2.4).
3. Gaussian modes with diagonal covariance var = 2P_tot². The beat-coupling and weight-window terms of MW11 Eq. 10 are neglected, as are ME07's "loose ends" (averaging the power over weight realisations, O(N^{-1/2})).
4. Independent, Poisson-distributed sightlines with no source clustering. MW11 §2.1: factor 1 + C_q n̄, negligible for n̄ ≲ 10⁻² Mpc⁻².
5. Uniform noise along a sightline, i.e. no intra-sightline noise variation (sky lines, masks). ME07 §II.B: "assuming pixels in the same spectrum have the same noise"; MW11 §2.4 σ_los heuristic.
6. SNR evaluated at one (z_q, λ_obs) per bin. There is no λ dependence across the 1040–1205 Å window or within the bin, and no spread in z_q.
7. Source density evaluated at a single z_source: midpoint rule for ∫dz_q dN/dz_q over a z_q window of width ≈ 0.15(1+z). ME07 treats dn_q/dm as constant; MW11 n̄ = N/A. This is not quantified here. The LBG tables span only z = 2.38–3.29 and are clamped outside.
8. Common forest length L for every source (ME07 L_q), with the forest range 1040–1205 Å. FR14 uses 985–1200 Å for BAO; ME07 uses 1041–1185 Å.
9. White pixel noise (MW11 Eqs. 3–4), not smoothed by the response. This is correct for white noise.
10. Continuum fitting is neglected. ME07 drops the first 2N_q,los k∥ modes; MW11 drops the first three line-of-sight modes and adds n̄⁻¹P_los^cont (§5). There is no k∥ cut in fishhighz, and ≈9% of the proxy forest-auto information lies at k∥ < 0.02 h/Mpc.
11. No Nyquist cut on k⊥. ME07 applies one; FR14 says it "definitely must be correct". MW11 instead lets the beat-coupling term cap the information. This matters for the QSO forest in bins 5–6 (k_Nyq = 0.13 and 0.08 h/Mpc), and it is not captured by the diagonal Gaussian covariance.
12. The aliasing P1D is analytic (PD2013) at z_eval, while densities and SNR are taken at z_source. The HCD/DLA, metal and damping-wing contributions to P_los are absent (MW11 §5, App. B).
13. Weights come from the auto only, and the same w is used in the crosses. This is optimal because the cross variance is also monotone in N_F (MW11 Eq. 34). The forest–forest cross correctly has no aliasing term (disjoint sightlines). The shared-sample terms (lya(qso)–qso, lya(lbg)–lbg) are neglected.
14. Weights are fixed at the fiducial model in the derivatives, and noise derivatives are absent (FR14 Eq. 11 convention).

## Doc mismatches

1. FOREST_WEIGHTING_DECISION.md §3 says "FR14 Eqs. 23–29 give the corresponding moment and weighting expressions". FR14 does not state that P_S includes aliasing ("signal power at some typical wavenumber"). Only ME07 does.
2. Decision §2 and §6 ("plausible treatment of intrinsic sightline fluctuations") do not note two things. First, the S→0 limit of the early form is exactly the MW11 minimum-variance weight. Second, at the chosen mode S L I1/B = 0.3–21, so the baseline operates far from that limit and its N_F is 2–22% above the minimum. The inverse-variance comparison (RESEARCH_BASELINE "D2", optional) is the decisive one, not a side check.
3. Decision §5 Step 2 and the LBG table say McDonald is "unavailable" in bin 2. On the production inputs McDonald converges (amplitude 0.906) at ≈300 updates. The 96-update cap and critical slowing (Λ ≈ 1.04) explain the failure. The bin 1/2 "failures" are in the amplitude, not in A or P_pix.
4. Decision §4 prescribes rtol = 1e-4, while the production INI uses 1e-5. This is disclosed in RESEARCH_BASELINE but not in the decision record.
5. `docs/methods/weights.md` calls early "the recommended research method". This is not supported by the numbers above. The statement "inverse_variance ... does not establish a multi-mode BAO optimum" is overly cautious: within fishhighz's covariance the information is monotone in N_F, and ivar is within ≤0.4% of the per-node optimum at all tested BAO modes.
6. `docs/methods/weights.md` calls I1/I2/I3 "prefix integrals". In the full-sum methods they are cumulative sums of the full-sample weights, and only the totals enter (also D1).
7. The `fishhighz/weights.py` and `kernels/weights.py` docstrings say the "Cumulative formulas ... implement McDonald & Eisenstein (2007)". The prefix (legacy) recurrence is not ME07's, since ME07 Eq. 14 integrates to m_max.
8. RESEARCH_BASELINE.md presents (2.4 deg⁻¹, 3.5×10⁻⁴ s/km) as "the representative mode" and defers mode sensitivity to D1. This is not ME07/FR14's (0.07, 0.5): k⊥ is 0.035 rather than 0.061, and S is about 2× larger. For the FKP-type baseline the mode controls a ~20% effect on N_F.
9. lyaforecast `weights.py`: the comment "weights include aliasing as signal" and the docstring "(McDonald & Eisenstein 2007)" are contradicted by the code. The `kt_w_deg = 2.4  # ~0.035 h/Mpc` value corresponds to k⊥ = k∥, not to μ = 0.5.
10. D1 F1 says "Real S and I1·L put it in the 0.1–1% range". This is refuted: the real-input N_F excess of early over ivar is 0.2–22%.

## Requested runs

All runs use the bundled `desi2_accuracy.ini` and the production numerics (step_scale 0.25, rtol 1e-5), with individual and joint σ(α∥, α⊥) per bin plus the combined result. Only the forest `ForestInput.weight_options` differ. They are injected programmatically after `prepare_survey`, because the INI parser only accepts `early_lyaforecast`.

1. R0: baseline, `early_lyaforecast` (reference).
2. R1: `method="inverse_variance"`, `alias=B` (the production auxiliary B = P1D(3.5×10⁻⁴)W² per bin and forest), `auxiliary_coordinates=None`. This is the decisive comparison.
3. R2: as R1 but with B_star at k∥ = 0.1 h/Mpc (q = 0.1/a_v s/km). This demonstrates the insensitivity to the reference mode.
4. R3: `mcdonald` with `max_updates=2000`. Report the amplitude and update count (LBG bin 2).
5. R4: `legacy`, iterations = 3 (lyaforecast weights on the fishhighz grid; equivalent to the full-compatibility estimator but with accuracy numerics).
6. R5: early with the auxiliary mode at ME07's comoving (0.07, 0.5): per bin, k_t = 0.0606·d_deg and k_p = 0.035/a_v.
7. R6 (optional, requires a code change): a k⊥ < π√n_2D/d_deg cut on all spectra involving lya(qso) in bins 4–6.

Expected outcome, from the proxies: R1 lowers the forest-auto σ by up to 10% (QSO bins 1–3, LBG bins 3–4). The joint 15×2pt change should be smaller and needs these runs.
