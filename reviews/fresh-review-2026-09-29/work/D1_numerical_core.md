# D1 audit: numerical core of fishhighz

Scope: `kernels/*`, `derivatives.py`, `fisher.py`, `covariance.py`, `grids.py`, `geometry.py`, `response.py`, `noise.py`, `weights.py`, `fields.py`, `_information.py`, `models/{templates,kaiser,biases,p1d,external}.py`; compared with `docs/methods/*.md` and `docs/research/RESEARCH_BASELINE.md`. Paths are relative to `lib/fishhighz/`; `fishhighz/` is the package. Scratch checks are in `/tmp/claude-81394/-global-cfs-cdirs-desicollab-users-acuceu-vega-dev/d013f31d-c80a-4f24-8aa6-2bb33b3bfb58/scratchpad/D1/` (`c1.py` to `c7.py`). None ran a survey forecast. No package file was modified.

## Summary

I found no major defect: no wrong factor, unit, ordering or derivative in the numerical core. Every item in the requested checklist reproduces either analytically or against an independent implementation (list below). Findings are one design-level suboptimality of the weight recurrence, which can reach the percent level, and a few latent or documentation issues.

## Findings

### F1. Full-sum recurrence weights do not minimise the aliasing plus pixel noise (design, minor; up to about 1-3 % in S+N in a toy)
- File: `fishhighz/kernels/full_sum_weights.py:31-52` (`update`); the production method is `early_lyaforecast` = `sum_historical`.
- Statement of the estimator. For weights w(m) the coefficients are A = I2/(L I1²) and P_pix = Δv I3/(L I1²). The total observed power at a mode is S + A·B + P_pix. Its dependence on w enters only through N_eff = A·B + P_pix ∝ Σρq w²(B + Δv σ²)/(Σρq w)². This is minimised by w ∝ 1/(B + Δv σ²), independent of S. That is exactly the seed (`seed`, full_sum_weights.py:14-18) and the documented `inverse_variance` option.
- The recurrence uses w ∝ 1/(σ² Δv + B + S·I1·L). The extra term S·I1·L is an FKP-like regulariser. It moves the weights away from the N_eff minimum whenever S·I1·L is not ≪ B.
- Evidence (`c4.py`, `c5.py`, toy: 80 magnitude nodes, B = 40 km/s, Δv = 60 km/s, σ² ∝ 10^(0.8m)). The table gives (S+N)/(S+N_seed):

| I1·L (deg⁻²) | S (deg² km/s) | historical | aliasing (`mcdonald`) |
| --- | --- | --- | --- |
| 30 | 1.2 | 1.0062 | 1.0015 |
| 30 | 0.3 | 1.0005 | 1.0242 |
| 10 | 1.2 | 1.0009 | 1.0197 |
| 100 | 1.2 | 1.0313 | 1.0206 |

  For I1·L ≪ B/S the recurrence collapses onto the seed, and the historical fixed point equals the seed to 1e-7. The seed weights give N_eff identical to the analytic optimum 1/(B+Δvσ²).
- Impact: σ(α) scales roughly as (S+N)/S, so the effect is 0.1-0.6 % at typical toy densities and about 3 % at the highest density tried. The bias is towards a slightly larger forest error, i.e. conservative. The `mcdonald` variant, which is self-consistent in A, is not better (up to 2.4 % worse), so the mismatch between the update alias term B/(I1 L) and the final A is not itself the cause. The docs already list the comparison with `inverse_variance` as the optional D2. I quantify it here only for the toy. Real S and I1·L (bin 1: about 30 forests deg⁻², S ~ 1 deg² km/s) put it in the 0.1-1 % range. Not verified on real inputs.

### F2. `k_max_<category>` cuts are not applied per pair in BAO mode (latent, minor)
- File: `fishhighz/survey_config.py:1257-1277` (outside the strict scope, found while checking the grid).
- The per-category `k_max_*` values only set the outer edge (max of the cuts) of a single shared grid. Every selected pair, and hence the joint covariance, is summed over all nodes up to the largest cut. `full_shape.py` has `_selected_nodes`, but BAO mode has no equivalent.
- The bundled INI sets none of these keys, so all cuts are 0.5 and the shipped baseline is unaffected. If a user sets e.g. `k_max_forest_forest < k_max_galaxy_galaxy`, the forest-forest information beyond its cut is silently included.

### F3. Volume element ignores curvature (latent, minor)
- File: `fishhighz/geometry.py:163-164`.
- V = Ω h³ ∫ c D_M²/H dz is exact for Ω_k = 0 only. For Ω_k ≠ 0 it needs the extra 1/√(1 + Ω_k H0² D_M²/c²). The docstring of `prepare_astropy_geometry` advertises support for curved FLRW ("uses comoving_transverse_distance"), which is not the correct volume there. The Planck18 baseline is flat (omk = 0), so there is no impact on the production numbers.

### F4. Stale docstring (doc, trivial)
- `fishhighz/grids.py:6-7`: "Initial Fisher calculations will hold covariance ... No Fisher calculation is implemented here." This is out of date; the Fisher code is in `fisher.py` and `forecast.py`.

## Modelling choices that change σ(α) but are documented (not defects)
- LSF: σ_v = c/(2√(2 ln 2) R) with R = 2500 is about 51 km/s; lyaforecast used σ = c/R (about 120 km/s). This is a deliberate, documented change in `RESEARCH_BASELINE.md` and `response.py`. It is a large lever on high-k∥ forest information. The value was not changed or audited here.
- The BAO-only model puts Q and the AP mapping on the wiggle component only. This is documented and equivalent to marginalising a free broadband shape. It is not the full-P AP.
- The weighting mode (2.4 deg⁻¹, 3.5e-4 s/km) is a single reference mode. k_t is in deg⁻¹ without a 2π, consistent with lyaforecast (k⊥ = 2.4/68 = 0.035 h/Mpc, matching the stated comment).
- The forest bias b_Lya = -0.1352 with fixed β_F = 1.45 and the galaxy f = f_CAMB(z_eval) are inputs, not numerics.

## Checks passed
1. a_v = H/(h(1+z)) [(km/s)/(Mpc/h)] and d_deg = h D_M π/180 (`geometry.py:165-166`) verified by dimensional analysis. Test values for z_eval = 2.35 are a_v = 105.1 and d_deg = 68.3 (Planck18; `c6.py`).
2. Forest noise conversion: (deg² km/s) → (Mpc/h)³ by d_deg²/a_v (`noise.py:82-84`). Derived independently: 1 deg² = d_deg² (Mpc/h)², 1 km/s = 1/a_v Mpc/h. Weighting-mode inverse a_v/d_deg² is correct (`weights.py`, `sample_auxiliary`).
3. Local galaxy density n = Σqρ (1+z)/c · a_v/d_deg² equals dN/dz/deg² · H/(c h³ D_M² (π/180)²) algebraically (`noise.py:32`). This is the correct number per (h/Mpc)³.
4. Forest density per velocity ρ = dN/dz (1+z_src)/c times L gives the forests per deg² crossing a pixel (Δln(1+z_q) = L/c). Correct.
5. q = kμ/a_v, in observed coordinates, for response and P1D. Response does not depend on AP (correct: instrument acts in observed coordinates).
6. Response: `np.sinc(x/π)` = sin x/x with x = qΔv/2 (Δv the full pixel width). Numerical check (`c1.py`) against sin x/x reproduces to 1e-16 for q = 0, 0.01, 0.02 s/km at Δv = 6 km/s. Gaussian exp(-q²σ²/2) correct.
7. Response applied once to the signal (W_iW_j) and to the Jacobian (`forecast.py:269`), squared once in aliasing, not applied to P_pix (`noise.py:83-84`), galaxies W = 1 (`InstrumentResponse(0,0)`).
8. Mode normalisation q_mode = k² w_k w_μ/(2π²), N = V·q_mode, μ ∈ [0,1]. Derivation: independent modes = ½ V ∫d³k/(2π)³ = V k² dk dμ/(4π²) on μ ∈ [0,1]. N_code = 2 N_ind, and Cov = (T T + T T)/N_code reduces to P²/N_ind for an auto. Brute-force FFT box count (`c1.py`, ratio 0.991 with 1/2 N_code) confirms N_code/2. Monte Carlo of the Wick formula for three correlated fields (`c2.py`) matches all pair combinations to 1-6 % (statistical). Fisher for ln A equals the Seo-Eisenstein integral V/(4π²)∫k²dk dμ to 2e-15.
9. Wick ordering: `fields.py` tables im, jn, in_, jm from (i,m),(j,n),(i,n),(j,m); kernel fills the upper triangle and mirrors it (`kernels/covariance.py`). Pair canonical order i ≤ j.
10. Gauss-Legendre grid: k-nodes edges + half-width (x+1) with weights half-width·w; μ = (x+1)/2 with weights w/2 (`grids.py:109-135`). ∑w_k = 0.49, ∑w_μ = 1, ∫μ² = 1/3 to 1e-16, ∑q_mode equals the analytic (0.5³ - 0.01³)/(6π²) exactly.
11. AP mapping k∥' = kμ/α∥, k⊥' = k√(1-μ²)/α⊥, Q = 1/(α∥α⊥²) (`kernels/kaiser.py:28-38`). This is the standard convention (P_obs(k) = P_fid(k/α)/(α∥α⊥²); ∫d³k P_obs is conserved). Remapped Kaiser factors use μ', damping uses the wiggle component's (k∥', k⊥') (the last loop iteration in `models/kaiser.py`). `KaiserModel` output agrees with an independent reimplementation (`c3.py`, max relative difference 0.0). The `_scales` bases satisfy α = √(α∥α⊥) and φ = α⊥/α∥ (also α_iso^3 = α∥α⊥²).
12. Cross damping: Σ_ij² as the arithmetic mean per axis makes D_ij = √(D_i D_j) exactly, as the docs claim.
13. Template h-rescaling k_fid = k h_t/h_fid, P_fid = P (h_fid/h_t)³ is the correct direction (P in (Mpc/h)³ scales as h⁻³ physical; k³P invariant). Numerical check with h_t = 0.6736, h_fid = 0.7: P_fid h_fid⁻³ = P_t h_t⁻³ to 1e-9 (`c3.py`).
14. Template spline domain: real template covers 1e-4 to 1147 h_fid/Mpc; the grid needs 0.01 to 0.5 with mapped k' = k/α (α = 1 ± 2.5e-4 in the production step). No node leaves the domain; out-of-domain would raise, not extrapolate.
15. Growth G = [σ8(z_eval)/σ8(z_ref)]² multiplies both components once (`_assemble`; `template.evaluate(..., growth=)` is only used for validation). Correct for σ8 ∝ D.
16. Central FD stencil: weights (-½, +½) with actual offsets, derivative = (f₊ - f₋)/(2h) (`derivatives.py:89`, `kernels/derivatives.py`). Forward stencil (-3f + 4f₁ - f₂)/(2h) correct. Step 1e-3 × step_scale 0.25 = 2.5e-4 (`survey_config.py`). Jacobian converged: σ(α∥) = 0.0205545, 0.0205538, 0.0205538 for step_scale 1, 0.25, 0.1.
17. End-to-end galaxy-only chain (real Planck18 template, astropy geometry, real grid 129 k-edges × 4 × 32 μ, QSO-like Kaiser, n̄ shot noise, central FD): σ(α∥), σ(α⊥) = (0.020554, 0.012988) agree to 1e-6 with an independent Seo-Eisenstein integral with its own derivatives (`c6.py`).
18. Fisher contraction: Cholesky forward-substitution `Σ(L⁻¹J)ᵀ(L⁻¹J)` equals Σ Jᵀ C⁻¹ J to 1.6e-15; the numba backend is bit-identical (`c7.py`). Normalisation D·L_R in `factor_covariance` is correct.
19. Sum over nodes includes quadrature weights only through N_modes (no double-counting); bins are independent, `fix_except(ap_i, at_i)` is the full per-bin block in BAO mode.
20. Weight kernels: I1, I2, I3, A = I2/(L I1²), P_pix = Δv I3/(L I1²) implemented as documented (`full_sum_weights.py:55-66`, `kernels/weights.py`). Direct sums, `_integrals` and `coefficients` agree to 1e-15 (`c4.py`).
21. Recurrence implemented = documented: `sum_historical` w = S'/(S' + σ²Δv/(I1 L)), S' = S + B/(I1 L); `sum_aliasing` S' = S + B I2/(I1² L); seed w₀ = (B/Δv)/(B/Δv + σ²); legacy cumulative form matches `weights.md`.
22. Convergence logic: `solve` returns the actually computed iterate at t = 2n (the doubled-count confirmation), never a substituted candidate; capped and failed states are labelled and rejected by `prepare_forest_weights`. Fixed-point residual at the returned state is 0 to 2e-16 in the toy; historical converges in 8 updates, aliasing in 20 (rtol 1e-6).
23. P1D (PD2013): identical to `lyaforecast/analytic_p1d_PD2013.py` to 7e-16 for z = 2.2, 3.0, 3.4; P(k0 = 0.009 s/km, z = 3) = πA/k0 = 22.34 km/s with A = 0.064, n = -2.55, α = -0.1, B = 3.55, β = -0.28, z0 = 3 (the published best-fit constants as I recall them; not re-verified against the paper here). The low-k floor sits at the exact maximum of the log-parabola (k_floor = k0 e^{5(2+n_z)}), so the plateau is continuous with zero slope.
24. Wick covariance and Fisher are ordered consistently (`selected_to_required` gather for J and factor blocks).

## Doc / code mismatches
- `grids.py:6-7` docstring (F4).
- `RESEARCH_BASELINE.md` states the default `rtol = 1e-4` and the S2/S3 state `1e-5`; the INI sets `weighting_rtol` (input policies) and `weight_rtol` (numerical) to 1e-5 and only the numerical one is read at `survey_config.py` (about line 1370), while the parser mirrors both. Consistent today; two keys for one quantity.
- `weights.md` describes `I1/I2/I3` as "prefixes"; in the production full-sum path these are cumulative sums of the final full-sample weights, not the prefix-consistent (limiting-magnitude) weights of the legacy recipe. Only the final entries enter A and P_pix.
- `kaiser.md` claim that cross damping is the geometric mean of the auto damping factors is correct (checked, item 12), but the phrase "mean of squared widths" and "geometric mean" describe the same object; no mismatch in code.
