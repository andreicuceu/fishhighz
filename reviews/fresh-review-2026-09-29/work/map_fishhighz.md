# Phase-1 map: fishhighz production path (Forecast("desi2_accuracy.ini").run()) — UNVERIFIED seed

Paths relative to lib/fishhighz/.

0. Entry
- Bundled INI fishhighz/data/desi2_accuracy.ini; defaults from fishhighz/accuracy.py:395-446 (INI_DEFAULTS). parse_survey_ini survey_config.py:227.
- Forecast.prepare (public.py:392-462): prepare_camb (public.py:253), load_template (:270), _build_readers (survey_config.py:1066), prepare_survey (:1126), prepare_bin per bin (forecast.py:71).
- Forecast.run (public.py:464-612): step_scale = 0.25 (:493) → run_forecast (forecast.py:318) → run_bin (:241). Per-bin reduction fix_except(ap_i, at_i) (:506) = principal submatrix. Individual spectra re-prepared with own covariance (_pair_spec :278, loop :525-551) — 78 extra preparations each re-solving weights.
- CLI: cli.py:491 --joint-only → run(individuals=False) raises in bao mode (public.py:489-490).

1. Cosmology/geometry
- z_eval = sqrt((1+zmin)(1+zmax))-1 (survey_config.py:928).
- CAMB (cosmology.py:294-435): σ8(z), f = fσ8/σ8 at z_eval, 2.406 (template), 2.3 (damping ref); h_fid; D_M = D_A(1+z).
- V = Ω h³ Σ w_z c D_M²/H (32 GL nodes) (geometry.py:115-191); a_v = H/(1+z)/h; d_deg = h D_M π/180.
- Template (models/templates.py:140-230, 278-337): k_fid = k h_t/h_fid, P_fid = P (h_fid/h_t)³; smooth = PKSB, wiggle = PK-PKSB; not-a-knot cubic in ln k; no extrapolation.
- G = (σ8(z_eval)/σ8(2.406))² on both components (survey_config.py:1285-1286, 1356; kernels/kaiser.py:68).

2. Survey
- Grid: k edges linspace(0.01,0.5,129), GL 4/interval, μ 32 on [0,1] (survey_config.py:1274-1283; grids.py:109). q_mode = k² w_k w_μ/(2π²) (grids.py:87).
- λ_obs = 1215.67(1+z_eval) (:1372); forest z_src = λ_obs/sqrt(1040·1205)-1 (:1388-1397); galaxies z_src = z_eval.
- Response: Δv = c·0.8/λ_obs; σ_v = c/(2500·2√(2ln2)); galaxies no response (:1373-1387; response.py:30-54).
- L_v = c ln(1205/1040) (:1426).
- Density (adapters/legacy_inputs.py:164-197; adapters/legacy_compat.py:61-148): n = raw·mask·target/Σ_{z>z_norm_min} raw_masked/(dz_first dm), RectBivariateSpline kx=ky=2; negatives and out-of-range floored to 1e-20 (legacy_compat.py:122). Galaxies density_magnitude_bounds=none (normalised over file mags, integrated over [16.1,26.75]).
- SNR (legacy_inputs.py:345-351; legacy_compat.py:173-256): gaussian σ=10 samples; σ² = 1/max(SNR√Δλ√(Nexp/Nexp_file),1e-10)²; bright clamp; out-of-range σ²=1e20.
- Magnitude quadrature (magnitude.py:11-58): breakpoints from density knots, spline roots, SNR mags; GL 16/segment; shared partition.
- Galaxy n̄ = Σ q dN/dzdm (1+z_eval)/c · a_v/d_deg² (noise.py:19-33; weights.py:35-51).

3. Forest weights
- weighting_method early_lyaforecast → sum_historical (kernels/full_sum_weights.py:11), called weights.py:301-325.
- Mode (weights.py:88-159) kt=2.4 deg⁻¹, kp=3.5e-4 s/km; k∥ = a_v kp, k⊥ = kt/d_deg; S = P3D_auto(k,μ) W² a_v/d_deg² (full Kaiser incl. BAO and damping); B = P1D(kp, z_eval) W².
- Seed w0 = (B/Δv)/(B/Δv+σ²) (full_sum_weights.py:14-18). Update (:31-52): w = S'/(S' + σ²Δv/(I1 L)), S' = S + B/(I1 L), I1 = Σρwq.
- Moments (kernels/weights.py:236-250): I1, I2, I3; A = I2/(I1² L); P_pixel = I3 Δv/(I1² L).
- Convergence (full_sum_weights.py:111-213): rtol 1e-5 (INI), min 3, stable 3, max 96, doubled-count confirmation; raises otherwise (weights.py:320).

4. Signal
- Biases (survey_config.py:1288-1299; models/biases.py:84-100): Lyα -0.1352, z_ref 2.33, α 2.9; QSO 3.54, 2.33, 1.44; LBG 3.3 const; LAE linear 1.76@2.5, 2.42@3.0 extrapolated (biases.py:56-73). β_F = 1.45; galaxies f_CAMB(z_eval).
- Kaiser (kernels/kaiser.py:337-342): forest b(1+βμ'²); galaxy b+fμ'².
- AP (kernels/kaiser.py:310-334; models/kaiser.py:249-297): wiggle only, k∥' = kμ/a∥, k⊥' = k√(1-μ²)/a⊥, Q = 1/(a∥ a⊥²); smooth unscaled.
- P_ij = G[F_i F_j P_sm + Q F'_i F'_j D_ij P_w(k')] (kaiser.py:355-375).
- Damping (survey_config.py:1305-1314; kaiser.py:345-352): Σ⊥ = 3.26 σ8(z_eval)/σ8(2.3)/√R, R=2 galaxies, 1 forests; Σ∥ = (1+f)Σ⊥; Σ²_ij = mean; D on wiggle coords.
- P1D (models/p1d.py:33-56) PD2013 at z_eval.
- Response W = sinc(qΔv/2π)·exp(-q²σ_v²/2), q = kμ/a_v; signal × W_iW_j (response.py:143-163).

5. Noise/covariance
- N_F = (A P1D(q) W² + P_pixel) d_deg²/a_v (noise.py:52-95); galaxies 1/n̄; diagonal only (noise.py:130-136).
- T_ij = W_iW_j P_ij + N_ij (forecast.py:165); C_AB = (T_im T_jn + T_in T_jm)/N_modes, N_modes = V q_mode (kernels/covariance.py:6-24).

6. Derivatives/Fisher
- Params ap_i, at_i per bin; step 2.5e-4; everything else fixed (survey_config.py:1206-1212).
- Central FD (derivatives.py:59-104).
- F = Σ (L⁻¹J)ᵀ(L⁻¹J), J = W_iW_j ∂P/∂θ (forecast.py:269-272; fisher.py:127-203).
- Bins summed (results.py:261). Bin 1 selection in INI: lya(qso)², lya(qso)×qso, qso².
- full_shape.py, bao_marginalized.py only for other [model] modes.

Production vs validation: production = top-level modules, models/*, kernels/*, accuracy.py, adapters/legacy_inputs.py, adapters/legacy_compat.py (lazy-imported by survey_config.py:1069-1070). Validation-only: fishhighz/validation/*, adapters/desi2_compatibility.py, adapters/lyaforecast.py. No production import of validation.

Tests by stage: 1 test_cosmology/geometry/templates; 2 test_survey_config/legacy_inputs/survey_composition/resources; 3 test_full_sum_weights/adaptive_weights/weights/weight_convergence/weight_limit/weight_range/weighting_w12; 4 test_kaiser/biases/p1d/response/models; 5 test_noise/covariance/grids/fields; 6 test_derivatives/fisher/results/forecast/public_forecast/desi2_accuracy_example. No production test pins real bundled-run σ values.

Flags (unverified): --joint-only crash; update aliasing B/(I1 L) vs final A = I2/(I1²L); Q on wiggles only; no shared-sample noise; galaxy normalisation ignores mag bounds; P1D/weights at z_eval vs density/SNR at z_src.
