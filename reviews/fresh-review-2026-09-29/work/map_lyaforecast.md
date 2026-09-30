# Phase-1 map: lyaforecast NewForecast path (DESI-2 15x2pt) — UNVERIFIED seed

Paths relative to lib/lyaforecast/lyaforecast/. Config: examples/desi2/lya_qso_lbg_lae_15x2pt.ini.

Setup: 5 tracers lya(qso), qso, lbg, lae, lya(lbg); 15 i<=j pairs (forecast_new.py:61-65). Per z bin: Covariance(zmin,zmax) -> compute_eff_density_and_noise -> model cache + observed cache (forecast_new.py:230-248); Fisher per spectrum + all 15 jointly (:255-309).

1. Cosmology
- z bins: 6 linear over [2.0,3.41]; zc = arithmetic mean (survey.py:78-82). Covariance uses z_mean = sqrt(lmin*lmax)/1215.67-1 (covariance.py:134-138) for conversions, volume, observed power. Fisher model cache evaluated at zc (forecast_new.py:235); damping growth indexed at z_centres.
- velocity_from_distance = 100 H(z)/H0/(1+z) (cosmoCAMB.py:95); velocity_from_wavelength = c/(λ_Lyα(1+z)) (:110); distance_from_degrees = D_A(1+z)(π/180)h (:140-141).
- Volume V = area·c·ln(lmax/lmin)·(dr/dθ)²/(dv/dr) at z_mean (covariance.py:210-242).
- P_lin CAMB at z_ref=2.3, recomputed each call, log grid 1000 pts [0.01,0.5], np.interp (cosmoCAMB.py:78-80). z evolution EdS ((1+z_ref)/(1+z))² (power_spectrum.py:75).
- growth_rate = fσ8/σ8; growth_factor_ratios σ8(z_i)/σ8(z_ref) (cosmoCAMB.py:45-57) used only in damping and discrete β.
- No wiggle split in cosmoCAMB; done in Fisher._get_p_pk.

2. Survey/spectrograph/tracers
- maglist = linspace(16.1,26.75,107) (survey.py:26-29).
- dn/dz/dm: counts/(dz dm), RectBivariateSpline kx=ky=2, no positivity clip (tracer.py:165-178, 216-229). Continuous tracers: mag cut then renormalise to target density over z>2.15 (:140-156). Discrete qso renorm z>2.15; lbg, lae over all z; config mag limits not applied to discrete except via file range (:79-80, 199-210).
- dn/dv = (dn/dz/dm)/(c/(1+z)) (weights.py:175-180); Lyα at z_qso = λc/sqrt(1040·1205)-1 (covariance.py:197-200); discrete at z_mean (weights.py:289).
- SNR tables: 12 files r=19.25..24.75, t=4000s, nexp=4, zq 2.0..4.75; gaussian_filter1d σ=10 along λ, RegularGridInterpolator (m,zq,λ) (spectrograph.py:102-132). lya(qso)×lya(lbg) spectrograph uses sqrt(SNR_qso SNR_lbg) (:98). m<19.25 capped; m>24.75 or out of range → noise 1e10 (:254-266).
- σ_N = 1/[SNR_Å·sqrt(Δλ_pix)·sqrt(nexp/4)] (spectrograph.py:269-277); Δλ = 0.8 Å; σ_N²(m) at one (zq, λc) per bin (weights.py:332-335).
- pix_kms = 0.8 Å·c/λc (covariance.py:172); L_q = c ln(1205/1040) (:230-233); N_pix = L_q/Δv (weights.py:270).
- res_kms = c/R ≈ 120 km/s (covariance.py:182-183) used as Gaussian σ in exp(-k²σ²/2) (spectrograph.py:301).

3. Forest weights (weights.py)
- Mode kt=2.4 deg⁻¹, kp=3.5e-4 s/km for bao, z-independent (:80-82). P3D_w = compute_p3d_kms_smooth(z_mean,kt,kp); P1D_w smoothed with kernel² (:105-109).
- Init: pix_var_1d = p1d_w/pix_kms; w = pix_var_1d/(pix_var_1d+noise_var) (:129-132).
- 3 iterations (:114-116): noise_power = σ_N²(m)/(I1(m)·N_pix) (:353-356); comment "weights include aliasing as signal" but signal_power = p3d_w; w = S/(S+noise_power) (:150-154).
- Integrals cumulative over magnitude (:199-251): int_1 = cumsum(dn·w·dm); int_2 = cumsum(dn·w²·dm); int_3 = cumsum(dn·w²·σ²·dm).
- P_w2D = I2/(I1² L_q); P_N_eff = I3·Δv/(I1² L_q) (covariance.py:292-294); only [-1] used later (:355, :359).
- Discrete compute_tracer_weights unused.

4. Signal
- Kaiser ∏ b_t(1+β_t μ²) (analytic_biases.py:116-121); no FoG/D_NL.
- Bias b = b_ref((1+z)/(1+z_ref))^α: Lyα -0.1352, 3.33, 2.9, β=1.45; QSO 3.54, 1.44; LBG const 3.3 (config); LAE linear 1.76@2.5→2.42@3.0 extrapolated (tracer.py:49-61). Discrete β = f(z)/b(z).
- P3D = P_lin(z_ref)·EdS·Kaiser (power_spectrum.py:176-183).
- W(k∥) = sinc(k∥Δv/2)·exp(-k∥²σ_res²/2) (spectrograph.py:298-303); W per Lyα leg (power_spectrum.py:148-152, 220-227); P1D gets W².
- P1D: PD2013 A=0.064,n=-2.55,α=-0.1,B=3.55,β=-0.28,k0=0.009,z0=3; flat below peak (analytic_p1d_PD2013.py:21-35).
- Damping (fisher.py:217-235): Σ⊥=3.26·σ8(z)/σ8(z_ref), Σ∥=(1+f)Σ⊥, recon /√R (R=2). Any pair name containing 'lya' unreconstructed; gal×gal reconstructed (:156-161).

5. Noise/covariance
- Lyα auto: P_tot = P3D_smooth + P_w2D[-1]·P1D_smooth(k∥) + P_N_eff[-1] then ×(dr/dθ)²/(dv/dr) (covariance.py:352-363).
- Discrete auto: P + 1/n, n = cumsum(dn/dv/dm dm)[-1] (:384-390).
- Crosses signal only (:415-418).
- Wick C_AB = (P_im P_jn + P_in P_jm)/N_modes (fisher.py:121-131). N_modes = V k² dk dμ/(2π²) (covariance.py:246-249); k linear 500 points [0.01,0.5]; μ 10 bin centres [0,1] (power_spectrum.py:43-51).

6. Fisher
- Params (α∥, α⊥) only (fisher.py:30, 54-73).
- Wiggles: 8th-order polynomial in normalised ln k fit to ln|P_model| (first 3 points weight 1e8), wiggle = P - sign·exp(poly), applied to full model incl. Kaiser, smoothing, EdS (_get_p_pk :182-198).
- Wiggles × damping (:156-161); backward FD dP/dlnk[1:] = ΔP/(dk/k[1:]), first 0 (:164-165), so includes damping-envelope derivative.
- F += outer([μ²,1-μ²])² Σ_k dPᵀC⁻¹dP (:58-71); dμ/dα neglected.
- Individual spectra use own C incl. autos (forecast_new.py:266-282); total uses 15×15.
- Bins reported independently.

Flags (unverified): prefix I1; aliasing missing in weighting signal; c/R as σ; either/or recon; zc vs z_mean; EdS vs CAMB growth; single z_q/λ per bin; polynomial de-wiggling on smoothed model; forest effectively m≤24.75; spline not clipped; discrete mag cuts ignored.
