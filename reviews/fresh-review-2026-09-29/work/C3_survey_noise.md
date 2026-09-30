# C3: survey inputs and noise model of fishhighz vs lyaforecast (scientific adequacy)

Scope: density tables, forest S/N model, forest length/continuum, P1D, shared-sample noise, bin-1 selection.
Method: code and paper-digest reading (A1 ME07, A2 MW11, A3 FR14, B_inventory, D1, D2) plus preparation-level
evaluations on the bundled inputs (no Forecast.run, no Fisher matrices). Scripts: scratchpad/C3/{dens,a1,a3,a3b,gal,nP,snr,snrlam,lda,lda2,cont}.py.
Density/weights numbers use the production readers (`prepare_survey`, `prepare_bin`, `prepare_forest_weights`,
`forest_noise`) with modified inputs; "N(0.14,0.6)" below is the total forest noise power at k=0.14 h/Mpc, mu=0.6.
D2 was read; bug-level items (M1 clamp, m2 first spacing, m3 floor) are taken as given and only their scientific size is added.

## Summary

The survey/noise inputs reproduce lyaforecast; the departures are choices, not errors. Findings that matter for sigma(alpha):
(i) the LBG/LAE tables have 5 z rows (cells 2.27-3.40) and vanish for r>24.25 (LBG, LAE) and r>23.5 (QSO). The r_max 25.75/26.25 variations are therefore null tests, and the r=24.75 SNR cutoff removes nothing.
(ii) The single-z_src local approximation for forests is biased by 4-24% in QSO-forest noise at high z (density falls steeply across 1040-1205 A), by up to +30% for lya(lbg) in bin 3, and is unsupported for lya(lbg) in bins 5-6, where the table ends at z=3.40 but z_src=3.40/3.66. Truncating at the table edge raises lya(lbg) noise by 1.4x (bin 5) and 6.8x (bin 6), or removes the forest (edge at 3.29).
(iii) Galaxy n-bar is shot-noise dominated (nP(0.14,0.6)=0.06-0.9), so first-spacing (+4.8% LBG, +2.6% LAE) and floor (+0.5%) normalisations act nearly linearly on sigma.
(iv) The QSO normalisation depends on a threshold sitting on a cell centre (z_norm_min=2.15): the count changes by 16% depending on a float comparison.
(v) Continuum-fit mode loss is neglected (toy: sigma_perp +8-20% for k_par<0.015-0.03).
(vi) Noise independence is right for lya(qso) x lya(lbg), defensible for lya x tracer if same-object pairs are excluded.
Pixel-width variations are expected to be null (sigma^2 Delta v is independent of Delta lambda).

## 1. Density normalisation

Bundled tables (checked): QSO z=0.05..4.75 (dz=0.1, 48 rows), r=16.1-24.9 (dm=0.2), counts nonzero only for r in [16.5,23.5]. LBG/LAE: z=2.38, 2.60, 2.83, 3.07, 3.29 (spacings 0.22, 0.23, 0.24, 0.22), r=21.75-26.75 (dm=0.1), counts nonzero only for r in [21.75,24.25]. Cutoffs are sharp (e.g. LAE row 2.38: 8.76 at r=24.25, exactly 0 at 24.35).

| Item | Verdict | Assessment |
|---|---|---|
| target_density, z_norm_min=2.15 (QSO) | questionable | Normalisation sums rows with z>2.15 strictly. The QSO row centre 2.15 is stored as 2.149999999999999467, so the row (cell 2.1-2.2, 16.9 of 106.9 normalised deg^-2) is excluded only by float representation. The effective target is "90 deg^-2 at z>2.2, extending to z=4.8". If "90" meant z>=2.1 the densities would be x0.842 (row count 710 of 4500 raw). The bin-integrated QSO counts in 2.0-3.41 are 40.8, 31.3, 21.2, 14.2, 9.0, 5.2 deg^-2 (sum 122). The physical definition of the target (redshift range, purity, redshift success) is undocumented in the repository (data README only states the file origin). FR14 applies an explicit 0.8 targeting-efficiency factor; no analogue here. |
| target for LBG/LAE | justified as a convention, questionable as a physical statement | Normalised over all rows (no z_norm_min); 380/430 = sum of counts. With true cell widths the integral is 398 (LBG) and 441 (LAE) deg^-2, see next row. |
| dz first spacing | LBG/LAE questionable; QSO justified | QSO uniform (exact). LBG/LAE rows divided by 0.22 although spacings are 0.22-0.24: total 398.2 (+4.8%) and 441.3 (+2.6%); rows 2.60/2.83/3.07 high by 2.3/6.8/4.5% vs mid-point widths (D2 m2; reproduced). All galaxies here have nP<1 (below), so sigma for LBG/LAE autos shifts by roughly the same percentage; lya(lbg) A scales as 1/rho. |
| density_magnitude_bounds survey vs none | justified, inert | Both table magnitude ranges lie inside [16.1,26.75]; both settings give identical densities. |
| Galaxy n-bar magnitude range | justified | Integration 16.1-26.75 covers all nonzero counts (QSO r<=23.5, LBG/LAE r<=24.25). Rectangle-rule vs quadrature differ by <0.6% (B row 13). |
| 2D quadratic spline (kx=ky=2, s=0) | questionable | Magnitude direction: ringing at the sharp cutoffs gives negative area of -0.44 to -0.60% of n (QSO, LBG, LAE, all z; my numbers -0.44..-0.52% QSO, -0.50..-0.60% LBG, -0.50..-0.60% LAE). Redshift direction for the 5-row LBG table: the interpolant peaks at 769 deg^-2 per unit z at z=2.93 while the nodal rows are 700 and 665; the lya(lbg) source density at z_src=2.894 is 759.7 (+8.5% vs the 2.83 row). Integral of the spline over 2.38-3.29 is 400.9 vs 382.3 (trapezoid of rows, +4.9%). LAE is benign (max 474 vs 468). Three-point or coarser tables are not constrained by a quadratic tensor spline; a monotone (PCHIP) or a cell-average-conserving interpolation is the scientifically defensible alternative. |
| Floor 1e-20 on negative interpolant | justified, small | Floor selects the positive part: n-bar +0.44-0.60% vs signed (lyaforecast keeps the negative part); I1 +0.04-0.8% (D2, largest for lya(lbg)). sigma effect below 0.3%. The floor is biased high relative to a conservative interpolant only by the same 0.5%. |
| Clamp outside z table | questionable | See items 2 and 7. LBG/LAE bins 1-2 (below 2.38) and lya(lbg) bins 5-6 (above 3.29) sit on clamped rows. |

## 2. Local-density approximation

Definition in code: galaxy n-bar = sum q rho(z_eval) a_v/d_deg^2 (one z); forests rho = dN/dz(z_src)(1+z_src)/c with z_src=lambda_obs(z_eval)/sqrt(1040*1205)-1; SNR at (z_src, lambda_obs(z_eval)). Verdict: questionable for forests (biases of 5-30% in noise), justified for QSO/LAE galaxies, questionable for LBG galaxies in bin 6 and bin 3.

Literature: ME07 takes all quantities at the central redshift of a wide bin (2.2<z<3.3) and states that a two-bin test lowers the combined error slightly, so single-z is the ME07 convention, tested by them at the ~percent level for a smooth luminosity function. FR14 digest gives no bin-integration prescription for the Lya block (LF tabulated per dz per deg^2). MW11 uses n-bar_eff at one representative z. None of them faces a 5-row histogram-like density.

Geometry of the forest average. A pixel at lambda_obs receives sources with z_s in [lambda/1205-1, lambda/1040-1]:
bin 2 [2.380, 2.916], bin 4 [2.855, 3.466], bin 6 [3.329, 4.016]; z_src (central) = 2.638, 3.149, 3.660.

Quantities from the production readers and weights (12-node Gauss-Legendre in z_s over the range, flattened (z_s,m) set through `prepare_forest_weights`, S and B fixed at the fiducial weighting mode). "N_src" = int dN/dz dz over the range / fiducial rho(z_src) L (the ratio requested: average over the range vs value at central z_src, within 0.3% from the (1+z) factor). "Noise" = N(0.14,0.6) ratio, with SNR also averaged over the range (fiducial = 1):

| bin | forest | N_src ratio | noise ratio (range-avg / central) |
|---|---|---|---|
| 1 | lya(qso) | 0.987 | 1.002 |
| 2 | lya(qso) | 1.044 | 0.957 |
| 3 | lya(qso) | 0.971 | 0.955 |
| 4 | lya(qso) | 1.041 | 0.895 |
| 5 | lya(qso) | 1.072 | 0.860 |
| 6 | lya(qso) | 1.269 | 0.763 |
| 2 | lya(lbg) | 1.233 | 0.965 |
| 3 | lya(lbg) | 0.802 | 1.306 |
| 4 | lya(lbg) | 0.858 | 1.069 |
| 5 | lya(lbg), table clamped at 3.29 row | 1.780 | 0.570 |
| 5 | lya(lbg), zero beyond 3.40 (cell edge) | 1.269 | 0.721 |
| 5 | lya(lbg), zero beyond 3.29 | 1.097 | 0.809 |
| 6 | lya(lbg), table clamped | 1.001 | 0.908 |
| 6 | lya(lbg), zero beyond 3.40 | 0.104 | 6.80 |
| 6 | lya(lbg), zero beyond 3.29 | 0 | no sources: forest absent |

Reading: for the QSO forest the fiducial is pessimistic (dN/dz convex and SNR at fixed lambda rises toward rest 1205 A, where sources with higher continuum dominate the noise-weighted average), up to 24% in noise in bin 6. For lya(lbg) the sign changes with the table structure (peaked at z=2.8-3.1); bin 3 fiducial is optimistic by 30% in noise because the quadratic spline peak coincides with z_src. The SNR tables vary strongly with z_s at fixed lambda (e.g. QSO r=22.25, lambda=5216: SNR/A 3.58, 1.37, 1.65 for z_s=3.33, 3.66, 4.02), so the source-redshift average of SNR is not represented by the central value.
Sensitivity of the ratio to the weights: rho x0.5 changes N(0.14,0.6) by x1.77-2.00 (A ~ 1/rho), consistent with the alias-dominated regime.

M1 (D2): LBG-forest source density in bins 5-6. What the table implies. Cells of the LBG/LAE tables (centres 2.38..3.29, half-widths ~0.11) span z=2.27-3.40. Nothing in the file says the sample is zero outside; the density beyond 3.40 is simply absent. The rows decline steeply (700, 665, 113 deg^-2 per unit z at 2.83, 3.07, 3.29), so the population behind the z=3.06-3.41 forest is at most the tail of the table. For bin 6 (lambda=5216 A) sources lie at z_s in [3.329,4.016]; only 3.329-3.40 (rest frame 1186-1205 A, about 10% of the forest path) is inside the table cell range. With the table cut at 3.40 the pixel is crossed by 8 LBG forests per deg^2 (N_src ratio 0.104); the clamp assigns 113.5 deg^-2 per unit z out to z_s=4.0 (78 per deg^2). For bin 5 (lambda=4931, z_s 3.09-3.74) the range-averaged clamped density is 1.78x the central value; cutting at 3.40 gives 1.27x, at 3.29 1.10x, so the clamp accounts for about 0.5 of the 1.78. In bin 6 the clamp makes lya(lbg) noise 6.8x too small relative to truncation at the cell edge (forest absent for a cut at 3.29); in bin 5 it is 1.4x too small relative to the cell-edge cut. Physically, whether LBGs exist at z>3.4 in the DESI-2 sample is a targeting question (u-dropout selection efficiency falls for z>3.3-3.5; g-dropouts at z~4 are a different sample). The current input cannot answer it, but constant continuation is the least conservative choice consistent with a decreasing table (last-row drop 5.9x from 3.07 to 3.29). Bins 5-6 lya(lbg) auto and its crosses depend on it. Verdict: unjustified as a default without a flag; report joint sigma with the density beyond the table set to zero (Requested run R1).

Galaxies. Bin-averaged (volume-weighted) over local n-bar at z_eval:

| tracer | bin1 | bin2 | bin3 | bin4 | bin5 | bin6 |
|---|---|---|---|---|---|---|
| QSO avg/local | 0.996 | 1.000 | 1.003 | 0.988 | 1.006 | 0.988 |
| LBG avg/local | (no table support) | 1.096 (0.959 support-limited) | 1.115 | 0.962 | 0.964 | 1.804 (1.767 support-limited) |
| LAE avg/local | (no table support) | 1.008 (0.870 support-limited) | 0.992 | 0.992 | 1.025 | 0.965 (0.927) |

Local n-bar (h/Mpc)^3 at z_eval: QSO 4.4e-5 ... 5.9e-6; LBG 8.3e-6, 8.2e-6, 4.9e-5, 1.8e-4, 1.8e-4, 3.0e-5; LAE 1.2e-4 ... 8.9e-5. Estimated nP at (0.14,0.6) with the fiducial biases (approx., growth from Planck18, P from the bundled template): QSO 0.33, 0.26, 0.19, 0.13, 0.09, 0.06; LBG 0.06, 0.06, 0.29, 0.93, 0.83, 0.12; LAE 0.17-0.28. All below or near unity (FR14 quotes nP_0.14,0.6=0.09 for DESI-1 QSO at z~2.5). The LBG bin-6 local value underestimates the bin-mean density by 1.8 because dN/dz falls from ~400 at z=3.18 to 113 at 3.29; that is conservative for bin 6.
For QSO and LAE galaxies the local approximation is accurate to <3%.

## 3. SNR model

| Item | Verdict | Assessment |
|---|---|---|
| sigma^2 = 1/(SNR^2 Delta lambda), SNR per A | justified | The tables are labelled "S/N per Ang for mean quasar with mean forest", so 1/SNR is the noise on delta_F if the noise is measured relative to F-bar C. Note the cancellation: the per-pixel weight and P_pixel=I3 Delta v/(I1^2 L) depend on sigma^2 Delta v, which is independent of Delta lambda at fixed SNR/A. The planned pixel 0.4/1.6 A variations are therefore null tests (only W changes; sinc at k=0.5, mu=1 is 0.996 (0.8 A), ~0.985 (1.6 A) in amplitude). Any larger change indicates a defect. |
| sqrt(N_exp/4) | questionable | Correct for sky/photon-limited noise at fixed exposure length; ignores read noise per exposure and assumes that a single exposure is the same length. The header EXPTIME=4000 is total time in the LBG generation command (`--total-exptime 4000 --nexp 4`) but the QSO command has only `--nexp 4`; the fishhighz reader message calls EXPTIME "per-exposure" and the annotated INI says "fixed header EXPTIME per exposure". If the QSO EXPTIME=4000 is per exposure, QSO and LBG tables differ in depth by 4x in time; must be confirmed against desi_quicklya.py defaults (not in the repository). |
| Gaussian smoothing sigma=10 samples along lambda | questionable | About 10 A (sample step 1 A); no smoothing along z_s, where the 0.25-spaced nodes vary non-monotonically (QSO r=19.25, lambda=3600: 5.8, 3.4, 4.0, 2.4, 1.4, 0.64 for z=2.0-3.25). The variation reflects the rest-frame position in the forest and the continuum, so it is physical, but linear interpolation between 0.25 nodes is coarse. The forest average of item 2 is a more faithful estimate than the central value. |
| Bright clamp r=19.25 | justified, inert | 1.6-3.3% of lya(qso) sources are brighter (rho-weighted). Un-clamping with SNR ~ flux (sky-limited) or ~ flux^0.5: A and P_pixel change by <5e-4 in every bin (weights saturate at w=1). |
| Faint cutoff sigma^2=1e20 for r>24.75 | justified, inert for the bundled data | LBG/LAE counts vanish for r>24.25 (LBG 100.0% of the counts have r<=24.25; spline-integrated density beyond 24.75 is <5e-4 of the total), QSO for r>23.5. The SNR table coverage (to 24.75) is adequate; the cutoff is not an artefact. Consequently the planned r_max 25.75/26.25 runs must reproduce the baseline within ~1e-3 (galaxy n-bar and forest I1 are unchanged; the quadrature partition changes slightly). |
| LBG forest quality at the faint end | questionable | rho-weighted median magnitude of LBG sources is 23.8; 58-63% have r>23.75 where SNR/A=0.27-0.45 (bins 1-6 at r=23.75) and 0.11-0.18 at r=24.25. The I1-weighted median magnitude is 23.5-23.8; the fraction of sources with w>0.5 is 2% (bin 2), 40%, 74%, 18%, 42% (bins 3-6). Forests at integrated S/N per forest ~1-3 cannot be continuum fitted individually; template continuum errors (ME07 and MW11 both add continuum error as noise or drop modes; MW11 Sec. 5: continuum error adds n^-1 P_los^cont) are not modelled. This is optimistic for LBGs, more than for QSOs. |
| SNR at one lambda per bin | questionable | Ratio <sigma^2 over the bin>/sigma^2(centre), from the SNR tables evaluated over lambda(z) with z_s from lambda: QSO (r=21-23) 1.17-1.19, 1.04, 0.76, 0.99-1.00, 1.04-1.05, 0.90-0.91 for bins 1-6; LBG 1.09, 1.01, 0.81, 1.10-1.11, 1.13-1.14, 1.06 (pure wavelength dependence, no density change). Bin 3 sits in a local dip of the table. Average over the bin of the result is ~1, so the joint sigma over six bins is less affected than single bins (by ~10-20% for bins 1, 3). |

## 4. Forest length, pixel, continuum

| Item | Verdict | Assessment |
|---|---|---|
| L=c ln(1205/1040)=44147 km/s (~420 h^-1 Mpc at a_v~105) | justified as length; questionable range | ME07: 1041-1185 A, L~330 h^-1 Mpc; FR14: 1041-1185 for broadband, 985-1200 (Lyb included) for BAO. Reaching 1205 A brings the pixel within 10 A of Lya emission (proximity zone, Lya wing, DLA wings, continuum error and QSO redshift-dependent emission-line residuals) and for LBGs into the Lya emission/absorption region of the galaxy. 1040-1185 changes L by 0.885 and z_src by lambda/sqrt(1040*1185): a listed run. |
| pixel width 0.8 A c/lambda | justified | 63 km/s (bin 1) to 46 km/s (bin 6); Nyquist q=0.05 s/km corresponds to k_par~5 h/Mpc, beyond k_max=0.5. |
| Continuum-fit mode loss | questionable (neglected) | ME07 drop the first 2 N_q,los discrete k_par modes ("makes no noticeable difference", cut of order 0.02 h/Mpc for their box); MW11 drop the first modes (mean flux/continuum). fishhighz keeps all k>=0.01 at all mu. A single forest has fundamental mode 2*pi/L=0.015 h/Mpc, and a mean+slope fit removes modes up to ~0.03. Toy (wiggle-only, Kaiser forest, template PK-PKSB, damping 3.26, mu uniform, fixed white noise level; N/P=0.3-10): sigma_par x1.000-1.001 and sigma_perp x1.08-1.10 for k_par<0.015; x1.004-1.009 and x1.19-1.22 for k_par<0.03; x1.02-1.04 and x1.41-1.45 for k_par<0.05. The transverse constraint comes from low mu, so the loss reaches 20-30% of alpha_perp information. This does not reproduce ME07's "no noticeable difference"; their cut is set by the survey box and may be smaller. A toy only; needs the forecast code with a k_par mask (R5). Also no k_perp Nyquist cut (ME07/FR14 impose one, ME07 finds it minor): fishhighz uses the aliasing term instead, which is the correct treatment. |

## 5. P1D

| Item | Verdict | Assessment |
|---|---|---|
| PD2013 fit at z_eval | justified with caveats | PD2013 is a BOSS DR9 measurement of the total forest 1D power (metals, DLAs included) over z=2.1-4.4, k~1e-3 to 2e-2 s/km (range from memory; check). Absorber z (z_eval) is the right redshift for P1D. Bin 1 (z_eval=2.115) is at the low-z edge. Evolution across a bin: ((1+z)/4)^3.55 changes by +-14% across the bin edges at fixed slope; second-order (curvature) error <1%. |
| low-k flattening | questionable but conservative w.r.t. the fit | Below k_floor (5.7e-4 s/km at z=2.35) P1D is held at the maximum of the fitted parabola in ln k. The BAO modes use q=k_par/a_v in 1e-4 to 5e-3 s/km, i.e. the low-k_par part of the BAO range (k_par<0.06 h/Mpc, q<5.7e-4 s/km) lies below the floor, where PD2013 has no data. The weighting mode (3.5e-4 s/km) is on the plateau. The plateau is physically reasonable (P1D flattens where it integrates the 3D power over k_perp) but untested. |
| PD2013 vs McDonald 2006 (ME07) | not quantified | The McDonald 2006 fit coefficients are not in the repository and I did not reproduce them from memory. Sensitivity instead: aliasing is 34-86% of lya(qso) forest noise at (0.14,0.6) (bins 1-6: 0.34, 0.46, 0.51, 0.63, 0.74, 0.86; at (0.2,0.9): 0.33-0.84) and 5-21% of lya(lbg) noise (0.08, 0.05, 0.08, 0.18, 0.21 for bins 2-6). A fractional P1D error e therefore gives about 0.3-0.9 e in QSO-forest noise (and through B in the weights); MW11: a factor 2 in P_los gives a factor 1.4 in n_eff. The P1D at the weighting mode is 16-64 km/s (z=2.1-3.3). A run with P1D x0.8/1.2 (R6) bounds the effect. |

## 6. Shared-sample noise

| Pair | Verdict | Assessment |
|---|---|---|
| lya(qso) x lya(lbg) | justified | Different sightline sets; sampling and pixel noise are independent; the common delta_F field is in the signal covariance through P_F cross terms, which is retained. |
| lya(qso) x qso and lya(lbg) x lbg | justified only with a condition | The same-object term (pixel on sightline i, object i) is the only source of a noise cross term N_Fq(k)=<delta_F(x_i) | object at x_i>, controlled by the linear 1D cross power P1D_Fq(k_par)=int d^2k_perp/(2pi)^2 b_F b_q(1+beta mu^2)(1+f mu^2/b_q)P_lin, which is not small: order 0.3-1 Mpc/h at k_par~0.05 (b_F b_q~0.7, P_lin(0.1)~500 (Mpc/h)^3 at z~2.6, k_perp range ~0.1), comparable to the auto P1D in comoving units (~0.3-0.6 Mpc/h). Estimators that discard pixel-object pairs on the same line of sight (my understanding of the standard Lya x QSO practice; to be verified in the analysis code) remove it exactly, at the cost of the r_perp=0 modes (negligible). The forecast neither states nor enforces this. The non-Gaussian couplings from common Poisson realisations (Cov(P_FF,P_qq) beyond 2 P_Fq^2) are also dropped. Sign of the bias in the joint sigma not established: B row 21 asserts that independence overestimates the gain; I did not verify. |
| lya x galaxy for different objects, e.g. lya(qso) x lbg | justified | Different objects. |

## 7. Bin-1 exclusion of LBG/LAE/lya(lbg)

Densities in bin 1 (z=2.0-2.235, z_eval=2.115):
- LBG and LAE tables have no cells below z=2.27; the local density is the clamped row-2.38 value: LBG 32.6, LAE 455.3 deg^-2 per unit z, giving N(bin 1)=7.7 (LBG) and 107 (LAE) deg^-2 in a range with zero table support. n-bar LBG 8.3e-6, LAE 1.16e-4 (h/Mpc)^3. Exclusion of the galaxy autos and crosses is justified: the input has no information there and the clamp would add a fictitious ~107 LAE deg^-2 (comparable to the bin-2 LAE count).
- lya(lbg) in bin 1 is different: z_src=2.383 lies inside the table (row 2.38: 32.4 per unit z), and at lambda=3787 the forest draws from z_s in [2.14,2.64] of which the table supports [2.27,2.64]; pixel-crossing count 38.6 deg^-2 vs 62.7 for QSO. Exclusion is a selection choice (bin-1 forest limited to the QSO sample; the SNR tables at 3650-3800 A drop by 40-50% for r=21-23 relative to bin-centre values), not forced by the data. The bin-1 lya(lbg) noise ratio (A B+P_pix)/S is 1.3 in bin 2 (comparable bin), so information gain could be non-negligible; size not estimated (R4).
- lyaforecast (B row 27): all 15 pairs in all bins with tables clamped to the row-2.38 values for z<2.38 (galaxy LBG/LAE at z_eval=2.115 and 2.350 both at row 2.38). In bin 1 lyaforecast therefore includes LBG/LAE at clamped densities (LAE 453 deg^-2 per unit z at z_mean=2.115) and lya(lbg) as in the table. fishhighz removes them. Direction: bin-1 sigma_fh > sigma_lf.

## Quantitative checks (summary of what was computed)

- Densities: normalised N per bin (2.0-3.41): QSO 40.8, 31.3, 21.2, 14.2, 9.0, 5.2; LBG 7.7, 8.4, 50.8, 156.3, 155.7, 48.1 (sum 427; 419 in bins 2-6); LAE 107.0, 107.8, 110.0, 95.1, 72.2, 76.6 (sum 569).
- Floor/signed: dN/dz floored vs signed at z_eval differ by 0.44-0.60% for all three tracers.
- Forest weights (production solver, rho as fiducial): A (deg^2) lya(qso) 0.0203, 0.0270, 0.0388, 0.0615, 0.114, 0.274; lya(lbg) bins 2-6: 0.0111, 0.00253, 0.00328, 0.0164, 0.0142. P_pixel: qso 0.63-2.9, lbg 1.4-3.7 deg^2 km/s.
- Aliasing fraction of forest noise at (0.14,0.6): see item 5.
- Bright clamp, floor, and faint cutoff: effects on A, P_pixel below 1e-3.

## Approximations and assumptions in the survey/noise inputs

1. Single evaluation redshift per bin (z_eval geometric): density, SNR, P1D, bias, growth, response. Local n-bar x volume.
2. Single source redshift and single observed wavelength for the forest SNR and density; no integral over source redshift or wavelength.
3. Sources treated as a complete uniform sample at the tabulated density: no fibre assignment/completeness, no targeting efficiency, purity, or redshift failure (FR14 has 0.8 and 0.73 factors).
4. LBG/LAE/QSO z coverage: tables clamp outside; LBG/LAE tables cover 2.27-3.40 only; QSO z 0.05-4.75.
5. Table magnitude edges sharp (r<=23.5 QSO, r<=24.25 LBG/LAE); no magnitude-dependent purity/efficiency.
6. Quadratic spline interpolation of coarse cell counts; first-spacing normalisation; 1e-20 floor.
7. SNR tables: one mean-QSO spectrum and one LBG template (index 2), mean forest; no diversity of continua, BAL, DLA masking losses, sky-line masks, or bad-pixel fractions. N_exp scaling as sqrt.
8. sigma_deltaF=1/SNR exactly: no continuum-fitting error, no mean-flux error, no flux-calibration error (MW11, ME07 discuss these).
9. Noise white in pixel, uncorrelated between pixels and sightlines (ME07 assumption); no sky-subtraction correlation.
10. Weights: single reference mode (k_t=2.4 deg^-1, q=3.5e-4 s/km); iterated FKP-like full-sample I1 (ME07/FR14 use P_S at k=0.07, mu=0.5).
11. Forest range 1040-1205 A (reaches Lya wings), L uniform, no proximity/wing exclusion.
12. P1D from PD2013 with plateau below 5.7e-4 s/km; aliasing uses the fiducial P1D times W^2; P_pixel unsmoothed.
13. Resolution: R read as lambda/FWHM (sigma_v=51 km/s at R=2500) instead of lyaforecast's sigma=c/R (120 km/s); documented.
14. No redshift-error damping for QSO/LBG/LAE: galaxy response W=1 (QSO/LBG redshift errors of a few 100 km/s give exp(-(q sigma_v)^2/2) power factors of 0.7-0.9 at k=0.2, mu=1, not included).
15. Shared-sample noise neglected: qso/lya(qso), lbg/lya(lbg); non-Gaussian couplings neglected.
16. Gaussian covariance, no survey window, area 5000 deg^2 (the same physical volume for all fields in a bin).
17. LAE bias linearly extrapolated outside 2.5-3.0 (1.10 at z=2.0; 2.80 at z=3.29; 1.56 at bin 2); LBG bias constant 3.3.
18. Bin 1 restricted to lya(qso), qso: 12 pairs excluded.
19. Galaxy density enters through local nP<1, hence sigma depends nearly linearly on the target normalisation (0.86-1.05 changes on LBG/LAE).

## Doc mismatches (documentation vs code/data)

1. docs/methods/survey.md states the strict readers "reject extrapolation, negative spline overshoot" and have "no floors, bright caps"; the native INI path uses `adapters/legacy_compat` (z clamp, floor 1e-20, bright clamp, sentinel 1e20, first spacing), documented only in INI_DEFAULTS (accuracy.py:36-87). The document should say which reader class the native path uses.
2. Annotated INI (examples/desi2_accuracy_annotated.ini) says NEXP scaling is "at fixed header EXPTIME per exposure"; the LBG headers were generated with `--total-exptime 4000 --nexp 4`; QSO headers have no exposure-time argument. Ambiguous.
3. SNR headers (`# using ...`) are stale relative to the per-file MAG: all QSO files state `--ab-magnitude r=24.75`, all LBG files `r=24.5`. data/README.md says headers are retained "unchanged" but does not warn that this line is not the file's magnitude.
4. data/README.md and INI comments do not state that LBG/LAE tables cover only z=2.27-3.40 (5 rows), that the QSO counts vanish for r>23.5 and LBG/LAE for r>24.25, and that r_max/density_magnitude_bounds are inert for the bundled tables.
5. Annotated INI describes SNR header columns as `SN(z=2.0) SN(z=2.5) ...`; actual grid is 0.25 in z (2.0-4.75).
6. `z_norm_min = 2.15` is documented as a strict threshold on raw rows; the meaning "z>2.2 in cell terms" and dependence on the float representation of the row centre 2.15 are not stated.
7. RESEARCH_BASELINE: "Correlations involving LBG, LAE, or Ly-alpha(LBG) are excluded" matches code; the rationale (no table support for LBG/LAE galaxies; lya(lbg) is supported) is not given.
8. survey.md `sample_forest_readers`: "Both readers use z_source; observed wavelength uses z_eval": consistent with code (checked).

## Requested runs (full forecasts needed)

R1. LBG/LAE/lya(lbg) density set to zero outside the table cell range [2.27,3.40] (toggle: `redshift_extension=zero` at the cell edge, and a variant with zero beyond the last node 3.29); report per-bin joint sigma(alpha_par, alpha_perp) and lya(lbg)-related individual spectra in bins 4-6; also bin 2 with support-limited LBG/LAE n-bar.
R2. Forest density and SNR averaged over z_s in [lambda/1205-1, lambda/1040-1] (my flattened (z_s,m) construction, lda.py), with clamped and with zero-beyond-3.40 tables; per-bin joint and individual sigma. Expected noise changes in the table above (QSO 0-24%, LBG bin 3 +31%).
R3. Normalisation variants: LBG/LAE with exact cell widths (x1/1.048, x1/1.026 in aggregate), QSO with z_norm_min=2.1 (x0.842 equivalent) and 2.2.
R4. Bin 1 with lya(lbg) and its crosses (table-supported) and, separately, with LBG/LAE clamped (lyaforecast behaviour), to quantify the cost of the exclusion.
R5. Continuum mask k_par>=0.015 and 0.03 h/Mpc (code option needed); toy predicts sigma_perp +8-10% and +19-22%, sigma_par <1%.
R6. P1D amplitude x0.8 and x1.2 (weights and noise), plus a variant without the plateau below 5.7e-4 s/km.
R7. Planned variations: r_max 25.75/26.25 (expected null within ~1e-3); densities x0.5/x2 (forest N(0.14,0.6) x1.77-2.0 / 0.5-0.6); N_exp 2/8 (informative, tests the sqrt scaling assumption and the EXPTIME semantics); R 2000/4000 (W only); pixel 0.4/1.6 A (expected null to <0.5%); forest range 1040-1185 (L x0.885) and 1050-1205 (L x0.934; z_src shifts).
