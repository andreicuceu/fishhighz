# E3: corrective runs requested by F1 (allocation 59090748, nid200232, 12 cases, all status ok)

Run definitions: /pscratch/sd/a/acuceu/fishhighz-review/harness/{cases3.py,toggles3.py,TOGGLES.md ("jobs3: F1 corrections")}. Raw results: /pscratch/sd/a/acuceu/fishhighz-review/runs/<case>/result.json. Tables below are the verbatim output of harness/analyse3.py.

## Coordinator reading of the numbers (joint, combined over bins, Δσ∥ / Δσ⊥ vs fh_baseline)
- k⊥ Nyquist cut (FR14), per-bin sightline density n_2D: +0.40 / +1.23 % (bin 6: +1.7 / +3.0 %); lya(qso) auto σ⊥ +2.1 %. With the more conservative MW11 n_eff for both forests: +1.29 / +5.52 % (bin 6 σ⊥ +10.3 %). The earlier fh_r6_nyq "null" is superseded.
- LBG-forest source density averaged over the 1040–1205 Å range and truncated at the table edge z=3.40 (lya(lbg) only): +2.12 / +1.51 % (bin 6 +21.7 / +16.8 %, lya(lbg) auto +11–12 %). Averaging with the production clamp: −0.55 / −0.39 %. The clamp therefore lowers the combined joint errors by ≈2.1–2.7 % relative to a truncated, range-averaged treatment.
- Bias in b·D terms (FR14): QSO b·D=1.2: +1.19 / +2.14 % (qso auto +16/+20 %); LAE b·D=0.89: −3.04 / −5.12 % (lae auto −22/−28 %). DESI DR1 Lyα b_F=−0.108, β_F=1.743 (as quoted in F1): +3.20 / +4.71 % (forest autos +15 / +30 %).
- Per-row cell-width density normalisation for LBG/LAE (removes the first-spacing inflation): +0.91 / +1.35 %.
- lyaforecast grid convergence: num_k 1000 → −0.04 / +0.33 %; num_k 2000, num_mu 40, 214 mag bins → +0.19 / +0.53 % vs lf_baseline. The lf baseline is converged to ≲0.5 % in the joint; the fh/lf ratios (0.942 / 0.965) change to ≈0.941 / 0.959 against the converged lf.

## 1. Joint, combined over bins (sigma_par, sigma_perp; Delta in %)

| case | baseline | sig_par base | sig_par | d_par % | sig_perp base | sig_perp | d_perp % |
|---|---|---|---|---|---|---|---|
| fh_nyq_qso | fh_baseline | 0.00719 | 0.00722 | +0.40 | 0.00502 | 0.00509 | +1.23 |
| fh_nyq_both | fh_baseline | 0.00719 | 0.00722 | +0.40 | 0.00502 | 0.00509 | +1.24 |
| fh_nyq_qso_mw11 | fh_baseline | 0.00719 | 0.00725 | +0.76 | 0.00502 | 0.00516 | +2.68 |
| fh_nyq_both_mw11 | fh_baseline | 0.00719 | 0.00728 | +1.29 | 0.00502 | 0.00530 | +5.52 |
| fh_lbgforest_zavg_trunc340 | fh_baseline | 0.00719 | 0.00734 | +2.12 | 0.00502 | 0.00510 | +1.51 |
| fh_lbgforest_zavg_clamp | fh_baseline | 0.00719 | 0.00715 | -0.55 | 0.00502 | 0.00500 | -0.39 |
| fh_qso_bD1p2 | fh_baseline | 0.00719 | 0.00728 | +1.19 | 0.00502 | 0.00513 | +2.14 |
| fh_lae_bD0p89 | fh_baseline | 0.00719 | 0.00697 | -3.04 | 0.00502 | 0.00477 | -5.12 |
| fh_bF_desidr1 | fh_baseline | 0.00719 | 0.00742 | +3.20 | 0.00502 | 0.00526 | +4.71 |
| fh_cellwidth_rows | fh_baseline | 0.00719 | 0.00726 | +0.91 | 0.00502 | 0.00509 | +1.35 |
| lf_conv_k1000 | lf_baseline | 0.00763 | 0.00763 | -0.04 | 0.00521 | 0.00522 | +0.33 |
| lf_conv_k2000 | lf_baseline | 0.00763 | 0.00765 | +0.19 | 0.00521 | 0.00523 | +0.53 |

## 2. Joint per bin (Delta sigma_par / sigma_perp in %)

| case | per-bin d_par / d_perp % |
|---|---|
| fh_nyq_qso | b1: +0.00 / +0.01<br>b2: +0.03 / +0.24<br>b3: +0.09 / +0.67<br>b4: +0.20 / +1.23<br>b5: +0.58 / +1.94<br>b6: +1.66 / +3.00 |
| fh_nyq_both | b1: +0.00 / +0.01<br>b2: +0.03 / +0.24<br>b3: +0.09 / +0.67<br>b4: +0.20 / +1.23<br>b5: +0.58 / +1.97<br>b6: +1.67 / +3.06 |
| fh_nyq_qso_mw11 | b1: +0.32 / +2.39<br>b2: +0.43 / +3.12<br>b3: +0.59 / +3.52<br>b4: +0.54 / +2.21<br>b5: +0.98 / +2.30<br>b6: +1.91 / +3.06 |
| fh_nyq_both_mw11 | b1: +0.32 / +2.39<br>b2: +1.43 / +7.99<br>b3: +0.88 / +5.76<br>b4: +0.73 / +3.65<br>b5: +1.77 / +5.71<br>b6: +3.07 / +10.27 |
| fh_lbgforest_zavg_trunc340 | b1: +0.00 / +0.00<br>b2: -0.21 / -0.25<br>b3: +2.32 / +2.63<br>b4: +0.47 / +0.51<br>b5: -1.90 / -1.68<br>b6: +21.67 / +16.75 |
| fh_lbgforest_zavg_clamp | b1: +0.00 / +0.00<br>b2: -0.21 / -0.25<br>b3: +2.32 / +2.63<br>b4: +0.39 / +0.41<br>b5: -3.64 / -3.30<br>b6: -1.81 / -1.60 |
| fh_qso_bD1p2 | b1: +2.85 / +5.54<br>b2: +1.94 / +3.92<br>b3: +1.21 / +2.37<br>b4: +0.64 / +1.14<br>b5: +0.74 / +1.22<br>b6: +1.15 / +2.14 |
| fh_lae_bD0p89 | b1: +0.00 / +0.00<br>b2: -7.45 / -12.79<br>b3: -4.51 / -7.83<br>b4: -1.96 / -3.21<br>b5: -1.54 / -2.37<br>b6: -2.24 / -3.84 |
| fh_bF_desidr1 | b1: +5.10 / +8.53<br>b2: +4.37 / +7.48<br>b3: +3.49 / +6.00<br>b4: +2.09 / +3.25<br>b5: +1.90 / +2.31<br>b6: +4.96 / +6.26 |
| fh_cellwidth_rows | b1: +0.00 / +0.00<br>b2: +0.19 / +0.23<br>b3: +0.80 / +1.03<br>b4: +1.71 / +2.44<br>b5: +1.42 / +1.94<br>b6: +0.00 / +0.00 |
| lf_conv_k1000 | b1: -0.18 / +0.19<br>b2: -0.12 / +0.26<br>b3: +0.06 / +0.43<br>b4: +0.01 / +0.36<br>b5: -0.11 / +0.27<br>b6: +0.06 / +0.46 |
| lf_conv_k2000 | b1: +0.05 / +0.39<br>b2: +0.12 / +0.47<br>b3: +0.28 / +0.63<br>b4: +0.23 / +0.56<br>b5: +0.11 / +0.46<br>b6: +0.27 / +0.65 |

## 3. Autos (and lya(qso)xqso), combined over bins: Delta sigma_par / sigma_perp in % (`-` = not computed or unavailable)

| case | lya(qso) auto | lya(lbg) auto | qso auto | lbg auto | lae auto | lya(qso)xqso |
|---|---|---|---|---|---|---|
| fh_nyq_qso | +0.46 / +2.08 | +0.00 / +0.00 | +0.00 / +0.00 | +0.00 / +0.00 | -0.00 / -0.00 | +0.61 / +2.74 |
| fh_nyq_both | +0.46 / +2.08 | +0.00 / +0.01 | +0.00 / +0.00 | +0.00 / +0.00 | -0.00 / -0.00 | +0.61 / +2.74 |
| fh_nyq_qso_mw11 | +1.48 / +7.72 | +0.00 / +0.00 | +0.00 / +0.00 | +0.00 / +0.00 | -0.00 / -0.00 | +1.86 / +9.92 |
| fh_nyq_both_mw11 | +1.48 / +7.72 | +1.11 / +6.80 | +0.00 / +0.00 | +0.00 / +0.00 | -0.00 / -0.00 | +1.86 / +9.92 |
| fh_lbgforest_zavg_trunc340 | +0.00 / +0.00 | +11.39 / +12.20 | +0.00 / +0.00 | +0.00 / +0.00 | +0.00 / +0.00 | +0.00 / +0.00 |
| fh_lbgforest_zavg_clamp | +0.00 / +0.00 | -3.42 / -2.50 | +0.00 / +0.00 | +0.00 / +0.00 | +0.00 / +0.00 | +0.00 / +0.00 |
| fh_qso_bD1p2 | +0.00 / +0.00 | +0.00 / +0.00 | +16.24 / +19.88 | +0.00 / +0.00 | +0.00 / +0.00 | +7.86 / +9.32 |
| fh_lae_bD0p89 | +0.00 / +0.00 | +0.00 / +0.00 | +0.00 / +0.00 | +0.00 / +0.00 | -22.06 / -28.35 | +0.00 / +0.00 |
| fh_bF_desidr1 | +15.20 / +29.64 | +15.14 / +29.71 | +0.00 / +0.00 | +0.00 / +0.00 | +0.00 / +0.00 | +7.07 / +14.54 |
| fh_cellwidth_rows | +0.00 / +0.00 | +1.60 / +2.24 | +0.00 / +0.00 | +2.75 / +3.18 | +1.94 / +2.15 | +0.00 / +0.00 |
| lf_conv_k1000 | -0.72 / -0.08 | +1.32 / +2.62 | -0.38 / +0.09 | -0.29 / +0.14 | -0.61 / -0.01 | -0.37 / +0.09 |
| lf_conv_k2000 | -0.50 / +0.11 | +1.53 / +2.79 | -0.19 / +0.26 | -0.09 / +0.33 | -0.43 / +0.16 | -0.16 / +0.28 |

## 4. Nyquist scales (from toggle_log; first record per bin and forest)

| case | z_eval | forest | n_2D [deg^-2] | n_eff MW11 | 1/A | k_Nyq(used) [h/Mpc] | cut applied |
|---|---|---|---|---|---|---|---|
| fh_nyq_qso | 2.1153 | lya(qso) | 63.54 | 20.56 | 49.18 | 0.389 | yes |
| fh_nyq_qso | 2.3504 | lya(qso) | 43.16 | 19.74 | 37.08 | 0.303 | yes |
| fh_nyq_qso | 2.3504 | lya(lbg) | 156.25 | 7.41 | 89.75 | 0.577 | no |
| fh_nyq_qso | 2.5856 | lya(qso) | 31.34 | 14.37 | 25.81 | 0.246 | yes |
| fh_nyq_qso | 2.5856 | lya(lbg) | 435.62 | 23.75 | 395.16 | 0.919 | no |
| fh_nyq_qso | 2.8207 | lya(qso) | 18.82 | 10.64 | 16.27 | 0.183 | yes |
| fh_nyq_qso | 2.8207 | lya(lbg) | 320.73 | 29.72 | 304.87 | 0.756 | no |
| fh_nyq_qso | 3.0558 | lya(qso) | 9.73 | 6.59 | 8.78 | 0.127 | yes |
| fh_nyq_qso | 3.0558 | lya(lbg) | 73.61 | 11.69 | 60.83 | 0.349 | no |
| fh_nyq_qso | 3.2909 | lya(qso) | 3.93 | 3.14 | 3.65 | 0.078 | yes |
| fh_nyq_qso | 3.2909 | lya(lbg) | 77.87 | 16.27 | 70.39 | 0.348 | no |
| fh_nyq_both | 2.1153 | lya(qso) | 63.54 | 20.56 | 49.18 | 0.389 | yes |
| fh_nyq_both | 2.3504 | lya(qso) | 43.16 | 19.74 | 37.08 | 0.303 | yes |
| fh_nyq_both | 2.3504 | lya(lbg) | 156.25 | 7.41 | 89.75 | 0.577 | yes |
| fh_nyq_both | 2.5856 | lya(qso) | 31.34 | 14.37 | 25.81 | 0.246 | yes |
| fh_nyq_both | 2.5856 | lya(lbg) | 435.62 | 23.75 | 395.16 | 0.919 | yes |
| fh_nyq_both | 2.8207 | lya(qso) | 18.82 | 10.64 | 16.27 | 0.183 | yes |
| fh_nyq_both | 2.8207 | lya(lbg) | 320.73 | 29.72 | 304.87 | 0.756 | yes |
| fh_nyq_both | 3.0558 | lya(qso) | 9.73 | 6.59 | 8.78 | 0.127 | yes |
| fh_nyq_both | 3.0558 | lya(lbg) | 73.61 | 11.69 | 60.83 | 0.349 | yes |
| fh_nyq_both | 3.2909 | lya(qso) | 3.93 | 3.14 | 3.65 | 0.078 | yes |
| fh_nyq_both | 3.2909 | lya(lbg) | 77.87 | 16.27 | 70.39 | 0.348 | yes |
| fh_nyq_qso_mw11 | 2.1153 | lya(qso) | 63.54 | 20.56 | 49.18 | 0.221 | yes |
| fh_nyq_qso_mw11 | 2.3504 | lya(qso) | 43.16 | 19.74 | 37.08 | 0.205 | yes |
| fh_nyq_qso_mw11 | 2.3504 | lya(lbg) | 156.25 | 7.41 | 89.75 | 0.126 | no |
| fh_nyq_qso_mw11 | 2.5856 | lya(qso) | 31.34 | 14.37 | 25.81 | 0.167 | yes |
| fh_nyq_qso_mw11 | 2.5856 | lya(lbg) | 435.62 | 23.75 | 395.16 | 0.215 | no |
| fh_nyq_qso_mw11 | 2.8207 | lya(qso) | 18.82 | 10.64 | 16.27 | 0.138 | yes |
| fh_nyq_qso_mw11 | 2.8207 | lya(lbg) | 320.73 | 29.72 | 304.87 | 0.230 | no |
| fh_nyq_qso_mw11 | 3.0558 | lya(qso) | 9.73 | 6.59 | 8.78 | 0.105 | yes |
| fh_nyq_qso_mw11 | 3.0558 | lya(lbg) | 73.61 | 11.69 | 60.83 | 0.139 | no |
| fh_nyq_qso_mw11 | 3.2909 | lya(qso) | 3.93 | 3.14 | 3.65 | 0.070 | yes |
| fh_nyq_qso_mw11 | 3.2909 | lya(lbg) | 77.87 | 16.27 | 70.39 | 0.159 | no |
| fh_nyq_both_mw11 | 2.1153 | lya(qso) | 63.54 | 20.56 | 49.18 | 0.221 | yes |
| fh_nyq_both_mw11 | 2.3504 | lya(qso) | 43.16 | 19.74 | 37.08 | 0.205 | yes |
| fh_nyq_both_mw11 | 2.3504 | lya(lbg) | 156.25 | 7.41 | 89.75 | 0.126 | yes |
| fh_nyq_both_mw11 | 2.5856 | lya(qso) | 31.34 | 14.37 | 25.81 | 0.167 | yes |
| fh_nyq_both_mw11 | 2.5856 | lya(lbg) | 435.62 | 23.75 | 395.16 | 0.215 | yes |
| fh_nyq_both_mw11 | 2.8207 | lya(qso) | 18.82 | 10.64 | 16.27 | 0.138 | yes |
| fh_nyq_both_mw11 | 2.8207 | lya(lbg) | 320.73 | 29.72 | 304.87 | 0.230 | yes |
| fh_nyq_both_mw11 | 3.0558 | lya(qso) | 9.73 | 6.59 | 8.78 | 0.105 | yes |
| fh_nyq_both_mw11 | 3.0558 | lya(lbg) | 73.61 | 11.69 | 60.83 | 0.139 | yes |
| fh_nyq_both_mw11 | 3.2909 | lya(qso) | 3.93 | 3.14 | 3.65 | 0.070 | yes |
| fh_nyq_both_mw11 | 3.2909 | lya(lbg) | 77.87 | 16.27 | 70.39 | 0.159 | yes |
(k_Nyq column: with `n_2D` for cases without `_mw11` and with MW11 n_eff for the `_mw11` cases.)

## 5. New vs related earlier cases (joint combined, Delta vs the fishhighz baseline, %)

| case | compared | d_par % | d_perp % | lya(qso) auto d_par % | lya(lbg) auto d_par % |
|---|---|---|---|---|---|
| fh_nyq_qso | fh_nyq_qso (new) | +0.40 | +1.23 | +0.46 | +0.00 |
| fh_nyq_qso | fh_r6_nyq | +0.00 | +0.00 | +0.00 | +0.00 |
| fh_nyq_both | fh_nyq_both (new) | +0.40 | +1.24 | +0.46 | +0.00 |
| fh_nyq_both | fh_r6_nyq | +0.00 | +0.00 | +0.00 | +0.00 |
| fh_nyq_qso_mw11 | fh_nyq_qso_mw11 (new) | +0.76 | +2.68 | +1.48 | +0.00 |
| fh_nyq_qso_mw11 | fh_r6_nyq | +0.00 | +0.00 | +0.00 | +0.00 |
| fh_nyq_both_mw11 | fh_nyq_both_mw11 (new) | +1.29 | +5.52 | +1.48 | +1.11 |
| fh_nyq_both_mw11 | fh_r6_nyq | +0.00 | +0.00 | +0.00 | +0.00 |
| fh_lbgforest_zavg_trunc340 | fh_lbgforest_zavg_trunc340 (new) | +2.12 | +1.51 | +0.00 | +11.39 |
| fh_lbgforest_zavg_trunc340 | fh_dzero_340 | +3.76 | +2.69 | -0.00 | +19.54 |
| fh_lbgforest_zavg_trunc340 | fh_zavg | -1.12 | -0.90 | -4.40 | -3.02 |
| fh_lbgforest_zavg_trunc340 | fh_lbgforest_zavg_clamp | -0.55 | -0.39 | +0.00 | -3.42 |
| fh_lbgforest_zavg_clamp | fh_lbgforest_zavg_clamp (new) | -0.55 | -0.39 | +0.00 | -3.42 |
| fh_lbgforest_zavg_clamp | fh_zavg | -1.12 | -0.90 | -4.40 | -3.02 |
| fh_qso_bD1p2 | fh_qso_bD1p2 (new) | +1.19 | +2.14 | +0.00 | +0.00 |
| fh_qso_bD1p2 | fh_qso_b12 | +5.40 | +9.40 | -0.00 | +0.00 |
| fh_lae_bD0p89 | fh_lae_bD0p89 (new) | -3.04 | -5.12 | +0.00 | +0.00 |
| fh_lae_bD0p89 | fh_lae_b089 | +6.09 | +9.45 | -0.00 | +0.00 |
| fh_bF_desidr1 | fh_bF_desidr1 (new) | +3.20 | +4.71 | +15.20 | +15.14 |
| fh_bF_desidr1 | fh_bF117_beta167 | +1.80 | +3.07 | +8.08 | +8.13 |
| fh_cellwidth_rows | fh_cellwidth_rows (new) | +0.91 | +1.35 | +0.00 | +1.60 |
| fh_cellwidth_rows | fh_cellwidth | +1.10 | +1.49 | -0.00 | +2.49 |
| lf_conv_k1000 | lf_conv_k1000 (new) | -0.04 | +0.33 | -0.72 | +1.32 |
| lf_conv_k2000 | lf_conv_k2000 (new) | +0.19 | +0.53 | -0.50 | +1.53 |

## 6. lyaforecast converged grid vs lf_baseline, fishhighz baselines relative to both (joint combined)

| case | ratio | d_par % | d_perp % |
|---|---|---|---|
| lf_conv_k1000 | lf_conv/lf_baseline | -0.04 | +0.33 |
| lf_conv_k1000 | fh_baseline / lf_conv_k1000  - 1 | -5.74 | -3.85 |
| lf_conv_k1000 | fh_baseline / lf_baseline  - 1 | -5.77 | -3.53 |
| lf_conv_k1000 | fh_cum6_lf / lf_conv_k1000  - 1 | +0.39 | +0.22 |
| lf_conv_k1000 | fh_cum6_lf / lf_baseline  - 1 | +0.35 | +0.55 |
| lf_conv_k1000 | fh_full_compat / lf_conv_k1000  - 1 | +0.04 | -0.33 |
| lf_conv_k1000 | fh_full_compat / lf_baseline  - 1 | +0.00 | -0.00 |
| lf_conv_k2000 | lf_conv/lf_baseline | +0.19 | +0.53 |
| lf_conv_k2000 | fh_baseline / lf_conv_k2000  - 1 | -5.95 | -4.04 |
| lf_conv_k2000 | fh_baseline / lf_baseline  - 1 | -5.77 | -3.53 |
| lf_conv_k2000 | fh_cum6_lf / lf_conv_k2000  - 1 | +0.17 | +0.02 |
| lf_conv_k2000 | fh_cum6_lf / lf_baseline  - 1 | +0.35 | +0.55 |
| lf_conv_k2000 | fh_full_compat / lf_conv_k2000  - 1 | -0.19 | -0.53 |
| lf_conv_k2000 | fh_full_compat / lf_baseline  - 1 | +0.00 | -0.00 |

| case | per-bin lf_conv vs lf_baseline d_par / d_perp % |
|---|---|
| lf_conv_k1000 | b1: -0.18 / +0.19<br>b2: -0.12 / +0.26<br>b3: +0.06 / +0.43<br>b4: +0.01 / +0.36<br>b5: -0.11 / +0.27<br>b6: +0.06 / +0.46 |
| lf_conv_k2000 | b1: +0.05 / +0.39<br>b2: +0.12 / +0.47<br>b3: +0.28 / +0.63<br>b4: +0.23 / +0.56<br>b5: +0.11 / +0.46<br>b6: +0.27 / +0.65 |

