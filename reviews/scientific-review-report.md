# Scientific review of FishHighz

Review date: 2026-09-24. This is an independent scientific and implementation assessment of the current working tree, including its local full-shape additions. It is a review, not a change of scientific prescription or acceptance of a forecast. Production code, inputs and defaults were preserved.

## Scientific conclusions

The Gaussian covariance, mode counting, mean-derivative Fisher contraction, full-sample forest moments, and nuisance marginalization are internally consistent in the tested calculations. The legacy magnitude-prefix recurrence is a real error relative to the published population-integral prescription, and the native FishHighz route avoids it. The adopted historical full-sum update is, however, a distinct approximation from both the aliasing-inclusive McDonald prescription and McQuinn–White sightline weighting. Its convergence establishes a numerical solution to that chosen recurrence, not optimality.

The most consequential remaining questions concern the information model: AP derivatives include changes of the damping envelope and Kaiser angle; mixed forest–galaxy reconstruction uses an assumed cross width; smooth full-shape information is undamped; and nuisance freedom, nonlinear physics and observational effects are limited. These are documented choices or physical approximations rather than newly discovered accidental numerical defects. Their validity must be assessed separately from the small numerical residuals below.

## Reviewed state and evidence

FishHighz HEAD: `d2e10d4b6b997de8df7712b76c5309f44874124d`. lyaforecast HEAD: `5abe8bcf8d12cc31d1f5ecd89c87e739c6a14b81`. The FishHighz working tree additionally modified `public.py`, `survey_config.py`, and two user documents, and contained untracked `full_shape.py` and its tests. The legacy working tree was clean. The [initial state](../.validation/scientific-review-20260924T164425Z/state/fishhighz.json), [legacy state](../.validation/scientific-review-20260924T164425Z/state/lyaforecast.json), diffs, untracked-file snapshots, and environment identify the actual code reviewed; a commit identifier alone does not describe it.

All 44 S2 input hashes matched their files. Eight of 61 recorded FishHighz module hashes differed. The inspected forecast, survey, and full-sum-kernel changes are docstring-only at executable-AST level; accuracy constants and magnitude integration have moved into shared production modules, with additional reader flexibility. The moved functions are not byte- or AST-identical, so numerical agreement is checked instead of inferred. Historical evidence is labeled throughout. The current full/fixed compatibility results independently reproduce saved S2 endpoints. The [source/input comparison](../.validation/scientific-review-20260924T164425Z/state/s2-evidence-applicability.json) and [semantic check](../.validation/scientific-review-20260924T164425Z/state/historical-semantic-check.json) retain the distinctions.

The full-shape common-cut series predates later category-cut/four-target additions. Its core covariance/Fisher/grid source agrees with the current files. The later LBG+LAE run matches the three changed full-shape source files. Neither saved series exercises unequal category cutoffs on real survey inputs; synthetic independent checks cover the interval algebra. This is a coverage limitation, not an unexecuted required forecast.

## Paper conventions and mathematical reconstruction

The retrieved versions are [McDonald & Eisenstein, arXiv:astro-ph/0607122v1](https://arxiv.org/pdf/astro-ph/0607122) (ME), [McQuinn & White, arXiv:1102.1752v3](https://arxiv.org/pdf/1102.1752) (MW), and [Font-Ribera et al., arXiv:1308.4164v2](https://arxiv.org/pdf/1308.4164) (FR). The ME arXiv revision is dated 2006; the journal reference is conventionally 2007. Equations and section numbers refer to these retrieved versions. ArXiv screenshot requests failed for relevant pages, so visual verification could not be completed; equations were checked against extracted text, dimensions, and independent algebra. This limitation does not turn a numerical reconstruction into a facsimile verification.

### Forest sampling and weights

With magnitude integration weights absorbed into $r_i=\rho_i q_i$,

$$
\begin{aligned}
I_1 &= \sum_i r_iw_i, & I_2 &= \sum_i r_iw_i^2, & I_3 &= \sum_i r_iv_iw_i^2, \\
A &= \frac{I_2}{LI_1^2}, & P_{\rm pixel} &= \frac{l_p I_3}{LI_1^2}.
\end{aligned}
$$

These are population moments, as in ME §II.B, equations 14–18 and FR §IV.B, equations 23–29. A magnitude-prefix $I_{1,i}$ defines a different nonlinear update, not a better quadrature. In the native angular/velocity convention, $\rho$ has units deg⁻² (km/s)⁻¹ mag⁻¹, $A$ is deg², and $P_{\rm pixel}$ is deg² km/s. The physical forest noise is

$$
 N_F(k,\mu)=\left[A P_{1D}(k\mu/a_v)W_F^2(k\mu/a_v)+P_{\rm pixel}\right]\frac{d_{\rm deg}^2}{a_v}.
$$

Here $a_v=H/[h(1+z)]$ and $d_{\rm deg}=hD_M\pi/180$. Instrumental response multiplies the 1D aliasing power once; additive pixel noise is not smoothed again.

The common seed is $B/(B+l_pv_i)$. The adopted update is

$$
 w_i' = \frac{PLI_1+B}{PLI_1+B+l_pv_i},
$$

whereas the implemented aliasing-inclusive alternative replaces $B$ in this expression by $BI_2/I_1$. ME specifies an aliasing-inclusive representative signal, but not a unique initialization/refresh/stopping algorithm. MW equations 11–13 instead yield sightline weights proportional to $[P_{\rm los}+P_{N,i}]^{-1}$, or normalized $\nu_i=P_{\rm los}/(P_{\rm los}+P_{N,i})$, under their estimator/covariance approximations. These are not interchangeable claims of optimality. The native representative weighting mode is (k_transverse,k_parallel)=(2.4 deg⁻¹,0.00035 s/km), so its comoving location varies with redshift; ME/FR discuss a representative comoving k=0.07 h/Mpc, mu=0.5. This is an additional declared approximation, not merely unit conversion. No new reference-mode sensitivity study was performed. The final $A,P_{\rm pixel}$ expressions remain common; the historical update does not replace final aliasing by $1/(LI_1)$.

At fixed intrinsic signal, an independent two-node Decimal calculation isolates the legacy bug: first-node weight 0.411764706 with full $I_1$, versus 0.333333333 with its prefix, a 19.05% difference. This illustrates the mechanism only; it is not the DESI-2 forecast impact. See the [forest derivation and controls](../.validation/scientific-review-20260924T164425Z/agents/forest-review.md).

### Signal, AP mapping and damping

FishHighz evaluates

$$
 P_{ij}=G\left[Q_s K_{i,s}K_{j,s}P_s(k_s)+Q_wK_{i,w}K_{j,w}D_{ij}(k_w,\mu_w)P_w(k_w)\right],
$$

with signed $P_w=PK-PKSB$, $P_s=PKSB$, forest $K_i=b_i(1+\beta_i\mu^2)$, galaxy $K_i=b_i+f\mu^2$, and fiducial $G=[\sigma_8(z)/\sigma_8(z_{\rm template})]^2$. Observed signal is $W_iW_jP_{ij}$.

Each component uses $k_\parallel'=k\mu/a_\parallel$, $k_\perp'=k\sqrt{1-\mu^2}/a_\perp$, and $Q=(a_\parallel a_\perp^2)^{-1}$, consistent with FR equation 19's fixed observed coordinates. In the `alpha_phi` basis, $\alpha=\sqrt{a_\parallel a_\perp}$, $\phi=a_\perp/a_\parallel$. This alpha is **not** the volume-isotropic $\alpha_{\rm iso}=(a_\parallel a_\perp^2)^{1/3}=\alpha\phi^{1/6}$.

For each axis, $\Sigma_{ij}^2=(\Sigma_i^2+\Sigma_j^2)/2$, and
$D_{ij}=\exp[-(k_\parallel'^2\Sigma_{ij,\parallel}^2+k_\perp'^2\Sigma_{ij,\perp}^2)/2]$.
This implies $D_{ij}=\sqrt{D_{ii}D_{jj}}$, coherent under a factorized propagator approximation for the oscillatory term. It does not derive the actual reconstructed galaxy–unreconstructed forest cross response from their auto widths. Even with fixed widths, AP remapping gives at the fiducial
$\partial_{a_\parallel}\ln D=\Sigma_\parallel^2k_\parallel^2$ and
$\partial_{a_\perp}\ln D=\Sigma_\perp^2k_\perp^2$.
The baseline widths are

$$
\begin{aligned}
\Sigma_\perp &= \frac{3.26\,\sigma_8(z)}{\sigma_{8,\mathrm{damping\,reference}}\sqrt{r}}\,h^{-1}\mathrm{Mpc}, \\
\Sigma_\parallel &= (1+f_{\mathrm{fid}})\Sigma_\perp,
\end{aligned}
$$

with $r=2$ for galaxies and $r=1$ for forests. Thus $r=2$ reduces width by $1/\sqrt{2}$, not the factor 0.5 used for FR's named 50% reconstruction. FR §IV.A.1 deliberately keeps damping and RSD outside isolated-BAO derivatives; its broadband prescription retains RSD derivatives but keeps information damping outside them. FishHighz therefore implements a different information prescription. Legacy back-differencing of the damped residual also includes a damping slope and is not the exact FR frozen-envelope prescription. See the [signal review](../.validation/scientific-review-20260924T164425Z/agents/signal-review.md).

### Covariance, information and full-shape parameters

For $0\leq\mu\leq1$, $N_q=V k^2w_kw_\mu/(2\pi^2)$,

$$
\begin{aligned}
T_{ij} &= W_iW_jP_{ij}+N_{ij}, \\
C_{(ij)(mn)} &= \frac{T_{im}T_{jn}+T_{in}T_{jm}}{N_q}, \\
F_{ab} &= \sum_q J_{q,a}^{T}C_q^{-1}J_{q,b}.
\end{aligned}
$$

The auto variance is $2T_{ii}^2/N_q$; no additional hemisphere factor belongs in Fisher information. These conventions match FR equations 11 and 20. Negative forest–galaxy power is retained. Unmeasured field spectra required by the selected covariance remain available. Individual-spectrum Fishers use their own variances and must not be summed to reconstruct the joint information.

Full shape introduces independent smooth/wiggle $(\alpha,\phi)$, bin-local growth, biases and common forest $\beta$; it has no direct legacy counterpart. The full active Fisher is inverted before selecting target covariance. Reported $f\sigma_{8,\rm fid}$ uses the complete covariance Jacobian, including off-diagonal entries. With free biases and fixed spectrum shape, writing $P_g=(b\sigma_8+f\sigma_8\mu^2)^2P_{\mathrm{shape}}$ shows that $\sigma_{8,\mathrm{fid}}$ can be a normalization convention for the identifiable f sigma8 product. Fixing sigma8 is therefore not, by itself, evidence of an artificially strong prior on that product; adding a separate amplitude can instead introduce a redundant degeneracy. The result remains conditional on shape, geometry, nuisance freedom and physical modeling, and is not an independently varied cosmological transfer-function forecast. Damping widths, weights, response, noise, grid and covariance remain fiducial during differentiation. Singular active constraints remain unavailable. Category cut intervals are disjoint and retain covariance closure. See the [covariance review](../.validation/scientific-review-20260924T164425Z/agents/covariance-review.md).

## Scientifically meaningful changes from lyaforecast

Verdicts concern scientific justification, not adoption. “Conditional” means sensible within the named model or input assumption, with the stated limitation. Entries marked “justified” need not produce smaller errors.

| Change or extension | Scientific effect and verdict | Principal source |
|---|---|---|
| Prefix normalization → full-population moments | Corrects the population-integral mismatch. **Justified.** | `kernels/full_sum_weights.py`; legacy `weights.py` |
| Intrinsic-only iterative signal → historical alias-containing update | Changes relative source weights; distinct from ME's full-moment alias signal and MW optimum. **Conditionally justified**, optimality unresolved. | `weights.py`, `kernels/full_sum_weights.py` |
| Three updates → confirmed adaptive stopping | Qualifies the selected recurrence, with explicit nonconvergence. **Justified numerically**; not a continuum theorem. | `kernels/full_sum_weights.py` |
| Pair-dependent preparation → a consistent state per observed forest | Maintains one noise/response state for each sampled field throughout covariance. **Justified.** | `forecast.py`, `survey.py`, `noise.py` |
| Rectangular magnitude integration → partitioned Gauss–Legendre | Improves integration while retaining declared support and interpolation. **Justified when refinement passes.** | `magnitude.py`, `survey_config.py` |
| Signed spline overshoot → named negative-density floor | Avoids negative physical source density but modifies interpolation support. **Conditionally justified**; a floor is not a measured selection function. | `adapters/legacy_compat.py` |
| Rectangular k and midpoint mu → composite Gauss–Legendre | Changes numerical measure, not the mode-count convention. **Justified when refinement passes.** | `grids.py` |
| Centre-volume approximation → integrated volume | More faithful homogeneous-bin geometry at common area. **Justified.** | `geometry.py` |
| Mixed arithmetic/geometric evaluations → consistent geometric redshift | Removes inconsistent fiducial evaluation; within-bin evolution still approximated. **Conditionally justified.** | `survey_config.py` |
| EdS power evolution → CAMB sigma8 ratio; exact bin growth rate | Uses the chosen cosmology consistently. **Justified**, conditional on fixed shape/scale-independent growth. | `cosmology.py`, `models/biases.py` |
| CAMB sampled spectrum → supplied Vega-format template | Changes spectrum interpolation and reference normalization. **Conditionally justified** by matching cosmology and template metadata. | `models/templates.py` |
| Polynomial peak extraction → signed supplied wiggle | Defines the feature explicitly and avoids response-dependent polynomial extraction. **Conditionally justified** by decomposition choice. | `models/templates.py`, `models/kaiser.py` |
| sigma=c/R → physical FWHM conversion | Corrects the meaning of resolving power when R=lambda/FWHM. **Justified.** | `response.py` |
| One selected forest response per pair → W_i W_j | Correctly represents different observed forest responses. **Justified**; equal-width baseline cannot test the distinction alone. | `response.py`, `forecast.py` |
| Undamped fiducial covariance → wiggle-damped covariance | Consistent with the chosen mean, but differs from FR information-damping treatment. **Conditionally justified.** | `models/kaiser.py`, `forecast.py` |
| Forest width for every mixed pair → mean squared auto widths | Changes effective mixed-pair reconstruction. **Conditionally justified** under factorized propagators; actual mixed-reconstruction physics remains **unresolved**. | `kernels/kaiser.py` |
| Projected backward slope → full wiggle inverse-AP derivative | Correct derivative of its defined model, but adds prefactor/angular/envelope information relative to phase-focused BAO. **Conditionally justified**; exact FR BAO equivalence is **unjustified**. | `models/kaiser.py`, `derivatives.py` |
| Rounded c → exact defined c | Small velocity/noise/response/geometry change. **Justified.** | `geometry.py` |
| Native selected bin 1 | Removes unsupported LBG/LAE-related observables before covariance/Fisher assembly. **Justified for this adopted sample**; historic all15 bin1 results are different selections. | native INI, `fields.py` |
| Generic full-shape targets, tied nuisance parameters and explicit priors | Extends the legacy two-parameter BAO calculation; marginalization and prior accounting are explicit. **Conditionally justified** by the intended inference. | `full_shape.py`, `parameters.py`, `results.py` |
| Category-specific cuts and four-target forest option | Retains covariance for active observables and makes parameter selection explicit. **Justified algebraically**; real mixed-cut validation remains limited. | `full_shape.py`, `survey_config.py` |
| Standalone native input preparation and serialized identities | Can change science if preparation differs; baseline comparison checks this boundary. **Justified when matched-input results agree.** | `public.py`, `survey_config.py` |

The P1D fit/low-q prescription, density target normalization, first-spacing raw-cell convention, source-centre approximation, SNR smoothing/clamps/sentinels, independent noise and Gaussian mode covariance are largely **inherited assumptions**, not new improvements. A new public interface does not itself improve them. FishHighz additionally exposes strict readers, alternative weights, arbitrary external models and supplied full noise; those capabilities do not automatically enter the native DESI-2 prescription.

## Numerical evidence

### Live legacy and compatibility comparison

The diagnostic prepares only the required selected spectra from current NewForecast components and calls the actual legacy Fisher implementation, including its derivative calculation. It independently contracts the FishHighz covariance on the same operands. It does not rerun the seven-case suite or silently retain all 15 bin-1 spectra. This targeted invocation does not test the unmodified all-bin public `new_run_forecast()` orchestration end to end.

The maximum legacy/full-compatibility Fisher residual is 1.4e-15. The full/fixed comparison with saved S2 Fisher matrices has maximum residual 6.02e-16. Selected results contain 78 individual and six joint constraints per profile. Elapsed time was 129.0 seconds with one numerical thread.

| Bin | Full sigma(ap) | Full sigma(at) | Fixed sigma(ap) | Fixed sigma(at) | Fixed/full radial change | Fixed/full transverse change |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 0.02472610 | 0.01804447 | 0.02514347 | 0.01835109 | +1.6880% | +1.6993% |
| 2 | 0.02046864 | 0.01462255 | 0.02006916 | 0.01429040 | -1.9517% | -2.2715% |
| 3 | 0.01799679 | 0.01237200 | 0.01735381 | 0.01188115 | -3.5728% | -3.9674% |
| 4 | 0.01501345 | 0.00966468 | 0.01468227 | 0.00943943 | -2.2059% | -2.3307% |
| 5 | 0.01635541 | 0.01079681 | 0.01605049 | 0.01061337 | -1.8643% | -1.6990% |
| 6 | 0.02313372 | 0.01741007 | 0.02158659 | 0.01642430 | -6.6878% | -5.6621% |

The fixed/full change groups restored full-sample normalization, changed iterative signal and confirmed stopping. It must not be reported as the isolated prefix-sum correction. The scalar fixed-signal control above isolates the bug itself.

### Current native accuracy and physical variations


The current six-bin native accuracy calculation reproduces the saved S2 accuracy Fisher matrices to at most 4.47e-9 relative for the joint constraints and 8.80e-9 for own-covariance individual spectra. All 78 individual and six joint two-parameter Fisher matrices have rank two. This is strong numerical continuity evidence, not scientific acceptance of the BAO signal or noise prescription. The 0.1% error-refinement criterion passes in all six directions tested.

| Bin | z interval | sigma(ap) | sigma(at) | rho(ap,at) |
| --- | --- | ---: | ---: | ---: |
| 1 | 2.000–2.235 | 0.02417269 | 0.01810703 | -0.414632 |
| 2 | 2.235–2.470 | 0.01929890 | 0.01409680 | -0.414161 |
| 3 | 2.470–2.705 | 0.01671429 | 0.01172291 | -0.410890 |
| 4 | 2.705–2.940 | 0.01424716 | 0.00935641 | -0.404950 |
| 5 | 2.940–3.175 | 0.01562206 | 0.01054263 | -0.404254 |
| 6 | 3.175–3.410 | 0.02076141 | 0.01616004 | -0.412447 |

Bin 1 includes only Lyα(QSO) auto, QSO auto and their cross. Bins 2–6 include all 15 spectra. The own-covariance individual Fisher matrix uses each selected spectrum's diagonal Wick variance, while the joint contraction retains the full within-bin covariance. The saved baseline NPZ files include observed Jacobian, all required total powers, mode counts, covariance and pair identities. Direct Wick reconstruction of bins 2 and 6 agrees with saved covariance at 1.2e-16 relative; direct contraction agrees with the native baseline Fisher at 1.6e-14. The small current/S2 difference is consistent with separate CAMB preparations and is not an algebraic contraction failure: current-vs-S2 total-power relative L2 differences are 4.74e-11 and 3.48e-10 in bins 2 and 6.

| One-at-a-time refinement | Maximum joint error change | Maximum individual error change | Maximum A change | Maximum P_pixel change |
| --- | ---: | ---: | ---: | ---: |
| AP step ×0.5 | 0.000229663% | 0.000208615% | 0 | 0 |
| k intervals 128→256 | 0.0000447496% | 0.0000699463% | 0 | 0 |
| μ order 32→64 | 0.00000133905% | 0.00000761850% | 0 | 0 |
| magnitude order 16→32 | 1.22e-13% | 1.02e-12% | 1.12e-12% | 7.11e-13% |
| redshift order 32→64 | 2.72e-10% | 2.72e-10% | 0 | 0 |
| weight rtol 1e-5→1e-6 | 1.12e-11% | 1.42e-10% | 1.13e-9% | 2.83e-10% |

Each forest-weight recurrence reports converged in these runs. The magnitude-order result is a finite refinement of the specified composite quadrature, not a proof of convergence to a continuum weighting limit. No failed dimension required another refinement.

Physical perturbations give substantially larger effects. The table reports percentage changes in joint sigma(ap), sigma(at) relative to the baseline; the raw records retain every individual error, correlation, rank, forest A and P_pixel, and actual numerical settings.

| Perturbation | Bin 1 | Bin 2 | Bin 6 |
| --- | ---: | ---: | ---: |
| All target densities ×0.8, QSO/LBG samples coherent | +10.20%, +14.20% | +9.68%, +13.68% | +12.26%, +15.59% |
| All target densities ×1.2, QSO/LBG samples coherent | -7.07%, -9.73% | -6.69%, -9.35% | -8.42%, -10.62% |
| LBG-forest sampled variance ×0.8 | — | -1.09%, -1.28% | -2.65%, -2.35% |
| LBG-forest sampled variance ×1.2 | — | +0.83%, +0.97% | +2.12%, +1.84% |
| LBG-forest pixel width ×1.2 with native SNR conversion | — | +0.00123%, +0.000324% | +0.00189%, +0.000405% |
| Observed k_max = 0.1 h/Mpc | — | +141.44%, +140.80% | +146.98%, +136.97% |
| Observed k_max = 0.2 h/Mpc | — | +11.95%, +13.93% | +14.00%, +14.56% |

Density and LBG-variance perturbations recompute forest weights and covariance. The LBG pixel-width experiment changes native velocity conversion and per-Å SNR variance; its small BAO effect is consistent with near cancellation in l_p v, while the nonlinear weight state and response were still recomputed. The lower k_max runs alter the fixed observed cuts and mode selection; they do not validate high-k signal physics in the baseline.

The fixed-covariance AP controls give the following percentage changes in joint sigma(ap), sigma(at). The AP prefactor is separated from the damping and RSD-angle terms.

| Controlled change | Bin 2 | Bin 6 |
| --- | ---: | ---: |
| Freeze damping derivative | +0.467%, +0.253% | +0.311%, +0.157% |
| Freeze RSD-angle derivative | +0.112%, -0.303% | +0.0557%, -0.345% |
| Freeze damping and RSD angle | +0.884%, -0.000662% | +0.520%, -0.162% |
| Freeze AP volume prefactor | +0.0959%, -1.371% | +0.0282%, -1.151% |
| Freeze damping, RSD angle and AP prefactor | +0.614%, -0.164% | +0.252%, -0.220% |
| Restore forest widths for mixed forest–galaxy wiggles | +4.130%, +2.294% | +3.108%, +1.377% |

The mixed-width case also updates fiducial total power and Wick covariance. It ranks above either isolated derivative term in these bins. Its Jacobian transformation uses the analytic first-order damping-ratio derivative about the fiducial AP point.

The operator controls use the **current native** baseline Jacobian, total power and covariance, with saved S2 template, field response and model fiducial values to reconstruct the isolated wiggle terms. This is a hybrid to avoid another CAMB run. The current/S2 relative L2 discrepancies in Jacobian and covariance are respectively (2.28e-10, 3.20e-10) in bin 2 and (2.40e-9, 1.74e-9) in bin 6, well below the reported control shifts. A representative mixed-pair five-point finite difference of the damping-plus-RSD-angle term agrees with the analytic subtraction to 5.09e-13 in bin 2 and 3.30e-13 in bin 6. These controls isolate adopted operator choices; they do not establish a physically preferred BAO estimator.

The full public six-bin `Forecast.run` gives 78 available individual constraints, six available joint constraints and 12 explicitly excluded bin-1 spectra. Relative to the direct saved-array calculation, its maximum Fisher discrepancy is 3.96e-16 and maximum BAO-error discrepancy is 4.44e-16. In a separate bin-2 area ×1.2 experiment, signal and noise arrays remain identical, the joint Fisher matrix scales by 1.2 to 3.12e-16 relative, and both errors scale by 0.912870929175277, equal to 1/sqrt(1.2). These checks establish public-route and volume normalization consistency within the model.

Reproducibility: [numerical driver](../.validation/scientific-review-20260924T164425Z/diagnostics/numerical_review.py), [operator controls](../.validation/scientific-review-20260924T164425Z/diagnostics/bao_operator_controls.py), [summary JSON](../.validation/scientific-review-20260924T164425Z/results/numerical/summary.json), [case settings and perturbation operands](../.validation/scientific-review-20260924T164425Z/results/numerical/case-settings.json), [public run](../.validation/scientific-review-20260924T164425Z/results/public-check.json), and [area check](../.validation/scientific-review-20260924T164425Z/results/area-check.json). These links are written for insertion into a FishHighz `docs/` review. Exact source hashes and per-bin input identities are structured in each case JSON. Effective INIs are saved for INI-expressible cases; the direct forest-variance multiplication is documented as a post-survey-assembly perturbation. All calculations used one numerical thread and actual bundled DESI-2 inputs.


The fresh native/fixed comparison, using the same selected spectra, gives the following joint error changes. This groups all accuracy-model and quadrature changes; it is not an isolated attribution.

| Bin | Native/fixed sigma(ap) change | Native/fixed sigma(at) change |
|---|---:|---:|
| 1 | -3.8610% | -1.3300% |
| 2 | -3.8380% | -1.3548% |
| 3 | -3.6852% | -1.3319% |
| 4 | -2.9635% | -0.8795% |
| 5 | -2.6693% | -0.6665% |
| 6 | -3.8227% | -1.6089% |

### Historical attribution, independently reconstructed

All 364 saved S4 records were reopened with `allow_pickle=False`. Independent scalar Wick blocks matched exactly; inversions reproduced saved errors to 6.67e-16. Twelve full fixed/accuracy endpoint Fisher contractions were reconstructed without FishHighz numerical routines, with maximum relative residual 4.10e-14. These validate saved operands and summaries, not a new calculation of 364 forecasts.

The largest isolated historical accuracy-side improvement is mixed-pair damping: introducing it into fixed compatibility changes joint radial errors by −4.105% to −2.376%, transverse by −2.383% to −0.954%. The grouped AP derivative change acts oppositely: +1.020% to +1.174% radial and +1.123% to +1.301% transverse. Fourier quadrature changes radial errors by −0.866% to −0.631%; other differences are usually smaller. These are conditional one-change effects, not additive shares. The historical endpoint accuracy/fixed changes are −3.861% to −2.669% radial and −1.609% to −0.666% transverse. The current endpoint comparison is reported separately above.

The [historical impact table](../docs/archive/reviews/s4-impact-table.md), [profile inventory](../docs/archive/reviews/s4-profile-inventory.md), and [fresh operand reconstruction](../.validation/scientific-review-20260924T164425Z/results/s4-reconstruction.json) retain all 15 attribution groups, reverse directions, interactions and residuals. A numerically favorable mixed-width choice remains physically unresolved.

### Full-shape evidence and independent checks

Independent scalar Wick, dense J-transpose-C-inverse-J, quadrature normalization, and category-interval controls passed. The dense synthetic Fisher residual was 4.44e-16. Saved common-cut results include 36 joint, 468 individual and 72 excluded records; all 36 combined Fisher arrays exactly reconstruct from bin blocks. The later LBG+LAE series provides five available joint and 15 available individual constraints with current extended source.

One weak individual QSO×LBG constraint at kmax=0.08 has a 1.13e-8 relative dense target-covariance reconstruction difference, above the nominal algebraic target. Its Fisher condition number is 2.58e9, and the direct inverse residual is 1.75e-8; this is explicitly a conditioning-limited check, not a claimed 1e-10 reconstruction. The synthetic well-conditioned checks meet the stated target. No pseudoinverse or regularization was introduced.

The existing test selections passed: 175 covariance/Fisher/full-shape tests, plus 382 signal/derivative/template/weights/noise/input/public-interface tests. The commands and logs are retained in the evidence directory. Some shell sessions emitted NERSC MUNGE warnings after otherwise successful commands; the commands exited zero and no numerical check depended on scheduler authentication.

The saved [full-shape cut comparison](../../../forecasts/full_shape/REPORT.md) gives, in bin 2, sigma(alpha_s)=0.01015 and sigma(f sigma8)=0.006948 at kmax=0.30, versus 0.007770 and 0.005911 at kmax=0.50 h/Mpc: reductions of about 23.4% and 14.9%. Wiggle-dilation errors change much less. Together with the new BAO cut comparison, this shows that broadband/growth information continues to increase after most of the adopted BAO information has accumulated. It does not validate that additional information against nonlinear or observational systematics.

## Bugs, mathematical mismatches and proposed dispositions

| ID | Finding and confidence | Consequence | Smallest justified disposition |
|---|---|---|---|
| B1 | **Confirmed legacy bug:** magnitude-prefix I1 is fed into a source-population recurrence; high confidence from source, papers and fixed-signal Decimal check. | Changes source weights and downstream aliasing/noise. Native FishHighz avoids it. | Retain literal legacy behavior only as a labeled reproduction control; do not restore it in native forecasts. |
| M1 | **Information-prescription mismatch:** AP-remapped damping and BAO RSD-angle derivatives differ from FR's treatment; high confidence. Both legacy and native differ in their own ways. | Additional model-dependent distance information; survey magnitude quantified above. | State this distinction in scientific use; if phase-only FR equivalence is intended, propose a separate frozen-envelope/angular derivative prescription for review. |
| A1 | **Mixed-pair reconstruction assumption:** geometric-mean damping is coherent under factorized propagators, but the cross-reconstruction physics is not established by the cited papers; high confidence in implementation, unresolved physical accuracy. | Historically the dominant isolated reduction in joint BAO errors. | Qualify results; physical validation or a conservative alternative is a separate decision, not an automatic bug fix. |
| A2 | **High-k full-shape assumption:** undamped smooth information with linear bias/RSD and no nonlinear nuisance model; high confidence. | Stronger high-k constraints are conditional model information, not validated DESI-2 performance. | State applicable model/cuts; assess nonlinear/observational systematics before scientific use of aggressive cuts. |
| A3 | **Growth interpretation:** f sigma8_fid is the identifiable amplitude coordinate for free biases and fixed shape; high confidence. | The product need not be artificially tightened by the normalization convention; shape, nuisance and physical assumptions still condition it. | Keep the conditioning explicit, including shared biases/beta, priors and geometry freedom. |

No new consequential FishHighz production bug was confirmed in the inspected and tested paths. This is a bounded review conclusion, not a proof that arbitrary extensions, external models or untested physical regimes are correct.

The evidence for each finding is localized as follows; detailed source line references are in the linked scientific subreports.

| Finding | Equation or assumption | Source and reproducible evidence |
|---|---|---|
| B1: population normalization | ME §II.B, Eqs. 13–19; FR §IV.B, Eqs. 23–29 | Legacy `weights.py:184–205,339–356`; FishHighz `kernels/full_sum_weights.py:14–65`; [fixed-signal Decimal control](../.validation/scientific-review-20260924T164425Z/agents/forest-scalar-check.py) and [results](../.validation/scientific-review-20260924T164425Z/agents/forest-scalar-results.json). The isolated six-bin error impact of the prefix correction alone was not quantified; the profile comparison also changes the signal and stopping. |
| M1: BAO derivative ownership | FR §IV.A.1, Eq. 16; §IV.B.1 | `kernels/kaiser.py:23–75`, `models/kaiser.py:247–290`; [operator controls](../.validation/scientific-review-20260924T164425Z/diagnostics/bao_operator_controls.py) and [bin-2 results](../.validation/scientific-review-20260924T164425Z/results/numerical/bin2-operator-controls.json). |
| A1: mixed damping | Factorized oscillatory propagators; FR auto-damping/reconstruction is insufficient to derive the actual mixed response | `kernels/kaiser.py`; [mixed-pair derivation](../.validation/scientific-review-20260924T164425Z/agents/adjudication.md), controls above and saved S4 attribution. The discrepancy from a true reconstructed cross spectrum remains unquantified. |
| A2: high-k full shape | Gaussian linear bias/RSD and fixed-template assumptions; FR §IV.A.2 broadband information damping | `models/kaiser.py`, `full_shape.py`; [signal review](../.validation/scientific-review-20260924T164425Z/agents/signal-review.md), saved cut comparison above. The forecast bias or error degradation from omitted physical effects is unquantified. |
| A3: growth amplitude | RSD amplitude coordinates; FR §IV.A.3 | `full_shape.py` and `parameters.py`; [dense covariance transformation](../.validation/scientific-review-20260924T164425Z/agents/covariance-check.py), [amplitude derivation](../.validation/scientific-review-20260924T164425Z/agents/adjudication.md). No correction to the coordinate transformation is indicated. |

## Scientifically relevant assumptions and approximations

| Assumption | BAO implications | Full-shape implications / evidence limit |
|---|---|---|
| Gaussian fields, independent Fourier cells, homogeneous common window | Standard Fisher approximation; finite geometry and mode coupling can alter precision. | Same limitation; no non-Gaussian covariance or window convolution. |
| Independent redshift bins, one z_eval per bin | Neglects cross-bin covariance and evolution within bins despite integrated volume. | Growth/bias/AP inference remains bin-centre approximation. |
| Shared volume and independent sampling noise | Disjoint sightline/target sampling is assumed; shared underlying signal remains correlated. | Real overlaps and continuum-projected noise need an observational covariance model. Zero cross noise is not zero signal covariance. |
| Effective central background source and fixed forest rest range | Approximates the source/absorber distribution; source redshift used for density and SNR. | Same inherited approximation, not an exact integration over each sightline length. |
| Adopted density/SNR interpolation, target normalization, floors and support | Can alter effective density and weighting; local sensible variations are not universal validation. | Same; the first-spacing raw-cell convention is not independent measurement of irregular cell edges. |
| Parametric P1D and low-q treatment | Sets aliasing and the representative weighting signal; not forced to equal a projection of the chosen P3D. | No automatic small-scale or cosmological consistency with an arbitrary external P3D. |
| Representative-mode, parameter-independent weights | Practical fixed fiducial estimator; convergence says nothing about optimum across modes. | Also fixed as physical parameters vary. D1/D2 comparisons remain deferred. |
| Fixed template/cosmology/growth shape and linear deterministic biases | Wiggle distances depend on decomposition and chosen nuisance treatment. | Conditional shape/growth constraints; no general cosmological transfer-function derivatives or scale-dependent growth model. |
| Wiggle-only damping, prescribed reconstruction, mixed widths | Information prescription differs from a conservative phase-only calculation. | Smooth information is not damped; the FR broadband information-loss prescription is absent. |
| Finite physical k cuts and numerical refinements | Refinement can establish stability of a specified integration. | It cannot establish validity of the linear model at k=0.3–0.5 h/Mpc. |
| Fixed covariance, noise and instrumental response in derivatives | Consistent with the fixed-covariance bandpower Fisher convention. | No extra covariance-derivative information; changing the data estimator is not silently differentiated. |
| Limited nuisance freedom and explicit priors only | Baseline has two AP targets and no unspecified broadband marginalization. | Bias/beta marginalization is real but does not cover nonlinear bias, template uncertainty or all astrophysical freedom. |
| Continuum fitting, survey-window mixing and contaminants not modeled | Precision may be optimistic or estimator-dependent, particularly low parallel modes. | Smooth-shape information is more exposed; metals, high-column-density absorbers and UV/thermal fluctuations are not covered by the simple signal. |
| No additional galaxy redshift-error/Fingers-of-God model in the native baseline | Instrumental forest smoothing does not account for these distinct effects. | High-k radial information remains conditional. |
| Gaussian parameter posterior, local derivatives, rank-based availability | Weak individual constraints need cautious interpretation. | A numerically finite inverse is not evidence of a Gaussian posterior; singular target sets remain unavailable. |

## Reproduction, limitations and review control

[Evidence directory](../.validation/scientific-review-20260924T164425Z/README.md). Scientific subreports: [forest](../.validation/scientific-review-20260924T164425Z/agents/forest-review.md), [signal](../.validation/scientific-review-20260924T164425Z/agents/signal-review.md), [covariance](../.validation/scientific-review-20260924T164425Z/agents/covariance-review.md), and [numerical](../.validation/scientific-review-20260924T164425Z/agents/numerical-review.md). Scripts retain exact operations and saved JSON/NPZ contain identities, matrices and statuses. Fresh runs use one numerical thread and current source imports. No full-shape real forecast, seven-case legacy suite, alternative-weighting/reference study, Slurm action, commit or production fix was performed.

The three scientific reviewers used GPT-6 Sol high. Numerical diagnostics used GPT-6 Sol medium. The coordinator controlled dispatch, reproduced legacy comparisons and historical attribution operands, checked reviewer mathematics and assembled this report. A bounded GPT-6 Astra high adjudication refined the growth-amplitude interpretation and the conditional propagator justification of mixed damping, without identifying a needed production correction; see [adjudication](../.validation/scientific-review-20260924T164425Z/agents/adjudication.md). Scientific acceptance and any change of prescription remain with the researcher.
