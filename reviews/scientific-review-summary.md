# FishHighz scientific review — summary

The review finds the tested FishHighz covariance, Fisher calculation and adopted full-sample forest weighting internally consistent. The important scientific qualifications concern the information attributed to BAO damping/reconstruction and the physical interpretation of high-k full-shape constraints. No new consequential FishHighz production bug was confirmed in the reviewed paths. The known legacy prefix-sum bug was independently reproduced and is avoided by the native forecast.

This assessment includes the current uncommitted full-shape implementation. It does not adopt a revised prescription. The [detailed report](scientific-review-report.md) contains equations, source references, numerical evidence, assumptions and proposed dispositions.

## Main conclusions

1. **The weight normalization correction is justified.** The papers use full-population moments, not a magnitude-dependent prefix in the iterative update. The native path uses the full moments and retains the correct final aliasing and pixel-noise coefficients. The historical `early_lyaforecast` update is a separate choice: it is not identical to the published aliasing-inclusive update, nor demonstrated to be globally optimal.
2. **The Gaussian covariance and information algebra pass independent checks.** Selected auto/cross spectra retain their mutual covariance, mode counts have the correct normalization, negative cross powers are preserved, and full-shape nuisance parameters are marginalized before reporting targets.
3. **The current BAO derivative is model-dependent.** It differentiates the AP prefactor, Kaiser angle and damping envelope as well as the wiggle scale. Font-Ribera et al. deliberately exclude damping/RSD changes from their isolated-BAO derivatives. FishHighz correctly differentiates its own declared model, but the results should not be described as an exact implementation of that paper's BAO information prescription.
4. **Mixed-pair reconstruction needs a physical qualification.** The mean squared auto-width rule is mathematically explicit; its treatment of a reconstructed galaxy field crossed with an unreconstructed forest is not independently validated by these papers. In the reconstructed historical attribution, this choice produces the largest isolated BAO error reduction.
5. **High-k full-shape forecasts remain conditional.** The smooth power carries undamped shape/RSD information, with a linear template and limited nuisance freedom. Numerically stable results at kmax=0.3–0.5 h/Mpc are not established DESI-2 performance predictions. With free biases, σ8_fid can simply normalize the identifiable fσ8 amplitude; its interpretation still depends on the fixed spectrum shape and physical/nuisance assumptions.

These distinctions follow from [McDonald & Eisenstein](https://arxiv.org/pdf/astro-ph/0607122), [McQuinn & White](https://arxiv.org/pdf/1102.1752), and [Font-Ribera et al.](https://arxiv.org/pdf/1308.4164), as specified equation by equation in the detailed report.

## Scientifically meaningful changes from lyaforecast

| Change | Assessment |
|---|---|
| Full-sample weight normalization, per-forest state and confirmed iterative stopping | **Justified** correction/implementation; finite convergence is not optimality. |
| Historical alias-containing weighting signal | **Conditionally justified**; distinct from the two papers' other weighting prescriptions. |
| CAMB growth, consistent evaluation redshift, integrated volume and physical resolving-power conversion | **Justified** within the selected cosmology, effective-redshift and instrument conventions. |
| Explicit signed template decomposition and revised interpolation | **Conditionally justified**; changes the feature definition and template source. |
| Composite Fourier/magnitude quadrature | **Justified** by the tested numerical refinements. |
| Negative-density flooring in the compatibility input adapter | **Conditionally justified**; changes spline treatment, not a direct empirical selection measurement. |
| Product of individual field responses rather than one response selected by pair | **Justified**, including unequal-response controls. |
| Wiggle-damped covariance and full inverse-AP BAO derivatives | Internally consistent; **conditionally justified** information prescription, differing from the papers' conservative treatment. |
| Mean squared cross damping widths | **Conditionally justified** as a product-propagator assumption; physical accuracy for mixed reconstruction remains unresolved. |
| Exact conversion constant, explicit bin-1 selection and reproducible native preparation | **Justified**; selection must be matched when comparing forecasts. |
| Independent smooth/wiggle dilations, growth/nuisance inference and category cuts | **New capabilities without a direct legacy counterpart**; algebra checked, interpretation conditional on the model. |

The detailed comparison separates these changes from inherited P1D, density-normalization, source-redshift, raw-cell-width, SNR-policy and independent-noise assumptions. It also distinguishes correcting the prefix bug from simultaneously changing the weighting signal and stopping rule.

## Numerical evidence

Fresh selected legacy/full-compatibility Fisher matrices agree to 1.4e-15. Both compatibility profiles reproduce saved S2 matrices to 6.1e-16. Fresh native accuracy agrees with the historical accuracy joint Fisher matrices to 4.5e-9, and the full six-bin numerical comparison is recorded separately from historical attribution.

The current native six-bin DESI-2 BAO accuracy forecast is numerically continuous with the saved S2 accuracy reference: maximum relative Fisher differences are 4.47e-9 jointly and 8.80e-9 for own-covariance spectra. The public six-bin interface reproduces the independent contraction to 3.96e-16 in Fisher matrices and 4.44e-16 in BAO errors. All 78 individual and six joint constraints are rank two; the prescribed 12 bin-1 spectra remain excluded. Direct Wick reconstruction and contraction pass at 1.2e-16 and 1.6e-14 respectively. Bin-2 area ×1.2 scales Fisher information by 1.2 and errors by 1/sqrt(1.2), while signal and noise remain identical.

Six independent numerical refinements have maximum BAO-error changes at or below 0.000230% across all joint and individual constraints, well under the 0.1% criterion. The strongest tested physical sensitivities are observed k_max and source densities: lowering k_max from 0.5 to 0.1 h/Mpc increases joint errors by 137–147% in bins 2 and 6; a coherent 20% density reduction increases joint errors by 9.7–15.6% in bins 1, 2 and 6. A 20% LBG forest-variance change produces roughly 0.8–2.7% joint changes; a 20% LBG pixel-width increase changes joint errors by at most 0.00189% with native SNR conversion. Restoring forest widths for mixed forest–galaxy wiggles changes bin-2 and bin-6 joint errors by (+4.13%, +2.29%) and (+3.11%, +1.38%) in (ap, at). Freezing the damping or RSD-angle AP derivative individually changes them by at most 0.47% and 0.35%, respectively.

The AP operator controls are explicitly hybrid local diagnostics: they use current native Jacobian, total power and covariance with saved S2 fiducial wiggle operands. Current/S2 Jacobian and covariance relative L2 differences are at most 2.40e-9 and 1.74e-9; a five-point finite difference checks the isolated analytic derivative to 5.1e-13. The mixed-width control updates fiducial total power and covariance but transforms the Jacobian analytically to first order. These results quantify prescription sensitivity and finite-setting stability; they do not validate reconstruction physics or nonlinear broadband information. Detailed tables, exact case operands, input hashes, status and Fisher arrays are in [numerical-section.md](scientific-review-report.md#current-native-accuracy-and-physical-variations) and [summary.json](../.validation/scientific-review-20260924T164425Z/results/numerical/summary.json).


The 364 saved attribution records were independently reconstructed without rerunning their forecasts. Historically, changing mixed-pair damping alone reduced joint radial BAO errors by 2.38–4.11%, partly opposed by the grouped AP derivative change, which increased radial errors by 1.02–1.17%. These conditional effects are not additive. Current derivative controls are reported separately.

Independent checks also covered 36 saved full-shape joint constraints and the subsequent current-source LBG+LAE calculation. A weak individual constraint has a conditioning-limited inverse reconstruction; it is explicitly qualified. Real unequal-category-cut forecasts were not rerun. The scientific test selections passed **557 tests**.

## Bugs and assumptions

The confirmed **legacy bug** is the use of prefix I1 values in a full-population weight recurrence. A fixed-signal two-node reconstruction changes the first weight by 19.05%; this is a mechanism demonstration, not the DESI-2 error shift. FishHighz retains that path only for labeled historical reproduction. No new consequential production defect was confirmed by this review.

The most relevant remaining assumptions are Gaussian independent Fourier cells and redshift bins; a common survey volume; independent sampling noise; central absorber/source approximations; prescribed P1D and density/SNR policies; fixed fiducial weights/noise/response; deterministic linear bias/RSD and fixed template shape; prescribed reconstruction/damping; and limited nuisance freedom. Continuum projection, realistic window coupling, nonlinear covariance/bias, redshift errors and contaminants are not comprehensively modeled. These omissions matter especially for broadband information.

Production code and defaults were unchanged. The deferred reference-mode and alternative-weighting studies remain deferred; no Slurm actions, commits or new full-shape forecasts were performed. Supporting [diagnostics and evidence](../.validation/scientific-review-20260924T164425Z/README.md) are linked from the report.

**Stopped for scientific review.** Numerical agreement does not constitute scientific acceptance or authorization to change a prescription.
