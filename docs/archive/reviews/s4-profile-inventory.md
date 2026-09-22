# S4 source inventory: fixed-compatibility versus revised accuracy

The endpoints change the power spectrum and BAO estimator as well as survey
integration and instrumental response. Both now use the same early-lyaforecast
full-sample recurrence. Consequently their remaining difference must not be
attributed to a change of weighting method. This inventory identifies the
physical and numerical dependencies; quantitative BAO effects belong to the
companion S4 controlled-comparison report. No forecasts were run for this
inventory. Source paths below are relative to this package, except explicit
`../lyaforecast/` paths. Line numbers refer to the inspected working tree.

| Quantity | Fixed-compatibility endpoint | Revised-accuracy endpoint | Attribution implication and source |
| --- | --- | --- | --- |
| Linear power shape and normalization | CAMB linear power sampled at 1000 logarithmic nodes and linearly interpolated, at its configured reference redshift | Supplied Vega FITS PK and PKSB, with wiggle PK−PKSB, evaluated through prepared template interpolation | This is an actual spectrum/decomposition change, not merely growth scaling. `../lyaforecast/lyaforecast/cosmoCAMB.py:59`, `power_spectrum.py:53`; `fishhighz/validation/accuracy.py:149`, `models/kaiser.py:276`. Separate spectrum-source and growth effects where feasible. |
| Power growth | EdS factor [(1+z_ref)/(1+z)]² | CAMB [sigma8(z)/sigma8(z_template)]² | `../lyaforecast/lyaforecast/power_spectrum.py:53`; `fishhighz/validation/accuracy.py:324`. Changes both covariance signal and auxiliary auto-P used to determine weights. |
| Mean versus covariance redshift | Mean/peak estimator at arithmetic bin centre; covariance, source queries, pixel widths and angular/velocity conversions at geometric centre sqrt[(1+z_lo)(1+z_hi)]−1 | Intrinsic mean, covariance, damping and exact CAMB quantities at geometric centre | `../lyaforecast/lyaforecast/forecast_new.py:215`, `:235`; `survey.py:66`; `covariance.py:136`; `fishhighz/validation/accuracy.py:313`, `:324`, `:417`. Source-query and geometry conventions already agree in their redshift choice: changing just the intrinsic mean is not a wholesale geometry-redshift change. |
| Galaxy RSD and forest biases | Galaxy beta=f(z)/b(z), with f linearly interpolated/extrapolated from the reference forecast's arithmetic-centre CAMB samples; forest beta=1.45; configured bias functions | Same bias prescriptions and forest beta; galaxies use common exact prepared CAMB f at geometric centre | `../lyaforecast/lyaforecast/analytic_biases.py:18`, `:85`; `fishhighz/validation/accuracy.py:317`, `:339`; `kernels/kaiser.py:37`. Neither endpoint assumes constant galaxy f(z_ref). Exact-f versus interpolated-f treatment remains distinct from EdS versus CAMB power growth. |
| Fiducial covariance damping | Linear full P has no BAO damping; damping is applied only to the extracted peak in the derivative | Wiggle damping also enters fiducial P and therefore covariance and auxiliary P | `../lyaforecast/lyaforecast/power_spectrum.py:156`; `fisher.py:155`; `fishhighz/models/kaiser.py:285`. Holding covariance fixed while changing derivatives does not measure the full damping effect. |
| BAO extraction and derivative | Degree-8 polynomial fitted to log(abs(observed P)+1e−12), first three nodes weighted by 1e8; backward k derivative of damped residual; first derivative node zero | Explicit template wiggle; symmetric numerical derivative of inverse AP mapping, AP prefactor, remapped RSD and remapped fixed-width damping; smooth component unscaled | `fishhighz/validation/numerics.py:138–182`; `models/kaiser.py:245–290`; `validation/accuracy.py:680`. These are different estimators. Equations below specify response and damping ownership. |
| Reconstruction and cross damping | Same configured reconstruction factor for every galaxy-only pair; any pair containing forest receives no reconstruction | Per-field galaxy reconstruction; forest unchanged; cross width² equals average of auto width² | `../lyaforecast/lyaforecast/fisher.py:155`, `:221`; `fishhighz/validation/accuracy.py:342`; `kernels/kaiser.py:46`. Forest–galaxy damping changes; galaxy autos and galaxy–galaxy pairs retain the reconstruction prescription, modulo redshift. Turning off all reconstruction is an optional sensitivity, not the endpoint difference. |
| Damping amplitude and growth rate | Sigma_perp=3.26 sigma8(z_arithmetic)/sigma8(z_ref)/sqrt(r_pair), Sigma_parallel=(1+f(z_arithmetic)) Sigma_perp | Same 3.26 reference prescription, evaluated at geometric centre, r per field | `../lyaforecast/lyaforecast/fisher.py:221–233`; `fishhighz/validation/accuracy.py:342–348`. Actual legacy forecasts pass bin index and use bin-dependent f; the code's index=None reference-f fallback is not the endpoint. |
| Resolving-power conversion | sigma_v=c/R | sigma_v=c/[2 sqrt(2 ln 2) R] | `../lyaforecast/lyaforecast/covariance.py:177`; `fishhighz/validation/accuracy.py:391`. Changes forest signal, aliasing B, and representative P/B and thus reconverged weights. Pixel top-hat remains the same prescription. |
| Response by pair versus field | Each pair uses a single selected forest pixel width for every forest factor, so forest–forest cross can use W_selected² | Product W_i W_j from independently specified field responses | `../lyaforecast/lyaforecast/covariance.py:146`, `power_spectrum.py:218`; `fishhighz/validation/accuracy.py:391`, `:689`. Covariance setup selects first continuous tracer, whereas weight setup can select the second; inspect captured pair inputs rather than infer ownership from labels. All six saved bins have identical pixel/resolution widths across their forest-associated pairs, so pair versus field response ownership is numerically equivalent here at common redshift/convention; this architectural distinction supplies no separate nonzero endpoint correction. |
| Magnitude measure | 107 equally spaced nodes, full rectangular endpoint weights | Composite Gauss–Legendre quadrature split at interpolation/policy boundaries; S2 final orders 16 | `fishhighz/validation/revised_compatibility.py:78`, `accuracy.py:189`, `:432`; `reviews/s2.md`. Includes galaxy density integration as well as forest moments. Recompute weights on each measure; keep support/policy fixed to isolate quadrature. |
| Negative density | Signed in-domain spline overshoot retained | Named extension replaces negative interpolants by 1e−20 | `../lyaforecast/lyaforecast/tracer.py:63`; `fishhighz/adapters/legacy_compat.py:110`. This change is distinct from quadrature and does not modify strict public defaults. |
| Density support/normalization | Survey-wide magnitude integration range; forest-specific magnitude cuts applied to raw table before spline preparation; target density normalization and z>2.15 convention where configured | Same retained conventions | `../lyaforecast/lyaforecast/tracer.py:127`; `fishhighz/validation/accuracy.py:283–299`, `:432`; `adapters/legacy_compat.py:83`. Do not invent a changed physical support merely from new quadrature nodes. Outside tabulated magnitude support both use the legacy 1e−20 floor. |
| Density cell widths | First raw redshift spacing used for all source cells | Explicitly retained legacy_first_spacing | `fishhighz/validation/accuracy.py:296`, `:661`. Unknown physical irregular-cell widths remain a shared assumption; ±10% width is a sensitivity, not an endpoint change. |
| SNR policy | Bright-end clamp; faint/redshift/wavelength out-of-domain large-noise sentinel; scaled SNR floor; exposure and pixel-width conversion | Retained per-population policy with explicit diagnostics | `../lyaforecast/lyaforecast/spectrograph.py` (`get_noise_rms`); `fishhighz/adapters/legacy_compat.py:153–252`, `validation/accuracy.py:479`. Removing bright/sentinel populations or varying floors is optional sensitivity. |
| Forest and galaxy noise | Fixed endpoint changes only forest-auto A and P_pixel relative to captured powers; cross total powers remain noise-free; galaxy auto has 1/n | Same forest physical noise form, per-field preparation; independent-sampling diagonal noise, galaxy 1/n | `fishhighz/validation/revised_compatibility.py:68–109`; `../lyaforecast/lyaforecast/covariance.py:335–418`; `fishhighz/validation/accuracy.py:566–598`. No new overlap or cross shot noise is inferred. Pair-specific legacy weight states on crosses do not directly add cross noise. |
| Volume | Area × c ln[(1+z_hi)/(1+z_lo)] × d_deg(z_geo)²/a_v(z_geo) | Integral of c D_M²/H over bin, with common survey area and h conversion | `../lyaforecast/lyaforecast/covariance.py:202`, `:237`; `fishhighz/geometry.py:163`. Isolated volume has sigma∝V^(−1/2); it must not silently change powers or Jacobians. |
| Fourier measure | Uniform k nodes with full endpoint rectangular weights; midpoint mu grid | Gauss–Legendre in fixed observed k range and mu∈[0,1] | `../lyaforecast/lyaforecast/power_spectrum.py:43`, `covariance.py:244`; `fishhighz/validation/accuracy.py:545`. Same observed cuts 0.01–0.5 h/Mpc and Gaussian mode-count normalization. Legacy peak fit/backward derivative also depends on its grid, so a node switch must state estimator treatment. |
| Conversion constant | c=299800 km/s | c=299792.458 km/s | `../lyaforecast/lyaforecast/cosmoCAMB.py:10`; `fishhighz/geometry.py:10`. Small actual conversion difference affecting velocity density, forest/pixel lengths, response and volume; identify explicitly or retain in a named conversion group. |
| Weight convergence tolerance | Confirmed adaptive rtol=1e−4 | S2 final rtol=1e−5 | `fishhighz/validation/profile_definitions.py:10`; `reviews/s2.md`. Same algorithm, numerical refinement only; S2 measured extremely small error changes under tightening. |

Saved S2 metadata confirms one unique (pixel width, resolution) tuple per bin across all forest-associated pairs: pixel widths 63.3298043347, 58.8849127896, 55.0233556274, 51.6373205842, 48.6440319735 and 45.9788792973 km/s, with resolution 119.92 km/s throughout. Thus response construction by pair versus field is a shared numerical behavior for this configuration, while FWHM conversion and mean-redshift response evaluation do change.

The reference S2 endpoints are immutable captured arrays for fixed-compatibility,
and a direct model calculation for accuracy. Current neighboring source explains
the prescriptions; any mismatch with captured metadata must be resolved in favor
of explicit saved operands for exact endpoint reconstruction.

## Derivative ownership

Write E_ij(k,mu) for the legacy peak extracted from the *already response-smoothed*
mean. Its implemented derivative is

    J_legacy = (mu², 1−mu²) D^-_lnk [ E_ij(k,mu) D_pair(k,mu) ],

where the backward operator is [X(k_n)−X(k_(n−1))]/(Delta k/k_n), with
first node zero, and

    D_pair = exp{−[Sigma_parallel² k² mu²
                  + Sigma_perp² k²(1−mu²)]/2}.

Thus response is embedded inside peak extraction and the subsequent k derivative.
It is not legitimate to replace this exactly by an unsmoothed peak derivative
multiplied by fixed response. Polynomial extraction itself need not commute with
response multiplication. Damping is inside this derivative as well.

Accuracy instead differentiates the wiggle contribution

    P^w_ij(alpha) = G Q K_i(mu') K_j(mu') P_w(k')
                   exp{−[Sigma_ij,parallel² k_parallel'²
                         + Sigma_ij,perp² k_perp'²]/2},
    k_parallel'=k mu/alpha_parallel,
    k_perp'=k sqrt(1−mu²)/alpha_perp,
    Q=1/(alpha_parallel alpha_perp²),
    J_accuracy = W_i(k,mu) W_j(k,mu) partial_alpha P^w_ij.

The widths are fixed inputs, but their damping function is differentiated through
the remapped coordinates. At alpha=1, partial_parallel ln D is
Sigma_parallel² k_parallel² and partial_perp ln D is
Sigma_perp² k_perp². The corresponding prefactor derivatives are −1 and −2;
mu' changes too. This is not the legacy mu²/(1−mu²) projection of a fixed-mu
k derivative. Instrumental response, survey weights and covariance are held fixed
in the accuracy derivative. No covariance-derivative Fisher information is added.
The global sign convention of a pure dilation derivative alone does not alter F;
the additional derivative terms and distinct extracted signal can alter F.

## Shared weighting and sampling assumptions

Both primary endpoints use full-sample moments
I1=sum rho_i q_i w_i, I2=sum rho_i q_i w_i²,
I3=sum rho_i q_i v_i w_i² and the early-lyaforecast recurrence

    w_i(new) = [P + B/(I1 L)] /
               [P + B/(I1 L) + v_i l_pixel/(I1 L)],
    A=I2/(I1² L),  P_pixel=l_pixel I3/(I1² L).

The common seed is B/(B+l_pixel v_i). Amplitude is retained in every nonlinear
update. Both use k_transverse=2.4 deg⁻¹ and k_parallel=0.00035 s/km;
P is the fiducial observed forest auto in angular/velocity units and
B=P1D_intrinsic W_field². Accuracy recomputes P and B from its own model and
response; it does not import compatibility weights. The physical forest noise is
[A P1D_intrinsic(q) W_field(q)²+P_pixel] d_deg²/a_v. Pixel and Poisson noise
are not instrument-smoothed. Sources: `fishhighz/kernels/full_sum_weights.py:15–74`,
`validation/profile_definitions.py:11`, `validation/accuracy.py:56`,
`validation/revised_compatibility.py:84`.

A fixed-comoving representative mode, McDonald recurrence, W12 inverse-variance
weights, floor variations, altered cell widths and population removals are
separately labeled sensitivities. They are not differences between these primary
endpoints. The P1D prescription and its low-k treatment, forest source-redshift
mapping and forest velocity-length definition are shared (apart from the small
c convention). Both retain independent bins, separate two-parameter BAO fits per
bin, selected-spectrum Wick covariance, bin-1 QSO-only selection, and all 15
spectra in bins 2–6. A joint result is not a sum of single-spectrum information.
