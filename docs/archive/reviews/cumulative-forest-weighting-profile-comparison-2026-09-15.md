# Cumulative forest-weight convergence

The prescription updates each magnitude sample using the **cumulative weighted density of brighter samples**:
\[
N_i^{(t)}=\frac{L}{\Delta v}\sum_{j\le i}r_jw_j^{(t)},\qquad
w_i^{(t+1)}=\frac{N_i^{(t)}}{N_i^{(t)}+v_i/S},
\]
where \(r_i=\rho_iq_i\), \(v_i\) is pixel variance, and \(S\) is the auxiliary forest auto-power.

Three distinctions matter:

- **Weight amplitude:** on every saved grid, the weights tend to zero. This explains the underflows.
- **Noise coefficients:** a common weight amplitude cancels from
  \[
  A=\frac{\sum r_iw_i^2}{L(\sum r_iw_i)^2},\qquad
  P_{\rm pixel}=\frac{\Delta v\sum r_iw_i^2v_i}{L(\sum r_iw_i)^2}.
  \]
  Therefore, vanishing weights do **not** establish divergent noise.
- **Scientific convergence:** the relative weight shape continues evolving, and its fixed-grid limit changes substantially with magnitude refinement. Holding the initial weight function fixed makes the tested magnitude integral converge; recomputing the cumulative weights restores the sensitivity.

In the saved 15×2pt study, iteration changes produced up to **5.0% changes in joint errors and 38.5% in individual-spectrum errors**, exceeding the \(0.1\%\) requirement. Step 13’s logarithmic representation addresses arithmetic range, but its revision-2 conclusion remains **unresolved continuum convergence**, pending independent review. Weight diagnosis (local-only path: `/global/cfs/cdirs/desicollab/users/acuceu/vega_dev/lib/fishhighz/reviews/step-12-r5.md:99`), Step 13 assessment (local-only path: `/global/cfs/cdirs/desicollab/users/acuceu/vega_dev/lib/fishhighz/reviews/step-13-weight-limit-r2.md:3`).

## Compatibility → accuracy: full scientific comparison

**Compatibility reproduces a specified finite legacy calculation. “Accuracy” names the profile being tested for convergence; it does not certify that convergence.**

In the last column:

- **Direct:** changes the cumulative iteration or its magnitude discretization.
- **Coupled:** changes inputs to the forest weighting/noise calculation.
- **Independent:** can be defined and tested without resolving that iteration. Its *forecast impact* can nevertheless depend on the unconverged forest covariance.

| Quantity or operation | Compatibility | Accuracy | Relation to weighting issue |
|---|---|---|---|
| **Underlying matter spectrum** | Upstream CAMB spectrum and its interpolation, referenced at \(z=2.3\). | Supplied Vega-format template; existing `PKSB` and signed `PK−PKSB`, with template interpolation. | **Coupled:** changes auxiliary \(S\), as well as forecast signal. |
| **Power-spectrum growth** | EdS scaling, proportional to \((1+z)^{-2}\). | CAMB \([\sigma_8(z)/\sigma_8(z_{\rm template})]^2\). | **Coupled:** changes \(S\). |
| **Evaluation redshifts** | Mean at arithmetic bin centre; covariance/noise at geometric centre. | Mean and covariance/noise consistently at geometric centre. CAMB is prepared at actual evaluation redshifts; galaxy \(f\) is evaluated there. Bias prescriptions are retained. | Mean-centre correction is **independent**; revised background/model values can change \(S\). |
| **BAO decomposition** | Degree-8 polynomial fit to log absolute observed power, retaining sign; \(10^{-12}\) regularizer and strongly weighted first three samples. | Uses the template’s supplied smooth/wiggle decomposition directly. | Derivative extraction is **independent**; the template model also enters \(S\). |
| **BAO damping in fiducial power** | Damping is applied to the extracted peak for the derivative; covariance uses the upstream undamped linear signal. | Wiggle damping is part of the intrinsic model used consistently for signal, covariance and derivatives; smooth power is undamped. | **Coupled:** forest auto-power also supplies \(S\). |
| **Cross-spectrum damping/reconstruction** | Any spectrum containing a forest uses unreconstructed forest damping; galaxy-only spectra use reconstruction factor 2. | Per-field widths: forest factor 1, galaxy factor 2; cross-width squares are the mean of the two auto-width squares. | **Independent:** changing forest–galaxy cross widths does not change the forest auto used for weighting. |
| **BAO parameter derivative** | Backward \(dP_{\rm peak}/d\ln k\), multiplied by \(\mu^2\) or \(1-\mu^2\); first node is zero. Instrument response is already inside the differentiated observed peak. | Numerical derivatives of the full wiggle-only inverse AP mapping, including \(Q=(a_\parallel a_\perp^2)^{-1}\), remapped RSD and damping. Smooth scaling is identity; observed-coordinate response stays fixed. | **Independent; derivative-step checks pass.** |
| **Fourier integration** | 500 endpoint-inclusive uniformly spaced \(k\) samples; ten midpoint \(\mu\) samples. | Composite Gauss–Legendre in \(k\), Gauss–Legendre in \(\mu\). Saved final runs use 128 \(k\) intervals × order 4, and \(\mu\) order 32. | **Independent; checks pass.** |
| **Survey volume** | Central conversion factors multiplied by logarithmic velocity depth and area. | Integrates \(dV/dz\) across the same bin; saved final volume quadrature order 32. | **Independent; checks pass.** |
| **Instrument resolution** | Treats \(c/R\) as Gaussian sigma. | Treats \(c/R\) as FWHM: \(\sigma_v=c/[R\sqrt{8\ln2}]\). | **Coupled:** changes auxiliary \(S\), initial-weight P1D term, aliasing and signal response. |
| **Response and weight ownership** | Pair-specific upstream preparation, including pair-specific forest inputs. | One response and one weight/noise preparation per observed field; forest weights use that field’s **auto**-P3D. Cross responses are \(W_iW_j\). | **Coupled:** establishes a consistent forest-noise prescription across spectra. |
| **Magnitude integration** | 107 endpoint-inclusive rectangular samples over \(16.1\le m\le26.75\). | Composite Gauss–Legendre on a fixed partition containing spline, support, SNR and density-zero boundaries. Orders 4–64 tested. | **Direct for forests:** refinement changes the prefix update. Galaxy-density integration does not have this iteration problem. |
| **Number of cumulative updates** | Exactly three. | Tests 3, 6, 12 and, when needed, 24. Latest saved 15×2pt primaries retain **six**, explicitly unconverged. | **Direct.** |
| **Negative interpolated density** | Literal upstream calculation retains negative spline overshoot. | Explicitly replaces negative interpolants with \(10^{-20}\); ordinary zero/positive in-domain values are retained. | **Coupled:** changes masses entering forest weights; also changes galaxy shot noise independently. |
| **Velocity conversion constant** | \(c=299800\ {\rm km\,s^{-1}}\). | \(c=299792.458\ {\rm km\,s^{-1}}\) in FishHighz conversions. | **Coupled, small:** enters density, pixel/forest lengths and response conversions. It is not the diagnosed cause of nonconvergence. |
| **Numerical admissibility and acceptance** | Reproduces captured arrays; audits field PSD separately and requires selected covariance to be positive definite. | Uses production positivity/PSD and arithmetic guards; explicitly tests Fisher/error stability under refinement, including individual spectra. | Guards **expose** underflow; stricter convergence tests reveal the scientific failure rather than causing it. |

These definitions are verified against the live accuracy recipe (local-only path: `/global/cfs/cdirs/desicollab/users/acuceu/vega_dev/lib/fishhighz/fishhighz/validation/accuracy.py`), compatibility calculation (local-only path: `/global/cfs/cdirs/desicollab/users/acuceu/vega_dev/lib/fishhighz/fishhighz/validation/compatibility.py`), and scientific comparison report (local-only path: `/global/cfs/cdirs/desicollab/users/acuceu/vega_dev/lib/fishhighz/reviews/step-12-r2.md:166`).

### Retained assumptions—not profile changes

Both retain the survey cases, fields, selected spectra, redshift bounds, area, magnitude bounds, exposure settings, intrinsic P1D prescription, auxiliary weighting coordinates, raw first-spacing density normalization, and legacy SNR clamps/sentinels. Pixel/Poisson noise remains unsmoothed; aliasing carries \(W^2\). Covariance and weights remain fixed during parameter differentiation.

The additional floor variations, bright/sentinel removal, and ±10% raw-cell-width variations are **sensitivity experiments**, not changes in the primary accuracy profile. They feed directly into forest weighting inputs, so their forest forecast effects also remain conditional on the finite iteration choice.

**Consequently, the issue does not invalidate the independently verified AP, damping, volume or Fourier-integration calculations. It prevents interpreting the resulting forest-containing forecast differences as converged uncertainty predictions.** Galaxy-only forecasts avoid this particular problem.

<oai-mem-citation>
<citation_entries>
MEMORY.md:109-117|note=[Located weighting conventions and historical validation context then checked live source]
</citation_entries>
<rollout_ids>
01a09b99-6021-77a2-9814-7b4c645c2fdb
</rollout_ids>
</oai-mem-citation>
