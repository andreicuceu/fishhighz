# FishHighz research baseline

This guide defines the recommended research prescription selected on 2026-09-21
from the S2--S4 calculation. Its recipe identity is
`early-lyaforecast-2026-09-18`, profile `accuracy`, method
`early_lyaforecast`, representative mode `(2.4 deg^-1, 0.00035 s/km)`, and
selection `bin1-qso-only; bins2-6-all`. Record this full identity with derived
results. Changing one of these choices creates a different recipe.

The prescription is a selected, numerically qualified set of forecast
assumptions. S4 found that the mixed forest--galaxy damping prescription is the
largest isolated contribution to its smaller joint BAO errors relative to fixed
compatibility. That attribution does not independently validate the
reconstruction physics or make smaller formal errors evidence of greater physical
accuracy.

## Observable, signal, and response

Each redshift bin is evaluated at its geometric centre

```text
1 + z_eval = sqrt[(1 + z_min)(1 + z_max)].
```

The volume is the supplied-area integral of `c D_M(z)^2/H(z)` over the bin. A
common volume is used for every field in a bin; bins are independent. The
Fourier measure is Gauss--Legendre quadrature within fixed observed-coordinate
`k` cuts and `mu in [0,1]`, with

```text
q_mode = k^2 w_k w_mu / (2 pi^2),    N_mode = V q_mode.
```

Here `k` is in `h_fid/Mpc`, power and volume are in `(Mpc/h_fid)^3`, and mode
counts are dimensionless. The velocity and angular conversion factors are

```text
a_v   = H(z_eval) / [h_fid (1 + z_eval)]       [(km/s)/(Mpc/h_fid)],
d_deg = D_M(z_eval) h_fid pi/180               [(Mpc/h_fid)/deg].
```

Quadrature nodes are not bandpower centres.

The built-in signal uses the supplied Vega decomposition
`P_nw = PKSB`, `P_w = PK - PKSB`. Signed components are preserved. Forest
Kaiser factors are `b_i(1 + beta_i mu^2)` and galaxy/QSO factors are
`b_i + f mu^2`. The template power growth is
`G=[sigma8(z_eval)/sigma8(z_ref)]^2`, prepared from the fiducial background at
the exact evaluation and template-reference redshifts. Instrument response is
owned by survey preparation and applied once to the observed mean:

```text
P_obs,ij = W_i W_j P_intrinsic,ij.
```

For resolving power `R`, the Gaussian velocity width is the physical FWHM
conversion

```text
sigma_v = c / [2 sqrt(2 ln 2) R].
```

Pixel top-hat and Gaussian response factors are per field. Noise is separate and
is not multiplied by this response.

For the real S2--S4 DESI-2 baseline, `f` is the exact CAMB fiducial growth rate
at `z_eval` and the forest beta is 1.45. For each field,

```text
Sigma_perp,i = 3.26 sigma8(z_eval)
               / [sigma8_damping_reference sqrt(r_i)]    [Mpc/h_fid],
Sigma_parallel,i = (1 + f) Sigma_perp,i,
```

with reconstruction factor `r_i=2` for galaxy/QSO fields and `r_i=1` for
forests. These are inputs of that validated recipe; the synthetic research
example uses explicit illustrative biases, beta, growth rate, and widths instead.

## Wiggle damping and AP derivatives

Only the wiggle component is damped. Each field has fixed parallel and transverse
auto widths. For a pair `ij`, on either axis,

```text
Sigma_ij^2 = (Sigma_i^2 + Sigma_j^2) / 2,
D_ij = exp[-(k_parallel^2 Sigma_parallel,ij^2
              + k_perp^2 Sigma_perp,ij^2) / 2].
```

This per-field crossed squared-width mean is retained for mixed forest--galaxy
pairs. Widths are fixed inputs, outside the parameter vector.

The validated BAO calculation assigns two parameters to each bin, `ap` and
`at`, and applies inverse AP mapping to the complete wiggle contribution:

```text
k_parallel' = k mu / ap
k_perp'     = k sqrt(1-mu^2) / at
Q           = 1 / (ap at^2).
```

The derivative therefore includes the AP prefactor, remapped template,
remapped Kaiser factors, and remapped fixed-width damping. The smooth component
is unscaled in this BAO-only calculation. The observed response, forest weights,
noise, integration grid, volume, and covariance stay fixed at their fiducial
values. FishHighz differentiates the selected mean only and adds no
covariance-derivative Fisher term.

The S2--S4 validation fit contains exactly the two AP parameters per bin. It adds
no broadband parameters, nuisance parameters, or priors, and does not marginalize
over an unspecified nuisance model. The public API is broader: callers may put
explicit nuisance parameters in `ParameterRegistry`, map them through
`BoundParameters`, and supply a prior precision matrix. Parameter roles and
bounds never create a prior or silently select a marginalized set.

## Forest weights and noise

For each forest field, supply ordered magnitude samples, positive quadrature
weights `q_i` in magnitudes, density
`rho_i` in `deg^-2 (km/s)^-1 mag^-1`, dimensionless pixel variance `v_i`,
forest velocity length `L` in km/s, and pixel width `l_pixel` in km/s. At the
representative mode, convert the angular/velocity components to the intrinsic
P3D query and evaluate the profile's own fiducial auto signal `P`:

```text
k_parallel = k_parallel_velocity a_v,
k_perp     = k_transverse_deg / d_deg,
P = P3D_intrinsic(k,mu) W_field(k_parallel_velocity)^2 a_v/d_deg^2
                                                        [deg^2 km/s],
B = P1D_intrinsic(k_parallel_velocity)
    W_field(k_parallel_velocity)^2                      [km/s].
```

The early-lyaforecast full-sample recurrence begins with

```text
w_i^(0) = B / (B + l_pixel v_i),
I1 = sum_i rho_i q_i w_i,
w_i^(n+1) = [P + B/(L I1)]
            / [P + B/(L I1) + l_pixel v_i/(L I1)].
```

Do not renormalize the weights between nonlinear updates. The public
`prepare_forest_weights(..., method="early_lyaforecast")` solver requires at
least three updates, three stable transitions, `rtol=1e-4`, doubled-count
confirmation, and a maximum of 96 updates by default. Failure is explicit; a
capped or otherwise unconfirmed last iterate is never substituted. S2 tightened
the final validation state to `rtol=1e-5` as a numerical refinement.

The returned coefficients use full-sample moments

```text
I2 = sum_i rho_i q_i w_i^2,
I3 = sum_i rho_i q_i v_i w_i^2,
A       = I2 / (L I1^2)                 [deg^2],
P_pixel = l_pixel I3 / (L I1^2)         [deg^2 km/s].
```

At a forecast mode `q=k mu/a_v` in `s/km`, physical forest noise is

```text
N_F(q) = [A P1D_intrinsic(q) W_field(q)^2 + P_pixel] d_deg^2/a_v.
```

Thus the representative `B` selects the weights; the final noise still uses its
own mode-dependent intrinsic P1D and response. Each forest has its own response,
P1D, density, variance and converged state. Galaxy auto noise is `1/n_bar`.
The selected baseline assumes independent sampling noise, with no inferred
forest--galaxy or inter-sample shot noise.

Strict public preparation accepts finite nonnegative density, variance and
quadrature inputs and rejects invalid domains and nonrepresentable arithmetic.
It performs no density floor, clamp, smoothing, renormalization, or sorting.
The `floor_negative` density policy used by the accuracy validation lives in
`fishhighz.adapters.legacy_compat`; that adapter either rejects or floors a
negative interpolant. Literal signed-density recurrence arithmetic is separate
and remains confined to the validation compatibility paths
`fishhighz.validation.compatibility_weights` and
`fishhighz.validation.revised_compatibility`. These validation-only behaviors do
not alter strict defaults.

## Covariance, selection, and returned metadata

For required field pairs,

```text
T_ij = W_i W_j P_ij + N_ij,
Cov[P_ij,P_mn] = (T_im T_jn + T_in T_jm) / N_mode.
```

The covariance keeps correlations among selected spectra at each Fourier cell.
It couples neither Fourier cells nor redshift bins. The joint Fisher matrix uses
this selected covariance directly; it is not the sum of individual-spectrum
Fisher matrices. In bin 1, future forecasts retain only the Ly-alpha(QSO) auto,
QSO auto, and their cross, including their full three-spectrum covariance.
Correlations involving LBG, LAE, or Ly-alpha(LBG) are excluded before
preparation. Bins 2--6 retain all 15 spectra.

`prepare_bin` returns immutable fiducial arrays for nodes, modes, responses,
pair-response products, signal, noise, total power and Cholesky factors. Its
diagnostics include the bin and field/pair identities, volume, node and pair
counts, response settings, noise convention, model/P1D call counts, galaxy
densities, forest provenance, and units. Each `ForestWeights` record includes
the method, inputs, weights, moments, `A`, `P_pixel`, and convergence record.
`run_bin` and `run_forecast` return the actual derivative-column and call
diagnostics plus `FisherResult`, whose data, prior, total information, rank,
condition, eigenvalues, scales and null directions remain explicit.

The compact [research example](examples/research_bao_forecast.py) constructs a
synthetic forest--galaxy auto/cross joint covariance through these public
methods. Install `fishhighz[templates]`, run it as a script, or import and call
`run()`. It is a self-contained API demonstration, not the DESI-2 validation.

## Profiles and numerical evidence

The three current validation profiles share recipe revision
`early-lyaforecast-2026-09-18` and the stated bin selection:

| Profile | Method | Scientific role |
| --- | --- | --- |
| `full-compatibility` | `legacy` | Literal captured estimator and inputs for historical reproduction. |
| `fixed-compatibility` | `early_lyaforecast` | Compatibility estimator with the selected full-sample weights, isolating weighting changes. |
| `accuracy` | `early_lyaforecast` | Recommended research prescription documented above. |

The legacy alias `compatibility` resolves to `full-compatibility` only in the
validation runner. These named profiles are validation recipes tied to captured
resources and evidence contracts; use explicit public methods to construct a new
forecast.

S2 found all 11 required forest preparations converged and all six accuracy bins
passed the declared finite refinements. The largest individual and joint BAO
error changes were 0.0011411% and 0.0012685%. This combines tested magnitude,
Fourier, volume, derivative-step, and stopping refinements; iteration convergence
alone would not establish finite-grid convergence or a continuum theorem.

The final S2/S3 evidence state used observed `k=[0.01,0.5] h_fid/Mpc`, 128
`k` intervals with Gauss--Legendre order 4 per interval, `mu` order 32, redshift
volume order 32, composite magnitude order 16 on the recorded partition, AP
finite-difference step `2.5e-4`, and final weight `rtol=1e-5`. The tighter
stopping tolerance is the retained qualified evidence state. The named recipe
identity and public adaptive solver default remain `rtol=1e-4`; callers must
record an override when reproducing the final S2/S3 state. The synthetic research
example uses smaller illustrative inputs and does not inherit the S2 numerical
qualification.

S3 found accuracy/fixed joint radial error changes of -2.6693% to -3.8610% and
transverse changes of -0.6665% to -1.6089%. S4 reconstructed the endpoints to
0.000129% in componentwise errors and found the mixed-pair damping change to be
the dominant isolated reduction, partially opposed by the full AP derivative
change. Its largest tested refinement change was 0.010783%. These results support
the stated numerical qualification and attribution within this model. They do
not constitute independent acceptance of reconstruction physics.

Reference-mode sensitivity is deferred to D1. D2 is the optional comparison of
the selected early-lyaforecast recurrence with the W12 inverse-variance method.
Future mode choices are not part of this baseline. General INI/CLI translation
and result serialization also remain outside the current research-ready Python
milestone.
