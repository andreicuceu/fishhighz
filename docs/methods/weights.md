# Forest weighting and sampling noise

`fishhighz.weights.prepare_forest_weights(field, geometry, response, *,
z_source, magnitudes, quadrature, rho, variance, length_velocity, method,
weights=None, iterations=None, signal=None, alias=None, auxiliary=None,
rtol=1e-4, min_updates=3, stable_steps=3, max_updates=96)` prepares
one immutable forest sample per field/bin. `field` is an `ObservedField`,
`geometry` a `BinGeometry`, and `response` the same `InstrumentResponse` used for
signal and final noise. Its full pixel width must be positive; Gaussian sigma
may be zero. Source redshift must exceed `geometry.z_eval`: densities/variances
are sampled at the source, while clustering, P1D and distance conversions use
`z_eval`. Field identity, redshifts, h convention, geometry conversion factors,
response and input arrays are retained. A fixed-count `iterations` value and an
adaptive outcome/counts/residuals are retained, but adaptive input controls must
be preserved separately in caller provenance. Noise evaluation rejects a
different field, evaluation geometry or response. Changing area alone does not
change this local noise preparation.

All magnitude arrays have shape `(n_magnitude,)`, including a one-node sample.
Magnitudes are strictly increasing, quadrature weights are explicitly positive
in magnitudes, and `rho` is already selected/normalized in
`deg^-2 (km/s)^-1 mag^-1`. `variance` is dimensionless delta-flux pixel **variance**,
not RMS. Forest length and response widths are km/s. No file counts, target
density, magnitude spacing, survey area, or source selection are inferred.
`density_per_velocity(dndzdm, *, z_source)` converts a normalized
`dN/(dz_source dm deg^2)` row by `(1+z_source)/299792.458` without normalization.
Legacy rectangular sums require constant `dm` for **every** node, including
both endpoints; their sum is not the node-span. Nonuniform positive quadrature
weights are equally valid.

Choose explicitly:

- `method="supplied", weights=w`: same-shaped nonnegative weights with positive
  weighted support. No auxiliary model calls or legacy settings are used.
  Uniform rescaling leaves final noise unchanged.
- `method="legacy", iterations=3, signal=S, alias=B`: initialize
  `w=B/(B+Delta_v*variance)` and apply exactly three simultaneous updates
  `w_new=S/(S+variance/(cumsum(rho*quadrature*w)*L_v/Delta_v))`.
  `iterations=0` returns initialization; any explicit nonnegative integer is
  supported. Zero-density samples get zero weight; positive-density samples
  with zero variance get one. This is the legacy cumulative recipe, not a
  convergence algorithm or a claim of global optimality. `S` and `B` must be
  strictly positive finite representable scalars.
- `method="early_lyaforecast"`: use the full-sample `sum_historical`
  recurrence. Its update signal is `S + B/(I1*L_v)`, with `I1` recomputed from
  the current weights. This is the recommended research method. With
  `iterations=None`, the adaptive solver requires at least three updates, three
  stable transitions, doubled-count confirmation, positive finite `rtol`
  (default `1e-4`) and a finite cap (default 96). Failure raises and no last
  iterate is substituted. An explicit nonnegative `iterations` value requests
  a labeled fixed-count diagnostic instead. The result retains the convergence
  outcome, counts and residuals; callers must retain nondefault input controls
  in their own provenance.
- `method="mcdonald"`: use the optional full-sample `sum_aliasing` recurrence,
  with update signal `S + B*I2/(I1**2*L_v)`. Adaptive and explicit fixed-count
  behavior is identical to `early_lyaforecast`, including explicit failure.
  This retained alternative is not the selected research baseline.
- `method="inverse_variance", alias=B_star`: compute
  `w=B_star/(B_star+Delta_v*variance)` once on positive-density nodes; store
  zero on zero-density nodes and exactly one for zero variance on support.
  Supply a strictly positive finite scalar `B_star` in km/s, already including
  the P1D response squared at your chosen fixed weighting mode. There is no
  default mode, automatic P1D sampling or 3D signal dependency. `weights`,
  `iterations` (even zero), `signal` and `auxiliary` must all be `None`.
  The returned method and alias record this choice, with no iteration history.
  It equals the legacy zero-iteration seed algebraically. For a nonnegative
  measure, common forest length/pixel width and independent sightline noise,
  it minimizes `Q=A*B_star+P_pixel` under the diagonal-covariance approximation
  at that weighting mode; this does not establish a multi-mode BAO optimum.

For the explicit inverse-variance option, pass these keywords directly in
`ForestInput.weight_options` and omit `auxiliary_coordinates`. Final noise
still uses its own mode-dependent **intrinsic** P1D and applies that mode's
response squared once to aliasing only. `B_star` does not replace this P1D;
pixel noise remains unsmoothed. Weights and coefficients, as well as final
fiducial noise, remain fixed during mean-model differentiation. Legacy
compatibility and all profile defaults are unchanged.

For actual callable integration, `sample_auxiliary(field, geometry, response,
prepared_p3d, theta_p3d, *, p1d_model, p1d_parameters, theta_p1d, k_t_deg,
k_p_velocity)` returns immutable samples usable as `auxiliary=sample` in place
of `signal/alias`. It routes only the matching forest auto through the original
field indices and binding, and calls independently bound P1D once. It computes
`k_parallel=a_v*k_p_velocity`, `k_transverse=k_t_deg/d_deg`,
`S=P3D_auto*W^2*a_v/d_deg^2` (deg² km/s) and `B=P1D*W^2` (km/s).
`S` is the caller's full-auto P3D prediction after the field response and
angular/velocity unit conversion shown above; it is not an unsmoothed intrinsic
scalar. A legacy method name containing “smooth” does not require a smooth-only
external provider. Explicit historical
auxiliary examples are `(2.4, 0.00035)` for BAO and `(7, 0.001)` for P1D, in
`deg^-1` and `s/km`. There is no default. Auxiliary modes may lie outside
forecast cuts, add no modes, and must remain inside provider/template domains.
The retained fiducial states have no provider object or mutable model cache.

The result exposes `weights`, `I1/I2/I3` prefix integrals, `A` in deg² and
`P_pixel` in deg² km/s, plus maximum absolute successive `weight_changes`.
For masses `r=rho*quadrature`, these are prefixes of `r*w`, `r*w²`, and
`r*w²*variance`, with `A=I2[-1]/(I1[-1]²*L_v)` and
`P_pixel=I3[-1]*Delta_v/(I1[-1]²*L_v)`. Zero prefixes are allowed without
optional divisions; unrepresentable arithmetic and empty weighted support fail.
Individual integral-product underflows are rejected before prefix accumulation,
even when another positive term would leave a nonzero total. Exact-zero products
remain valid. Extreme weight normalizations can therefore be rejected while an
equivalent representable normalization succeeds; exposed inputs are never rescaled.

`fishhighz.noise.forest_noise(weights, field, geometry, response, k, mu, p1d)`
accepts paired observed `(n_node,)` arrays and intrinsic nonnegative P1D samples
in km/s at **q=k*mu/a_v**. Obtain these through `evaluate_p1d` with its own
callable/binding/state; replacing P3D never chooses P1D. The immutable result has
`aliasing=A*P1D*W(q)^2*d_deg^2/a_v`,
`pixel=P_pixel*d_deg^2/a_v`, and `total`, all `(Mpc/h_fid)^3`.
Pixel noise is unsmoothed; signed sinc lobes are squared for aliasing.
`galaxy_noise(n_bar)` gives `1/n_bar` for explicit positive comoving number
density `(h_fid/Mpc)^3`. `local_galaxy_density(dndzdm, quadrature, geometry)`
uses `sum(quadrature*dndzdm)*(1+z_eval)/c*a_v/d_deg²`; this is a local-density
approximation, not a count integrated over the bin. Supply the chosen `n_bar`
directly for a count/volume convention.

`prepare_noise(selection, n_node, *, diagonal=None,
independent_sampling=None, full=None)` returns immutable
`(n_node, n_required)` noise in `selection.required_pairs` order. Generated
noise requires `diagonal={active_field_id: node_array, ...}` with exact active
coverage and explicit `independent_sampling=True`. Unused fields need no noise
or model calls. Full packed `full=N` is a **replacement**, not an addition;
combining it with generated settings is rejected. Signed crosses, zero and
singular matrices are supported, with the existing normalized PSD tolerance
checked independently of signal. No overlap is inferred from physical labels.

Prepare each field's weights/noise once and share them across pairs. Add noise
once using `combine_observed_power(W_i*W_j*P_ij, N)` for covariance; the default
mean excludes known noise. Freeze weights, response, P1D/noise, modes and factors
while differentiating intrinsic P3D. A different P1D with fixed supplied weights
changes aliasing only; explicitly regenerating legacy weights can also change
pixel noise. No automatic P1D/noise or covariance derivatives are added.

Run `python examples/weighted_noise.py` for two distinct forests, a galaxy,
negative signal crosses, explicit independent sampling, an independent matrix
covariance/amplitude Fisher oracle and fixed-state checks. Inputs are synthetic;
this is not a physical DESI-2 forecast. Raw readers and overlap-derived noise
remain outside these array APIs. The cumulative formulas retain
lyaforecast GPLv3/McDonald & Eisenstein (2007) provenance in the source.

