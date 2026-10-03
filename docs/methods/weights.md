# Forest weighting and sampling noise

The weighting prescriptions, stopping rule and the central and integrated source
modes are described in the [forest-weighting note](forest-weighting.md). This
page gives the call contracts.

`fishhighz.weights.prepare_forest_weights(field, geometry, response, *,
z_source, magnitudes, quadrature, rho, variance, length_velocity, method,
signal=None, alias=None, auxiliary=None, rtol=1e-4, min_updates=3,
stable_steps=3, max_updates=96)` prepares one immutable forest sample per
field/bin. `field` is an `ObservedField`, `geometry` a `BinGeometry`, and
`response` the same `InstrumentResponse` used for signal and final noise. Its
full pixel width must be positive; Gaussian sigma may be zero. Source redshift
must exceed `geometry.z_eval`: densities/variances are sampled at the source,
while clustering, P1D and distance conversions use `z_eval`. Field identity,
redshifts, h convention, geometry conversion factors, response, input arrays and
the convergence record are retained. Noise evaluation rejects a different field,
evaluation geometry or response. Changing area alone does not change this local
noise preparation.

All magnitude arrays have shape `(n_magnitude,)`, including a one-node sample.
Magnitudes are strictly increasing, quadrature weights are positive in
magnitudes, and `rho` is already selected/normalized in
`deg^-2 (km/s)^-1 mag^-1`. `variance` is the dimensionless delta-flux pixel
**variance**, not RMS. Forest length and response widths are km/s. No file
counts, target density, magnitude spacing, survey area or source selection are
inferred. `density_per_velocity(dndzdm, *, z_source)` converts a normalized
`dN/(dz_source dm deg^2)` row by `(1+z_source)/299792.458` without normalization.
Rectangular sums require constant `dm` at **every** node, including both
endpoints; nonuniform positive quadrature weights are equally valid.

`method` is `"early_lyaforecast"` (`sum_historical` recurrence, the adopted
prescription) or `"mcdonald"` (`sum_aliasing` recurrence). The reference scalars
`signal` (S, deg² km/s) and `alias` (B, km/s) are strictly positive, finite and
include the response, or are taken from an `auxiliary` sample. The iteration is
adaptive; non-convergence within `max_updates` raises `ValueError`, and no last
iterate is substituted. The result exposes `weights`, the full-sample integrals
`I1`, `I2`, `I3` (`r = rho*quadrature`: sums of `r*w`, `r*w**2`,
`r*w**2*variance`), `A = I2/(I1**2 L_v)` in deg², `P_pixel = I3 Delta_v/(I1**2
L_v)` in deg² km/s, and the convergence record. Zero-density nodes receive
zero weight. Unrepresentable arithmetic and empty weighted support fail; products
that underflow individually are rejected, and inputs are never rescaled.

For actual callable integration, `sample_auxiliary(field, geometry, response,
prepared_p3d, theta_p3d, *, p1d_model, p1d_parameters, theta_p1d, k_t_deg,
k_p_velocity)` returns immutable samples usable as `auxiliary=sample` in place
of `signal/alias`. It routes only the matching forest auto through the original
field indices and binding, and calls independently bound P1D once. It computes
`k_parallel=a_v*k_p_velocity`, `k_transverse=k_t_deg/d_deg`,
`S=P3D_auto*W^2*a_v/d_deg^2` (deg² km/s) and `B=P1D*W^2` (km/s).
`S` is the caller's full-auto P3D prediction after the field response and
angular/velocity unit conversion shown above; it is not an unsmoothed intrinsic
scalar. Explicit auxiliary examples are `(2.4, 0.00035)` for BAO and `(7, 0.001)` for P1D, in
`deg^-1` and `s/km`. There is no default. Auxiliary modes may lie outside
forecast cuts, add no modes, and must remain inside provider/template domains.
The retained fiducial states have no provider object or mutable model cache.

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
while differentiating intrinsic P3D. A different P1D with fixed weights
changes aliasing only; regenerating the weights can also change pixel noise. No automatic P1D/noise or covariance derivatives are added.

Run `python examples/weighted_noise.py` for two distinct forests, a galaxy,
negative signal crosses, explicit independent sampling, an independent matrix
covariance/amplitude Fisher oracle and fixed-state checks. Inputs are synthetic;
this is not a physical DESI-2 forecast. Raw readers and overlap-derived noise
remain outside these array APIs. The recurrences follow the early lyaforecast
code (GPLv3) and McDonald & Eisenstein (2007).

## Integrated forest sources

`prepare_integrated_forest_weights(field, geometry, response, source, *, method,
...)` generalizes the single-source-redshift weights to a forest field whose
background sources and pixels are integrated inside a redshift bin (see the
[survey page](survey.md) for the geometry). `source` is an
`IntegratedForestSource`; only `method='early_lyaforecast'` is supported, with the
same `rtol`, `min_updates`, `stable_steps`, `max_updates`, `signal`, `alias` and
`auxiliary` options and the same convergence requirement.

The quadrature measure of pixel node `p` and magnitude node `j` is
`mu_pj = (1+z_q,p) w_y,p c w_u,p dn/(dz dm)(z_q,p, m_j) w_m,j / L_bin` in deg^-2,
with `L_bin = c(u2-u1)`. The unchanged full-sum kernels run on the flattened
arrays with `density := mu`, `quadrature := 1` and `length := 1`, so the central
moments `L*sum(rho q w)` become

- `N1 = sum(mu w)`, replacing `L*I1[-1]`,
- `N2 = sum(mu w^2)` and `N3 = sum(mu w^2 v)`, replacing `L*I2[-1]`, `L*I3[-1]`
  (`N1`, `N2` in deg^-2 and `N3` in deg^-2 times the dimensionless variance `v`),

and the recurrences use `S = P + B/N1` with per-pixel noise `v*Delta_v/N1`, where
`Delta_v` is the bin's fixed velocity pixel width. The coefficients are
`A = N2/N1^2` (deg^2) and `P_pixel = Delta_v N3/N1^2` (deg^2 km/s), so
`forest_noise` consumes `IntegratedForestWeights` exactly as `ForestWeights`. The
weights have shape `(n_pixel, n_magnitude)`; `z_eff = sum(z_pix mu w)/N1` is a
diagnostic and does not enter the forecast.

Exact limit. If `(1+z_q) dn/dz_q` is independent of `z_q` (constant `dn/dv_q`)
and the variance depends on magnitude only, the integrals over `y` and `ln(lambda)`
of the overlap factor give `N1 = K (b-a)` for the untruncated window, with
`K = (1+z) dn/dz`, `[a,b] = ln[lambda_r,min, lambda_r,max]`, which equals
`L*I1` of the central weights with `rho = K/c`. All weights, `A` and `P_pixel`
then coincide with the central results (measured relative differences 6e-14;
asserted at 1e-12 in the tests). The pixel width convention below is what
separates the two modes otherwise.

Pixel width. The S/N tables are per angstrom. The pixel variance is queried at
`Delta_lambda(lambda) = Delta_v lambda/c`, with `Delta_v` the same fixed velocity
width as the response (converted at the evaluation wavelength). Then
`v*Delta_v = c/(lambda SNR_A^2 N_exp/N_file)` is the noise power per unit velocity
and does not depend on the angstrom width; at `lambda_c` it equals the central
value.

If the density covers the window but every pixel lies outside the S/N table, the
`1e20` sentinel variance of `ForestSNR` gives `N1 ~ 0` and a huge noise instead
of an error, as in the central mode; this cannot occur for the DESI-2 tables. The
error for absent density coverage is raised at the source construction.
