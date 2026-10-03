# Prepared surveys and raw inputs

`fishhighz.survey.BinSpec(id, geometry, grid, p3d, responses, *, ...)` supplies
explicit existing objects. `p3d` owns the only pair selection and field order.
Use `forests={id: ForestInput(...)}`, `galaxies={id: n_bar}` and
`independent_sampling=True`, or `full_noise=N` as a complete replacement with
none of those generated inputs. Every field requires an `InstrumentResponse`,
including identity settings for galaxies. Only active fields require noise data.

```python
from fishhighz.forecast import prepare_bin, run_bin, run_forecast

prepared = [prepare_bin(spec) for spec in bin_specs]
one = run_bin(prepared[0], batch_size=64, step_scale=1.0)
combined = run_forecast(prepared, batch_size=64, prior_fisher=prior)
refined = run_forecast(prepared, step_scale=0.5, prior_fisher=other_prior)
```

Bins share the **same** `ParameterRegistry` object. Explicit bindings tie shared
physics or give distinct nuisance IDs; unused global columns stay zero. Bin IDs
must be unique, redshift interiors must not overlap, and h_fid and reused field
identities must agree. Unequal widths, gaps and non-midpoint z_eval values are
supported: each geometry integrates its own bounds and area. Raw source-table
cells never define forecast-bin bounds or evaluation redshifts. Caller bin order is preserved. `BinRun.result` is a
zero-prior `FisherResult`; `ForecastRun.combined` adds the supplied prior once.
Singular bin information is valid; no automatic inversion or marginalization
occurs. Inspect `result.diagnostics` or explicitly request marginalized errors.

Preparation freezes intrinsic fiducial powers, response products, known noise,
mode counts and Cholesky factors. Prepared arrays are immutable owned snapshots:
`k`, `mu`, `modes` are `(node,)`; `power`, `noise`, `products`, `total` are
`(node, required_pair)`; `response` is `(node, field)`; `factors` are
`(node, selected_pair, selected_pair)`. Noise is absent from the default mean.
Repeated runs differentiate that mean and reuse every fixed survey quantity.
Numerical step overrides, `step_scale`, and `numerical` follow the existing
`evaluate_derivatives` contract. Errors retain bin, provider/field, global node
slice and stencil context. No cuts or steps are adjusted automatically.

`batch_size=None` uses all nodes; a positive integer bounds transient Jacobians
with consecutive C-order slices. Full factors remain in memory, requiring
O(n_node*n_selected²) storage. This is not disk streaming or a performance claim.
Small accumulation roundoff can differ with batch size. `BinRun.calls` records
actual derivative provider calls aggregated across batches, including repeated
fiducial calls; `columns` records methods/stencils. `PreparedBin.diagnostics`
records preparation calls, units, IDs, volume, response/noise conventions and
source provenance; `weights` retains per-field weighting diagnostics. Model
callables must be deterministic and must not be externally mutated between
preparation and runs. Changing the fiducial/model requires fresh preparation;
callable internals are not deep-copied, hashed, or serialized.

Install `fishhighz[survey]` for the lazy SciPy readers in
`fishhighz.adapters.legacy_inputs`; normalized arrays and external models need
only NumPy. Standalone reader signatures are:

- `DensityReader(path, *, semantics='cell_count_per_deg2', target_density=...,
  z_norm_min=..., magnitude_bounds=None, redshift_widths=None, width_policy=None,
  label='density')`. The semantics argument
  is required. Set target and threshold explicitly to `None` for no normalization
  and whole-grid selection. Inclusive raw magnitude masks precede the strict
  `z > z_norm_min` raw count sum. Divide each row by Delta_z[i]*dm exactly once;
  later magnitude quadrature is separate. Supply `redshift_widths` as a positive
  (n_z,) vector aligned with reconstructed sorted redshifts for physical cell
  measures, including irregular centres. Widths are owned and immutable; no
  physical widths are inferred from centres or forecast bins. Without widths,
  the default (`width_policy=None` or `'uniform'`) requires uniform redshifts.
  Explicit `width_policy='legacy_first_spacing'` uses z[1]-z[0] for every row,
  including irregular tables, solely as the labeled legacy conversion. It is
  never selected automatically. Widths and compatibility policies conflict.
  Provenance retains axes, effective widths, their ordering and policy label. `query(z, magnitudes)` returns ordered differential
  density. `local_galaxy_density(geometry, magnitudes, quadrature)` explicitly
  chooses the local z_eval approximation, without an area multiplier.
- `SNRReader(paths, *, smoothing, label='SNR')` parses header magnitudes and
  metadata, sorts files, and requires identical grids. `smoothing='legacy'` uses
  sigma=10 **sample indices**, reflect boundaries and truncate=4 on wavelength
  only; `'none'` explicitly disables smoothing. Linear interpolation follows.
  `query(*, z_source, magnitudes, wavelength)` returns SNR per Angstrom.
  `variance(..., pixel_width_angstrom, exposure_count, exposure_time=None)` returns
  dimensionless delta-flux variance `1/(SNR²*Delta_lambda*N_exp/N_exp_file)`.
  Explicit exposure time must match file EXPTIME.
- `sample_forest_readers(density, snr, geometry, response, *, z_source,
  magnitudes, pixel_width_angstrom, exposure_count, exposure_time=None)` returns
  immutable `rho`, `variance`, and plain `provenance`. Both readers use z_source;
  observed wavelength uses z_eval. It verifies the pixel velocity width using
  c=299792.458 km/s. Pass these arrays in `ForestInput.weight_options` with explicit
  quadrature, forest velocity length, and the Step 10 weighting settings.
  `ForestInput` also requires an independent P1D callable, binding and state.
  Optional `auxiliary_coordinates=(k_t_deg,k_p_velocity)` selects one auto/P1D
  query for `legacy`, `early_lyaforecast`, or `mcdonald` weights; supplied and
  inverse-variance weights bypass auxiliary calls.

Readers retain resolved paths, SHA-256 hashes, owned tables, domains, units,
normalization/exposure settings and interpolation provenance. They reject
extrapolation, negative spline overshoot, irregular magnitude grids, inconsistent
SNR grids and unrepresentable arithmetic. Uniform-axis checks allow only
64*eps64*max(1,max(abs(axis))) absolute coordinate roundoff. There are no floors,
bright caps, population averaging or automatic resource lookup. Nonuniform raw
redshifts without explicit widths or legacy_first_spacing fail with an actionable
error. With the default `interpolation='piecewise_constant'` the cells must also
tile each axis: a gap or overlap between `centre + width/2` and the next lower
edge larger than the same roundoff tolerance (for example rounded centres under
`legacy_first_spacing`) raises `ValueError`; contiguous explicit nonuniform
widths and the spline are accepted. Normalization remains a raw count sum before cell-width division; its
measure is distinct from the later magnitude integral.
Live readers never enter the prepared forecast or derivative loops. Adapted
interpolation conventions retain lyaforecast GPLv3/source provenance in the module.

Run `python examples/survey_forecast.py` for generated raw fixtures, two distinct
forests plus a galaxy, two unequal independent bins with explicit non-midpoint
evaluation redshifts, shared target/distinct nuisances,
repeated runs and explicit prior/error diagnostics. Its outputs are synthetic,
not DESI-2 forecasts. The native `Forecast` interface described in the [API reference](../api/forecast.md) provides
INI configuration, a CLI and JSON/NPZ output; unchanged lyaforecast INIs are not
translated.

Input snapshots preserve non-object ndarray dtypes and values until scientific
validation. Complex, Boolean, string and object scientific arrays are rejected
consistently with the direct weighting API; no imaginary parts are discarded.
Accepted prepared scientific arrays are immutable float64. Plain metadata arrays
retain their dtype (including indices and Boolean flags). Object arrays cannot be
frozen as immutable byte-backed values and are rejected; use typed arrays or
plain nested metadata instead.

## Integrated forest-source mode

The default central mode samples density and S/N of each forest field at a single
source redshift per bin, `z_s = lambda_Lya(1+z_eval)/sqrt(lambda_r,min lambda_r,max) - 1`.
`[input policies] forest_source_integration = integrated` (see the
[INI reference](../user/ini.md#integrated-forest-sources)) replaces it by an
integral over the background sources and forest pixels contributing to the bin.
Galaxy fields, `z_eval`, the response, the auxiliary reference mode and the
magnitude partition (nodes and weights) are unchanged.

Geometry. With `u = ln(lambda)`, `y = ln(1+z_q)`, the bin slice
`[u1, u2] = ln[lambda_Lya(1+z_min), lambda_Lya(1+z_max)]` (the actual bin bounds,
not `z_eval`), `L_bin = c(u2-u1)` and the forest `[a,b] = ln[lambda_r,min,
lambda_r,max]`, a source at `y` contributes pixels `u in [max(u1, y+a),
min(u2, y+b)]`. The source window is `y in [u1-b, u2-a]` intersected with
`[ln(1+min_zq_forest), ln(1+max_zq_forest)]`. The measure in deg^-2 is

`mu = (1+z_q) c w_y w_u (dn/dz_q dm) w_m / L_bin`,

so that the sums of the central mode, `L*I1` etc., become `N1 = sum(mu w)`, `N2`,
`N3` ([weights](weights.md#integrated-forest-sources)). For an untruncated window
and constant `dn/dv_q`, `N1` equals `(b-a) K` and the central `L*I1`.

Quadrature. Composite Gauss-Legendre in `y` with panel boundaries at the window
ends, the overlap kinks `u1-a` and `u2-b`, all density-table redshift cell edges
and all S/N-table source-redshift nodes inside the window (boundaries closer than
64 eps are merged); the overlap length is piecewise linear and integrated exactly,
and the piecewise-constant density and piecewise-linear S/N are smooth inside each
panel. Within each overlap, `forest_lambda_panels` equal wavelength panels each
carry a `forest_lambda_order` rule; `forest_zq_order` is the order per `y` panel.
The geometry identity `sum(w_y overlap) = (b-a)(u2-u1)` holds to roundoff.
Refinement of the orders is tested by the order ladder in
`scripts/check_integrated_forest_refinement.py`, which also compares N1-N3 with a
high-resolution midpoint rule; the 0.1 per cent BAO-error target applies.

Inputs. Density is the piecewise-constant table (the spline is rejected) queried
on the `y` nodes by `LegacyDensity.sample_grid` or `DensityReader.query_grid`;
the pixel variance is queried on (z_q, lambda, magnitude) by
`LegacySNR.variance_grid` or `SNRReader.variance_grid` (legacy sentinel, bright
clamp and floor, or strict domain errors, as for the scalar queries, with fallback
counts recorded in the provenance). Queries are chunked to at most 2**20 points.
The pixel width in angstrom is `Delta_v lambda/c` with the fixed velocity width of
the bin response.

No coverage. If no source redshift in the window has density above the legacy
floor (for example the density table ends below the window), preparation raises
`ValueError` naming the field and bin; no floor is applied. If the density covers
the window but all pixels fall outside the S/N table, the legacy adapter's sentinel
variance gives `N1 ~ 0` instead of an error, as in the central legacy path.

Cost. The weight solve runs on `n_pixel * n_magnitude` nodes (about 11 million for
the default orders of one DESI-2 bin), roughly 1.5 minutes per field and bin on one
login-node thread; reduce the orders for quick tests.
