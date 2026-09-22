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
error. Normalization remains a raw count sum before cell-width division; its
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

