# Native forecast and configuration

`Forecast(source=None)` parses the bundled native INI when no source is given.
A supplied source is a native INI path or supported parser input. Construction
parses configuration; `prepare()` returns a cached `PreparedForecast`, and
`run()` returns `SurveyResult`. Explicit background, template, and reader objects
or factories replace the corresponding default preparation stages. Factories
receive `(config)`, `(config, background)`, and `(config)`, respectively.

`run(batch_size=2048, step_scale=None, numerical=False)` retains joint and
individual constraints. Individual spectra use their own covariance; they are
not blocks of an inverse joint covariance. With `step_scale=None`, the configured
AP step is divided by 0.001 to obtain the derivative multiplier. Prepared arrays
and fiducial weights remain fixed during differentiation.

`SurveyResult.individual`, `.joint`, and `.excluded` preserve explicit constraint
records. Each `SpectrumConstraint` records bin bounds/redshift, pair, parameter
IDs, status, optional `FisherResult`, dimensionless dilation errors and their
correlation. A missing or singular constraint is not a zero error.
`save(path)` creates a new directory containing JSON identities/statuses and NPZ
matrices; an existing directory raises `FileExistsError`.

Malformed native configuration raises `ValueError`; unsupported schemas raise
`UnsupportedSchemaError`. Missing files and optional dependencies report their
underlying file/import errors. Scientific domain and shape errors are not
silently repaired. Bundled resource paths obtained with `bundled_path` or
`bundled_paths` are context-managed and valid within their context.

The default CAMB background requires the `camb` extra, FITS templates the
`templates` extra, and raw survey readers the `survey` extra. Install the required extras together for native survey preparation. CAMB distances use Mpc, H uses km/s/Mpc, and growth
rates and sigma8 are dimensionless; the fixed `h_fid` defines forecast units.

## `fishhighz.public`

```{eval-rst}
.. autoclass:: fishhighz.public.Forecast
   :members: prepare, run
   :undoc-members:
```

```{eval-rst}
.. autoclass:: fishhighz.public.PreparedForecast
   :members: prepared_bins, bin_specs
   :undoc-members:
```

```{eval-rst}
.. autoclass:: fishhighz.public.SpectrumConstraint
   :members: result, available
   :undoc-members:
```

```{eval-rst}
.. autoclass:: fishhighz.public.SurveyResult
   :members: individual_spectra, joint_constraints, constraints, fisher_results, save
   :undoc-members:
```

## `fishhighz.survey_config`

```{eval-rst}
.. autoclass:: fishhighz.survey_config.FieldConfig
```

```{eval-rst}
.. autoclass:: fishhighz.survey_config.BinConfig
```

```{eval-rst}
.. autoclass:: fishhighz.survey_config.SurveyConfig
   :members: observed_fields, provenance
   :undoc-members:
```

```{eval-rst}
.. autoclass:: fishhighz.survey_config.PreparedSurvey
   :members: bin_specs
   :undoc-members:
```

```{eval-rst}
.. autoclass:: fishhighz.survey_config.UnsupportedSchemaError
```

```{eval-rst}
.. autofunction:: fishhighz.survey_config.parse_ini
```

```{eval-rst}
.. autofunction:: fishhighz.survey_config.parse_survey_ini
```

```{eval-rst}
.. autofunction:: fishhighz.survey_config.prepare_ini
```

```{eval-rst}
.. autofunction:: fishhighz.survey_config.prepare_survey
```

## `fishhighz.resources`

```{eval-rst}
.. autofunction:: fishhighz.resources.bundled_resource
```

```{eval-rst}
.. autofunction:: fishhighz.resources.bundled_path
```

```{eval-rst}
.. autofunction:: fishhighz.resources.bundled_paths
```

```{eval-rst}
.. autofunction:: fishhighz.resources.resolve_input_path
```

## `fishhighz.cosmology`

```{eval-rst}
.. autoclass:: fishhighz.cosmology.CAMBBackground
   :members: z_bins, H, D_M, sigma8_zbins, growth_rate_zbins, f, H_values, D_M_values, f_values, h_fid, sigma8, results, z_to_index, index, sigma8_at, growth_rate_at, growth_rate, hubble_at, transverse_distance_at, hubble_parameter, comoving_radial_distance, transverse_comoving_distance
   :undoc-members:
```

```{eval-rst}
.. autofunction:: fishhighz.cosmology.prepare_camb
```

```{eval-rst}
.. autofunction:: fishhighz.cosmology.prepare_camb_background
```

