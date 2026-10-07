# Geometry, response, weighting, and noise

See [geometry and response](../methods/geometry.md) and
[forest weighting](../methods/weights.md) (prescription: [forest-weighting note](../methods/forest-weighting.md)) for formulas and worked examples.
Response amplitudes have shape `(node, field)`; their pair products multiply
intrinsic signal once. Noise powers have shape `(node, required_pair)` and are
added after signal response. Geometry uses a fixed fiducial h and one common
volume per bin. Velocity widths use km/s, `a_v` converts comoving length to
velocity, P1D uses km/s, and P3D/noise use `(Mpc/h_fid)^3`.

Prepared weighting quantities belong to a specific field, geometry, and response.
Supplied full noise replaces generated noise; independent sampling must otherwise
be explicit. Invalid or nonfinite inputs, mismatched contexts, and unsupported
noise combinations raise `ValueError`. Astropy geometry is optional through the
`cosmology` extra; the supplied-background path does not need it.

## `fishhighz.geometry`

```{eval-rst}
.. autoclass:: fishhighz.geometry.BinGeometry
```

```{eval-rst}
.. autofunction:: fishhighz.geometry.prepare_geometry
```

```{eval-rst}
.. autofunction:: fishhighz.geometry.prepare_astropy_geometry
```

```{eval-rst}
.. autofunction:: fishhighz.geometry.mode_counts
```

```{eval-rst}
.. autofunction:: fishhighz.geometry.wavenumber_comoving_to_velocity
```

```{eval-rst}
.. autofunction:: fishhighz.geometry.wavenumber_velocity_to_comoving
```

```{eval-rst}
.. autofunction:: fishhighz.geometry.p1d_velocity_to_comoving
```

```{eval-rst}
.. autofunction:: fishhighz.geometry.p1d_comoving_to_velocity
```

```{eval-rst}
.. autofunction:: fishhighz.geometry.width_velocity_to_comoving
```

```{eval-rst}
.. autofunction:: fishhighz.geometry.width_comoving_to_velocity
```

## `fishhighz.response`

```{eval-rst}
.. autoclass:: fishhighz.response.InstrumentResponse
```

```{eval-rst}
.. autofunction:: fishhighz.response.pixel_width_angstrom_to_velocity
```

```{eval-rst}
.. autofunction:: fishhighz.response.gaussian_sigma_angstrom_to_velocity
```

```{eval-rst}
.. autofunction:: fishhighz.response.gaussian_fwhm_velocity_to_sigma
```

```{eval-rst}
.. autofunction:: fishhighz.response.resolving_power_fwhm_to_sigma
```

```{eval-rst}
.. autofunction:: fishhighz.response.velocity_response
```

```{eval-rst}
.. autofunction:: fishhighz.response.prepare_response
```

```{eval-rst}
.. autofunction:: fishhighz.response.pair_response
```

## `fishhighz.weights`

```{eval-rst}
.. autoclass:: fishhighz.weights.AuxiliarySamples
```

```{eval-rst}
.. autoclass:: fishhighz.weights.ForestWeights
   :members: validate_context
   :undoc-members:
```

```{eval-rst}
.. autofunction:: fishhighz.weights.density_per_velocity
```

```{eval-rst}
.. autofunction:: fishhighz.weights.sample_auxiliary
```

```{eval-rst}
.. autofunction:: fishhighz.weights.prepare_forest_weights
```

```{eval-rst}
.. autoclass:: fishhighz.weights.IntegratedForestWeights
```

```{eval-rst}
.. autofunction:: fishhighz.weights.prepare_integrated_forest_weights
```

## `fishhighz.noise`

```{eval-rst}
.. autoclass:: fishhighz.noise.ForestNoise
```

```{eval-rst}
.. autofunction:: fishhighz.noise.local_galaxy_density
```

```{eval-rst}
.. autofunction:: fishhighz.noise.galaxy_noise
```

```{eval-rst}
.. autofunction:: fishhighz.noise.forest_noise
```

```{eval-rst}
.. autofunction:: fishhighz.noise.prepare_noise
```


## `fishhighz.survey_tables`

Source-density and forest S/N tables. `DensityTable` and `SNRTable` read and
validate the raw tables and raise outside their domains; `SourceDensity` and
`ForestSNR` apply the boundary conventions of the native INI path (source
density $10^{-20}$ outside the cells, S/N clamped at the brightest tabulated
magnitude and a $10^{20}$ variance outside the table).

```{eval-rst}
.. autoclass:: fishhighz.survey_tables.CellHistogram2D
```

```{eval-rst}
.. autoclass:: fishhighz.survey_tables.DensityTable
```

```{eval-rst}
.. autoclass:: fishhighz.survey_tables.SNRTable
```

```{eval-rst}
.. autoclass:: fishhighz.survey_tables.SourceDensity
```

```{eval-rst}
.. autoclass:: fishhighz.survey_tables.ForestSNR
```

```{eval-rst}
.. autofunction:: fishhighz.survey_tables.sample_forest_tables
```

```{eval-rst}
.. autofunction:: fishhighz.survey_tables.density_discontinuity_edges
```
