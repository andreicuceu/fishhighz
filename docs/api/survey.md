# Geometry, response, weighting, and noise

See [geometry and response](../methods/geometry.md) and
[forest weighting](../methods/weights.md) for formulas and worked examples.
Response amplitudes have shape `(node, field)`; their pair products multiply
intrinsic signal once. Noise powers have shape `(node, required_pair)` and are
added after signal response. Geometry uses a fixed fiducial h and one common
volume per bin. Velocity widths use km/s, `a_v` converts comoving length to
velocity, P1D uses km/s, and P3D/noise use `(Mpc/h_fid)^3`.

Prepared weighting quantities belong to a specific field, geometry, and response.
Supplied full noise replaces generated noise; independent sampling must otherwise
be explicit. Invalid or nonfinite inputs, mismatched contexts, and unsupported
noise combinations raise `ValueError`. Astropy geometry is optional through the
`cosmology` extra; the supplied-background path requires only NumPy.

## `fishhighz.geometry`

```{autoclass} fishhighz.geometry.BinGeometry
```

```{autofunction} fishhighz.geometry.prepare_geometry
```

```{autofunction} fishhighz.geometry.prepare_astropy_geometry
```

```{autofunction} fishhighz.geometry.mode_counts
```

```{autofunction} fishhighz.geometry.wavenumber_comoving_to_velocity
```

```{autofunction} fishhighz.geometry.wavenumber_velocity_to_comoving
```

```{autofunction} fishhighz.geometry.p1d_velocity_to_comoving
```

```{autofunction} fishhighz.geometry.p1d_comoving_to_velocity
```

```{autofunction} fishhighz.geometry.width_velocity_to_comoving
```

```{autofunction} fishhighz.geometry.width_comoving_to_velocity
```

## `fishhighz.response`

```{autoclass} fishhighz.response.InstrumentResponse
```

```{autofunction} fishhighz.response.pixel_width_angstrom_to_velocity
```

```{autofunction} fishhighz.response.gaussian_sigma_angstrom_to_velocity
```

```{autofunction} fishhighz.response.gaussian_fwhm_velocity_to_sigma
```

```{autofunction} fishhighz.response.resolving_power_fwhm_to_sigma
```

```{autofunction} fishhighz.response.legacy_resolving_power_to_sigma
```

```{autofunction} fishhighz.response.velocity_response
```

```{autofunction} fishhighz.response.prepare_response
```

```{autofunction} fishhighz.response.pair_response
```

## `fishhighz.weights`

```{autoclass} fishhighz.weights.AuxiliarySamples
```

```{autoclass} fishhighz.weights.ForestWeights
:members: validate_context
:undoc-members:
```

```{autofunction} fishhighz.weights.density_per_velocity
```

```{autofunction} fishhighz.weights.sample_auxiliary
```

```{autofunction} fishhighz.weights.prepare_forest_weights
```

## `fishhighz.noise`

```{autoclass} fishhighz.noise.ForestNoise
```

```{autofunction} fishhighz.noise.local_galaxy_density
```

```{autofunction} fishhighz.noise.galaxy_noise
```

```{autofunction} fishhighz.noise.forest_noise
```

```{autofunction} fishhighz.noise.prepare_noise
```

