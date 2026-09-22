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
.. autofunction:: fishhighz.response.legacy_resolving_power_to_sigma
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

