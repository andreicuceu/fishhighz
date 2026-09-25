# Native INI reference

The current schema is `fishhighz-native-survey`, version `1`. It implements the
BAO accuracy prescription by default, with independent AP parameters per redshift
bin, and an explicit full-shape mode.
Unknown sections/options, duplicate definitions and `[DEFAULT]` values are
rejected. Sections and string values are case sensitive; option names are case
insensitive. Use full-line comments: trailing comments are not stripped. Numeric
lists accept commas or whitespace except where the option specifies otherwise.
Legacy lyaforecast INIs are not translated. No environment-variable
interpolation or URL/download syntax is supported.

Download the {download}`compact native recipe <../../fishhighz/data/desi2_accuracy.ini>`
or the {download}`complete annotated recipe <../../examples/desi2_accuracy_annotated.ini>`.
The section-by-section reference below includes the maintained annotated INI
directly, so the displayed options and downloadable file have one source.
Values described as defaults apply **with** `[prescription] name=accuracy`;
example survey, cosmology and field values are not implicit defaults.

## Required sections and defaults

Always supply `[schema]`, `[cosmology]`, `[survey]`, `[fields]`, one `[field ID]`
for each listed ID and one `[pairs bin N]` per bin. `[prescription]` is optional.
With `name=accuracy`, the complete `[model]`, `[input policies]` and `[numerical]`
sections may be omitted. The prescription also defaults `band=r`,
`lya_rest_angstrom=1215.67`, `evaluation_redshift=geometric_1plusz` and derives
`num_z_bins` from `z_edges` when omitted. Without the prescription, all displayed
keys in the model, policies, numerical, cosmology and survey sections are
required, subject to the weight-tolerance alias rule. Field requirements remain
as stated below. Explicit settings override defaults only within supported
values; fixed policy labels declare the implemented recipe rather than selecting
arbitrary alternative algorithms.

`[numerical] weight_rtol` and `[input policies] weighting_rtol` are aliases:
setting only one updates both; two explicit values must agree numerically.
Their native default is `1e-5`. The separate generic `accuracy.STOPPING` constant
uses `1e-4`; it is not the default for this native INI interface.

Overrides are recorded in saved effective settings. The prescription revision
identifies its starting defaults and does not establish independent numerical
qualification of altered inputs.

## `[schema]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 41-45
```

## `[prescription]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 46-51
```

## `[cosmology]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 52-72
```

## `[survey]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 73-102
```

## `[fields]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 103-110
```

## `[field lya(qso)]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 111-181
```

## `[field qso]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 182-193
```

## `[field lbg]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 194-205
```

## `[field lae]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 206-217
```

## `[field lya(lbg)]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 218-233
```

## `[model]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 234-265
```

## BAO with marginalized tracer nuisances

Add `[model] mode=bao_marginalized` to a compact accuracy INI. The resolved
model uses `parameterization=alpha_iso_phi`, `parameter_names=alpha_iso,phi`,
`smooth_scaling=identity`, `wiggle_scaling=alpha_iso_phi`,
`growth_rate=fixed_fiducial_f`, `biases=marginalized` and
`forest_bindings=shared`. An expanded INI must set these labels explicitly.
The wiggle mapping is
`a_parallel=alpha_iso*phi**(-2/3)` and
`a_perp=alpha_iso*phi**(1/3)`; both smooth dilations equal one.

Each bin retains the two dilation targets. The galaxy Kaiser factors use a
single fixed `f_fid(z)` and independent QSO, LBG and LAE biases. Both forest
samples share one Lyα bias and beta. Only nuisances represented in the selected
spectra remain active; there are seven parameters for a full 15-spectrum bin.
No nuisance prior is imposed. Template normalization, growth, damping widths,
weights, response, noise and covariance remain fixed during mean derivatives.

`[numerical] scale_step` and `relative_step` set dilation and relative nuisance
steps, respectively, with defaults `0.00025` and `0.001`.
`Forecast.run(step_scale=0.5)` halves these steps. A blank `selected` value
excludes a bin and preserves its original index. `Forecast.run(individuals=False)`
or CLI `--joint-only` skips individual spectra. Category cut controls and template
coverage checks follow the full-shape mode described below.

## Full-shape mode

Add `[model] mode=full_shape` to the compact accuracy recipe. This fills the
following supported full-shape defaults; existing explicit BAO scaling options
must be replaced when adapting an expanded BAO INI:

```ini
[model]
mode = full_shape
parameterization = alpha_phi
parameter_names = alpha_w, phi_w, alpha_s, phi_s, f
smooth_scaling = alpha_phi
wiggle_scaling = alpha_phi
growth = camb_sigma8_ratio
growth_rate = free_f
biases = marginalized
forest_bindings = shared
reported_growth = f_sigma8_fid

[numerical]
k_min = 0.01
k_max = 0.20
scale_step = 0.00025
relative_step = 0.001
```

`growth=camb_sigma8_ratio` fixes the template power normalization to the fiducial
background. The separate `growth_rate=free_f` fits the galaxy Kaiser growth rate;
results report it as `f*sigma8_fid(z)`. The wiggle and smooth components use
independent `alpha_phi` mappings, with `a_parallel=alpha/sqrt(phi)` and
`a_perp=alpha*sqrt(phi)`, including separate angular and Fourier-volume factors.

For isotropic dilation coordinates, set `parameterization`, `smooth_scaling`,
and `wiggle_scaling` to `alpha_iso_phi`, and set `parameter_names` to
`alpha_iso_w, phi_w, alpha_iso_s, phi_s, f`. A forest-only
`target_set=dilation_only` omits `f` from that list. Each component then has
`a_parallel=alpha_iso*phi**(-2/3)` and `a_perp=alpha_iso*phi**(1/3)`.
The two bases are distinct: `alpha_iso=alpha*phi**(1/6)`. The INI basis is
saved with the result, so historical `alpha_phi` uncertainties retain their
original interpretation.

Each bin marginalizes independent QSO, LBG and LAE biases, one Lyα bias shared
by both forest samples and one shared Lyα beta, without external priors. Absent
tracer nuisances are removed; all five targets remain. The bundled selection
therefore has 8 active parameters in bin 1 and 10 in bins 2–6. All parameters
are independent between bins. Damping widths, response, weights, noise and
covariance remain fiducial during mean differentiation.

`scale_step` is the absolute step for each dilation parameter; `relative_step`
multiplies the absolute nonzero fiducial growth or nuisance value, with an
absolute fallback at zero. `Forecast.run(step_scale=0.5)` halves all steps.
These initial steps require numerical qualification for each analysis.
Full-shape-only options are rejected in BAO mode; omitting `mode` preserves BAO.

Optional `[numerical]` values `k_max_galaxy_galaxy`,
`k_max_galaxy_forest`, and `k_max_forest_forest` set observed-coordinate
cutoffs by spectrum category; each defaults to `k_max`. All selected spectra
within a category share its cutoff. The integration is divided at distinct
cutoffs, with each interval using the covariance closure of its active
spectra. Interval Fisher matrices are added before nuisance marginalization.
`Forecast.run(individuals=False)` computes joint bin constraints without
the selected individual-spectrum calculations; the default computes both.

For forest-only selections, `[model] target_set=dilation_only` keeps four
dilation targets and fixes the Kaiser growth rate to its fiducial value.
Forest bias and beta remain marginalized. The default `target_set=full`
keeps all five targets. A blank `selected` value in a full-shape `[pairs bin N]`
section excludes that bin, retains its original bin index in results, and
creates an explicit joint record with status `excluded`.

Edit `k_min` and the applicable `k_max` values to change Fourier cuts. They remain fixed
observed-coordinate limits; quadrature nodes, weights and mode counts follow
any finite `0 < k_min < k_max`. Both mapped components must remain inside the
template domain at every derivative stencil evaluation. Unsupported coverage
raises an error, without clipping or changing cuts. `k_intervals`, `k_order`
and `mu_order` control radial and angular refinement. Saved settings record
requested limits and effective quadrature; see the [results guide](results.md).

## `[input policies]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 266-330
```

## `[numerical]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 331-359
```

## `[pairs bin 1]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 360-372
```

## `[pairs bin 2]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 373-375
```

## `[pairs bin 3]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 376-378
```

## `[pairs bin 4]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 379-381
```

## `[pairs bin 5]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 382-384
```

## `[pairs bin 6]`

```{literalinclude} ../../examples/desi2_accuracy_annotated.ini
:language: ini
:lines: 385-386
```

## Compact installed recipe

The bundled recipe uses the named prescription to keep survey choices visible
while taking the defaults documented above. Package resource paths make it
independent of the checkout directory.

```{literalinclude} ../../fishhighz/data/desi2_accuracy.ini
:language: ini
```

For physical interpretation see the
[research baseline](../research/RESEARCH_BASELINE.md); for running and saving
see the [Python guide](python.md) and [results guide](results.md).
