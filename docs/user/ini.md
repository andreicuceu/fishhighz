# Native INI reference

The current schema is `fishhighz-native-survey`, version `1`. It implements the
BAO accuracy prescription with independent AP parameters per redshift bin.
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
