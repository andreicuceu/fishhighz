# Python interface

## Construction, preparation and calculation

```python
from fishhighz import Forecast

forecast = Forecast("my-survey.ini")
config = forecast.config
prepared = forecast.prepare()
result = forecast.run(batch_size=2048)
```

Construction reads and validates the native INI and records its identity/hash.
`Forecast()` loads the bundled recipe. `Forecast("desi2_accuracy.ini")` uses a
local file when one exists and otherwise selects the bundled recipe; other
missing paths raise an error. Prefer `Forecast()` when package defaults are
intended.

`prepare()` prepares the background, signed template and input readers, assembles
native `BinSpec` records, then prepares fixed grids, forest weights, noise and
covariance. The returned `PreparedForecast` exposes `config`, `survey`,
`background`, `template`, `bins` (also `prepared_bins`), `bin_specs` and
`input_identity`. Successful preparation is cached on this `Forecast` instance;
repeated `prepare()` calls return it. Construct a new instance for changed inputs.

`run()` uses that prepared state and calculates each bin's joint result plus
individual spectra with their own covariance. Repeated calls recalculate Fisher
constraints. Its arguments are Python options, not INI keys:

| Argument | Default | Meaning |
| --- | --- | --- |
| `batch_size` | `2048` | Integration batch size. |
| `step_scale` | `None` | Multiplier of parameter finite-difference steps; `None` uses INI `ap_step / 0.001`. |
| `numerical` | `False` | Force numerical model derivatives when true, rather than using available analytic derivatives. |

The native registry has fiducial `ap=at=1` in each bin and reference parameter
steps `0.001`. Thus the default effective AP step is the INI `ap_step`.
Weights, instrumental response and covariance are fixed at the fiducial model;
Fisher information comes from derivatives of the predicted mean.

## Inject prepared inputs

Prepared objects or factories can replace file-backed preparation:

```python
forecast = Forecast(
    "my-survey.ini",
    background=background,
    template=template,
    readers=readers,
)
```

Here `background`, `template` and `readers` are caller-supplied objects. The
alternative keyword arguments are `background_factory(config)`,
`template_factory(config, background)` and `readers_factory(config)`. Factories
are called during preparation only when the corresponding object is absent;
explicit objects take precedence. Defaults prepare every omitted component.
Injection does not bypass native schema or consistency validation.

The background supplies `h_fid`, `hubble_parameter(z)`,
`transverse_comoving_distance(z)`, template-growth/damping reference redshifts,
`sigma8_damping_reference`, and sigma8/growth at the requested redshifts. The
accepted callable pairs are `sigma8_at`/`growth_rate_at` or `sigma8`/`growth_rate`;
exact matching entries in `z_bins`, `sigma8_zbins`, `growth_rate_zbins` are also
supported. Template `z_ref` must match configured `template_redshift`, template
`h_fid` must match the background, and the background's reference redshifts must
match the configuration. The template must implement the spectrum evaluation
contract used by `KaiserModel`, not merely expose metadata.

The reader mapping must cover exactly the configured field IDs. Each value is a
mapping with `density` and, for forests, `snr`. Density adapters expose `sample`
or `query`; forest SNR adapters expose `sample` or `variance`. They must also
provide the input structure needed to partition magnitude integration. Using
the supplied reader classes retains those conventions; arbitrary arrays are
not substitutes for reader objects. Input provenance labels injected objects
and factories separately and retains known hashes without reopening their files.

## Lower-level scientific API

`parse_survey_ini(path)` validates without running CAMB. With an already prepared
background and template, `prepare_survey(config, background=..., template=...,
readers=...)` or `prepare_ini(path, ...)` assembles native survey records without
running CAMB. These functions return `PreparedSurvey`; fixed covariance
preparation belongs to `prepare_bin` or `Forecast.prepare()`.

For models, selections, parameters or noise beyond the native schema, assemble
`BinSpec` objects and use `prepare_bin` and `run_forecast` directly. The
[research baseline](../research/RESEARCH_BASELINE.md) explains these choices in
physical terms. Its synthetic example is available as a
{download}`Python script <../../examples/research_bao_forecast.py>` and requires
`fishhighz[templates]`. Historical compatibility examples additionally require
the optional lyaforecast reference installation; they are separate from the
native interface.
