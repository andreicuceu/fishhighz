# Linear-power templates

External-model users can install `fishhighz` with NumPy alone. To prepare local
Vega-format FITS templates, install `fishhighz[templates]`; Astropy reads FITS
and SciPy prepares coefficients. Development acceptance uses `.[dev,templates]`
and executes every template test (missing extras are not silently skipped).
Neither extra is imported by package/module import or repeated evaluation.
Missing preparation dependencies give an installation message when requested.

Obtain a decomposed template externally, for example a file produced by Vega's
`bin/make_template.py`. FishHighz reads that format without running Vega, CAMB or
the generator; it bundles no cosmological template assets. Use a local file whose
provenance and scale coverage are appropriate for your application:

```python
import numpy as np
from fishhighz.models.templates import load_template

# Replace this path with your externally prepared Vega-format file.
template = load_template("/path/to/linear_power.fits", h_fid=0.67)
k_min, k_max = template.domain  # converted h_fid/Mpc; not forecast cuts
k = np.geomspace(k_min, k_max, 25)
components = template.evaluate(k)  # z_ref, G=1; (query,2)
smooth, wiggle = components.T  # template.component_names == ("smooth", "wiggle")
dp_dk = template.evaluate(k, derivative=1)
# Illustrative caller-supplied growth ratio 0.9, squared exactly once:
at_other_z = template.evaluate(k, z=template.z_ref + 0.2, growth=0.9**2)
```

The loader selects the unique named `PK` binary table, wherever it occurs in the
file. Scalar real columns `K`, `PK`, `PKSB` must have at least four finite samples,
with strictly positive increasing K. Extra columns/extensions are ignored.
`PKSB` is the supplied smooth component and `PK-PKSB` is the wiggle component;
negative values in **either** are preserved. No positivity constraint, logarithm
of amplitudes, sideband fit, smoothing or new decomposition is applied.

The table header's `ZREF` and `H0/100` specify z_ref and h_template. If absent,
supply `z_ref=...` and/or `h_template=...` explicitly to `load_template`.
Primary-header values do not substitute for missing table metadata. Supplied and
header values must agree within `64*eps64*max(1,abs(header),abs(supplied))`, comparing
h after H0/100; consistent table values take precedence. ZREF must be finite and
nonnegative and H0/h values finite and positive. Original scalar header metadata
(including H0, F_ZREF, SIGMA8 and cosmology) remains immutable provenance; it never
implies a growth function or a forecast parameter. COMMENT/HISTORY cards are not
retained. `source_path` and `source_sha256` identify the exact single byte snapshot
parsed, even if the file changes later. Prepared objects retain no file handles.

Absent TUNIT cards mean the documented Vega convention: k in h_template/Mpc,
power in (Mpc/h_template)^3. The bounded unit-spelling table accepts the following
case-sensitive forms, ignoring spaces and treating `**` like `^`:

| Column | Accepted explicit units |
| --- | --- |
| K | `h/Mpc`, `h Mpc^-1`, `h Mpc-1` |
| PK, PKSB | `(Mpc/h)^3`, `(Mpc/h)3`, `Mpc^3/h^3`, `Mpc3/h3`, `Mpc^3 h^-3`, `Mpc3 h-3` |

Other spellings, including physical `1/Mpc` or `Mpc^3`, are rejected. This is
not a general unit conversion library. Explicit h_fid is always required, and
preparation applies `k_fid=k_file*h_template/h_fid` and
`P_fid=P_file*(h_fid/h_template)^3` once. The same factor scales smooth and wiggle.
Source `k_file`, `pk_file`, `pksb_file` and both h values remain available to audit
this conversion. It changes units for a fixed cosmology; it is not AP scaling.
Unrepresentable conversions or collapsed logarithmic knots fail.

For in-memory input, use
`prepare_template(k, pk, pksb, z_ref=..., h_template=..., h_fid=..., metadata=None)`.
It follows the same source units and numerical validation and only needs SciPy.
Its optional metadata is copied scalar/string provenance; required numerical
inputs are explicit. Both factories return `PowerTemplate`. All source/converted
samples and prepared arrays are owned read-only float64. Component order is
`("smooth", "wiggle")`; converted `full` retains the original full-power samples.

Interpolation is a **not-a-knot cubic spline in ln(k_fid)** for both signed linear
amplitudes. `log_k` stores knots; `coefficients` stores descending cubic powers
of `ln(k)-log_k[interval]`, shaped `(4,interval,component)`. SciPy constructs these
once. NumPy interval lookup and Horner kernels evaluate them; no SciPy object is
retained or invoked. Interpolating both components and summing reproduces full
power interpolation up to floating-point arithmetic. `derivative=1` includes the
1/k chain factor and returns dP/dk in `(Mpc/h_fid)^4`, not dP/dln(k).

Queries are positive finite nonempty 1D arrays. Order/repetitions, read-only and
noncontiguous arrays and consecutive slices are supported; outputs are owned
C-order float64 `(query,2)`. Exact converted endpoints are accepted, with the
upper endpoint in the last interval. Even a one-ULP excursion fails before
logarithms with requested/available domain context. No extrapolation, clipping
or added modes occurs; later transformed queries require caller-provided coverage.

At exactly z_ref, G defaults to 1. At another nonnegative finite z, supply the
positive finite **power factor** `growth=G=[D(z)/D(z_ref)]^2`. Both powers and
k derivatives are multiplied by G once. Explicit G at z_ref must equal 1 within
`64*eps64*max(1,abs(G))`; arbitrary free amplitude changes are not this operation.
Reference metadata/coefficients stay fixed. This layer supplies no growth solver,
Kaiser factors, BAO/AP transformations, damping, survey physics or forecast cuts.
Analytic toy sampling convergence and local format checks do not establish
physical cosmology or survey-forecast accuracy.

