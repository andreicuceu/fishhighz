# Kaiser signal and BAO dilation

Run `python examples/builtin_forecast.py` after installing `.[dev,templates]`
(or use `fishhighz[templates]` for template preparation). This standalone example
constructs a synthetic oscillatory template and compares BAO-only dilation and
common AP plus galaxy f forecasts. It supplies fixed responses and diagonal
noise explicitly, checks grid/template h_fid equality, and factors covariance
once per mode. Nuisance priors are explicit toy external constraints. No survey
noise preparation or physical forecast-accuracy claim is implied.

`KaiserModel` and `Scaling` live in `fishhighz.models.kaiser`. Construct one model
per redshift bin using a prepared `PowerTemplate` and the original ordered
`ObservedField` sequence. All listed fields are supported by this model; every
one needs a bias and both fixed widths. The signature is:

```python
model = KaiserModel(
    template,
    fields,
    biases={"F": "bf", "g": "bg"},
    betas={"F": "beta"},  # exactly the forest IDs; no galaxy beta
    widths={"F": (4.0, 2.0), "g": (3.0, 1.0)},  # parallel, transverse
    f="rate",  # one shared galaxy f; omit for forest-only models
    local_names=("bf", "bg", "beta", "rate", "w_ap", "w_at"),
    wiggle=Scaling("ap_at", ap="w_ap", at="w_at"),
    # omitted smooth is fixed identity: BAO mode
)
```

Every bias, forest beta, shared galaxy f, or scaling coordinate accepts a fixed
finite number or an explicitly named local slot string. `local_names` defines
exact theta_local order. Each declared slot must occur exactly once in these
settings: unknown, missing, unused and duplicate assignments fail. Use distinct
local slots mapped to one global ID through existing `BoundParameters` to tie
fields or smooth/wiggle blocks; neither field labels nor names imply sharing.
There is no additional registry or expression language. In particular, a width
or G cannot be a free slot. Explicitly choose all finite biases/betas/f or caller
registry bounds; no positivity prior is inferred for them.

Wrap the callable with `P3DProvider(label, model, bound_parameters, owned_pairs)`
and `PreparedP3D` exactly as for an external provider. Pairs use original integer
field indices, with requested order preserved, including reversed pairs and
covariance-required unselected spectra. The callable signature is unchanged:
`model(theta_local, z, k, mu, pairs) -> (node,pair)`. Use the Step 06 numerical
derivative path with explicit steps; a production analytic Jacobian is not added.

Each scaling block chooses exactly one basis and its two coordinates:

| Public basis | Coordinates | a_parallel | a_perp | Vega name |
| --- | --- | --- | --- | --- |
| `ap_at` | `ap`, `at` | ap | at | ap_at |
| `alpha_phi` | `alpha`, `phi` | alpha/sqrt(phi) | alpha*sqrt(phi) | phi_alpha |
| `alpha_iso_epsilon` | `alpha_iso`, `epsilon` | alpha_iso*(1+epsilon)^2 | alpha_iso/(1+epsilon) | aiso_epsilon |

Positive scales/alpha/phi and epsilon>-1 are required, with finite representable
derived factors and Q. `aiso_aap` is unsupported. Smooth and wiggle blocks may
use separate bases, fixed values, partial freedom or independent local slots.
To implement common AP scaling, tie distinct slots in both same-basis blocks to
the same global IDs. No nonlinear cross-basis ties are inferred. Dilations are
template-coordinate parameters; interpreting them as BAO distance/sound-horizon
ratios or geometric AP ratios belongs to the caller, without inferred cosmology.

For each component independently, Fourier mapping uses
`k_parallel=k*mu/ap`, `k_perp=k*sqrt(1-mu²)/at`, their radial magnitude and
`mu_component=k_parallel/k_component`, plus `Q=1/(ap*at²)`. Both the template and
Kaiser angular factor use these component coordinates. Exact identity scaling
preserves the supplied radial coordinates, including template endpoints.

Forest factors are `b*(1+beta*mu_component²)`; galaxy factors are
`b+f*mu_component²`, including b=0 without division. Signed pair products preserve
negative forest crosses. Wiggle damping is
`exp[-(k_parallel,w²*Sigma_parallel,ij²+k_perp,w²*Sigma_perp,ij²)/2]`, with each
pair width squared the mean of its two field widths squared. Autos recover their
widths and cross damping is the geometric mean of the two auto factors. Widths
are explicit nonnegative constants in Mpc/h_fid; zero gives unit damping and a
large finite exponent may underflow to zero. Smooth power is never damped.
Changing f does not change widths or G; varying wiggle coordinates changes D.

The output is `G*(Q_s*B_i,s*B_j,s*P_s + Q_w*B_i,w*B_j,w*D_ij*P_w)` in
`(Mpc/h_fid)^3`, before instrument response or noise. G defaults to 1 only at
z_ref; otherwise prepare with explicit fixed `z=...` and `growth=G`, following
Step 07's squared-growth-ratio convention. Calls at a different z fail. The
model evaluates the reference template and multiplies by G once; retained F_ZREF
and SIGMA8 metadata do not set f or infer evolution. Observed-total PSD validation
remains in covariance preparation; separate scalings/signed components need not
produce a physical total for arbitrary parameters, and are never repaired here.

Provide template coverage for **both mapped components at every fiducial and
finite-difference state**. Domain failures report component, mapped range,
available domain and scales; the derivative engine adds provider/stencil context.
No clipping, extrapolation, dropped nodes, changed steps or automatic regridding
occurs. The model has no second set of forecast cuts or modes. Observed nodes,
weights, volume and k_min/k_max remain fixed, and callers must use the template's
h_fid convention without silent coordinate conversion during derivatives.

Preparation owns fixed settings/widths and integer free-slot maps. Repeated
calls use only NumPy, evaluate each template/factor batch once per component,
reuse results across pairs, and hold no mutable parameter-dependent cache.
Read-only/noncontiguous paired nodes and consecutive slices are supported.
Lazy optional-import checks, independent scalar/product-rule tests and bounded
individual Vega-method comparisons cover the implementation; they do not
validate a complete Vega pipeline, survey model or JIT performance.

