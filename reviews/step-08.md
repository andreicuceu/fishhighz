# Step 08 handoff — revision 1

Implemented `IMPLEMENTATION_STEP.md`, revision 1 (2026-09-12). Ready for user
and independent review; acceptance and progression are not claimed. Step 09
has not started. Base commit remains
`0d69786a06d5d676564a51fad14a7156951c2228`, with preexisting Steps 03–07 still
uncommitted and preserved.

## Interfaces and scientific behavior

Added `fishhighz/models/kaiser.py`, `fishhighz/kernels/kaiser.py`,
`tests/test_kaiser.py`, `examples/builtin_forecast.py` and this report. Updated
README.md. No earlier scientific source/test, dependency declaration or planning
file changed. NumPy remains the only unconditional runtime dependency; the
existing optional Astropy/SciPy `templates` extra handles template preparation.

`Scaling(basis, **coordinates)` accepts exactly one of:

| Basis | Coordinates | ap | at | Vega convention |
| --- | --- | --- | --- | --- |
| ap_at | ap, at | ap | at | ap_at |
| alpha_phi | alpha, phi | alpha/sqrt(phi) | alpha*sqrt(phi) | phi_alpha |
| alpha_iso_epsilon | alpha_iso, epsilon | alpha_iso*(1+epsilon)^2 | alpha_iso/(1+epsilon) | aiso_epsilon |

`KaiserModel(template, fields, *, biases, betas, widths, f=None, local_names=(),
smooth=None, wiggle=None, z=None, growth=None)` prepares one bin's intrinsic
signal. `biases` covers all listed field IDs; `betas` covers exactly forest IDs.
Every supported field has explicit fixed `(parallel,transverse)` widths.
One shared f is required with galaxies and forbidden for forest-only models.
Each bias/beta/f/scaling setting is a fixed finite number or named local slot;
local_names gives the exact ordered argument vector. Each slot must be assigned
exactly once. Unknown, unused, duplicate and wrong-kind settings fail. Ties use
distinct local names mapped to one existing global ID with `BoundParameters`;
physical labels never infer ties. Widths/G cannot be free slots. Omitted scaling
blocks are fixed identity; separate bases, partial freedom and same-basis ties
are supported. No nonlinear cross-basis tie or aiso_aap is added.

The callable `model(theta_local,z,k,mu,pairs)` returns owned float64
`(node,pair)` power, in requested original field-index order. It works unchanged
inside P3DProvider/PreparedP3D for the full required covariance closure. Public
validation covers finite real local states, exact vector length, paired nonempty
k>0 and mu in [0,1], integer field indices and the fixed bin redshift. Finite
signed biases/betas/f have no implicit prior; galaxy b=0 remains valid.

Each component independently uses inverse Fourier mapping, its own transformed
mu and `Q=1/(ap*at^2)`. Derived scales and Q must be positive, finite and
representable; Q uses log-space exponentiation to avoid overflowing the volume
product. Exact identity/isotropic radial mapping avoids rejecting valid template
endpoints through angular roundoff. Both component domains are checked at every
parameter state, including derivative perturbations. Errors identify the
component, mapped range, available domain and scaling state; Step 06 adds
provider/parameter/stencil context. No clipping, extrapolation, regridding,
changed cuts, omitted nodes or step shrinking occurs.

Forest factors are `b*(1+beta*mu_component^2)` and galaxy factors are
`b+f*mu_component^2`, without division by bias. Signed pair products are used.
Wiggle damping is Gaussian in its transformed parallel/transverse coordinates,
with pair width squared equal to the mean of the two squared field widths.
Smooth power is not damped. Cross damping equals the geometric mean of the
auto factors at the same coordinates. Zero widths give unit damping, including
very large finite coordinates; a large finite exponent may underflow to zero.
The kernel multiplies coordinates by effective widths before squaring, avoiding
spurious infinity-times-zero. Unrepresentable squared widths/nonfinite exponents
fail. No f-dependent width prescription exists.

Final power is `G*(Q_s*B_i,s*B_j,s*P_s + Q_w*B_i,w*B_j,w*D_ij*P_w)` in
`(Mpc/h_fid)^3`, before response/noise. Template evaluation is at z_ref, and fixed
G multiplies the result once. Preparation reuses Step 07's z/G validation:
nonreference z requires explicit G; calls at another z fail. h_fid is the
unchanged template convention. No F_ZREF/SIGMA8 inference, growth solver,
geometry, noise preparation or covariance derivative is introduced.

## Numerical flow and ownership

Preparation owns read-only widths, their squares, fixed values, forest masks,
and integer destination/source slot maps. Kernels are `_resolve`, `_scales`,
`_coordinates`, `_field_factors`, `_damping` and `_assemble`; all take numeric
arrays/scalars/indices and contain no metadata, I/O or model dispatch.
Templates and field factors are evaluated once per component/query batch and
reused across requested pairs. Tests count two template calls for a three-pair
request. No per-node interpolation loop, node-by-knot matrix, coefficient
rebuilding or mutable parameter-dependent cache exists. A→B→A calls, read-only
noncontiguous inputs and concatenated slices preserve values and preparation.
Changing f leaves stored widths and G unchanged.

Production derivatives use the accepted Step 06 finite-difference engine with
explicit steps. Tied local slots move together under global perturbations;
there is no separate built-in differencing implementation or production analytic
Jacobian. Independent test oracles cover product-rule b/beta/f responses,
zero-bias crosses, component isotropic dilation including dD/dalpha, and the
identity basis chain matrices `[[1,-1/2],[1,1/2]]` and `[[1,2],[1,-1]]`.
Fixed-covariance Fisher matrices obey the corresponding basis transformation.

## Acceptance evidence

Evidence directory: `.validation/step08-r1-20260912T204258/` (`$RUN` below).

- **60 focused tests passed in 1.59 s.**
- Complete quick runner: **468 passed in 22.66 s**, Ruff lint and all 61
  formatting checks passed. The previous 408 tests remain retained.
- **15 bounded individual Vega-method comparisons passed**, maximum absolute
  discrepancy `8.881784197001252e-16`.
- Fresh installed-wheel built-in evaluation/derivative/Fisher probe and both
  runnable examples passed from `/tmp` with `-I`.
- Final report formatting, Git whitespace and preservation/hash checks passed.

Tests include all 15 original-order pairs of five fields (two forests, three
galaxies), negative forest bias and zero galaxy bias; reordered/reversed/subset
pairs; equivalent nonidentity bases and independent smooth/wiggle angular
mapping; isolated Fourier Q; fixed/partial/free/tied settings; unequal widths,
axis and strong-damping limits; unit conversion with G!=1; validation failures;
identity endpoints and both component domains; derivative-boundary failure and
successful padded-template evaluation with unchanged grid/cuts/weights/modes.

Scalar assembly tests use rtol around `3e-14` with `2e-13` absolute allowance
for signed polynomial spline arithmetic. Product-rule numerical derivatives use
`rtol=3e-11, atol=3e-11`; damped isotropic dilation uses
`rtol=3e-8, atol=2e-8` at h=2e-5. Basis Jacobians use
`rtol=2e-7, atol=2e-8` and transformed Fisher matrices
`rtol=3e-8, atol=2e-9` at h=1e-5. These finite-difference tolerances account for
second-order truncation and subtraction; no reference numerical tolerance changed.

The first focused run passed 53 tests and failed one newly introduced angular
roundoff assertion: equivalent-basis conversion differed by `1.6653e-16`
absolute (`5.2284e-16` relative) against rtol=5e-16. Its new-test threshold was
adjusted to 1e-15 to cover the multiple root/mapping operations. No model formula
was changed to force agreement. Formatting also removed initial multi-statement
lines. An additional robustness check motivated coordinate-times-width damping
arithmetic, preserving exact unit damping for zero widths even at k=1e200.
All final checks pass; initial/final logs are retained.

## Synthetic forecasts and convergence

`examples/builtin_forecast.py` uses a fixed 600-knot toy template over [0.01,0.7],
smooth `100/(1+2k)` plus `8*sin(110k)*exp(-(k/0.4)^2)`, and a fixed small grid
covering [0.04,0.25]. It explicitly checks grid/template h_fid=0.7. A forest/galaxy
selected subset `[gg,Fg]` requires unselected FF for covariance. Supplied fixed
response factors and positive diagonal noise are used; covariance is constructed
and factored once per mode, then reused across all derivative steps/convergence.
Snapshots verify fixed grid arrays, responses, noise and factors.

BAO mode fixes smooth coordinates and f; common AP+f ties both component blocks
and varies shared galaxy f. Explicit toy nuisance prior sigmas are bg=0.5,
bf=0.2, beta=0.5; bounds/steps are never priors. At derivative step scales
1,1/2,1/4 (base absolute steps 0.001), results are approximately:

| Mode | Target errors at finest step | Max Jacobian changes: h→h/2, h/2→h/4 | Max Fisher changes |
| --- | --- | --- | --- |
| BAO | ap 0.15176947, at 0.08821526 | 0.04308982, 0.01077242 | 0.01309456, 0.00327862 |
| AP+f | ap 0.11563680, at 0.07452078, f 0.27286151 | 0.04298211, 0.01074549 | 0.01488057, 0.00373089 |

Changes decrease by about four. The explicit h/2 versus h/4 derivative study
passes `atol=1e-5, rtol=2e-3`; coarse/fine Fisher agrees at
`rtol=2e-3, atol=1e-7`. Total ranks with stated nuisance priors are 5 and 6.
These are step-size checks on one fixed interpolant, not a new template-sampling
study; Step 07 separately assessed interpolation. Independent polynomial oracles
supply derivative accuracy checks. Exact installed results are in
`wheel-probe.json`. A separate forest-only regression has an exact zero global
f column, rank deficiency and infinite conditional f error with no prior hiding
it; joint marginalized errors correctly fail.

## Bounded reference comparison

`reference.py` hashes and parses the two read-only source files and extracts only
ScaleParameters.ap_at/phi_alpha/aiso_epsilon and
PowerSpectrum.compute_kaiser/compute_peak_nl using AST. It imports no Vega
package/dependency graph and runs no full model or fit. Fresh minimal power
objects prevent approximate cache comparisons from hiding changes. Both Gaussian
widths are always provided. Galaxy comparisons use beta=f/b only at nonzero b;
zero-bias behavior is independently tested above.

Three nonidentity scale settings across three bases, three Kaiser cases and
three Gaussian cases pass. `reference.json` retains exact inputs, source hashes
and errors. These checks validate only individual conventions, not Fourier Q,
inverse mapping, per-field width combinations or separate-component assembly;
those have independent tests. No new real FITS check, CAMB or lyaforecast run
was required, and existing template/reference evidence remains untouched.

## Commands, versions and exact wheel

All numerical commands used `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1` on the login node. From the component root:

```bash
.venv/bin/python -m pytest tests/test_kaiser.py -q > "$RUN/focused-final2.log" 2>&1
.venv/bin/python "$RUN/reference.py" > "$RUN/reference.log" 2>&1
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh > "$RUN/checks.log" 2>&1
.venv/bin/python -m build --wheel --outdir "$RUN/dist" > "$RUN/build.log" 2>&1
.venv/bin/python -m venv "$RUN/wheel-env"
"$RUN/wheel-env/bin/python" -I - <<'PY'
import importlib.util
assert importlib.util.find_spec('fishhighz') is None
PY
"$RUN/wheel-env/bin/python" -m pip install \
  "$RUN/dist/fishhighz-0.1.0.dev0-py3-none-any.whl[templates]" \
  > "$RUN/wheel-install.log" 2>&1
```

The first build approval request was interrupted by the workspace-path update;
it produced no artifact. The resumed request built one exact wheel in the fresh
directory. `fresh-environment.log` records no preexisting FishHighz. Approved
network access supplied isolated build/install dependencies; no shared environment
was modified. NERSC MUNGE socket messages appeared after successful development
subprocesses; exit statuses/checks passed and no Slurm action occurred.

With absolute paths, the probe ran from `/tmp`:

```bash
"$RUN/wheel-env/bin/python" -I "$RUN/probe.py" \
  "$RUN/dist/fishhighz-0.1.0.dev0-py3-none-any.whl" "$ROOT" "$RUN" \
  > "$RUN/wheel-probe.log"
```

It checks quiet/lazy imports, NumPy-only unconditional dependency metadata and
the retained template extra, known constant-model power 20, bias derivative 20
and Fisher 800, plus both copied standalone examples. All **21 package Python
modules** match source/wheel/installed bytes; METADATA/WHEEL also match. Both
example copies match their source. Origins are under the fresh environment's
site-packages, outside source imports; exact paths are in `wheel-probe.json`.
No second base-only environment was required for this step.

Exact wheel: `$RUN/dist/fishhighz-0.1.0.dev0-py3-none-any.whl`.
SHA-256: `a2224a6179d55578e22f95f28a1bbd37b8f0c48339c5cb0b1938a1f1f248de9f`.
Python 3.13.15, NumPy 2.5.3, SciPy 1.18.1, Astropy 8.0.1; development tools:
pytest 9.1.1, Ruff 0.16.7, build 1.6.1, pip 26.2.1. Package version remains
0.1.0.dev0. See development/probe version records.

## Preservation and limits

`before.json`/`preservation.json` confirm README.md is the only changed preexisting
maintained file. `source-sha256.txt`, `step.diff`, `git-status.txt`, `tracked.diff`
and `base-commit.txt` identify this step and the prior dirty tree. The report's
hash is external to avoid self-reference. No earlier scientific API, test, review,
planning document, dependency declaration or reference evidence changed.

No real lyaforecast quick capture or full suite was required or run. Synthetic
and individual-method evidence does not establish physical survey/cosmology
accuracy, realistic noise, a full Vega pipeline, autodiff, JIT performance or a
Python-version matrix. Arbitrary signed decompositions and separate scalings can
produce unphysical totals; Step 04 validates them without model repair. No commit,
push, planning revision, delegated agent or Step 09 work occurred. Stop for user
and independent review.

Changed implementation/documentation hashes:

```text
8c46fe3d60b74f3473f77cb858d5a0c7c09ecb8ea87cf8708832ba70f7b96b3b  fishhighz/models/kaiser.py
6e4f3c682a2068b10dca8b37d355fed7df0efa2ad00f4d57094a01fc8ae1cdac  fishhighz/kernels/kaiser.py
cf6207cafed2bb9e3dbeaac30a6313f87eb68e429417a6bae7fd4d709050f7df  tests/test_kaiser.py
0e11f63d05142b91df7b4e5d08068674011e38e65048f3245ff90376272889ec  examples/builtin_forecast.py
a9b19518a58272e85acff9c96ac620104e60fb15c00c55c1fd4cee1c895ab425  README.md
```
