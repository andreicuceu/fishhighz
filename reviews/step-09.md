# Step 09 handoff — revision 1

Implemented `IMPLEMENTATION_STEP.md`, revision 1 (2026-09-13), following the
user's explicit dispatch. Ready for user and independent review; acceptance and
progression are not claimed. Step 10 has not started. Base commit remains
`0d69786a06d5d676564a51fad14a7156951c2228`.

## Scope and preservation

Added `fishhighz/geometry.py`, `fishhighz/response.py`,
`fishhighz/kernels/response.py`, `fishhighz/models/p1d.py`, five dedicated test
files (`test_geometry.py`, `test_response.py`, `test_p1d.py`,
`test_survey_composition.py`, `test_survey_imports.py`),
`examples/survey_primitives.py`, and this report. Updated README.md and added the
optional `cosmology = ["astropy", "scipy"]` extra to pyproject.toml. NumPy remains
the sole unconditional dependency; the templates extra is retained.

Evidence directory: `.validation/step09-r1-20260913-implementation/`, called
`$RUN` below. `before.json` and `git-before.txt` snapshot the pre-dispatch dirty
state before source edits. Of 62 preexisting maintained files, 60 are unchanged;
only README.md and pyproject.toml changed. `preservation.json` verifies this.
Both old examples, all accepted scientific modules/tests, prior reports,
AGENTS.md, and IMPLEMENTATION_STEP.md are unchanged. The pre-planning snapshots
under `.validation/step09-plan-r1-20260913T002414Z/` were not overwritten.
Existing Steps 03–08 remain substantially uncommitted/untracked. No sibling
repository, reference artifact, planning/governance document, or license changed.

## Public behavior

`prepare_geometry(z_min, z_max, *, z_eval, area_deg2, h_fid, hubble,
transverse_distance, z_order)` samples independent ordinary H/D_M callables in
read-only 1D batches, at the quadrature nodes and the explicit evaluation
redshift. H is km/s/Mpc, transverse comoving D_M is Mpc. Positive interior
samples, explicit bin/evaluation/area/h/order constraints, exact output shapes,
and representable conversions/integrands/volumes are checked. The full-sky
public area endpoint is recognized before radians conversion can overshoot by
one ulp. No background consistency, flatness, interpolator, cosmology, or
quadrature convergence is inferred.

`BinGeometry` retains immutable bytes-backed arrays, scalar metadata, integrated
`Omega*h_fid^3*sum(w_z*c*D_M^2/H)`, evaluation H/D_M, `a_v`, `d_deg`, and
c=299792.458 km/s. It contains no background callable or cosmology object. Array
writeability cannot be re-enabled. Explicit bidirectional wavenumber/P1D/width
conversion helpers use fixed a_v; `mode_counts(geometry, grid)` checks exact
h_fid agreement and positive representable V*q_mode without changing any grid
array or introducing another mode factor.

`prepare_astropy_geometry(cosmology, z_min, z_max, *, z_eval, area_deg2, h_fid,
z_order)` lazily loads optional dependencies, requires a caller-created FLRW,
and converts H and transverse comoving distance using explicit Astropy units.
It delegates to the same core integrator, retains no Quantities, and does not
replace h_fid with cosmology.h. Missing optional support fails at adapter use
with a `fishhighz[cosmology]` install message.

`InstrumentResponse(pixel_width_velocity, gaussian_sigma_velocity)` requires
explicit nonnegative full pixel and Gaussian one-sigma widths in km/s.
`prepare_response(fields, k, mu, *, a_v, settings)` returns immutable
W(node,field), with a complete field-ID mapping and original field order.
`velocity_response(q, *, pixel_width_velocity, gaussian_sigma_velocity)` gives
the corresponding (node,) transfer. Observed q=k*mu/a_v is fixed; sinc lobes stay
signed, zero q and zero widths give exact identity, and strong Gaussian
attenuation may underflow to zero. Nonfinite/unrepresentable arguments fail.
The numeric kernel receives only arrays; name/unit validation stays on the host.

`pair_response(W, selection)` gathers W_i*W_j using original required-pair
indices. These products multiply intrinsic signal and Jacobians before selecting
means. Supplied noise gets no automatic response. Explicit unit-bearing helpers
separate wavelength full-pixel/sigma inputs, Gaussian FWHM, ordinary R_FWHM, and
labeled legacy c/R-as-sigma. README documents the narrow-width approximation and
explicit Ly-alpha wavelength 1215.67*(1+z_eval) Angstrom.

`default_p1d(theta_local, z, k_parallel_velocity)` retains the exact PD2013
constants and redshift-dependent stationary floor, in intrinsic km/s units.
`p1d_floor(z)` exposes that floor. The callable rejects nonempty local vectors,
invalid redshifts/shapes/NaNs/negative k and unrepresentable positive outputs;
it supports zero k, read-only arrays, slices, and preserved input order.
The default connects through `BoundParameters(registry, (), {})`; an external
P1D has its own explicit binding and identical downstream W²/a_v treatment.
No P3D call obtains P1D. Direct z>-1 mathematical validity does not change the
existing evaluator's nonnegative-z requirement or imply empirical fit accuracy.
The adapted routine identifies its lyaforecast/PD2013 provenance and GPLv3 terms;
both repository LICENSE texts were inspected, and no licensing change was made.

## Acceptance results

- Focused final checks: **87 passed in 1.76 s** (`focused-final2.log`).
- Complete quick runner: **555 passed in 20.99 s**, Ruff lint passed, all
  **73 files** formatted (`checks-final.log`). All 468 previous tests remain.
- **33 bounded legacy-method comparisons passed** (`reference-final.log`,
  `reference.json`), using matched backgrounds/constants where appropriate.
- Installed final-wheel geometry, response, P1D, amplitude Fisher, curved Astropy,
  all three examples, and blocked-optionals subprocess passed (`wheel-probe.log`,
  `wheel-probe.json`, `blocked-optional.log`). All **25 package Python modules**
  match source/wheel/installed bytes; METADATA/WHEEL match too.
- Final Git whitespace, preservation, source hashes and report checks passed;
  detailed manifests and dirty state are saved alongside the report evidence.

Ordinary tests cover multiple redshifts/h values, analytic constant-H shell
volume, h scaling, reciprocal conversions and invariant k*P1D, custom and
Gauss–Legendre modes, full-sky/zero-lower-bin endpoints, overflow/underflow
rejection, owned immutable preparation, stopped background calls, and optional
flat/open/closed cosmologies with h_fid unequal to cosmology.h. Response tests
include independent scalar sin(x)/x and Gaussian oracles, zero/small arguments,
a sinc zero/negative lobe, very strong attenuation, explicit galaxy identity,
unequal forest settings, wavelength/velocity/FWHM equivalence, pair order and
noncontiguous slices. P1D tests use the original scalar power-law expression at
z=2,3,4, the pivot, each floor and both sides, positivity, plateau, slicing,
independent substitution and actual first/second derivative behavior.

Conversion/response/P1D comparisons use rtol at most 5e-14 for ordinary scalar
cases. The sinc zero uses atol=3e-16 to cover normalized-sinc argument roundoff.
P1D finite differences above the floor use rtol=3e-9 at relative step 1e-5;
join offsets 1e-2,1e-3,1e-4 give quadratically decreasing power departures and
linearly decreasing dimensionless forward slopes (~9.90e-4, 9.99e-5, 1.00e-5).
Right curvature approaches -0.2 P/k_floor² while left curvature is zero; the
finite-offset check allows 0.4% truncation at offset 1e-3. No old tolerance changed.

Volume convergence on [2,3], with explicit orders 4/8/16, is recorded in
`convergence.json`. Astropy tests compare to its independent comoving-volume
shell, including curvature, rather than repeating the production quadrature.

| Background | Relative error, order 4 | Order 8 | Order 16 |
| --- | --- | --- | --- |
| Einstein–de Sitter analytic shell | 1.1143e-9 | 8.8818e-16 | 6.6613e-16 |
| FlatLambdaCDM, H0=68, Om0=0.3 | 2.1146e-9 | 4.4409e-16 | 7.7716e-16 |
| LambdaCDM, Ok0=+0.15 | 1.3465e-9 | 4.7740e-15 | 5.2180e-15 |
| LambdaCDM, Ok0=-0.15 | 3.4767e-9 | 4.4409e-16 | 4.4409e-16 |

The finest rules meet 1e-10 (EdS) and 1e-9 (Astropy) gates; order-8/16 differences
are at roundoff. Volume convergence is separate from fixed Fourier quadrature.

The new example has two distinct forests and a galaxy, unequal instrument
responses, signed forest/galaxy cross-power, a selected subset requiring full
covariance closure, and positive fixed synthetic diagonal noise. An independent
field-matrix oracle checks covariance and analytic amplitude information.
Covariance is factored once and reused through three numerical derivative step
sizes, with snapshots of geometry/grid/W/modes/factors. The amplitude Fisher is
approximately **30334.17490694**, with error **0.00574161287**. These are synthetic
normalization checks, not physical forecast predictions. Separate one-field tests
verify Var(P)=2*T²/N and area doubling: modes/Fisher double, covariance halves,
and errors scale by 1/sqrt(2). Permutations/slices retain observable order and
information. P3D parameter changes leave P1D and fixed preparation untouched.

## Bounded reference interpretation

`reference.py` hashes live source and extracts only named functions/methods via
AST into minimal fake objects. `reference.json` records the reference Git
revision/status, absolute source paths/hashes, all extracted entry points,
inputs/results/errors, constants and limitations. No CAMB initialization, SNR
file, survey construction, Vega model, or actual legacy forecast was executed.
The largest absolute comparison discrepancy is 1.863e-9 for a large volume;
all direct comparisons meet rtol=5e-13.

Checks cover standalone `P1D_z_kms_PD2013`,
`PowerSpectrum.compute_p1d_kms/compute_p1d_hmpc`,
`Spectrograph.smooth_kernel_kms`, `CosmoCamb`'s four conversion methods, and
`Covariance._get_redshift_depth/_get_survey_volume/_get_num_modes/_get_pix_kms/
_get_res_kms`. Nonzero response probes match intrinsic/smoothed/converted P1D.
The reference velocity smoothing has zero-mode 0/0; the new analytic limit is
one. The reference comoving wrapper clamps k*pixel to 1e-5, producing a tiny
zero-mode attenuation; the new response has no such clamp. These differences
are disclosed in `reference.json`, not hidden by modifying production physics.

With identical EdS H/D_M and the actual legacy geometric redshift centre, the
integrated versus legacy volume differences after matching c are:

| Bin | Relative integrated/legacy - 1, same c |
| --- | --- |
| [2,3] | -0.0028118605500 |
| [2.395,2.405] | -2.9749747e-7 |
| [2.39995,2.40005] | -3.1468605e-11 |

The production/reference light-speed ratio 299792.458/299800 is tracked
separately. Broad-bin equality is not asserted; the thin-bin limit converges.
Direct velocity sigma is unchanged by c; the labeled R helper explicitly differs
from the reference by the light-speed ratio. No saved end-result bundle is used
as evidence for these new primitives.

## Commands, environment and exact artifact

Commands run from the package root unless noted. Numerical commands used
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` on the login node.
`$ROOT` denotes the absolute component root; `$RUN` is the absolute evidence path.

```bash
.venv/bin/python -m pip install -e '.[dev,templates,cosmology]'
.venv/bin/python -m pytest tests/test_geometry.py tests/test_response.py \
  tests/test_p1d.py tests/test_survey_composition.py tests/test_survey_imports.py -q
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh
.venv/bin/python "$RUN/reference.py"
.venv/bin/python "$RUN/convergence.py"
.venv/bin/python -m build --wheel --outdir "$RUN/dist-final"
.venv/bin/python -m venv "$RUN/wheel-env"
"$RUN/wheel-env/bin/python" -I -c \
  "import importlib.util; assert importlib.util.find_spec('fishhighz') is None"
```

The one fresh environment initially installed the first exact wheel from `dist/`
with `[templates,cosmology]`, downloading its optional dependencies. A late
full-sky boundary regression then required a corrective second artifact in
`dist-final/`; both build logs/artifacts are preserved. The final wheel replaced
the first in that same isolated environment offline:

```bash
"$RUN/wheel-env/bin/python" -m pip install --no-index --no-deps --force-reinstall \
  "$RUN/dist-final/fishhighz-0.1.0.dev0-py3-none-any.whl[templates,cosmology]"
# From /tmp, using absolute paths:
"$RUN/wheel-env/bin/python" -I "$RUN/probe.py" \
  "$RUN/dist-final/fishhighz-0.1.0.dev0-py3-none-any.whl" "$ROOT" "$RUN"
git diff --check
```

This corrective rebuild is the only deviation from the single-build intention;
there is one isolated environment, no wheel/sdist/version matrix. Initial absence
of FishHighz is in `fresh-environment.log`; final byte/metadata identity prevents
a stale same-version install from passing. Example copies are byte-compared and
executed from /tmp with -I. Installed module origins are under the fresh
`wheel-env/lib/python3.13/site-packages/fishhighz/`; the complete resolved paths
are in `wheel-probe.json`. Quiet/lazy imports and optional metadata are checked
separately from a new subprocess actively blocking Astropy/SciPy while running
the entire supplied-background/external composition example. Adapter use then
fails with the actionable extra message.

Final wheel SHA-256:
`3f6802f8f745be8555ba5ba34300d83c841bca4d0d3322321ebbba9fe8b6777f`.
Python 3.13.15, FishHighz 0.1.0.dev0, NumPy 2.5.3, Astropy 8.0.1,
SciPy 1.18.1, pytest 9.1.1, Ruff 0.16.7, build 1.6.1, pip 26.2.1.
`versions.json` and the installed probe record versions independently.

## Failures, repairs and limits

The initial snapshot command used unavailable `python`; it was rerun with
`.venv/bin/python` before any source edit. Early new tests incorrectly reversed
extreme values for inverse conversions, and the sinc-zero oracle needed 3e-16
absolute allowance for normalized-sinc arithmetic. Both test-only corrections
are retained in logs. Ruff found two assigned lambdas in a new test; they became
named functions. Initial focused tests then passed, and the first quick run had
554 passing tests. A later full-sky test found the real one-ulp endpoint rejection,
which was fixed in geometry and covered by regression. Its overflow test fixture
also needed larger k to actually overflow modes. Final focused and quick checks
pass after these repairs. No accepted scientific tolerances or modules changed.

The initial no-build-isolation editable install lacked setuptools; the isolated
retry could not resolve PyPI under the sandbox. Approved network access completed
the development install, builds and fresh optional dependency install. Logs are
retained as `dev-install*.log`, `build*.log`, `wheel-install*.log`. Only the local
environments were modified. Harmless NERSC MUNGE socket messages appeared after
successful development checks; commands returned success, and no Slurm action
occurred.

No real legacy capture or full suite was requested or run under the quick-only
policy; this is not a failed/outstanding acceptance gate. No source-density/SNR
loading, forest weights, effective density, noise estimation, overlap model,
survey parser, growth solver, intrinsic-model change, serialization, orchestration,
JIT, or performance study is included. No commit, push, delegation, planning
revision or Step 10 implementation occurred. Stop for user and independent review.
