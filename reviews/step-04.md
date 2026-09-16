# Step 04 handoff — revision 1

Implemented `IMPLEMENTATION_STEP.md`, revision 1 (2026-09-11). Ready for user
review; acceptance and progression are not claimed. Step 05 has not started.
Base commit remains `0d69786a06d5d676564a51fad14a7156951c2228`, with the accepted
Step 03 implementation and planning changes already uncommitted at dispatch.

## Interfaces and numerical contract

Added `fishhighz/covariance.py`, `fishhighz/kernels/covariance.py`, the minimal
kernel initializer, and `tests/test_covariance.py`. Updated README.md with a
cross-only covariance example. No existing Python source/test, dependency,
planning document, reference tool, or historical report was edited.

- `combine_observed_power(signal, noise)` adds supplied observed clustering
  (already containing field responses) and noise in exactly matching
  `(n_node,n_required)` arrays. No broadcasting, smoothing, noise generation,
  off-diagonal assumptions, or mean subtraction occurs. Components may have
  signed entries; physical validation applies to the total at the next boundary.
- `gaussian_covariance(total_power, mode_counts, selection)` accepts a prepared
  `PairSelection`, total power in its `required_pairs` order, and positive finite
  `(n_node,)` counts, including fractional values. It returns owned C-contiguous
  float64 `(n_node,n_selected,n_selected)` blocks in selected-pair order.
- `_gaussian_covariance_kernel(total_power, mode_counts, im, jn, in_, jm, out)`
  is a private numeric-only function. It overwrites preallocated float64 `out`
  and returns None. Inputs must already have valid shapes, dtypes, and lookup
  indices; output must not alias inputs. It computes one triangle and mirrors it.

The kernel evaluates `(T_im*T_jn + T_in*T_jm)/mode_counts` using the accepted
four lookup tables. For an auto spectrum, variance is `2*T_ii**2/mode_counts`.
With `mode_counts = V_fid*grid.q_mode`, both conjugate hemispheres are already
counted on mu in [0,1]; no additional weight, volume, or factor of two is applied.
Input power uses `(Mpc/h_fid)^3`, covariance `(Mpc/h_fid)^6`, and counts are
dimensionless. Nodes retain Step 03's C-order flattening, with mu fastest. Blocks
never couple nodes or redshift bins. Inputs and selection arrays are unchanged.

Fixed fiducial covariance is implemented from supplied arrays: there is no theta
argument or automatic covariance/weight update. Known noise contributes to total
power even when it is subtracted from the mean. Survey preparation and Fisher
calculations remain unimplemented; covariance changes alone add no information.

## Physical validation and bounded memory

Validation explicitly verifies that required pairs contain every canonical pair
among active selected fields. Wholly unused fields do not enter the dense check.
One active field matrix is reconstructed at a time. Negative auto power is
rejected exactly. An exactly zero auto requires an exactly zero row/column,
including arbitrarily small nonzero cross power; otherwise the node and field
are reported. Zero fields remain in the returned covariance, without a noise
floor or automatic removal.

For positive autos, divide each axis by its auto-power square root to form the
normalized correlation matrix R. This avoids forming potentially unrepresentable
products of auto powers and prevents large amplitudes from hiding invalid weak
correlations. Let n be the number of positive-auto fields and eps64 the float64
machine epsilon. The exact dimensionless tolerances are:

- `abs(R_ij) <= 1 + 64*eps64*n`;
- `lambda_min(R) >= -64*eps64*n*max(1,max(abs(eigvalsh(R))))`.

Both checks permit only roundoff-sized violations. Full eigensystem validation
rejects indefinite multi-field matrices even when all pairwise bounds pass.
Diagnostics identify the first failing node, involved field IDs, and the
correlation magnitude or minimum eigenvalue and allowed tolerance. Valid
semidefinite totals, zero fields, and all-zero matrices pass without clipping,
jitter, Cholesky, inverses, or regularization. The selected covariance need not be
singular just because total field power is singular; solvability is later work.

Beyond its output, the kernel allocates one node-vector and uses views plus
in-place ufunc outputs. Inspection confirms no full covariance-shaped temporary,
name/configuration work, validation, or eigensolve inside the kernel. Public
validation normalizes input arrays once; physical checks and output-finiteness
checks operate one node at a time. Independent node slices concatenate to the
whole-batch result. Nonfinite addition or covariance arithmetic raises ValueError
with location/context. Intermediate overflow is an error even if another
arithmetic ordering could have avoided it; no units or amplitudes are changed.

## Checks and artifacts

Commands ran from the FishHighz root except for the installed probe. The existing
local development environment contained Python 3.13.15, NumPy 2.5.3, pytest
9.1.1, Ruff 0.16.7, build 1.6.1, and pip 26.2.1. NumPy remains the sole runtime
scientific dependency. Full version metadata is in
`.validation/step04/dev-versions.json`.

| Command | Actual result |
| --- | --- |
| `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest tests/test_covariance.py` | Initial 39 focused tests passed in 1.09 s; `focused.log`. |
| `PATH="$PWD/.venv/bin:$PATH" scripts/check.sh` | Final **186 tests passed in 25.67 s**, including 40 new covariance tests and all 146 accepted tests. Ruff lint passed; 33 files already formatted. `checks.log`. |
| `git diff --check` | Passed. |
| `.venv/bin/python -m build --wheel --outdir .validation/step04/dist` | One wheel built in a fresh output directory with isolated dependencies; `build.log`. |
| `.venv/bin/python -m venv .validation/step04/wheel-env` | New isolated installation environment created. |
| `.validation/step04/wheel-env/bin/python -m pip install .validation/step04/dist/fishhighz-0.1.0.dev0-py3-none-any.whl` | Installed that exact artifact and declared NumPy 2.5.3; `wheel-install.log`. |

All logs above are under `.validation/step04/`. The runner sets the three
numerical thread limits to one. Build/install used approved network-enabled
execution because sandbox PyPI DNS had failed in the previous step. No failing
test or unresolved check remains in this assignment.

Tests cover analytic one/two-field matrices (negative and zero cross powers),
noise-only inputs, explicit cross noise, fractional counts, scaling identities,
grid/volume normalization, five fields/all 15 spectra, permuted selected principal
submatrices, unused fields, and column errors. They cover singular and zero
limits, full three-field indefiniteness, valid/invalid disparate scales (including
cross-only autos of 1e-300 and 1e300), and preservation of roundoff-sized values
without regularization. Deterministic `L@L.T` cases match an independent dense
fourth-moment Wick tensor and satisfy symmetry/PSD within stated test roundoff.
Shape/type/finiteness/overflow errors, owned contiguous outputs, noncontiguous
inputs, unchanged input/selection data, kernel parity, and sliced-batch equality
are exercised. The retained reference-tool tests remain synthetic.

### Fresh installed-wheel evidence

Exact wheel: `.validation/step04/dist/fishhighz-0.1.0.dev0-py3-none-any.whl`.
SHA-256: `0453ee2432d9fda20e99ca4a23a890847ede10738ed0696d765ba7848daf53e5`.
No earlier same-version environment was reused. The probe ran from `/tmp` with
`-I` and all three thread limits equal to one. Its invocation was equivalent to
this command with ROOT set to the absolute FishHighz workspace path:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$ROOT/.validation/step04/wheel-env/bin/python" -I \
  "$ROOT/.validation/step04/probe.py" \
  "$ROOT/.validation/step04/dist/fishhighz-0.1.0.dev0-py3-none-any.whl" \
  "$ROOT/.validation/step04/readme-example.py"
```

`wheel-probe.json` records exact absolute paths and results. Both
`fishhighz.covariance` and `fishhighz.kernels.covariance` resolved beneath
`.validation/step04/wheel-env/lib/python3.13/site-packages/fishhighz/`.
All ten installed Python modules matched the selected wheel bytes. Quiet package
import without loading NumPy or neighboring packages passed, as did the README
example, explicit analytic two-field matrix, direct kernel parity, singular
input, and invalid-power diagnostics. Probe Python/NumPy versions matched the
development environment. The source home path resolves to CFS in artifact paths;
writes stayed within the authorized FishHighz package.

## Git identity and review limits

`before.json` captured preexisting source/test/report and governance hashes;
`preservation.json` confirms README.md was the only changed preexisting file in
that snapshot. Step 03's uncommitted work was preserved. `git-status.txt` and
`tracked.diff` retain the full final Git state. The inline hashes below identify
this assignment's implementation, tests, and README; `source-sha256.txt` also
includes this handoff and the current governance inputs. The report's own hash
is kept separately to avoid self-reference.

No real lyaforecast quick capture or full suite ran; neither is required here,
and a full run was not requested. Existing reference evidence was preserved.
Validation is synthetic, on Python 3.13, with no claim of survey-forecast
convergence, invertibility, performance measurement, or JIT speedup. No commit,
push, planning revision, or Step 05 implementation was performed. Stop for user
review and independently requested review before progression.

```text
f2dbfd11da26911963c84abca5d5fe5a1927e91651238cb9a7c2c97a69f5fbee  README.md
95be1ca5803b2e5b2535ea3e354459ed6ee449ef97fadfca0e55d5c97341e180  fishhighz/covariance.py
c406667a62606e02e368e0c2901b6db1fd3d29a574bd5695979ac30b8e9a0301  fishhighz/kernels/__init__.py
325e919bcc96bccd19ab6798db588924895a2a25d6dac3eb297b1ebf695194f4  fishhighz/kernels/covariance.py
cc266ced0043f84112f49ea0e85c7f6688d10aa4c62e4337fd4819c1c2f070ce  tests/test_covariance.py
```
