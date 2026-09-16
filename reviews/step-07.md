# Step 07 handoff — revision 1

Implemented `IMPLEMENTATION_STEP.md`, revision 1 (2026-09-12). Ready for user
and independent review; acceptance and progression are not claimed. Step 08
has not started. Base commit remains
`0d69786a06d5d676564a51fad14a7156951c2228`; preexisting uncommitted Steps 03–06,
planning documents, reviews and reference evidence are preserved.

## Delivered interfaces and dependency choice

Added `fishhighz/models/templates.py`, `fishhighz/kernels/templates.py` and
`tests/test_templates.py`. Updated README.md and pyproject.toml; no earlier
scientific source/test or planning/governance file changed. NumPy remains the
only unconditional runtime dependency. The confirmed `templates` extra contains
Astropy and SciPy. Development instructions use `.[dev,templates]`; tests import
the extras directly, so absent dependencies cannot silently skip acceptance.

| Interface | Contract |
| --- | --- |
| `prepare_template(k, pk, pksb, *, z_ref, h_template, h_fid, metadata=None)` | Validated in-memory preparation from explicit source units and metadata; returns `PowerTemplate`. Requires SciPy only for coefficient construction. |
| `load_template(path, *, h_fid, z_ref=None, h_template=None)` | Reads a local Vega-format PK table, resolves required metadata, and prepares the same numerical record. Astropy is imported only here. |
| `PowerTemplate.evaluate(k, *, derivative=0, z=None, growth=None)` | Owned C-order float64 `(query,2)` components, labeled `("smooth","wiggle")`; derivative 0 gives power, 1 gives dP/dk. |
| `PowerTemplate.domain` | Closed converted k interval in h_fid/Mpc; not a forecast scale selection. |
| `_intervals(log_knots, log_query)` | Array-only searchsorted interval lookup; last endpoint uses the last interval. |
| `_evaluate(log_knots, coefficients, k_query, derivative)` | Array-only Horner evaluation/differentiation in log k, including the derivative's 1/k chain factor. |

`PowerTemplate` is a prepared record; callers use the validated factories.
`k_file`, `pk_file`, `pksb_file`, converted `k`, `full`, `components`, `log_k`
and `coefficients` are owned read-only float64 arrays. Coefficients have shape
`(4,interval,component)` in descending powers of `ln(k)-log_k[interval]`.
Metadata is a copied immutable mapping of ordinary scalar/string values. Source
path/hash, z_ref and both h conventions remain available. No FITS handle, memory
map, Astropy object or SciPy spline object remains in the record.

The implementation follows the format written by the read-only
`../vega/bin/make_template.py`; no generator/scientific code or cosmological
asset was copied, changed, executed or bundled. No Vega/CAMB runtime import,
sideband fit or template regeneration occurs.

## FITS, units and amplitude

The loader hashes and parses the same single byte snapshot. It finds the unique
named `PK` binary table independent of extension position and rejects ambiguous,
missing or wrong-type tables. K/PK/PKSB must be scalar real finite matching
columns with at least four nodes. K must be positive and strictly increasing;
no sorting, duplicate removal or resampling occurs. Unrelated extensions/columns
are ignored. Negative PKSB samples are legal; `P_nw=PKSB` and `P_w=PK-PKSB`
are prepared once and retain their signs. Full power is not treated as wiggles.

Required ZREF and H0 come from the PK table, not the primary header. Missing
values require explicit z_ref/h_template. Provided/header values must agree
within `64*eps64*max(1,abs(header),abs(provided))`, comparing h after H0/100;
consistent table values take precedence. z_ref>=0 and h/H0>0 must be finite.
Original scalar header metadata, including F_ZREF/SIGMA8/cosmology, is retained
without scientific inference. COMMENT/HISTORY are omitted from the provenance
mapping. The in-memory constructor takes required numerical metadata explicitly;
its optional mapping is provenance only.

Absent TUNIT cards mean the documented Vega h-scaled convention. The bounded,
case-sensitive spelling table ignores spaces and normalizes `**` to `^`:

- K: `h/Mpc`, `h Mpc^-1`, `h Mpc-1`.
- PK/PKSB: `(Mpc/h)^3`, `(Mpc/h)3`, `Mpc^3/h^3`, `Mpc3/h3`,
  `Mpc^3 h^-3`, `Mpc3 h-3`.

Other spellings, including physical `1/Mpc` and `Mpc^3`, fail with column/unit
context. This deliberately does not implement a general unit converter.
Explicit h_fid applies `k_fid=k_file*h_template/h_fid` and
`P_fid=P_file*(h_fid/h_template)^3` once to the signed components. Source samples
remain auditable. This is a unit change for one cosmology, not AP rescaling.
Nonfinite/unrepresentable conversion, vanished nonzero powers, non-increasing
converted knots or indistinguishable log knots fail.

Evaluation defaults to exactly z_ref and G=1. Another finite nonnegative redshift
requires explicit positive finite `growth=G=[D(z)/D(z_ref)]^2`. G scales powers
and derivatives exactly once. At z_ref, explicit G must agree with 1 at the same
64-eps roundoff tolerance. Optional cosmology/growth-rate metadata never supplies
a growth function. Metadata and coefficients stay unchanged across requests.

## Interpolation, domain and numerical evidence

Both signed linear amplitudes use the same explicitly not-a-knot cubic spline
in ln(k_fid). SciPy prepares its coefficients once. Repeated NumPy evaluation
uses interval indices and a few query-by-component arrays, not a dense
query-by-knot matrix. No preparation library or I/O is invoked. dP/dk has units
`(Mpc/h_fid)^4`, including 1/k, while power has `(Mpc/h_fid)^3`.

Unsorted/repeated positive finite 1D queries and read-only/noncontiguous slices
retain their order. Exact converted endpoints are accepted. A one-ULP excursion
fails before logarithms and reports requested/available domains; no extrapolation,
clipping, padding or new Fourier modes occur. Coefficient/evaluation overflow
raises contextual errors. There is no compilation or speedup claim.

Analytic irregular-log-k cubic tests cover signs, exact endpoints, knots,
one-ULP knot neighbors, off-knot values and derivative continuity. Tolerances
are `rtol=2e-13`, `atol=2e-14` for power and `atol=3e-13` for dP/dk, covering
float64 spline preparation/Horner arithmetic and the small-k chain factor.
An independent SciPy off-knot check covers nonuniform log spacing, signed
components and component-sum equality with a separately prepared full-power
spline. Physical-unit checks compare power/h^3 and derivatives/h^4 between
h_template=0.7 and h_fid=0.5.

The oscillatory toy is smooth `1000/(1+5k)` plus wiggle
`50*exp(-(k/0.4)^2)*sin(110k)` over [0.02,0.5]. Analytic values/derivatives are
compared at 1,501 off-knot queries in [0.020013,0.499987]. Maximum absolute
errors (meaningful through wiggle zero crossings) were:

| Knots | Smooth power | Wiggle power | Smooth dP/dk | Wiggle dP/dk |
| --- | --- | --- | --- | --- |
| 48 | 7.9354e-5 | 15.8686 | 1.4701e-2 | 2991.30 |
| 96 | 4.7453e-6 | 0.884331 | 1.7741e-3 | 168.793 |
| 192 | 2.8928e-7 | 0.146098 | 2.1716e-4 | 114.506 |
| 384 | 1.7860e-8 | 0.0127082 | 2.6732e-5 | 19.1581 |
| 768 | 1.0942e-9 | 0.000827743 | 3.2853e-6 | 2.48694 |

The first focused run had 75 passing tests and one inadequate sampling-study
fixture: 48/96/192 knots did not satisfy its derivative target. The final
regression uses 48/192/768 knots, requires decreasing errors, explicitly rejects
coarse-grid adequacy, and tightens final limits to smooth/wiggle power
`2e-9/1e-3` and derivative `4e-6/3`. No numerical implementation or established
tolerance was weakened. The full measured sequence is `sampling-errors.json`.
This is template sampling convergence, not IntegrationGrid or parameter-step
convergence, and not a physical forecast accuracy claim.

## Acceptance checks and read-only references

Evidence directory: `.validation/step07-r1-20260912T194343/` (`$RUN` below).
Focused template tests: **76 passed in 1.80 s**. The complete quick runner
passed **408 tests in 19.87 s**, Ruff lint and 55 formatting checks, with no skips.
The earlier 332 tests remain retained. Coverage includes FITS schema/types,
snapshots, units/metadata failures, explicit growth, signed components,
polynomial/oscillatory accuracy, strict domains, owned arrays, slices, lazy
imports and missing-extra diagnostics. Evaluation after closing/changing the
source file and disabling FITS/SciPy calls passes.

Both specified local Vega reference files were loaded read-only with h_fid=0.6736:

| File | z_ref | Nodes | Negative PKSB samples | SHA-256 |
| --- | --- | --- | --- | --- |
| Planck18_z_2.406.fits | 2.406 | 814 | 121 | b4a73103e1105b7f9cdb59bbe8133ff0f752b05bc27580dd526f80b9c2fdc26a |
| DESI-2024_z_2.33.fits | 2.33 | 814 | 121 | e3413abafaa419856fbfefab553c4132185285c4c95150d0cc3cd19d56080cdd |

Both have H0=67.36 and converted domain [0.0001,1152.499999999999]. Hashes match
the plan and remain unchanged after reading. Knot reconstruction of full power
has maximum absolute error `1.1368683772161603e-13`; all 813 off-knot midpoint
power/derivative evaluations are finite. See `reference.py`, `reference.json`
and `reference.log`. This wide domain is not a validated forecast scale range.
No Vega import, CAMB run or reference forecast occurred.

## Commands, versions and fresh wheel

All numerical commands use `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1` on the login node. From the component root:

```bash
.venv/bin/python -m pip install -e '.[dev,templates]' > "$RUN/dev-install.log" 2>&1
.venv/bin/python -m pytest tests/test_templates.py -q > "$RUN/focused-final.log" 2>&1
.venv/bin/python "$RUN/reference.py" > "$RUN/reference.json"
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh > "$RUN/checks.log" 2>&1
.venv/bin/python -m build --wheel --outdir "$RUN/dist" > "$RUN/build.log" 2>&1
.venv/bin/python -m venv "$RUN/wheel-env"
"$RUN/wheel-env/bin/python" -m pip install numpy==2.5.3 > "$RUN/numpy-install.log" 2>&1
"$RUN/wheel-env/bin/python" -I - <<'PY'
import importlib.util
for name in ('fishhighz', 'scipy', 'astropy'):
    assert importlib.util.find_spec(name) is None
PY
"$RUN/wheel-env/bin/python" -m pip install --no-index --no-deps \
  "$RUN/dist/fishhighz-0.1.0.dev0-py3-none-any.whl" > "$RUN/wheel-install.log" 2>&1
```

The fresh-environment assertion is logged in `fresh-environment.log`. With
absolute `$RUN`/`$ROOT`, the following ran from `/tmp` with `-I`:

```bash
"$RUN/wheel-env/bin/python" -I "$RUN/base-probe.py" \
  "$RUN/dist/fishhighz-0.1.0.dev0-py3-none-any.whl" "$ROOT" \
  "$RUN/external_forecast.py" > "$RUN/base-probe.json"
```

The base probe verifies Astropy/SciPy are actually absent, quiet lazy module
imports, actionable missing-extra errors, tied external derivatives/Fisher,
the Step 06 duplicate-stencil repair and retained standalone forecast example.
It passed before installing extras for the same exact artifact:

```bash
"$RUN/wheel-env/bin/python" -m pip install \
  "$RUN/dist/fishhighz-0.1.0.dev0-py3-none-any.whl[templates]" \
  > "$RUN/templates-install.log" 2>&1
# Again from /tmp, with absolute paths:
"$RUN/wheel-env/bin/python" -I "$RUN/templates-probe.py" \
  "$RUN/dist/fishhighz-0.1.0.dev0-py3-none-any.whl" "$ROOT" "$RUN" \
  > "$RUN/templates-probe.log"
```

The installed generated-FITS probe passed with maximum cubic power error
`1.2212453270876722e-15` and derivative error `1.4210854715202004e-13`. It checks
h-unit conversion, explicit G, one-ULP domain rejection, slicing and evaluation
with preparation-library calls disabled. All **19 Python package modules**
match source/wheel/installed bytes both before and after extras installation;
METADATA/WHEEL and optional-dependency metadata match. Origins are under the
fresh environment's site-packages, not the checkout. Home and resolved CFS paths
identify the same authorized package. See both probe JSON records.

One wheel was built in the fresh directory:
`$RUN/dist/fishhighz-0.1.0.dev0-py3-none-any.whl`.
SHA-256: `1af6778f6b68439001e3665f1b09eb69d0237b80db7fbbcb0000dc4232b72d57`.
Development and installed template environments use Python 3.13.15, NumPy 2.5.3,
SciPy 1.18.1, Astropy 8.0.1, pyerfa 2.0.1.5 and astropy-iers-data
0.2026.9.7.0.56.14. Development tools: pytest 9.1.1, Ruff 0.16.7, build 1.6.1,
pip 26.2.1. See `dev-versions.json`. Approved network access supplied build/optional
dependencies; no shared environment was modified.

NERSC MUNGE socket messages appeared after successful development subprocesses;
exit codes and pytest/check outcomes remained successful. Reference stdout was
preserved in `reference.log`, with its JSON payload separated into `reference.json`.
No scheduler action was requested. The initial sampling failure described above
was corrected; no required acceptance failure remains.

## Preservation and limits

`before.json` and `preservation.json` confirm README.md and pyproject.toml are
the only changed preexisting maintained files. All earlier scientific APIs,
tests, reports and governance files remain byte-identical. New files and final
hashes are in `source-sha256.txt`; `step.diff` isolates these changes from the
prior dirty tree. `git-status.txt`, `tracked.diff` and `base-commit.txt` record
Git state. The report hash is stored externally to avoid self-reference.
`git diff --check` and the final report formatting check pass.

No real lyaforecast quick capture or full suite was required or run; full mode
was not requested. No scientific reference evidence was replaced. This step
validates template format, unit/amplitude handling and spline behavior, not
physical cosmology, survey convergence, Kaiser/BAO physics, a Python-version
matrix, autodiff or JIT performance. No commit, push, planning revision or
Step 08 work was performed. Stop for user and independent review.

Changed implementation/documentation hashes:

```text
3afc1943ddece810da0153af15b1d32056ab48d1cf20ae782e2d195f4376763a  fishhighz/models/templates.py
89c4fb990ca646e30cd972bc1d36ebd71372fde70f9e402f821c935e5de96f9f  fishhighz/kernels/templates.py
e2149539c3e35ed7ddb57e2f4c72f79403cf1aa1d19bd540359b57f2d0f54fa6  tests/test_templates.py
9dd8597404972d524c1be818a762578c003be8f62d47260d691c2d1699e4a16d  README.md
48fe9cc190ba4f800e9d9c3f04601b0ba235373c479e161568a337be4359e9ab  pyproject.toml
```
