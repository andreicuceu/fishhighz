# Step 09 revision 1 independent review

Date: 2026-09-13.
Outcome: review passed; no required changes. The user reported no comments.
Ready for the user's acceptance/progression decision; Step 10 is not planned.
No production source, ordinary test, example, dependency, or implementation
handoff was changed by this review. No repair revision is needed.

## Assessment

Reviewed the exact revision-1 assignment, handoff, all four new source modules,
five new test files, example, documentation/dependency changes, and reference/
installation evidence. The implementation satisfies the bounded assignment.

Geometry integrates the transverse-distance volume element with an explicit
redshift rule and common bin/evaluation settings. The optional Astropy adapter
uses explicit units and preserves h_fid independently of cosmology.h. Owned
immutable arrays retain no background callable/cosmology state. Mode counts
reuse V*q_mode without changing the fixed Fourier grid or normalization.

Instrument response is a signed field transfer in observed coordinates, with
analytic zero-mode identity and explicit pixel/sigma/FWHM conventions. Ordinary
resolving power and labeled legacy compatibility remain distinct. Required-pair
products multiply mean powers and Jacobians once; supplied noise is separate.
PD2013 constants and its moving stationary floor are retained in an independent
zero-local-parameter callable. Response and comoving conversion are explicit
consumer operations. No forest-noise preparation or orchestration was added.

## Fresh review validation

All numerical commands used OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, and
MKL_NUM_THREADS=1 on the login node. Artifacts are under
`.validation/step09-review-r1/`, separate from implementation evidence.

From the package root:

```bash
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh \
  > .validation/step09-review-r1/quick.log 2>&1
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python .validation/step09-review-r1/verify.py \
  > .validation/step09-review-r1/verification.log 2>&1
git diff --check
```

- **555 tests passed in 21.58 seconds**, no skips; Ruff lint and formatting passed
  (74 files reported formatted by the current runner). This is a fresh review
  run, distinct from the handoff's 20.99-second/73-file result.
- **48 supplied-background volume cases** use independently integrated monomial
  volume elements with varied bins, H/D_M, area, h_fid, and quadrature orders.
  Maximum relative error is **5.00e-15**, below rtol=8e-13.
- **272 P1D points** at eight redshifts from -0.9 to 9 compare to a 60-digit
  Decimal implementation of the original power-law expression. Samples include
  zero, each independently calculated floor and its neighborhood, and a broad
  logarithmic k range. Maximum relative error is **2.66e-15**, below 5e-13.
  Negative redshift probes test the direct callable's mathematical contract;
  they do not change the existing evaluator's nonnegative-redshift boundary.
- **40 observed nodes and three fields** compare signed pixel/Gaussian responses
  with independent scalar sin(x)/x expressions. Maximum absolute error is
  **1.25e-16**, within rtol=5e-13, atol=3e-16; negative lobes are exercised.
- Complete-pair amplitude Fisher agrees with the independent field-matrix trace
  identity, with reported relative difference **0.0**. A selected subset with
  additional covariance closure also matches independently assembled scalar
  field-matrix covariances. These checks use supplied positive noise.
- The initial reviewer script chose pixel widths too small to reach a negative
  sinc lobe at its own geometry/grid, failing its fixture assertion. The review
  fixture widths were increased to 8000/13000 km/s solely to exercise that
  mathematical limit. The original log is retained as verification-initial.log.
  No production code or acceptance tolerance changed.

The ordinary tests additionally reran EdS and flat/open/closed Astropy volume
convergence, mode/area scaling, conversion inverses, P1D join limits, identity/
strong-response boundaries, immutable preparation, and fixed derivative inputs.
Synthetic response widths and backgrounds are numerical tests, not instrument
recommendations or physical survey-accuracy claims.

## Source and artifact identity

Before status-document updates, all **13** handoff manifest hashes matched live
files. The implementation's **62-file** pre-dispatch snapshot confirms exactly
README.md and pyproject.toml changed among those files; the other 60 match.
The review's own **79-file** starting snapshot and exact assignment copy are
preserved in before.json and assignment-r1.md. Historical reports/evidence are
unchanged. A HEAD diff alone would mix Steps 03–09, so these manifests identify
the actual reviewed implementation.

Final wheel SHA-256:
`3f6802f8f745be8555ba5ba34300d83c841bca4d0d3322321ebbba9fe8b6777f`.

The reported initial installation was followed by a corrective final-wheel
reinstall into the same isolated environment. This is disclosed and acceptable:
explicit final-artifact and installed-payload comparison verifies the code being
reviewed. The reviewer did not create another environment or rebuild the wheel.

After reading probe.py, the reviewer reran it from /tmp using the existing
implementation wheel-env/bin/python with `-I`, the exact `dist-final/` wheel,
the absolute source root, and the new reviewer output directory. All three
example copies in that directory were byte-compared with source. Original probe
outputs were not overwritten.

The rerun passed quiet/lazy imports, optional dependency metadata, **25-module**
source/wheel/installed byte identity, METADATA/WHEEL identity, geometry/P1D/
response and amplitude-Fisher oracles, curved Astropy volume, all three examples,
and a separate subprocess actively blocking optional imports. Module origins
resolve inside the isolated environment's site-packages. Results are in
wheel-probe.log, wheel-probe.json, and blocked-optional.log. This is a fresh
execution of the existing installation's probe, not a new fresh installation.
Python 3.13.15, NumPy 2.5.3, Astropy 8.0.1, SciPy 1.18.1.

## Reference evidence and limits

Inspected the bounded AST-based legacy comparison script and its **33** saved
comparisons; all **five** referenced source hashes still match. The script uses
minimal fake backgrounds, isolates the 299792.458/299800 light-speed difference,
documents legacy zero-mode/small-argument response behavior, and compares broad
and narrow bins without forcing integrated and bin-centre volumes to agree.
Those saved comparisons were verified, not rerun. Fresh analytic tests and
installed probes above supply independent evidence for the current FishHighz.

No real legacy quick capture, full forecast suite, CAMB initialization, SNR-file
loading, or Slurm action was required or performed. Full execution remains
explicit-request only and is not an outstanding acceptance gate. NERSC MUNGE
socket messages followed successful quick checks; commands returned zero.

The source/quick/numerical/installation checks leave no required findings.
Status documents record review passed while preserving revision 1's scientific
assignment. Await the user's next instruction; no commit, push, delegation,
or Step 10 preparation occurred.
