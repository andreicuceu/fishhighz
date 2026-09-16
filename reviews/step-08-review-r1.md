# Step 08 revision 1 review

Date: 2026-09-12.
Outcome: review passed; no required changes. Ready for progression when the user
requests Step 09. No implementation source or ordinary tests were changed by the
reviewer.

## Assessment

Reviewed the assignment, handoff, host model and kernels, new tests, README,
synthetic example, and saved reference/packaging evidence. The implementation
satisfies the bounded intrinsic Kaiser/BAO assignment.

Fixed/free settings and explicit global ties reuse the accepted callable and
binding interfaces. Forest b/beta and shared galaxy f remain distinct, including
zero galaxy bias. Both components receive their own inverse Fourier mapping,
angular factors, and volume normalization. Gaussian damping uses the wiggle
coordinates and mean squared widths; varying f changes neither widths nor G.
G is applied once. Both mapped domains are validated without moving observed
cuts or derivative steps.

Preparation reuses the accepted template, evaluations batch template/factor work
across pairs, and there is no mutable parameter-dependent cache. Numerical
derivatives and covariance/Fisher calculations reuse the accepted machinery.
The example keeps response, noise, grid, and covariance factors fixed and
discloses its synthetic nuisance priors. No survey accuracy or production
analytic-derivative claim is made. No repair is required.

## Independent validation

From the package root:

```bash
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python .validation/step08-review-r1/verify.py \
  > .validation/step08-review-r1/verification.log 2>&1
git diff --check
```

- **468 quick tests passed in 28.46 seconds**, including 60 new tests, with no
  skips. Ruff lint and formatting passed (62 files already formatted).
- Eighteen independent polynomial-template cases cover all nine combinations of
  smooth/wiggle bases at nonidentity anisotropic settings, each with two
  unit/growth/damping setups. Each evaluates five fields and all 15 pairs at
  23 nodes, including axis limits and zero galaxy bias.
- The oracle evaluates known polynomial amplitudes directly without production
  interpolation or kernels. Maximum absolute power error was `1.14e-13`, with
  comparison tolerance `rtol=2e-12, atol=3e-12`.
- Complex-step differentiation of that independent analytic formula supplies
  12 derivative columns per case. Production finite differences use explicit
  steps of `1e-5`. Across **216 columns**, maximum absolute error was `1.64e-7`;
  maximum error divided by `max(1,abs(reference derivative))` was `1.29e-8`,
  below the declared `3e-7` tolerance. This exercises anisotropic angular,
  volume, and damping derivatives away from identity.
- Pair-column reordering and concatenated node slices agree exactly. Reversing
  field indices agrees within roundoff. An initial reviewer assertion wrongly
  demanded bitwise equality for reversed fields despite changed association of
  Q*B_i*B_j; its largest discrepancy was `2.85e-14` absolute (`4.35e-16`
  relative). That assertion now uses `rtol=2e-15, atol=3e-14`. No production
  code or established test tolerance changed. The initial log is preserved.

The verifier/results are under `.validation/step08-review-r1/`. Numerical
commands used the three single-thread limits. NERSC MUNGE socket messages
appeared after successful quick-test subprocesses; checks exited zero.

## Source, artifact, and reference evidence

All five inline handoff hashes match current files. The exact wheel hash is
`a2224a6179d55578e22f95f28a1bbd37b8f0c48339c5cb0b1938a1f1f248de9f`, and all
21 package Python modules match source. The 59-file pre-dispatch snapshot
confirms README.md was the only changed preexisting maintained file. These
checks preceded reviewer status/documentation updates. Earlier scientific code,
tests, reports, dependencies, and governance were preserved.

Inspected the bounded Vega script and its 15 saved passing cases. Both reference
source hashes still match. It extracts only the intended methods, supplies both
widths, and uses fresh objects to avoid the approximate cache comparison. This
review verified saved evidence and source identity; it did not rerun these
comparisons or a full Vega model.

The reviewer reran the probe in `.validation/step08-r1-20260912T204258/`, using
its existing isolated wheel-env/bin/python with -I from /tmp and all three
thread limits. Arguments were the exact wheel, absolute source root, and the
new `.validation/step08-review-r1/` output directory. Both examples were copied
there after confirming byte identity; original evidence was not overwritten.

The probe passed quiet/lazy imports, dependency metadata, source/wheel/installed
module and metadata identity, the constant-model power/derivative/Fisher oracle,
and both built-in BAO/AP+f and retained external examples. Imports resolved under
the isolated environment's site-packages. Results are in wheel-probe.log and
wheel-probe.json in the reviewer directory. This reruns the reported fresh
installation's probe; it does not claim a second fresh installation. Versions:
Python 3.13.15, NumPy 2.5.3, SciPy 1.18.1, Astropy 8.0.1.

## Limits and progression

No real lyaforecast quick capture or full suite was required or run. Full mode
remains restricted to explicit user requests. These checks do not establish
physical survey/cosmology accuracy, realistic noise, template-sampling convergence
for a cosmology, or compiled performance.

Step 08 revision 1 and its handoff satisfy the assignment. No repair revision is
requested. Status documents record the review; await the user's instruction
before preparing Step 09. No commit or push was performed.
