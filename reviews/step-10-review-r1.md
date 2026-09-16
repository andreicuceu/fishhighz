# Step 10 revision 1 independent review

Date: 2026-09-13.
Outcome: **changes required**, one medium-priority finding (R1).
The user reported no comments. Step 10 is not accepted or ready for progression.
[Revision 2 instructions](../IMPLEMENTATION_STEP.md) retain the original scope
and add the bounded repair below. No production code was changed by this review.

## R1 — P2: mixed underflow silently loses a pixel-noise contribution

Locations: [kernels/weights.py:25](../fishhighz/kernels/weights.py#L25) and
[weights.py:282](../fishhighz/weights.py#L282).
The left-associated product `r*w*w*variance` can underflow before the variance
multiplier restores its mathematical magnitude. Underflow is ignored, and the
host guard only detects loss when the **entire final** P_pixel is zero. Another
positive contribution lets an incorrect positive coefficient pass validation.

Reproduction uses two ordered magnitudes, `rho=quadrature=[1,1]`, supplied
`weights=[1,1e-200]`, `variance=[1e-150,1e300]`, and `L_v=Delta_v=1`.
Geometry is the ordinary synthetic background used by the review, with
z_eval=2.4 and z_source=3; its values do not affect this coefficient.

- Returned I3 is `[1e-150,1e-150]` and P_pixel is `1e-150`.
- A 100-digit Decimal oracle gives I3 approximately `[1e-150,1e-100]`
  and P_pixel approximately `1e-100`. All cumulative outputs in this ordering
  are representable as finite float64 values.
- Multiplying both supplied weights by `1e100` returns the correct P_pixel
  `1e-100`. Thus the accepted results violate uniform-weight normalization
  invariance by a factor of `1e50`.
- Reversing the sample rows also reproduces the error. The same four probes
  reproduce it against the exact installed wheel, outside the checkout.

These are deliberately extreme numerical inputs, not realistic DESI-2 densities
or variances. The ordinary/reference results below are unaffected. Nevertheless,
the explicit assignment requires unrepresentable arithmetic to be diagnosed,
and the supplied-weight API currently returns silently wrong known noise. This
is a numerical validation defect, not a proposed change to weighting science.

Close R1 by computing the affected contributions safely or rejecting their
unrepresentable intermediate arithmetic clearly. Detect individual lost terms;
checking only final coefficient positivity is insufficient. Preserve legitimate
zero density/weights/variance. Add independent regressions for both row orders,
the rescaled control, and zero controls, without weakening tolerances or adding
floors/regularization. The repair assignment gives exact acceptance requirements.

## Fresh validation and source assessment

Read the exact revision-1 assignment, implementation handoff, all three new
scientific modules, three new test files, example, README, and relevant existing
geometry/response/model/covariance contracts. Per-field auto routing, independent
P1D, source/evaluation units, cumulative simultaneous iterations, unsmoothed
pixel/Poisson noise, explicit independence, full replacement, and independent
noise PSD validation follow the agreed conventions. No additional required
scientific or architectural finding was identified.

New evidence is separate under `.validation/step10-review-r1/`. Commands used
OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, and MKL_NUM_THREADS=1 on the login node.

```bash
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh > .validation/step10-review-r1/quick.log 2>&1
.venv/bin/python .validation/step10-review-r1/verify.py
.venv/bin/python .validation/step10-review-r1/range_regression.py source
.venv/bin/python .validation/step10-review-r1/reference.py
git diff --check
```

- **628 tests passed in 31.92 seconds**, no skips; Ruff lint and formatting
  passed (83 files reported formatted). This rerun does not close R1, which the
  ordinary suite lacks. MUNGE socket messages accompanied successful checks;
  no Slurm action occurred.
- **256 independent Decimal weighting cases** (64 varied samples at 0/1/3/6
  updates), including zero density/variance, check weights, all prefixes, A,
  P_pixel, and scalar response/noise. Maximum relative error across these and
  the independent matrix check was **1.56e-15**, within rtol=5e-13, atol=0.
- A four-node, three-field complete-pair Fisher calculation with signed PSD
  noise matches the independent field-matrix trace identity; amplitude Fisher
  is **24857.631425126114**. Invalid noise is rejected independently of signal.
- **R1 regression fails as expected on both live source and installed wheel**:
  two unscaled cases fail; two rescaled controls pass. Exact outputs and
  assertion failures are preserved in range-source.json/log and
  range-wheel.json/log. These are actual open failures, not passing acceptance.
- The first general reviewer fixture combined very small pixels with a fixed
  Gaussian width and large negative-lobe test coordinates, triggering the
  intended noise representability error. Its log is preserved as
  verification-initial.log. The fixture now uses sigma=0.05*pixel to test ordinary
  algebra separately from R1. No package code or tolerance changed.

## Reference and installed evidence

All **88** final maintained-file hashes, **33** implementation evidence hashes,
and **34** reference source/resource hashes matched before governance edits.
The implementation's 80-file pre-dispatch manifest confirms README.md was the
only changed preexisting file; the other 79 remain identical. The review's own
88-file starting snapshot, Git status, exact original assignment, and governance
copies are preserved. Base commit remains
`0d69786a06d5d676564a51fad14a7156951c2228`; earlier dirty/untracked work is retained.

Read the bounded AST reference script and relevant live legacy methods, then
reran an identical script copy with outputs in the review directory: **81
comparisons passed**, covering synthetic arrays and the five real populations
at one representative bin point. No CAMB, NewForecast, Fisher forecast, or
capture ran. This is a fresh bounded input/intermediate comparison, not a real
DESI-2 forecast. The handoff's synthetic background/S/B, four-magnitude subgrids,
normalization policies, rejected negative spline rows, and light-speed/response
differences remain explicitly disclosed. Original evidence was not overwritten.

After inspecting probe.py, reran the existing installation's probe from /tmp
with `-I`, absolute source path, exact wheel, and the new review output directory.
All **28 module** source/wheel/installed byte comparisons, METADATA/WHEEL,
optional metadata, scalar weight/noise and Fisher checks, four byte-matched
example copies, Astropy geometry, and blocked optional-import subprocess passed.
Separately ran range_regression.py with this same isolated interpreter; R1
also fails there. These executions reuse the reported installation; the reviewer
did not build or install a new wheel.

Exact wheel SHA-256:
`e379ddee37694daabb1acc9220d32d6ba1f9f7ec6df1a593f08cdd61e47a2cd2`.
Python 3.13.15, NumPy 2.5.3, Astropy 8.0.1, SciPy 1.18.1. Probe results/origins
are in wheel-probe.json/log and blocked-optional.log.

Design, roadmap, package guidance, and the current assignment now record the
repair state. Historical handoffs/reviews and implementation artifacts remain
unchanged. No source, ordinary tests, example, dependency, or sibling checkout
was edited. No full forecast suite, commit, push, dispatch, or Step 11 planning
occurred. Full execution was not requested and is not an outstanding gate.
