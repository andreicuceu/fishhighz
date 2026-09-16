# Step 10 revision 2 independent review

Date: 2026-09-13.
Outcome: **review passed; R1 resolved; no further required changes**.
The user reported no comments. Acceptance and progression remain with the user;
no Step 11 planning or new repair revision is needed.

## Repair assessment

Reviewed the exact revision-2 assignment, current handoff, archived revision-1
handoff, prior independent finding, live source changes, new regression tests,
README clarification, and installation/preservation evidence.

The repair uses the explicitly permitted rejection option. The three individual
products in `kernels/weights.py:_integrals` execute under a local NumPy inexact
underflow trap. The existing host boundary converts failures into ValueError
with field identity and the underlying arithmetic diagnosis. This prevents
another positive contribution from concealing a lost I3 term. Multiplication
order, prefixes, coefficient formulas, and the legacy iteration remain unchanged.

The original unscaled R1 fixtures now fail clearly in both sample orders.
Their representable rescalings succeed with P_pixel=1e-100. Exact-zero products
remain supported. The README explains the numerical-range rejection; the API
does not promise to recover results from every extreme weight normalization.
No high-precision runtime dependency, implicit renormalization, or scientific
convention change was introduced.

## Fresh review validation

Evidence is in `.validation/step10-review-r2/`, separate from both implementation
runs and the original failing review. Numerical execution used single-thread
OMP/OPENBLAS/MKL limits on the login node.

```bash
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh
.venv/bin/python .validation/step10-review-r2/range_regression.py source
.venv/bin/python .validation/step10-review-r2/verify.py
git diff --check
```

- **638 tests passed in 36.52 seconds**, no skips; Ruff lint and formatting
  passed (86 files). The original 628 tests are retained unchanged.
- The original byte-identical four-case R1 regression passed on source and
  the exact installed revision-2 wheel: two diagnosed underflow rejections and
  two successful rescaled controls, with zero coefficient discrepancy.
- Independent checks exercised **12 operation-boundary rejections**, covering
  all three products, both sample orders, and caller underflow policies of
  ignore/raise. Every failure retained field/underflow context, and the caller's
  NumPy arithmetic policy was restored after the call.
- **Six exact-zero controls** and an **exact subnormal control** passed. The
  latter retains first-prefix I1=2^-1030 and I2=2^-1060, demonstrating that the
  trap distinguishes exact subnormal products from inexact underflow.
- **96 supplied-weight Decimal comparisons**, including 48 normalization pairs,
  passed for all prefixes and both coefficients at rtol=5e-13, atol=0. Maximum
  prefix relative error across the accepted independent fixtures was 2.23e-16.
- The existing isolated installation's probe was rerun from /tmp with `-I`,
  the exact wheel path, absolute source root, and the new review output directory.
  All 28 module source/wheel/installed byte comparisons, METADATA/WHEEL checks,
  scalar weight/noise and PSD checks, amplitude Fisher, four byte-matched example
  copies, Astropy geometry, and blocked-optional-import subprocess passed.
  This was a new probe execution of the reported installation, not a new build
  or fresh installation by the reviewer. The R1 replay used that same interpreter.

Logs and JSON results are quick.log, range-source/wheel.log/json,
verification.log/json, wheel-probe.log/json, and blocked-optional.log.
All review checks passed without fixture or production corrections. MUNGE socket
messages appeared alongside successful quick checks; no scheduler action occurred.

## Identity, preservation, and limits

Before governance updates, all **91** final maintained-file hashes and **30**
implementation evidence hashes matched live artifacts. The pre-dispatch manifest
confirms exactly four of 89 preexisting files changed: README.md, the two
weighting modules, and reviews/step-10.md. The other 85 are byte-identical.
The archived reviews/step-10-r1.md matches the original handoff hash and its
implementation archive. A direct comparison against the revision-1 wheel
confirms only the two reported production modules changed; the diff is saved.
The review's own 91-file manifest, original assignment and governance snapshots
protect existing dirty/untracked work. Base commit remains
`0d69786a06d5d676564a51fad14a7156951c2228`.

Exact reviewed wheel SHA-256:
`670a143f19fd9a507e6c0ecb07c43ba0c8a5bb988bdb2bf28393525fbc6b04ac`.
The probe records installed origins and Python 3.13.15, NumPy 2.5.3,
Astropy 8.0.1, and SciPy 1.18.1. NumPy remains the sole unconditional dependency.

All **34** legacy source/resource hashes still match. The **81** bounded
reference comparisons from revision 1 remain historical evidence, not new runs.
Reusing them is appropriate for this validation-only repair: reference inputs,
formulas, and ordinary multiplication order are unchanged. This review makes
no new physical DESI-2 forecast claim. No real quick capture or full forecast
suite was requested or run; full execution is not an outstanding gate.

Only review evidence, this report, and status/governance documents were written.
Production code, ordinary tests, examples, dependencies, implementation handoffs,
earlier review evidence, and sibling repositories were preserved. Revision 2's
repair/scientific instructions remain intact, with its status marked passed.
No commit, push, dispatch, acceptance, or progression occurred.
