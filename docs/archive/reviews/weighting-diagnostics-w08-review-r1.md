# Scientific conclusion and impact of proposed changes

**W08 revision 1 requires one bounded representability correction.** The
implemented `method="inverse_variance", alias=B_star` reproduces the intended
fixed reference for ordinary and saved inputs. Independent checks confirm the
exact rational coefficients, equality to independently supplied weights and the
legacy zero-iteration seed, immutable iteration-free metadata, correct
mode-dependent P1D/response ownership, and fixed weights/noise during survey
mean differentiation. No legacy/supplied compatibility regression or ordinary
noise regression was found.

One positive instrumental-power product can nevertheless underflow silently.
Because `l_p*v` is evaluated under `under="ignore"`, positive finite inputs can
produce zero before addition to B_star, yielding an incorrect accepted weight
even when the final A and P_pixel remain representable. This does not change the
reported rational or saved W06 coefficients, but it contradicts W08's retained
range-safety contract and overstates the justified input domain. The smallest
resolution is to reject inexact underflow of positive `l_p*v` in the new branch,
with one reproducer and one exact-zero-variance control. Do not alter the legacy
branch, `_iterate`, `_integrals`, noise, survey or profile defaults.

No substantive user feedback accompanied the request; its feedback field was
left as the placeholder. The exact r1 assignment is archived and W08 revision 2
now specifies only this repair. This review is not user acceptance, profile
adoption or permission to advance.

## Exact implementation and retained scientific behavior

The live hashes match the r1 handoff: `fishhighz/weights.py` is
`36d4d489...1a16`, the unchanged weighting kernel is `3867c9ee...e66f`, the
new focused test is `64814534...c38`, and the saved report is
`9b64dd91...666`. `IMPLEMENTATION_STEP.md` remains unchanged at
`3675cd86...c310`; package Step 13 was not reviewed or modified.

The preserved pre-W08 copy shows that only the new method branch and docstring
were added to `weights.py`. The supplied and legacy branches are otherwise
byte-identical. Source inspection establishes:

- `weights.py:268–275` rejects weights, iterations, signal and auxiliary for
  the new method, validates positive finite B_star, and records no iterations;
- `weights.py:291–296` constructs the fixed vector directly and never invokes
  `_iterate` or a provider;
- `noise.py:52–91` retains
  `A*P1D(q)*W(q)^2 + P_pixel`, with the common 3D conversion once and no response
  on pixel noise;
- `forecast.py:98–138` evaluates final P1D/noise once during preparation, while
  `forecast.py:238–278` differentiates only the mean using fixed prepared state.

The returned r1 object records `method="inverse_variance"`, `alias=B_star`,
`signal=iterations=auxiliary=None`, an empty immutable `weight_changes`, owned
immutable arrays and the existing field/geometry/response context. Zero density
stores zero weight, and zero variance on positive support stores exactly one.
The fixed B_star is not substituted for final P1D and is not reoptimized at
each Fourier mode.

## Independent coefficients, noise and survey ownership

The focused invocation passed **109 tests in 3.21 s**:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 \
  .venv/bin/python -B -m pytest -q \
  tests/test_inverse_variance_weights.py tests/test_weights.py \
  tests/test_weight_range.py tests/test_weight_models.py tests/test_noise.py \
  tests/test_forecast.py::test_generated_spies
```

The independent review script imports neither the implementation check script
nor its results. Exact Fraction arithmetic with q=(1/2,2), rho=(2,1/2),
v=(1,4), and L=l_p=B_star=1 gives r=(1,1), weights=(1/2,1/5), moments
(7/10,29/100,41/100), A=29/49, P_pixel=41/49 and Q=10/7. The public method
agrees to at most `1.88e-16` relative error; independently supplied weights and
legacy zero iterations agree exactly. Traps on `_iterate`, auxiliary sampling,
P3D and P1D provider functions remain uncalled.

At q=(0.003,0.02) s/km with intrinsic P1D=(4,9) km/s, B_star=2 km/s and a
nontrivial sinc/Gaussian response, independent scalar aliasing, unsmoothed pixel
noise and totals agree exactly. The supplied-weight route gives identical noise.

An independently constructed tiny forest bin gives one P3D and one P1D call per
preparation for both inverse-variance and supplied weights. The weights,
coefficients, noise and Fisher matrix agree exactly. After two preparations and
two two-parameter central-difference runs, P3D has 12 calls (one preparation and
five mean calls per method) while P1D remains at two. Runtime traps verify that
`run_bin` does not repeat weight, response, noise or P1D preparation.

The saved coefficient check selects only accuracy bin 0, order-32
`lya(qso)` from the assigned report after verifying its exact hash, identity,
partition and 5,184 nodes. Scalar 80-digit Decimal arithmetic forms rho_v,
nu and all moments directly from exact saved binary values. The public method
returns:

| Quantity | Public method | Independent Decimal |
| --- | ---: | ---: |
| A [deg^2] | 0.029419845663300446 | 0.029419845663300400 |
| P_pixel [deg^2 km/s] | 0.31435492218062416 | 0.31435492218062315 |
| Q [deg^2 km/s] | 0.7955679523681858 | 0.79556795236818397 |

The maximum coefficient discrepancy is `3.18e-15`; agreement with W06's printed
reference is `1.54e-15`, both below rtol=5e-12 with zero absolute tolerance.
This is a coefficient-only synthetic geometry check, not a geometry or forecast
replay. The inherited density/SNR policies were neither reapplied nor physically
validated.

## Required finding W08-R1: positive-product underflow bypass

At `weights.py:281–295`, the new denominator is evaluated inside a surrounding
`under="ignore"` state. The independent reproducer uses

```text
tiny = nextafter(0, 1)
l_p = B_star = tiny
q = (1, 1), rho = (1/8, 1/8), v = (0, 1/2), L = 1.
```

The positive product `tiny*(1/2)` underflows to zero. Revision 1 accepts
weights=(1,1), A=4 and P_pixel=tiny; therefore neither the positive-weight check
nor the final coefficient checks detect the lost instrumental contribution.
The mathematical second weight depends on that positive term, so accepting one
is not a faithful representable evaluation of the stated estimator.

Possible impact is confined to the admitted extreme float64 boundary; it does
not affect the ordinary examples or the saved DESI-2 coefficient sample. It is
still consequential to the justified range claim because W08 explicitly retains
rejection of unrepresentable positive arithmetic. The minimal resolving check is
the reproducer above: revision 2 must raise contextual `ValueError`, while an
otherwise identical exact-zero-variance input must continue to return weight one.
A local raised-underflow multiplication is sufficient; no range extension,
rescaling, logarithmic evaluation or algorithm change is warranted.

## Evidence, execution and disposition

Independent evidence is in
`.validation/forest-weight-diagnostics/w08-review-r1-independent/check.py` and
`result.json`. The successful bounded invocation completed in **0.323 s**:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 30 \
  .venv/bin/python -B \
  .validation/forest-weight-diagnostics/w08-review-r1-independent/check.py
```

An initial review-script invocation failed before saved coefficient preparation
because its synthetic geometry used z_eval=2.5 above the saved
z_source=2.3830. The reviewer corrected the geometry to the assigned saved
z_eval=2.1152848986890427; this was a review-fixture error, not a production or
scientific failure. Final Ruff lint/format and `git diff --check` pass.

The exact r1 assignment is preserved as
`reviews/weighting-diagnostics-w08-instructions-r1.md` with SHA-256
`a1906a78...ad2`. `WEIGHTING_DIAGNOSTIC_STEP.md` is now revision 2 and carries
only W08-R1. No production repair was implemented by the reviewer. W07's BAO
comparisons remain historical evidence; W08 produced no new survey result.
No real forecast, broad suite, wheel, environment, profile change, physical
input-policy change, package Step 13 work, subsequent step, Slurm action, agent
dispatch, commit or push was performed. Stop for user review and approval.
