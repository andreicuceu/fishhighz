# Step 13 revision 1 independent review

## Outcome

**Requires revision.** The finite-recurrence reformulation and the numerical
fixed-grid results are well supported. The submitted evidence establishes strong
grid dependence over the saved orders and substantial iteration drift. It does
not yet establish the report's unconditional absence of a finite continuum limit
for the full nonlinear recurrence on all twelve real populations. Two additional
required corrections concern source/attempt binding and diagnostic arithmetic.

Reviewed against the exact revision-1 assignment, preserved in
[step-13-instructions-r1.md](step-13-instructions-r1.md). The implementation and
scientific reports remain unchanged. New checks, original document snapshots,
source identity and deliberately corrupted **copies** are in
step13-review-r1 (local-only path: `../.validation/step13-review-r1/`). Step 12 R2 remains closed;
its R3 and scientific acceptance remain open. This review introduces Step 13
finding identifiers S13-R1/R2/R3 to distinguish them from Step 12 findings.

No production code, tests or implementation scripts were edited. No new survey
forecast, raw-reader interpolation, external model evaluation, installation,
Slurm operation, agent dispatch, commit or push occurred.

## 1. Scientific results that survive independent checks

### 1.1 Amplitude decay and noise convergence are different questions

For positive supported samples, the live recurrence is

\[
 w_i^{t+1}=\frac{N_i^t}{N_i^t+v_i/S},\qquad
 N_i^t=\frac{L}{\Delta_v}\sum_{j\leq i}r_jw_j^t.
\]

Writing \(w=su\) and retaining both \(s\) and \(u\) is algebraically equivalent
to this recurrence. The implementation uses logarithms for both. Its moments
cancel the common amplitude before evaluation:

\[
 A=\frac{\sum r_i u_i^2}{L(\sum r_i u_i)^2},\qquad
 P_{\rm pixel}=\frac{\Delta_v\sum r_i u_i^2v_i}{L(\sum r_i u_i)^2}.
\]

The units remain deg² and deg² km/s, respectively. Renormalizing the original
weights before the nonlinear update would change its denominator and therefore
the estimator; this prototype retains that amplitude dependence.

Independent long-double direct updates reproduce A and P_pixel through 24
updates on **all 60 saved bin/population/order combinations**, with maximum
relative discrepancies **1.04e-13 and 1.21e-13**. A separate 160-digit Decimal
recurrence carried a three-cell example far beyond float64 weight underflow;
prototype discrepancies were **3.42e-14 in A**, **8.22e-15 in P_pixel**, and
**5.17e-13 in relative log shape**. Its log-amplitude discrepancy was 9.10e-13
absolute. These checks independently connect the diagnostic to the recurrence;
the implementation's existing 80/160-digit test compares Decimal with Decimal
and does not itself make that connection.

Thus the old arithmetic failures need not imply that the coefficient ratios
are undefined. The unchanged production integrals still reject 56 of 200 saved
weight states. The diagnostic preserves all 144 available and 96 unavailable
historical coefficient attempts. This numerical improvement does not by itself
establish a converged estimator or authorize production adoption.

### 1.2 Fixed-grid decay and limiting shapes

For positive variance, the zero-weight derivative is lower triangular with
\(d_i=(L/\Delta_v)Sr_i/v_i\). The saved spectral radii are
**2.99338e-5 to 0.407247644**, with unique dominant entries. Since
\(0\leq F(w)\leq Mw\), \(w_t\leq M^t w_0\to0\) at each fixed grid. This
comparison handles nonnormal transients; the spectral radius is not a one-step
norm-contraction factor and is not a shape-convergence tolerance.

The fixed-grid eigenvector calculation is independently reproducible without a
dense eigensolver or the implementation's forward eigenvector recurrence. Let
\(p_i=r_i x_i/\sum r_jx_j\), \(C_i=\sum_{j\leq i}p_j\), and
\(\lambda=\max d_i\). The eigenvector equation gives

\[
 p_i=(d_i/\lambda)C_i,\qquad
 C_{i-1}=(1-d_i/\lambda)C_i,\qquad C_n=1.
\]

Starting at the last cell and working backwards to the unique dominant cell
therefore gives the limiting probability distribution. All 60 distributions
match the saved values to **1.06e-14 absolute**; independently summing
\(A=\sum p_i^2/(Lr_i)\) and
\(P_{\rm pixel}=\Delta_v\sum v_i p_i^2/(Lr_i)\) reproduces the fixed-grid
coefficients to **4.01e-13 and 6.39e-13 relative**.

For these subcritical positive finite grids, the nonlinear perturbations are
summable: if \(\epsilon_t=\|Mw_t\|_\infty\), then
\(Mw_t/(1+\epsilon_t)\leq F(w_t)\leq Mw_t\), and
\(\sum_t\epsilon_t<\infty\). Applying these bounds from a sufficiently late
starting time supports the unique-dominant normalized limit. The report should
state this argument and quantify a bound wherever it uses a finite-time
approximation. The check above validates the eigenvector coefficients, not a
claim that 1024 updates have reached them.

Indeed, at order 64, the asymptotic P_pixel is still **38.5–43.0% larger** than
the 1024-update value for QSO forests, and **68.7–101.0% larger** for LBG forests.
These ratios use the finite-iteration value as denominator. Small absolute
weights have concealed substantial continuing changes in their relative shape.

### 1.3 Refinement dependence is large and physically consequential

The independently verified fixed-grid limits give these representative values:

| Population, bin 0 | Order 16 | Order 32 | Order 64 |
| --- | ---: | ---: | ---: |
| QSO-forest A [deg²] | 0.0374662 | 0.0572858 | 0.0971897 |
| QSO-forest P_pixel [deg² km/s] | 12.6340 | 24.7726 | 49.0883 |
| LBG-forest A [deg²] | 2.34876 | 4.39865 | 8.26859 |
| LBG-forest P_pixel [deg² km/s] | 7671.75 | 14849.5 | 29171.9 |

Every one of the twelve populations shows more than 80% growth in P_pixel on
both 16→32 and 32→64 refinement. The three-order empirical exponent is
0.9615–0.9798. The effective-measure exponent is negative throughout,
−0.9079 to −0.2671. This is compelling evidence for concentration over the
sampled refinements, and clearly fails the proposed 1e-4 coefficient screen.
No normalization of the common amplitude removes these shape changes.

The concentration statistic has a useful direct interpretation:
\(R_{\rm eff}=[\sum_i p_i^2/r_i]^{-1}=1/(LA)\). It measures the effective
sightline measure used by the normalized weights, rather than their amplitude.
The pixel coefficient additionally weights this concentration by variance.
Its exponent can therefore differ from the exponent of A.

The position of the largest diagonal entry is **not** the location at which the
limiting sightline distribution concentrates. For bin-0/order-64 QSO forests,
those values are m=19.3488 and mean m=23.4428 (width 0.1441); for LBG forests
they are m=22.5244 and mean m=24.3517 (width 0.06862). The stable dominant-cell
coordinate quoted in the handoff is insufficient to establish concentration at
a stable physical feature. Use the distribution, its width and its tail measure.

At fixed geometry and P1D/response,
\[
 N_F=\left[A P_{1D}W^2+P_{\rm pixel}\right]d_{\rm deg}^2/a_v.
\]
Both terms are nonnegative. Increasing P_pixel raises an additive noise floor
that smoothing of the alias term cannot cancel. This explains why a converged
finite-grid ratio is insufficient for a grid-independent forest forecast.
It does not determine the numerical changes in joint or individual-spectrum
BAO errors; those require a separately authorized forecast with an adopted
prescription. Step 12's six failed accuracy bins remain failed.

## 2. Required findings

### S13-R1 — high: finite refinement screens are promoted to a continuum theorem

Locations: diagnose_weight_limit.py:421 (local-only path: `../scripts/diagnose_weight_limit.py:421`),
diagnose_weight_limit.py:610 (local-only path: `../scripts/diagnose_weight_limit.py:610`), and
scientific report:125 (local-only path: `step-13-weight-limit-r1.md:125`), especially its concluding
claims at lines 196–208.

The negative verdict uses two relative coefficient changes greater than 0.8,
shrinking effective measure and a small dominant-coordinate displacement. These
are empirical screens. The report supplies the constant-coefficient **linearized**
Volterra counterexample, but no continuum concentration argument whose
assumptions are verified for the full nonlinear saved-input problem. There is
no analysis establishing exchange of iteration/refinement limits or delimiting
which ordering the negative result actually excludes. An absolute coefficient
change also does not encode the direction of change.

The review's logical counterexample supplies the finite convergent sequence
\(A_n=P_n=n/(1+n/10000)\), with effective measure \(1/A_n\) and coordinate
\(22+n^{-2}\), at n=16,32,64. The classifier returns
`no_suitable_finite_grid_independent_limit`, although both coefficients tend to
10000. This is a test of the **decision rule**, not a claimed DESI recurrence
or new physical source model. It demonstrates that the implemented premises
cannot establish their reported conclusion.

A possible derivation for the iteration-first limit follows from the independent
backward relation above:
\[
 C_i=\prod_{j>i}(1-d_j/\lambda)
 \leq\exp[-\sum_{j>i}d_j/\lambda].
\]
For a consistently refined non-atomic measure, if max(r_i/v_i)→0 while each
fixed terminal interval retains positive integral of dR/v, the distribution
concentrates towards the faint support boundary. If its tail measure R_tail→0,
Cauchy–Schwarz gives \(LA\geq p_{\rm tail}^2/R_{\rm tail}\); a positive lower
bound on variance in that tail also bounds P_pixel from below. These are
**conditional sufficient assumptions**, not a verified theorem for every
saved survey population. The reconstruction of the continuous measure, support,
variance bounds, quadrature family and nonlinear limit must be documented.

**Closure:** establish a controlled result with explicitly checked assumptions
and the intended limit ordering, or classify the continuum question as
unresolved while retaining the strong measured refinement dependence. A finite
candidate also needs the originally required remainder/concentration argument;
passing screens alone cannot yield `finite equivalent limit supported`.
Do not select a new estimator or alter convergence thresholds to close this.

### S13-R2 — high: finalization does not bind source identity or attempted scope

Locations: weight_limit.py:724 (local-only path: `../fishhighz/validation/weight_limit.py:724`),
diagnose_weight_limit.py:547 (local-only path: `../scripts/diagnose_weight_limit.py:547`), and
diagnose_weight_limit.py:577 (local-only path: `../scripts/diagnose_weight_limit.py:577`).

The controller validates source files independently, then discards the source
objects during finalization. Its expected bin/field/order is copied from the
report being checked. The batch input hash is an internal checksum; it is not
compared with the corresponding authoritative source arrays/scalars.

Confirmed probes on separate copies:

- Swapping the two bin-0/order-4 population labels, preserving all arrays and
  the complete sixty-identity inventory, **passes full finalization** and emits
  the global negative conclusion.
- A self-consistent trajectory from masses multiplied by 1.01, retaining the old
  claimed provenance, passes batch validation.
- Removing the 1024 checkpoint from a copied order-4 batch, updating its own
  attempts, cap, shapes and hashes, and removing its optional fixed-mesh screen
  **passes full finalization** while its summary still declares the full
  checkpoint list through 1024 (`omitted-checkpoint-probe.json`). The controller
  must bind checkpoints to declared attempted scope, not reconstruct that scope
  from the surviving arrays.

The independent review compared every original batch's actual masses, variance,
magnitudes and four scalars with the appropriate saved source inputs: all 60
match. Thus this is a demonstrated validator defect, **not evidence that the
submitted unmodified trajectories contain swapped or fabricated inputs**.
Step 12's separate per-spectrum validator remains unchanged and R2 stays closed.

**Closure:** reconstruct expected source arrays/scalars and identity from the
verified immutable context; require comparison per population/quantity, bind
source metadata and original attempts/caps, validate derived historical and
summary fields, and reject all mutations through both isolated and complete
bundle routes. A correctly recorded failed/capped attempt must remain inspectable
without being counted as completed or converged. Preserve all original evidence.

Documentation correction: the claimed twelve order-64 source/profile matches
are actually **four order-32 matches** (both populations in bins 0–1) and
**eight order-64 matches** (bins 2–5). All sixty inputs do match their diagnosis
sources; these are distinct provenance comparisons.

### S13-R3 — medium: a representable diagnostic still fails through a scalar ratio

Locations: weight_limit.py:141 (local-only path: `../fishhighz/validation/weight_limit.py:141`),
plus the corresponding spectrum/eigenvector expressions at lines 315, 363,
378 and 384.

`log(length / pixel)` forms the potentially unrepresentable ratio before taking
its logarithm. For the finite positive one-cell inputs
r=1e300, v=1e-300, L=1e-300, Delta_v=1e300, S=B=1, the exact values are
w0=1/2, d=1, w1=1/3 and **A=P_pixel=1**. The diagnostic instead raises at
`log(0)` because L/Delta_v underflows. This is an artificial range control,
not a realistic survey input; it disproves the current unrestricted finite-
positive domain claim without affecting the saved DESI numerical conclusions.

**Closure:** retain scalar products/ratios in logarithmic form too, or explicitly
validate and document a narrower numerical domain before arithmetic. Add this
control and an overflow counterpart with a Decimal oracle; retain the production
path and its guards. Add an actual prototype-versus-Decimal comparison past
underflow, rather than only Decimal precision against itself.

## 3. New verification and preserved evidence

- `PATH="$PWD/.venv/bin:$PATH" scripts/check.sh`: **1303 passed, 25 skipped,
  187.73 s**; Ruff lint/format passed. Nonfatal site MUNGE messages occurred after
  pytest; the command exited zero. The suite ran with one-thread settings.
- Existing exact installed wheel, isolated imports outside the checkout:
  **121 affected regressions passed**, plus **all six examples**. This is a
  rerun of the existing installation, not a new build/install.
- All **53 source/wheel/installed Python modules** and METADATA/WHEEL bytes
  match. All **52 r6 modules are unchanged**; the sole added package module is
  `validation/weight_limit.py`. All five handoff source-hash entries match.
- Wheel SHA256:
  `01956a930d288cb25b0f4b8d9eb3bedb246453a085ad301ecb11b91e440f3d93`.
  Original summary SHA256:
  `94f069593115e521e0787789195037708fccc9506371a5df954c3894288c6aee`.
- Independent backward-distribution coefficients, extended-precision finite
  trajectories, Decimal controls and source comparisons are recorded in
  `independent-results.json`; full-finalization corruption probes have their own
  explicitly labeled directories. No corrupted copy is scientific evidence.
- The first numerical-check script ended at the deliberately induced range
  exception before serializing its results. A second invocation catches and
  records that expected failure; both logs are preserved. No production repair
  was applied to make the check pass.

The implementation's original 840-checkpoint run and 0.812 s maximum trajectory
time remain historical implementation evidence. The quick review does not
repeat a real forecast or infer new Fisher convergence. The inherited NumPy
matrix-speed target and explicit-Numba/SciPy-blocked limitation remain unchanged.

## 4. Documents and next action

[Revision 2](../notes/IMPLEMENTATION_STEP.md) retains this same diagnostic step and
adds concrete closure tests for S13-R1/R2/R3. The design, roadmap and package
standing status are synchronized. The scientific handoff and revision-1 evidence
are preserved as historical submissions, including their overstrong conclusion.
README's statement that Step 13 is deferred should be corrected by the implementer
when the revised diagnostic documentation is written.

No scientific prescription must be selected merely to perform these repairs.
The user controls revision-2 approval/dispatch, feedback, scientific acceptance
and progression. Do not advance to Step 14 or infer authorization for a real run.
