# Scientific conclusion and impact of proposed changes

**W05 revision 1 passes independent scientific review within its assigned
scope. No scientifically consequential correction is needed.** The explicit
reference nu=B/(B+l_p v) is algebraically identical to the legacy seed, including
amplitude. The positive-measure minimum of Q is correctly derived. Independent
arithmetic confirms that the original prefix t=3 and saved W04 fixed/refreshed
t=24 vectors have different shapes and respectively 14.3618%, 12.8358% and
19.1309% higher Q than nu on the unchanged signed grid. Refreshed Q exceeds
fixed Q by 5.57897%, despite its smaller aliasing contribution.

The signed-grid ranking establishes no physical optimum or forecast improvement.
W04's finite-count stability remains historical evidence, not a minimum or
continuum-convergence proof. No repair revision or scientific choice is required.
**Awaiting user review and acceptance; this pass grants no progression or
production adoption.** Package Step 13 and its acceptance status are unchanged.

## Assignment and inspected evidence

Reviewed the [W05 report](weighting-diagnostics-w05-r1.md) against the
[archived exact assignment](weighting-diagnostics-w05-instructions-r1.md).
The request supplied no substantive feedback; its feedback placeholder remained.
Read workspace/package AGENTS.md, design section 0, main roadmap handover/progress
register, diagnostic roadmap, W04 report/review and W02's response audit.
Inspected the complete W05 script, W04 seed/moment expressions, live legacy
weights.py:104–133 and :184–253, covariance.py:287–294 and :352–363, and the
small FishHighz weighting kernel. No production module was imported or changed.

[McQuinn & White, section 2.2, equations 11–13](https://arxiv.org/pdf/1102.1752)
uses inverse total sightline power, neglecting subdominant off-diagonal covariance;
the factorized weights and independent sightline noise are stated with equation 2.
Its arbitrary normalization B is distinct from this diagnostic's forest P1D B.
The report correctly maps individual instrumental power to l_p v, retains the
saved response in B/P0, and includes aliasing exactly once in Q. It makes no
claim that the legacy recurrence is the published inverse-variance estimator.

## Independent derivation and positive control

For r_i>=0 and D_i=B+l_p v_i>0, Cauchy–Schwarz applied to
sqrt(r_i D_i) w_i and sqrt(r_i/D_i) gives

```text
I1² <= [sum r_i D_i w_i²] K,       K=sum r_i/D_i,
Q >= 1/(L K).
```

Equality requires w_i=c/D_i on positive support, with nonzero c; zero-measure
nodes are immaterial. Choosing c=B gives nu, I1=B K,
n_eff=L B K, and Q[nu]=1/(L K)=B/n_eff. The report's completed-square identity
is also correct. Signed r invalidates the nonnegativity argument, even when
K and the measured coefficients happen to be positive. For nonnegative physical
power, minimizing Q minimizes the diagonal auto-power variance proportional to
(P0+Q)² at fixed mode count under the stated approximation only.

Dimensions agree: I1/I2/I3 are deg^-2 (km/s)^-1, K is deg^-2 (km/s)^-2,
A is deg², n_eff is deg^-2, and Q and P_pixel are deg² km/s.
Multiplying numerator and denominator of the legacy seed by l_p proves its
identity with nu; neither expression depends on P0 or an iteration count.

Independent exact Fraction arithmetic for the assigned two-bin example gives
nu=(1/2,1/5), K=n_eff=7/10, moments (7/10,29/100,41/100),
A=29/49, P_pixel=41/49 and Q=10/7. Uniform Q=7/4 is larger.
The direct numerator and effective-density expressions agree exactly, and
scaling nu by 3 preserves A, P_pixel and Q exactly.

## Independent saved-vector arithmetic

The review script imports neither implementation helpers nor recurrence code.
It converts exact binary inputs and saved weights to 80-digit Decimal values,
forms direct scalar sums, and independently computes nu in Decimal. All four
vectors' moments and coefficients reproduce the report and W05 arrays. Prefix
A/P_pixel also match the authoritative row's final coefficients; both t=24
moment/coefficient sets match the W04 NPZ and summary. Seed, normalized-shape
differences, all requested fractional contrasts, direct-numerator Q versus
A B+P_pixel, and both effective-density identities pass. Relative tolerance is
5e-12 with zero absolute tolerance and exact zero handling; fractional and
shape differences use an absolute 5e-12 bound. Maximum relative discrepancy
is **8.861e-16**. No tolerance was relaxed.

| Vector | A B [deg² km/s] | P_pixel [deg² km/s] | Q [deg² km/s] |
| --- | ---: | ---: | ---: |
| nu | 0.481006065095 | 0.314093507492 | 0.795099572586 |
| Prefix t=3 | 0.361263005838 | 0.548026906250 | 0.909289912088 |
| Fixed P0 t=24 | 0.356264911394 | 0.540892274852 | 0.897157186246 |
| Refreshed t=24 | 0.338140036238 | 0.609069323592 | 0.947209359829 |

Aliasing dominates only for nu; pixel noise dominates for the other vectors.
The independently evaluated refreshed/fixed Q contrast is
0.05578974827451095. All moments, D, K and coefficients are finite and positive;
K=2.84882279513e-5 and formal n_eff=20.5422664189 deg^-2.
The high-precision comparison excludes roundoff/cancellation as an explanation
for the percent-level ordering; it does not establish physical admissibility.

All 107 magnitudes, signed densities and variances match the authoritative
compatibility-bin-0 row and W04/W05 arrays exactly: z=[2.0,2.235], endpoints
16.1 and 26.75, nine negative densities. dm, L, l_p, P0 and B match both
summaries exactly. All ten source/input hashes recorded in W05 matched live
files before status synchronization, including the authoritative report hash
3e5200d3e5b1bab3e7570c61fc02651caa907ae30816036b0b1b323098a87ed6.
The exact assignment was archived before updating its status. The review script
uses that identical archive for subsequent hash checks.

## Execution and disposition

New evidence: `.validation/forest-weight-diagnostics/w05-review-r1-independent/`
contains `scalar_check.py` and `result.json`. From the package root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 30 .venv/bin/python -B .validation/forest-weight-diagnostics/w05-review-r1-independent/scalar_check.py
.venv/bin/ruff check scripts/compare_inverse_variance_forest_weights.py .validation/forest-weight-diagnostics/w05-review-r1-independent/scalar_check.py
```

The numerical check passed in 0.047 s internally, 1.18 s command wall time,
with one process and numerical thread limits of one. An earlier command used
a duplicated relative output path and failed before creating/running the script;
no numerical timeout or numerical failure occurred. Ruff passed after formatting
and import sorting of the review-only script. No implementation script was edited. After substituting the identical instruction
archive path and formatting the review script, the same bounded check passed
again in 0.035 s internally; Ruff and tracked `git diff --check` also passed.

Only review evidence/report, an exact instruction archive and diagnostic-status
text were added or changed. The W05 scientific assignment remains revision 1;
no repair or subsequent step is prepared. Historical reports, inputs, responses,
production code and unrelated dirty/untracked work are preserved. Package
IMPLEMENTATION_STEP.md remains at SHA-256
3675cd861fc7909d464d9f7e60f7725f1a9e2f5cb3fa8a9f9849f3a30552c310.
No trajectories, forecasts, broad suites, W01 repair, additional densities,
modes, populations, bins, magnitude refinement, Slurm, dispatch, commits or
pushes were performed. Stop for user review.
