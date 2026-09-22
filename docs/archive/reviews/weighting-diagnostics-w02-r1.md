# Forest weighting diagnostic W02 revision 1

## Scientific conclusion

The live lyaforecast and FishHighz legacy paths implement a triangular,
magnitude-prefix recurrence. A common full-sample integral gives a different
update even for two positive bins. Their final aliasing and pixel-noise
coefficients nevertheless have the same algebraic form as the published
coefficients. Agreement of those final formulas does not imply agreement of
the weights used in them.

The live lyaforecast weighting signal is smoothed intrinsic P3D, with no
sampling-aliasing contribution, despite its comment to the contrary. FishHighz's
auxiliary sampler likewise adds no aliasing. The historical C++ recurrence
remains unverified. The paper wording does not completely specify how an
aliasing-inclusive representative signal is initialized or refreshed.

The three requested analytic checks pass, as do independent float64 evaluations
and a tiny check of the live FishHighz array kernel. A scalar integral alone
does not guarantee a nonzero fixed point. These statements establish equation
differences and toy limits only: they establish neither the origin nor the size
of a real-survey forecast discrepancy, continuum convergence, or optimality.

Status: implementation submitted for independent/user review. The user's current
approval supersedes the assignment's awaiting-approval text; planning files were
not changed. W01 r3 remains archived and unassigned; package Step 13 is unchanged.

## Literature inspected and comparison

Primary sources read directly through the arXiv PDF text:

- [McDonald & Eisenstein](https://arxiv.org/pdf/astro-ph/0607122), here ME07:
  PDF marked **astro-ph/0607122v1, 7 July 2006**, 18 pages; printed title-page
  date August 20, 2018. Section II.B, printed pp. 4–5, equations (13)–(19).
- [Font-Ribera et al.](https://arxiv.org/pdf/1308.4164), here FR14:
  **1308.4164v2, 15 May 2014**, 41 pages; title-page date May 16, 2014.
  Section IV.B, printed pp. 9–10, equations (23)–(29).

Journal versions were not inspected. PDF screenshot requests failed with cache
misses; the equations and adjacent paragraphs above were available as PDF text.
The differing ME07 title-page date is recorded rather than interpreted as a
different arXiv revision. No optional McQuinn–White comparison was needed.

All local references below describe the live checkout. `LF` denotes
`../lyaforecast/lyaforecast/`; `FH` denotes `fishhighz/`.

| Quantity | ME07, pp. 4–5 | FR14, pp. 9–10 | Executable code |
| --- | --- | --- | --- |
| Magnitude integral | (14), (15), (18): bounds −∞ to sample m_max; scalar I1 | (26)–(28): observed luminosity function, no displayed bounds; scalar I1 | LF weights.py:184–253 uses cumsum with common dm. FH kernels/weights.py:16 uses prefixes; :31–33 retains cumulative integrals. |
| Update and normalization | (19): S/(S+PN), unity low-noise normalization; following paragraph PN=v lp/(I1 L) | (29) and following paragraph: same PN and implicit weight dependence | LF weights.py:135–156, :255–273, :339–358 uses a different N_i for every i. FH kernels/weights.py:9–20 matches this on positive support. No inter-update normalization. |
| Initial iterate | No explicit seed in this discussion | No explicit seed in this discussion | LF weights.py:122–133 and FH kernels/weights.py:12 use B/(B+lp v). |
| Representative S | Following (19): constant representative mode k=0.07 h/Mpc, mu=0.5; explicitly includes aliasing | Following (29): same representative k, mu; aliasing inclusion not explicit, refers back to ME07 | LF weights.py:104–109, :149–154 fixes intrinsic smoothed P3D. FH weights.py:89–152 samples intrinsic auto/P1D once. |
| Final coefficients | (16)–(18): A=I2/(L I1²), P_pixel=lp I3/(L I1²); (13): P3D+A B+P_pixel | (23)–(28): same decomposition and coefficients | LF covariance.py:287–294, :352–363 uses final elements. FH kernels/weights.py:23–37 uses final elements; noise.py:52–85 applies response/conversion. |

FR14's omitted bounds do not specify a running upper bound tied to the magnitude
being weighted. Its less explicit signal wording does not demonstrate a change
of physical prescription. A literal intrinsic-signal reading and an
aliasing-inclusive reading inherited from ME07 remain distinct interpretations.

## Explicit code recurrences, measure, and units

Let q_i be the magnitude quadrature, rho_i=dn/dm per volume,
r_i=rho_i q_i, v_i the dimensionless fractional-flux pixel variance, L the forest
length, and lp the pixel width. In the code's angular/velocity convention:

| Quantity | Units |
| --- | --- |
| rho_i | deg⁻² (km/s)⁻¹ mag⁻¹ |
| q_i, r_i | mag; deg⁻² (km/s)⁻¹ |
| w_i, v_i, L/lp | dimensionless |
| L, lp, B (1D power) | km/s |
| I1, I2, I3, N_i | deg⁻² (km/s)⁻¹ |
| S, v_i/N_i, P_pixel | deg² km/s |
| A | deg² |

In comoving coordinates, replace deg² by transverse length squared and km/s by
radial length: r and N are inverse volume, B is length, A is area, and S and
P_pixel are volume. In particular **v_i lp is 1D noise power**, whereas
v_i/N_i is a 3D quantity. Neither is generally the final P_pixel.

The live legacy measure is q_i=dm=maglist[1]−maglist[0] for every node, including
both endpoints. LF survey.py:26–29 constructs a finite linspace from the configured
bright to faint limits; this is a rectangular node sum, not an exact integral
from −∞ or trapezoidal endpoint weighting. LF weights.py:158–182 evaluates the
source density at z_q and divides dn/dz/dm by c/(1+z_q). This source-velocity
measure must not silently be replaced by a density at the forest redshift.
FH weights.py:33–48 provides the corresponding explicit conversion and :264–272
forms masses rho*q; quadrature is caller-supplied.

For fixed S and B the initialization and live update are

```text
w_i^(0) = B / (B + lp v_i)
I1,i^(t) = sum_(j<=i) r_j w_j^(t)
N_i^(t) = (L/lp) I1,i^(t)
w_i^(t+1) = S N_i^(t) / (S N_i^(t) + v_i).
```

The full-sample comparison replaces I1,i by I1=sum_j r_j w_j and N_i by a
single N for all magnitudes. For a separately solved sample ending at m_a,
every weight within that sample would share I1^(a). In the live triangular
recurrence, each weight instead uses its own prefix. The cumsum comment at LF
weights.py:203 describes magnitude-limit plotting, but the executable call
through get_np_eff_lya feeds those prefixes into the update. The arrays therefore
do more than report a family of limiting-magnitude forecasts.

After the chosen number of updates, the full-sample coefficients are

```text
I1 = sum_i r_i w_i
I2 = sum_i r_i w_i^2
I3 = sum_i r_i v_i w_i^2
A = I2 / (L I1^2)
P_pixel = lp I3 / (L I1^2).
```

LF executes exactly three updates (weights.py:114–119); FH requires an explicit
count. LF constructs coefficient arrays, then selects [-1] in the forest auto
total power. FH returns scalar final coefficients. Multiplying all weights by
c>0 leaves both coefficients invariant, but multiplies N by c and changes the
next fixed-S update. Initialization fixes an amplitude, not just relative weights.

### Signal, response, and aliasing ownership

LF weights.py:151 says “weights include aliasing as signal”; the following
assignment is simply signal_power=self._p3d_w. Tracing that variable resolves the
comment/code disagreement: power_spectrum.py:109–154 converts intrinsic P3D and
applies one field transfer per forest leg; :156–183 returns linear matter power
times analytic bias. analytic_biases.py:92–121 contains Kaiser factors, with no
sampling density or aliasing term. P1D is independently evaluated and smoothed
at power_spectrum.py:79–107. Thus B enters initialization but is not added to S
during these updates. This is executable behavior, not an inference from names.

The LF BAO mode is (kt,kp)=(2.4 deg⁻¹,0.00035 s/km), and its P1D option uses
(7,0.001), at weights.py:73–88. The implied comoving k,mu depend on redshift;
these constants are not identically the papers' fixed comoving mode. The response
is the sinc pixel transfer times Gaussian resolution transfer
(spectrograph.py:279–304), squared for forest auto power and P1D. The suffix
`smooth` here denotes instrumental smoothing.

FH sample_auxiliary (weights.py:89–152) uses the caller's intrinsic auto P3D and
independent P1D, the specified auxiliary coordinates, and W². It converts S
by a_v/d_deg², where a_v=dv/dchi and d_deg=dchi_perp/dtheta_deg. The kernel holds
the supplied S and B fixed. A caller could supply a different fixed S; the kernel
does not construct or refresh a weight-dependent aliasing contribution.

An **interpretation**, not an implemented rule, is
S[w]=P3D_obs(k_w,mu_w)+A[w] B_obs(k_parallel,w), with both powers evaluated at
the representative mode using consistent response and units. Recomputing this
after each update couples S to I1 and I2; fixing it at some reference weights
defines a different recurrence. The inspected paper discussions do not resolve
that algorithmic choice or the seed. No aliasing-inclusive recurrence was run.

Final code noise uses A times the response-smoothed P1D plus unsmoothed P_pixel.
Multiplication by d_deg²/a_v converts angular/velocity 3D power into comoving
volume units (LF covariance.py:354–363; FH noise.py:52–85). Final aliasing thus
exists even though it is absent from the live weighting signal.

FH host validation requires nonnegative density and variance, positive quadrature,
some positive density support, and positive S/B for legacy iteration
(weights.py:213–290). Zero-density bins get zero iterative weights. LF's literal
arithmetic has no corresponding positivity guard. Historical W01 signed inputs
must not be identified with these positive-support checks; nothing was floored.

## Hand calculation and independent tiny checks

All examples below use dimensionless units and at most two bins.

### Two integration domains

For r=(1,1), v=(1,4), L=lp=S=B=1, initialization gives w=(1/2,1/5).
The prefixes are (1/2,7/10), whereas the full integral is 7/10.
Hence one update gives

```text
prefix: ( (1/2)/(1+1/2), (7/10)/(4+7/10) ) = (1/3, 7/47)
full:   ( (7/10)/(1+7/10), (7/10)/(4+7/10) ) = (7/17, 7/47).
```

The bright-bin difference is 4/51; the faint-bin values agree for this first
update. They need not agree after subsequent updates, since the bright-bin
weight contributes to the next full integral. No later two-bin updates were run.

### Amplitude normalization

Initially I1=7/10, I2=29/100, I3=41/100, giving A=29/49 and P_pixel=41/49.
Under w→2w these integrals become 7/5, 29/25, 41/25, leaving both ratios
unchanged. But the next prefix update becomes (1/2,7/27) and the full update
becomes (7/12,7/27). Normalizing between updates therefore changes this recurrence.

### One-bin fixed points

With L=lp=S=v=1, f(w)=rw/(1+rw). Its roots satisfy
w[rw+1−r]=0, so w*=0 or 1−1/r. The latter exists as a positive physical weight
only for r>1. f'(w)=r/(1+rw)² gives the following exact values:

| r | Roots | Slopes at listed roots | Physical behavior |
| --- | --- | --- | --- |
| 1/2 | 0, −1 | 1/2, 2 | Zero attracts; negative root is inadmissible. |
| 1 | 0 (coincident roots) | 1 | Linear stability is marginal. |
| 2 | 0, 1/2 | 2, 1/2 | Zero repels positive perturbations; positive root attracts. |

At r=1, exactly 1/w_(t+1)=1/w_t+1 for w_t>0, giving
w_t=w_0/(1+t w_0). Zero is attracting from the physical positive side, with
algebraic rather than geometric decay. This follows analytically, not from a
finite trajectory. The zero-weight boundary makes the original coefficient ratios
0/0; a limiting ratio is a separate question, not a nonzero-weight solution.

The standalone script uses Fraction identities, separately evaluated float64
equations, one two-bin call to the live _iterate kernel, and _integrals on the
initial and doubled weights. Fixed-point substitutions are checked exactly and
numerically; derivatives and the marginal reciprocal identity are checked
exactly. Float tolerances are rtol=1e-12, atol=1e-14. All passed.

## Bounded historical source search

The search was confined to known local lib/lyaforecast and lib/fishhighz filenames,
lyaforecast README/source documentation, existing Git history for C/C++ file
extensions, and C++-related commit messages. It ended within five minutes, with
no explicit public C++ source link found, so zero public source links were followed.
No shared-root traversal, clone, compilation, historical execution, or recapture
was performed. The absent lib/lyaforecast/docs directory was recorded and not
searched elsewhere.

No .cpp/.cc/.cxx/.h files or matching historical paths were returned. Two existing
lyaforecast commits refer to comparisons with C++:
`c904b912107946f33b951a8dd0f1917b3f1aa0a7` (2018-10-05) and
`8d5286f51adb357faa37832215dff620126a3efd` (2017-06-28).
Their inspected statistics list Python/notebook changes, not C++ routines.
The initial README contains only the project description. These are evidence of
comparison activity, not evidence of the original recurrence or its equivalence.
Availability elsewhere remains unknown; the search is complete for W02's bounds.

## Execution and source identity

Only this report and `scripts/check_forest_weight_equations.py` were added.
No separate pytest file or artifact framework was necessary.

FishHighz HEAD: `0d69786a06d5d676564a51fad14a7156951c2228`.
Initial tracked modifications: AGENTS.md, IMPLEMENTATION_STEP.md, README.md,
pyproject.toml; extensive untracked implementation, scripts, tests, and historical
reports were already present, including all three FH source files cited below.
lyaforecast HEAD: `5abe8bcf8d12cc31d1f5ecd89c87e739c6a14b81`, initially clean.
Live source was read directly rather than inferred from HEAD.

SHA-256 values identifying the untracked scientific source inspected:

```text
3867c9eee6d43f22fca878b16c9f38bafc65a55061355d79828fcb135cc5e66f  fishhighz/kernels/weights.py
fdfcef33eed0955729f8b4b0086d25c4360f5bda4ed999850dfd9e0f8658aabf  fishhighz/weights.py
9cb580b1cb52860d3b65572bb18dfb707795b30a9646ff92778b930b429ef7b6  fishhighz/noise.py
```

IMPLEMENTATION_STEP.md SHA-256 before writing:
`3675cd861fc7909d464d9f7e60f7725f1a9e2f5cb3fa8a9f9849f3a30552c310`.

Actual numerical/lint commands from lib/fishhighz:

```bash
.venv/bin/ruff format scripts/check_forest_weight_equations.py
.venv/bin/ruff check scripts/check_forest_weight_equations.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python scripts/check_forest_weight_equations.py
```

Ruff formatted one new file and reported all checks passed. The numerical command
exited 0 with Python 3.13.15 / NumPy 2.5.3, printing all roots/slopes and the exact
fractions above, followed by PASS. The command wrapper also emitted two MUNGE
socket errors after the successful script output; no scheduler command was issued.
The existing .venv interpreter points to the NERSC Python 3.13 installation;
no environment or dependency was changed.

Final checks also passed:

```bash
.venv/bin/ruff format --check scripts/check_forest_weight_equations.py
.venv/bin/ruff check scripts/check_forest_weight_equations.py
sha256sum fishhighz/kernels/weights.py fishhighz/weights.py fishhighz/noise.py IMPLEMENTATION_STEP.md
git diff --check
git status --short -- reviews/weighting-diagnostics-w02-r1.md scripts/check_forest_weight_equations.py
git -C ../lyaforecast status --short
```

The four hashes matched those recorded above. The two deliverables are new
untracked files, lyaforecast remains clean, and `git diff --check` emitted no
errors (it checks tracked changes, not these untracked deliverables).

Decisive read-only/source-search commands from the workspace root included:

```bash
git -C lib/fishhighz status --short
git -C lib/lyaforecast status --short
git -C lib/fishhighz rev-parse HEAD
git -C lib/lyaforecast rev-parse HEAD
rg --files lib/lyaforecast lib/fishhighz | rg '\.(cpp|cc|cxx|h)$|README|HISTORY'
rg -n -i 'c\+\+|original|github.com|forecast.*code' lib/lyaforecast/README.md lib/lyaforecast/docs
git -C lib/lyaforecast log --all --oneline -- '*.cpp' '*.cc' '*.cxx' '*.h'
git -C lib/lyaforecast log --all --oneline --regexp-ignore-case --grep='C++\|original\|McDonald'
git -C lib/lyaforecast log --all --oneline -- README.md
git -C lib/lyaforecast show 4b42a35:README.md
git -C lib/lyaforecast show c904b91 --stat
git -C lib/lyaforecast show 8d5286f --stat
```

An initial equivalent `rg --files` invocation with `-g` filename filters was
blocked by the traversal hook before execution; the successful command above
retained the same two bounded checkout roots. Source reads used `nl -ba`,
`sed -n`, and bounded `rg -n` against the files cited in this report.

## Historical evidence, uncertainties, and review boundary

The pre-existing [W01 r2 review](weighting-diagnostics-w01-review-r2.md) supports
finite-update sensitivity on its signed 107-node legacy example. That evidence
was read as history, not rerun or repaired here. None of its survey/Fisher results
are new W02 measurements. W01-R4/R5 remain recorded and deferred.

Remaining uncertainties are the historical C++ implementation, the papers'
unspecified seed and aliasing-signal refresh convention, and the survey-scale
consequence of any equation choice. The positive toy checks do not resolve the
signed-input behavior or the failed forest convergence criteria.

The audit motivates a controlled scientific comparison, subject to the user's
decision: **for otherwise identical inputs and response, how much do a common
sample integral and the treatment of aliasing in S separately change the prepared
noise coefficients?** No such comparison is executed or assigned here, and no
weighting prescription is adopted. Stop for independent/user review.
