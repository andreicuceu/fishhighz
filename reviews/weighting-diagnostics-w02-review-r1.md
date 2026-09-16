# Scientific conclusion and impact of proposed changes

**W02 revision 1 passes independent scientific review within its stated scope.
No scientifically consequential changes are needed.** The live implementations
feed magnitude-prefix integrals back into the weights; replacing these by a
common sample integral changes the update even with two positive bins. Their
weighting signal contains no sampling-aliasing contribution. Agreement of the
final noise-coefficient formulas therefore does not establish agreement of the
weight preparation. The report correctly leaves the historical C++ recurrence,
the aliasing-signal iteration convention, and the survey-scale consequences
unresolved. No repair revision or new weighting prescription is proposed.

Date: 2026-09-15. Reviewed assignment: W02 r1, now preserved in
[the exact instruction archive](weighting-diagnostics-w02-instructions-r1.md).
Reviewed deliverable: [implementation report](weighting-diagnostics-w02-r1.md)
and `scripts/check_forest_weight_equations.py`. No substantive user feedback was
supplied; the prompt retained its feedback placeholder.

**Awaiting user review and acceptance.** This pass neither accepts the forecasts
nor authorizes progression. Package Step 13 and its acceptance status are unchanged.

## Direct literature and source assessment

I read the equations and adjacent discussion directly in the available arXiv
PDF text, rather than relying on the implementation's literature summary.
The inspected versions are astro-ph/0607122v1 (7 July 2006; the PDF title page
also prints August 20, 2018) and 1308.4164v2 (15 May 2014). Journal versions were
not inspected.

- **Integration and normalization:** ME07 section II.B, printed pp. 4–5,
  equations (14)–(19), defines one sample integral with upper bound m_max and
  normalizes the weight to unity in the low-noise limit. Its effective pixel
  density is common to that sample. The adjacent discussion includes aliasing
  in the representative signal. It does not specify an initial iterate or
  enough algorithmic detail to reconstruct how that signal is initialized or
  refreshed. Treating a weight-dependent signal as an interpretation, rather
  than a verified historical algorithm, is justified.
  [McDonald & Eisenstein](https://arxiv.org/pdf/astro-ph/0607122).
- **Later formulation:** FR14 section IV.B, printed pp. 9–10, equations
  (23)–(29), retains scalar integrals over the observed luminosity function and
  the same coefficient structure. Undisplayed integration bounds do not imply
  a magnitude-dependent upper bound. Its signal wording does not establish a
  deliberate removal of aliasing from the earlier prescription.
  [Font-Ribera et al.](https://arxiv.org/pdf/1308.4164).

Here LF denotes `../lyaforecast/lyaforecast/`, and FH denotes `fishhighz/`,
relative to the package root. The decisive live calls are:

| Question | Directly inspected source | Assessment |
| --- | --- | --- |
| Do prefixes enter the update? | LF `weights.py:136–156`, `:184–206`, `:255–273`, `:339–358`; FH `kernels/weights.py:9–20` | Yes. I1 is an array of prefixes and enters the corresponding magnitude's denominator. These are not separately solved limiting-magnitude samples. |
| What fixes the amplitude and update count? | LF `weights.py:104–133`; FH `kernels/weights.py:12–19` | Initialization is B/(B+lp v); LF makes three updates. FH takes an explicit count. Neither rescales between updates. |
| Is aliasing hidden in S? | LF `power_spectrum.py:109–183`, `analytic_biases.py:92–121`; FH `weights.py:89–152` | The LF path returns intrinsic matter power times bias and instrumental transfer factors. The FH sampler applies response/conversion to the supplied intrinsic auto. Neither adds aliasing. The LF comment at `weights.py:151` contradicts its executable path. |
| How is final noise formed? | LF `covariance.py:287–294`, `:352–363`; FH `kernels/weights.py:23–37`, `noise.py:52–85` | Final full-sample I1/I2/I3 form A and P_pixel. Aliasing receives smoothed P1D; pixel noise is unsmoothed. |
| What measure and domain are used? | LF `survey.py:26–29`, `weights.py:158–204`; FH `weights.py:33–48`, `:213–290` | LF uses a finite linspace with a full dm at each node and source-redshift velocity conversion. FH requires explicit positive quadrature and nonnegative density/variance, with positive weighted support. |

The dimensional mapping is consistent: I1 has units deg^-2 (km/s)^-1;
N=(L/lp)I1 has the same units. Thus v/N and P_pixel are deg^2 km/s,
whereas lp v and B are km/s and A is deg^2. The comoving power conversion
is d_deg^2/a_v, inverse to the auxiliary S conversion. LF's transfer in
`spectrograph.py:279–304` is a sinc times Gaussian, with two forest legs for
auto-power. No additional response or noise term is hidden in these calls.
The report properly restricts equivalence to the supported inputs and does not
apply positive-measure conclusions to W01's signed legacy example.

## Independent minimal checks and results

I independently evaluated exact rational expressions in an inline Python
process, then compared one two-bin update through each live implementation.
I did not import or run the implementer's diagnostic functions. Its short
script was inspected; rerunning it would add no necessary scientific evidence.

For r=(1,1), v=(1,4), L=lp=S=B=1:

| Quantity | Independent result |
| --- | --- |
| Initial weights | (1/2, 1/5) |
| One prefix update | (1/3, 7/47) |
| One full-sample update | (7/17, 7/47) |
| Full minus prefix, bright bin | 4/51 |
| A and P_pixel before update | (29/49, 41/49) |
| A and P_pixel after doubling initial weights | (29/49, 41/49) |
| Prefix update from doubled weights | (1/2, 7/27) |
| Full update from doubled weights | (7/12, 7/27) |

For the one-bin map f(w)=rw/(1+rw), subtraction gives
f(w)-w=w[(r-1)-rw]/(1+rw). Its fixed points are zero and 1-1/r;
the latter is physical only for r>1. The slopes are r and 1/r.
Exact substitutions and slope evaluations pass for r=1/2, 1, 2.
At r=1, the identity 1/f(w)=1/w+1 yields
w_t=w_0/(1+t w_0) for positive w_0. This proves algebraic attraction
to zero without extrapolating a finite trajectory. The original coefficient
ratios at exactly zero remain undefined. This verifies the report's distinction
between a nonzero fixed point and a possible limiting noise ratio.

Execution used the existing interpreter, from `lib/fishhighz`:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -B - <<'PY'
from fractions import Fraction as Q
import runpy
import numpy as np
w = [Q(1,2), Q(1,5)]
v = [1,4]
def update(x, prefix):
    return [sum(x[:i+1] if prefix else x)/(sum(x[:i+1] if prefix else x)+v[i]) for i in range(2)]
def coefficients(x):
    return (sum(t*t for t in x)/sum(x)**2, sum(v[i]*x[i]**2 for i in range(2))/sum(x)**2)
assert update(w,True)==[Q(1,3),Q(7,47)]
assert update(w,False)==[Q(7,17),Q(7,47)]
assert coefficients(w)==coefficients([2*t for t in w])==(Q(29,49),Q(41,49))
assert update([2*t for t in w],True)==[Q(1,2),Q(7,27)]
assert update([2*t for t in w],False)==[Q(7,12),Q(7,27)]
for r in [Q(1,2),Q(1),Q(2)]:
    for x in [Q(0),1-1/r]:
        assert r*x/(1+r*x)==x
        assert r/(1+r*x)**2==(r if x==0 else 1/r)
for x in w:
    assert 1/(x/(1+x))==1/x+1
kernel=runpy.run_path('fishhighz/kernels/weights.py')
actual,_=kernel['_iterate'](np.ones(2),np.array(v,dtype=float),1,1,1,1,1)
np.testing.assert_allclose(actual,[float(x) for x in update(w,True)],rtol=1e-12,atol=1e-14)
LF=runpy.run_path('../lyaforecast/lyaforecast/weights.py')['Weights']
obj=LF.__new__(LF)
obj._p3d_w=1
obj._forest_length=obj._pix_kms=1
obj.maglist=np.array([0.,1.])
obj._zq=2
obj._lya_tracer=object()
obj._get_dn_dkmsdm=lambda *args: np.ones(2)
obj._get_pix_var_m=lambda: np.array(v,dtype=float)
np.testing.assert_allclose(obj._compute_weights_lya(np.array([float(x) for x in w])),actual,rtol=1e-12,atol=1e-14)
print('PASS: independent exact arithmetic; one-bin roots/slopes and marginal identity; one two-bin update through each live implementation.')
print('prefix',update(w,True),'full',update(w,False),'coefficients',coefficients(w))
PY
```

All assertions passed, exit status 0, approximately 1.8 seconds wall time.
The wrapper emitted two MUNGE socket messages after the successful numerical
output; no scheduler command was issued. The LF object uses supplied toy density
and variance, bypassing all survey/model/input preparation.

## Historical source, identity, and retained uncertainty

A bounded read of existing LF history confirmed that the two cited commits
`c904b912107946f33b951a8dd0f1917b3f1aa0a7` and
`8d5286f51adb357faa37832215dff620126a3efd` discuss comparison with C++ but list
Python/notebook changes. `git log --all --oneline -- '*.cpp' '*.cc' '*.cxx' '*.h'`
returned no paths. This supports the report's limited conclusion: the historical
recurrence has not been verified. It does not prove that source is unavailable
elsewhere. No broader historical search is needed for that explicitly uncertain
conclusion.

FishHighz HEAD is `0d69786a06d5d676564a51fad14a7156951c2228`;
LF HEAD is `5abe8bcf8d12cc31d1f5ecd89c87e739c6a14b81`, with a clean LF checkout.
FishHighz already had tracked changes to AGENTS.md, IMPLEMENTATION_STEP.md,
README.md and pyproject.toml, and extensive untracked source/reports/tests.
The three scientific-source SHA-256 values in the implementation report match
live files. The separate IMPLEMENTATION_STEP.md hash also matches
`3675cd861fc7909d464d9f7e60f7725f1a9e2f5cb3fa8a9f9849f3a30552c310`.

W01's reviewed finite-update sensitivity is historical context only. Its
W01-R4/R5 findings remain deferred. W02 neither measures a real-survey effect
nor establishes continuum convergence or optimality. The proposed controlled
comparison is scientifically motivated as a question for the user, without
selecting its signal convention or authorizing execution.

## Review disposition and document changes

No consequential finding requires revision of the scientific assignment.
The original W02 r1 instructions were archived byte for byte before updating
its status. The diagnostic roadmap/register, design section 0, main roadmap's
diagnostic handover paragraph and package AGENTS.md now record this review and
pending user acceptance. No package Step 13 text or acceptance state was changed.
Historical reports, implementation code and the implementer's W02 script/report
were preserved. Only this review, the instruction archive and those diagnostic
status entries were written.

No production edit, survey calculation, Fisher forecast, broad test suite,
validation infrastructure, W01 repair, Slurm action, agent dispatch, commit,
push or subsequent-step preparation was performed. Stop for user review.
