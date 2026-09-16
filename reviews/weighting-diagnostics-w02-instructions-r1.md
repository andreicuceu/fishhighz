# Forest weighting diagnostic W02, revision 1: equations and source audit

Status: instructions prepared at the user's request; awaiting user approval and
separate implementation-agent dispatch. W01's finite-update scientific result is
independently supported. Its remaining evidence-validation findings are preserved,
not resolved or prerequisites for W02. Package Step 13 and scientific acceptance
of the forecasts are unchanged.
Created: 2026-09-15.

## 1. Scientific question and completion condition

What weighting prescription is actually specified by McDonald & Eisenstein
(2007) and Font-Ribera et al. (2014), and which parts agree with or differ from
the live lyaforecast and FishHighz recurrences?

Resolve the magnitude integration domain, scalar versus prefix-dependent I1,
normalization, the signal used in the nonlinear update, and the distinction
between weight preparation and final aliasing/pixel-noise coefficients. Determine
what can be established about the original C++ implementation from available
source. Separate published equations, executable code, comments, and inference.

Completion means a source-supported equation comparison and the few analytic or
tiny numerical checks needed to substantiate it. An explicitly unresolved
historical implementation or paper ambiguity is a valid outcome. This audit
does not have to identify the cause of the real-survey convergence failure.

## 2. Read only the relevant context and source

Read workspace/package AGENTS.md, design section 0, the main roadmap handover
and progress register, and the [diagnostic roadmap](../../FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md).
Read [W01's review](reviews/weighting-diagnostics-w01-review-r2.md), especially its
scientific result, and its implementation report as needed. The archived
[W01 revision-3 instructions](reviews/weighting-diagnostics-w01-instructions-r3.md)
are historical; do not execute their evidence-repair assignment.
`IMPLEMENTATION_STEP.md` belongs to package Step 13 and is context only.

Inspect these live source paths, following only dependencies needed to resolve
an equation, unit, or caller convention:

- `../lyaforecast/lyaforecast/weights.py`: initialization, iteration, I1/I2/I3,
  effective pixel density and magnitude-dependent noise.
- `../lyaforecast/lyaforecast/covariance.py`: conversion of integrals into
  aliasing/pixel noise and selection of the final magnitude element.
- `../lyaforecast/lyaforecast/power_spectrum.py` and `survey.py` only as needed
  to establish whether evaluated power already contains aliasing, the response,
  and the magnitude measure. Verify actual filenames before following them.
- `fishhighz/kernels/weights.py`, with `fishhighz/weights.py` and `noise.py` only
  where host validation, support or units affect the comparison.

Start with Git status and preserve dirty/untracked files. File/line references
must describe the live checkout, not an assumed GitHub version.

Primary literature:

1. [McDonald & Eisenstein (2007)](https://arxiv.org/pdf/astro-ph/0607122),
   section II.B, equations (13)–(19) and the surrounding weighting discussion.
2. [Font-Ribera et al. (2014)](https://arxiv.org/pdf/1308.4164), section IV.B,
   equations (23)–(29) and the paragraph defining PS and PN.

Record the version actually read; distinguish the available arXiv text from a
journal version not inspected. Recheck equations directly rather than treating
our earlier literature summary as authoritative. The later
[McQuinn & White prescription](https://arxiv.org/pdf/1102.1752), equations
(11)–(13), is optional context for interpreting initialization; evaluating or
adopting it as a survey alternative is outside W02.

For historical C++ source, make one bounded search of known local forecast
checkout filenames and their existing documentation/history. Follow at most two
explicit public source links if found; cap this search at five minutes. Do not
scan shared filesystem roots, clone unrelated repositories, contact authors or
compile/run historical code. If unavailable, give the searched locations and
state that code ancestry does not establish the original recurrence. If found,
record source identity and inspect only the relevant weight/integral routines.

## 3. Required equation comparison

Use one notation and give its mapping to each source. For a discrete magnitude
measure, let r_i = (dn/dm)_i q_i, v_i = sigma_N,i^2, forest length L, pixel width
l_p, S the 3D weighting signal and B the 1D power at the weighting mode. Keep
transverse area and radial-length units explicit, including the code's angular
and velocity convention. Do not equate 1D noise v_i l_p with 3D noise.

Write down and check against the source:

- Initialization and its normalization.
- Prefix I1,i = sum_{j<=i} r_j w_j versus full-sample I1 = sum_j r_j w_j.
- N_i = (L/l_p) I1,i or the common N = (L/l_p) I1; the corresponding update
  w_i' = S N_i / (S N_i + v_i), or with N for all magnitudes.
- Full-sample I2 = sum_i r_i w_i^2 and I3 = sum_i r_i v_i w_i^2;
  A = I2/(L I1^2), P_pixel = l_p I3/(L I1^2).
- Where S is held fixed, where B enters, and whether an aliasing term is present
  in S. If an aliasing-inclusive interpretation gives S = P3D + A[w] B,
  identify it as such, specify the mode/response and its weight dependence,
  and distinguish it from a fixed S. Do not assume this is an implemented rule.

Explain whether cumulative arrays merely report forecasts for different limiting
magnitudes, or feed each prefix back into that same magnitude's weight. A family
of samples with separately solved weights is not automatically equivalent to
one triangular recurrence. Inspect the power provider as well as comments before
claiming that aliasing is absent or present.

Do not infer a change of physical prescription solely from the less explicit
2014 wording about S. Preserve alternatives when the source does not resolve
one. Document positive-support restrictions in FishHighz separately from the
literal signed legacy inputs in W01; do not floor or otherwise modify them.

## 4. Minimal independent scientific checks

Use at most one short standalone diagnostic script with ordinary assertions;
a separate pytest file is optional only if it adds an independent check. No
production code changes or reusable solver/validation framework are required.
Analytic calculations can be written in the report. Numerical checks should
complete in seconds, with at most three magnitude bins and a few updates.

1. **Distinguish the two integration meanings.** In dimensionless toy units use
   r=(1,1), v=(1,4), L=l_p=S=B=1 and w=(1/2,1/5). Derive one prefix and one
   full-sample update by hand or exact rational arithmetic, then compare any
   numerical implementation against those independent values. The bright-bin
   updates must differ and the final-bin updates must agree for this first
   update. This tests an equation difference, not a real-survey forecast effect.
2. **Check normalization where it matters.** On that positive example, multiply
   all initial weights by a fixed positive factor, e.g. 2. Verify that A and
   P_pixel are unchanged while the next fixed-S nonlinear update generally
   changes. Do not normalize between iterations and call the recurrence unchanged.
3. **Check that a scalar integral alone does not guarantee a nonzero fixed point.**
   For one bin, L=l_p=S=v=1 and r=1/2, 1, 2, derive the fixed points of
   w' = r w/(r w+1), their physical admissibility and local stability. Explain
   the marginal r=1 case without extrapolating a finite trajectory. A short
   substitution/derivative check suffices; no long iteration experiment is needed.

Require exact agreement for rational identities, or rtol=1e-12 and atol=1e-14
for these order-unity float64 examples. A failed analytic identity or unit
mapping must be resolved or explicitly restrict the conclusion. Do not add
large parameter sweeps, precision ladders, or artificial failure-path tests
unless a concrete numerical ambiguity in these checks requires one.

These checks support only the audited equations and their stated toy limits.
They do not prove continuum convergence, survey-wide optimality, or that either
source discrepancy explains W01. No saved-survey recurrence, magnitude refinement,
physical input preparation, Fisher sum or new forecast is authorized in W02.

## 5. Minimal implementation and handoff

Allowed new deliverables are `reviews/weighting-diagnostics-w02-r1.md` and, if
needed, `scripts/check_forest_weight_equations.py` plus one focused test file.
Keep small outputs in a new `.validation/forest-weight-diagnostics/w02-r1-<UTC>/`
directory if useful. Reuse the existing interpreter. Use one process and set
OMP_NUM_THREADS, OPENBLAS_NUM_THREADS and MKL_NUM_THREADS to 1 for numerical runs.
Do not modify production code, W01 tooling/tests/evidence or historical reports.

The report must contain:

- A short scientific conclusion with its assumptions and remaining uncertainty.
- A compact paper/code comparison table covering the items in section 3, with
  paper version/page/equation and live source/function/line references.
- The explicit recurrences, units, analytic examples and actual check results.
- Historical C++ availability and the limits of any inference about its behavior.
- Exact commands and the relevant source revisions/dirty state; record hashes
  only for the few local files that support the comparison if needed to identify
  uncommitted source. Distinguish earlier W01 evidence from new W02 checks.
- Whether the evidence justifies a subsequent controlled comparison, phrased as
  a scientific question for the user, not instructions for another numbered step.

Run the tiny script/assertions, any genuinely necessary focused test, and Ruff
on newly written Python only. No full test suite, new wheel/environment,
packaging check, JSON schema, artifact replay system, provenance mutation suite,
or general evidence-validator repair is required. Missing historical source does
not justify expanding scope. Stop for independent/user review.

## 6. Independent review: scientific consequence determines further work

Review only the source, derivations and checks needed to answer W02. Read the
cited equations and decisive source calls directly; independently reproduce the
small examples rather than only rerunning the implementer's functions. No broad
regression suite or W01 evidence-hardening work belongs in this review.

Propose further changes only when they are likely to affect the scientific
conclusion or its justified scope. Examples include a wrong integration domain,
unit conversion, hidden aliasing contribution, invalid fixed-point argument,
incorrect source identity, or an unsupported claim about historical code.
A hypothetical validator weakness, formatting preference, larger test inventory,
or stronger provenance machinery alone is not a reason for another revision.
Use a narrower conclusion or an explicit uncertainty when that answers W02
without further implementation. Do not turn optional hardening into a to-do list.

For every proposed correction, state the evidence, the smallest necessary
change/check, and how it could change the conclusion (or its applicability).
The review must lead with a short **Scientific conclusion and impact of proposed
changes** summary. If no scientifically consequential changes are needed, say
so explicitly. Write `reviews/weighting-diagnostics-w02-review-r1.md` and stop.
Revise this same assignment only for consequential findings or explicit user
feedback, preserving the prior instructions. Do not silently fix code.

The user retains approval, dispatch, acceptance and progression. Do not select a
production prescription, weaken scientific tolerances, run a full forecast,
execute package Step 13, add future numbered steps, dispatch agents, commit or push.
