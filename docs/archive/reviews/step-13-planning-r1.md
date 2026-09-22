# Step 13 revision 1: planning assessment of the normalized weighting limit

The user selected investigation of a normalized/asymptotic formulation of the
existing cumulative rule. This authorizes Step 13 planning and progression of
the investigation; it does not select another estimator, waive convergence, or
scientifically accept the Step 12 accuracy forecasts.

The live implementation is `fishhighz/kernels/weights.py` and its caller
`fishhighz/weights.py`. Put r_i=rho_i*q_i, v_i=variance_i,
alpha=L/Delta_v, a_i=v_i/S, and B=alias. On positive density support it uses

    w_i(0) = 1/(1 + Delta_v*v_i/B)
    N_i(t) = alpha * sum(j<=i, r_j*w_j(t))
    w_i(t+1) = N_i(t)/(N_i(t)+a_i).

The last expression assumes positive denominator; zero-variance and zero-density
branches must retain the actual caller's semantics. The noise coefficients are

    A = sum(r*w*w) / (L * sum(r*w)^2)
    P_pixel = Delta_v * sum(r*w*w*v) / (L * sum(r*w)^2).

They depend on weight shape, not common amplitude. A is deg^2; P_pixel is
deg^2 km/s. Their invariance under w -> c*w does not make the nonlinear update
invariant: its denominator retains the absolute amplitude through N.

## Exact finite-iteration reformulation to investigate

For positive variance, write w=s*u, with max(u)=1 and h_i=alpha*sum(j<=i,r_j*u_j).
Then b_i=h_i/(a_i+s*h_i), b_max=max(b), and

    s_next = s*b_max
    u_next = b/b_max.

This is an algebraic identity for the existing finite recurrence, not a new
weighting rule. A numerically stable realization must retain s, preferably via
log(s), and must also retain or bound small relative components of u. Merely
resetting the weight maximum to one after each original update changes the
recurrence. A log-domain formulation or another demonstrably equivalent method
belongs first in a diagnostic prototype. Existing production guards must remain.

## Why iteration convergence is insufficient

Near zero amplitude on a fixed grid, the linearized map is
M_ij=alpha*r_j/a_i for j<=i. Its spectrum is its diagonal, but the limiting
normalized shape also depends on repeated/near-equal eigenvalues, accessibility
from the initial weights and nonnormal transients. Spectral radius below one
proves decay under the established comparison bound; it does not establish
finite noise coefficients under grid refinement.

The distinction is visible in an analytic control. On a continuous unit interval
with constant density, variance and positive initial weight, the linearized
prefix operator is Kf(x)=c*integral(0..x,f(y)dy). Starting with f=1 gives a
normalized shape u_t(x)=x^t. If the integrated density is R, then

    A_t = (t+1)^2 / ((2*t+1)*L*R)
    P_pixel,t = Delta_v*v*A_t.

A_t grows approximately as t/(2*L*R), even though the normalized shape remains
bounded. This is a control for the linearized continuum operator, not a proof
about the full nonlinear DESI-2 recurrence. It shows why a finite useful limit
must be established rather than assumed. The step must examine iteration-first,
refinement-first and joint limits, with a declared continuous magnitude measure.

Planning-only checks in `.validation/step13-plan-r1/analytic_controls.py` verify
15 scalar nonlinear recurrences against a 90-digit closed form and seven
linearized-continuum integrals against Gauss-Legendre quadrature. They pass;
no model, forest input preparation or forecast was evaluated. These controls
inform the assignment and do not resolve R3 for the real samples.

## Recommended investigation sequence

1. Reproduce and document the finite recurrence, noise ratios and units.
2. Derive the exact amplitude/shape evolution and analytic counterexamples.
3. Build a small independent diagnostic evaluator and demonstrate finite-step
   equivalence, including through the original underflow boundary.
4. Analyze discrete eigenstructure and continuum concentration before expensive
   iteration scans; distinguish evidence for a limit from apparent plateaus.
5. Apply bounded saved-array diagnostics to the two forest populations in all
   six 15x2pt bins; compare coefficients and concentration measures across the
   existing magnitude orders, without regenerating the survey forecast.
6. Report a supported finite limit, demonstrated lack of a suitable limit, or
   a precisely bounded unresolved result. Production adoption and real forecast
   verification follow only after user review of the scientific conclusion.

The [Step 13 instructions](../notes/IMPLEMENTATION_STEP.md) define the detailed tests,
limits and handoff. The previous active instructions are preserved exactly in
[step-12-instructions-r6.md](step-12-instructions-r6.md). The former prospective
performance/documentation steps move to 14/15; no detailed work is added to them.
