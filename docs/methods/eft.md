# Full EFT external provider

`fishhighz.models.eft.EFTModelAdapter` wraps the current prepared
`eftkernels.EFTModel` API for forest and selection-free galaxy 3x2pt. Install
eftkernels separately (or the optional `eft` extra); importing FishHighz itself
does not import it. The native INI model remains the existing built-in route.

The adapter explicitly selects `master_integrals="analytic"`: the bundled
scalar master tables use a frequency spacing different from the current
FFTLog grid. Analytic masters are evaluated at the actual prepared exponents.
The legacy table option and its frozen numerical answers remain untouched.
Independent smooth halo SPT tests pass at 1e-4 of linear power; remaining
forest selection-operator differences and high-k validity are documented in
the study report, and its forest/joint forecasts remain exploratory.

The adapter supplies the ordinary external provider signature with paired
`k, mu` nodes and `(node, pair)` output. Fields are F=0 and G=1. It uses
dictionary parsing for the halo mapping, and changes growth with
`update_fiducial`, keeping the linear spectrum and sigma8 fixed. `Drat=1`.

Smooth and wiggle components have independent `alpha_s, phi_s` and
`alpha_w, phi_w` coordinates. For either component,
`q_parallel=alpha/sqrt(phi)`, `q_perpendicular=alpha*sqrt(phi)`, and the power
Jacobian is `alpha**-3*phi**-.5`. This is not the `alpha_iso_phi` convention.
The wiggle loops are the difference of full and smooth loop evaluations;
generating loops from wiggle input alone would lose mixed terms.

The loop order includes every implemented forest bias, an independent galaxy
`b2, bG2, bGamma3`, independent tracer `cs, cs2, cs4`, and mixed `cx`.
The optional `cs6` angular extension and an additional galaxy k^4 velocity
counterterm are not included in the example study. Pair stochastic residuals
are independent polynomial coefficients, separate from known observational
noise; the legacy half-sum cross prescription is disabled. The example leaves
these residuals free under its bias/counterterm prior studies.

The tree option uses the EFT tree matter components with fixed fiducial IR
damping, multiplied by the exact varying K1 kernels. Growth information from
an assumed damping width is excluded. Consequently forest-only tree spectra
retain the exact f/eta-bias degeneracy. This differs from varying the package's
tree-only IR damping with f, and is an explicit forecasting convention.

`examples/eft_3x2_study.py` is a reproducible idealized local study. It writes
the complete means, Jacobians, Fisher matrices, SVD ranks/null directions,
dimensionless higher-EFT principal modes, priors and numerical settings.
Run with single-threaded BLAS; install matplotlib to build the companion
technical report. The study's sample is not a validated DESI forecast. High-k
model validity and remaining independent FFTLog loop residuals must be assessed
from the delivered report; regression agreement alone does not settle them.
