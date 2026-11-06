# Full EFT external provider

`fishhighz.models.eft.EFTModelAdapter` wraps the current prepared
`eftkernels.EFTModel` API for forest and selection-free galaxy 3x2pt. Install
eftkernels separately (or the optional `eft` extra); importing FishHighz itself
does not import it. The native INI model remains the existing built-in route.

This route requires the companion [eftkernels PR #1](https://github.com/andreicuceu/eftkernels/pull/1)
branch, including the local analytic-only API changes that remove the
`master_integrals` keyword. Scalar masters are evaluated analytically at the
runtime FFTLog exponents, with supported resolutions `nmax=128`, `256`, and
`512`. No master-integral tables or selection option are needed.

The EFT branch also corrects eight forest P13 factors. These corrections
propagate directly into the adapter's loop spectra; earlier forest forecasts
are not unchanged numerical references. See eftkernels'
`documents/forest_p13_correction.md` and `documents/scalar_master_integrals.md`,
with the scalar-master and halo SPT checks in `tests/test_master_integrals.py`
and `tests/test_analytic_masters.py`. Adapter regression agreement does not
establish FFTLog convergence or physical validity at high k; the forest and
joint forecasts remain exploratory.

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

## Reproducing the study and report

Run the following from the FishHighz checkout, with the companion eftkernels
checkout at `../eftkernels` and its reviewed API described above installed.
The commands below were checked against eftkernels commit
`72dbe1e05a6071527963b87225bc98b9f49669ae`. The validation generator is part
of that source checkout, rather than the installed wheel. It computes the
independent SPT quadratures; no precomputed loop-sum JSON is required.

```bash
python -m pip install -e ../eftkernels -e . matplotlib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
python examples/eft_3x2_study.py --output results/main
python examples/eft_3x2_study.py --output results/refined --eft-mu 128 --output-k 840
python examples/eft_3x2_study.py --output results/step_half --step-scale 0.5
python examples/eft_3x2_study.py --output results/quadrature --k-order 20 --mu-order 32
python examples/eft_3x2_study.py --output results/fft128 --nmax 128
python ../eftkernels/tests/validate_halo_integrals.py
python examples/eft_3x2_report.py --results results --output report
```

On PowerShell, replace the `export` line with
`$env:OMP_NUM_THREADS=$env:OPENBLAS_NUM_THREADS=$env:MKL_NUM_THREADS=$env:NUMBA_NUM_THREADS='1'`.
The main run uses radial/angular quadrature orders 14/20, 64 EFT angular
samples, 420 output wavenumbers and FFTLog N=256. Each variant changes only
the specified setting: interpolation resolution, derivative step, integration
quadrature, or FFTLog resolution, respectively. These are numerical comparisons,
not certification of convergence or validity at a given kmax.

The halo validator writes `results/halo_integrals.json` and
`results/halo_loop_sums.json` relative to the working directory. The latter is
the report input. Each study writes `arrays.npz`, `forecasts.json` and
`metadata.json`; regenerate older outputs because the report now also uses the
saved whitened derivatives. The report checks all prerequisites before writing
plots, LaTeX tables, `results/forecast_summary.csv` and
`results/diagnostics.json`. It can be launched from another directory with
explicit `--results` and `--output` paths; its study import is relative to the
script. Singular higher-EFT-fixed fits record rank and null directions, with
unidentifiable errors represented by JSON null and no regularization.

The adapter caches the fiducial model set and reuses a separate working set
for growth perturbations, avoiding all loop preparation on return to the
fiducial growth rate. It does not share eftkernels' private preparation state
between tracer models. AP coordinates and the volume factor use the same
Kaiser helpers as the built-in model. This changes the last-bit AP arithmetic
from the initial adapter, including exact isotropic boundary handling; it
does not change the alpha/phi convention or numerical tolerances.

Templates must supply `k`, `pk_full`, `pk_nw`, `growth_rate`, `rd_mpc_h` and
`z_fiducial`. A FITS input to `read_template` must contain `PKNW`; the adapter
does not silently replace it with `pk_smooth`. Paired nodes must have positive
k and mu on [0, 1]; negative mu is rejected, since the real-valued provider
does not represent odd-mu contributions.
