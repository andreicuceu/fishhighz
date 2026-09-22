# Step 2 native cosmology and bias implementation handoff

Status: implementation complete and ready for independent/user review. This
handoff does not record scientific acceptance or authorize Step 3.

## Scientific result and scope

FishHighz now has a NumPy-only importable cosmology boundary with lazy CAMB
preparation. `prepare_camb` uses the bundled byte-preserved
`camb_configs/Planck18.ini` by default and prepares one CAMB result surface for
the complete exact requested redshift set. The resulting immutable arrays retain
the caller's order
and an explicit redshift-to-index map for

```text
H(z) [km/s/Mpc], D_M(z) [Mpc], sigma8(z), f(z)=fsigma8(z)/sigma8(z).
```

The template-growth normalization and damping-reference normalization are
separate named redshifts and sigma8 quantities. Geometry callables use the
retained CAMB result surface at arbitrary quadrature nodes; growth lookup is
exact-redshift only and fails closed for unprepared values. CAMB is an optional
`fishhighz[camb]` dependency with an actionable missing-extra diagnostic.

The preparation requests the complete unique set in strictly decreasing
redshift order in one CAMB solve. It validates the returned
`Params.Transfer.PK_redshifts` surface against that exact request, then maps
`sigma8` and `fsigma8` by the validated redshift keys back to caller order.
There is no unlabelled array reversal. The example routes volume preparation
through the native `hubble_parameter` and transverse-comoving-distance methods,
so the public D_M contract is used rather than a radial-distance attribute.
Supplying both template-redshift keyword aliases is rejected.

The native bias utilities reproduce the existing lyaforecast analytic
constants/powers, fixed Ly-alpha beta, non-Ly-alpha `f/b` beta, and registered
linear tabulated bias with endpoint-slope extrapolation. The accuracy example
now uses these utilities and the native cosmology preparation; it no longer
constructs `lyaforecast.CosmoCamb` or `lyaforecast.PowerSpectrum`. Existing
Kaiser, signed template, damping, AP derivative, covariance, weight and Fisher
recipes are unchanged.

`fishhighz/validation/accuracy.py` and its historical evidence machinery still
retain their `CosmoCamb`/`PowerSpectrum` imports. They are validation-only
reproduction paths and were deliberately left unchanged; replacing them would
broaden this bounded native-example step and alter historical evidence scope.

## Changed files

- `fishhighz/cosmology.py`: lazy CAMB preparation and immutable background.
- `fishhighz/models/biases.py`, `fishhighz/models/__init__.py`: analytic and
  tabulated bias utilities.
- `fishhighz/__init__.py`, `pyproject.toml`: NumPy-only public exports and the
  optional `camb` extra.
- `examples/desi2_accuracy.py`: native cosmology/bias integration with all
  survey preparation and S2--S4 controls retained.
- `tests/test_cosmology.py`, `tests/test_biases.py`: deterministic controls for
  ordering, exact-redshift failure, optional imports, formulae and extrapolation.

The reviewed Step-1 resource/package changes and their files under `fishhighz/data/`
were preserved without byte changes.

## Checks

Focused command:

```bash
PATH="$PWD/.venv/bin:$PATH" pytest -q \
  tests/test_import.py tests/test_forecast_imports.py tests/test_cosmology.py \
  tests/test_biases.py tests/test_desi2_accuracy_example.py
```

Result: 15 passed. Ruff check and format check passed for every changed Python
source/test file; `git diff --check` passed. The fake-CAMB surface checks the
single bulk request, caller-order remapping, returned-redshift metadata
mismatch, H/D_M extraction, f extraction, immutability, exact-redshift failure,
alias rejection and missing-extra behavior.

The available CAMB 2.0.1 interpreter was used for two bounded one-thread
comparison attempts with the byte-identical upstream and bundled Planck18.ini
(SHA-256 `45a04472fb946a2306b0c9668081922e6b0bddd11dfe28f7634aac34d6db9199`).
The original per-redshift implementation was stopped as materially unsuitable
for login-node use. After the bulk implementation, the native solve completed
in 59.4 s for the six geometric evaluation redshifts plus template `z=2.406`
and damping `z=2.3`; the subsequent legacy comparison was stopped under the
same login-node runtime/resource policy before completion. No live numerical
equivalence result is therefore certified. The one-bulk-solve design is retained
as the bounded performance correction; a completed live CAMB comparison remains
an explicit review limitation.

The bounded existing suite reached 1463 passed and 25 skipped. Three failures
were outside this step: two historical synthetic independent-accuracy fixture
failures in `tests/test_step12_performance.py` and one optional Matplotlib
failure in `tests/test_three_profile_plots.py`; none involved the changed native
cosmology, bias, or example code.

No real forecast, accepted scientific equivalence claim, Slurm action, commit,
push, or roadmap advancement was performed. Scientific equivalence is limited
to the tested analytic/tabulated formulae and the fake-CAMB bulk-order contract;
the live CAMB/version comparison and forecast-level reproduction remain for
review.
