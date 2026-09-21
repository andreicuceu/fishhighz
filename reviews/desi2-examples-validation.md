# DESI Run-2 example validation

## Scientific result

All three examples completed the selected six-bin `lya_qso_lbg_lae_15x2pt`
forecast once.  Each produced 78 available individual-spectrum constraints and
six available joint constraints, with the 12 excluded bin-1 spectra recorded
separately.  No selected result was unavailable or non-finite.  Direct inversion
of every saved 2x2 Fisher matrix reproduced its saved uncertainties and
correlation.

The joint results are `(sigma_ap, sigma_at, correlation)`:

| Bin | Full compatibility | Fixed compatibility | Accuracy |
|---:|---|---|---|
| 1 | (0.02472610, 0.01804447, -0.432481) | (0.02514347, 0.01835109, -0.432006) | (0.02417269, 0.01810703, -0.414632) |
| 2 | (0.02046864, 0.01462255, -0.430652) | (0.02006916, 0.01429040, -0.431184) | (0.01929890, 0.01409680, -0.414161) |
| 3 | (0.01799679, 0.01237200, -0.426455) | (0.01735381, 0.01188115, -0.427555) | (0.01671429, 0.01172291, -0.410890) |
| 4 | (0.01501345, 0.00966468, -0.420230) | (0.01468227, 0.00943943, -0.420975) | (0.01424716, 0.00935641, -0.404950) |
| 5 | (0.01635541, 0.01079681, -0.419483) | (0.01605049, 0.01061337, -0.420211) | (0.01562206, 0.01054263, -0.404254) |
| 6 | (0.02313372, 0.01741007, -0.426903) | (0.02158659, 0.01642430, -0.429338) | (0.02076141, 0.01616004, -0.412447) |

## Reproduction against S2/S3 evidence

The independent checker compared `settings.json` and `results.npz` with
`.validation/s3/three-profile-tables.json`, including the raw Fisher matrices.
Maximum relative differences in the reported `(sigma_ap, sigma_at,
correlation)` were `9.23e-16`, `8.70e-16`, and `4.43e-9` for full, fixed, and
accuracy respectively.  Maximum Fisher-element relative differences were
`7.14e-16`, `7.39e-16`, and `8.82e-9`.

The first accuracy comparison used `rtol=1e-9` and failed on a `2.30e-9`
relative transverse-error difference.  The final checker uses `rtol=1e-8` for
both saved Fisher and summary comparisons.  This is far below the S2 0.1%
finite-refinement target, but the tolerance change is recorded rather than
treated as an exact reproduction.  A bounded metadata comparison found exact
volume and representative `B`; maximum relative differences were `2.98e-9` in
growth rate and sigma8, `5.97e-9` in template power growth and representative
`P`, `2.98e-9` in damping widths, `1.04e-9` in `A`, and `2.21e-9` in
`P_pixel`.  The small numerical difference is not assigned to a single cause.

Accuracy also reproduced the S2 magnitude-node counts
`(2592, 2592, 2592, 2624, 2688, 2640)` and converged update counts: QSO forest
`(22, 20, 18, 16, 14, 12)` and LBG forest `(not selected, 36, 36, 30, 28, 26)`.

Machine-readable comparisons are in
`.validation/desi2-examples/comparison.json`; per-profile outputs are under
`.validation/desi2-examples/{full-compatibility,fixed-compatibility,accuracy}/`.

## Execution and environment

The Sol (medium) validation agent, rather than the coordinator, executed the
examples sequentially on the login node with `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, and `MKL_NUM_THREADS=1`.  The environment was Python
3.13.0, NumPy 2.4.6, SciPy 1.15.3, Astropy 8.0.1, CAMB 2.0.1, and the default
NumPy FishHighz Fisher backend.  Source revisions were FishHighz
`0cfcf2a0e7ae25c424a05604b64835f2c1bcde40`, lyaforecast
`5abe8bcf8d12cc31d1f5ecd89c87e739c6a14b81`, and Vega
`01f1a9784767e43e68c89ee4464c2bdcb15420c0` with pre-existing dirty FishHighz
work preserved.

Commands, run exactly from `lib/fishhighz`, were:

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=.:../lyaforecast /usr/bin/time -p .validation/s1-s3-env/bin/python examples/desi2_full_compatibility.py --reference ../lyaforecast --output .validation/desi2-examples/full-compatibility
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=.:../lyaforecast /usr/bin/time -p .validation/s1-s3-env/bin/python examples/desi2_fixed_compatibility.py --reference ../lyaforecast --output .validation/desi2-examples/fixed-compatibility
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=.:../lyaforecast /usr/bin/time -v .validation/s1-s3-env/bin/python examples/desi2_accuracy.py --reference ../lyaforecast --template ../vega/vega/models/Planck18/Planck18_z_2.406.fits --output .validation/desi2-examples/accuracy
```

Wall times were 124.27 s, 123.96 s, and 177.72 s.  Accuracy used 494,228 kB
peak resident memory and no swap.

The shared result helper acquired optional per-record `parameter_ids` after the
full-compatibility run completed.  Therefore that retained output has an
unambiguous global `parameter_order=[ap, at]` but no repeated per-record IDs;
fixed and accuracy each have IDs on all 84 records.  This metadata-only change
did not alter numeric arrays, and the final serialization tests cover it.  The
completed full calculation was not rerun.

## Focused checks

The final combined command was:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=. .validation/s1-s3-env/bin/python -m pytest -q tests/test_compatibility_bao.py tests/test_compatibility_weights.py tests/test_full_sum_weights.py tests/test_desi2_compatibility_examples.py tests/test_desi2_accuracy_example.py tests/test_research_bao_example.py
```

Result: **41 passed in 1.22 s**.  `ruff check` passed for the three examples,
shared helper, compatibility adapter, and two new test modules; `ruff format
--check` reported all seven files already formatted.  Both short help commands
and the accuracy help command passed from the live source checkout.

This validates numerical reproduction and the implemented API interactions.  It
does not add scientific acceptance beyond the retained S2--S4 evidence.
