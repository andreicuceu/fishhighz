# Step 5 live validation: native DESI-2 reproduction

This is a bounded live validation of the standalone FishHighz implementation
after explicit authorization for the native-versus-reference CAMB comparison,
the selected six-bin DESI-2 run, and the 78+6 result comparison. It is
implementation evidence only; it is not scientific acceptance of the forecast
or adoption of any cosmological/weighting profile.

## Execution and provenance

All commands ran serially on login node `login07`, without Slurm, with
`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, and `MKL_NUM_THREADS=1`. The
interpreter was
`.validation/s1-s3-env/bin/python` (Python 3.13.0, NumPy 2.4.6, SciPy 1.15.3,
Astropy 8.0.1, CAMB 2.0.1, lyaforecast 0.1.0, FishHighz 0.1.0.dev0).
The complete environment and controls are in
`.validation/standalone-live-20260922T152624Z/environment.txt` and
`runtime-controls.txt`; the final numerical commands are in `commands.txt`,
and the complete CAMB attempt/output-copy chronology is in
`camb-attempts.json`.

The exclusive evidence directory is
`.validation/standalone-live-20260922T152624Z/`. The native and reference
Planck18 INIs were copied respectively from `fishhighz/data/camb_configs` and
`../lyaforecast/lyaforecast/resources/camb_configs`. Both have SHA-256
`45a04472fb946a2306b0c9668081922e6b0bddd11dfe28f7634aac34d6db9199`.

## Native/reference CAMB comparison

The native solve used one bulk `prepare_camb` call with the six geometric
`z_eval` values, template-growth redshift 2.406, and damping-reference
redshift 2.3. The reference solve reproduced the historical construction:
`CosmoCamb(Planck18.ini, z_ref=2.3,
z_centres=sorted({six z_eval, 2.406}))`. Reference growth and sigma8 were
matched by explicit redshift keys, not array position. The six geometric
redshifts were:

```text
2.1152848986890427, 2.3504402695765223, 2.5855752676523194,
2.8206936542989154, 3.0557983184571693, 3.2908915157575356
```

The final machine-readable result is `camb-comparison-corrected.json`:

| quantity | maximum relative difference |
| --- | ---: |
| H | 0.0 |
| D_M | 0.0 |
| sigma8 | 6.246761000721432e-11 |
| growth rate f | 3.334302704998877e-09 |

The native one-bulk solve took 61.31406985601643 s and the reference
template-plus-bins solve took 114.41983850701945 s. The `/usr/bin/time -p`
wall time for the combined comparison was 176.74 s. Raw stdout, stderr, and
timing are retained as `camb-final.*`.

Two earlier harness results are retained, but explicitly superseded: the first
reversed the damping/template reference roles, and the second used unsorted
reference bin inputs with positional mapping. Neither is used as evidence for
the final result.

## Public six-bin native forecast

The authorized public calculation executed exactly
`Forecast("desi2_accuracy.ini").run()` from the FishHighz checkout. It returned
78 selected individual constraints, 6 joint constraints, and 12 explicitly
excluded records. The run's internal elapsed time was 131.26374126703013 s;
`/usr/bin/time -p` reported 131.76 s wall time. The complete native result is
retained in `native-result/settings.json`, `native-result/results.npz`, and the
expanded finite-value record dump `native-raw.json`. Standard output, standard
error, and timing are `native-forecast.*`.

## Strict comparison with saved accuracy results

`compare_accuracy.py` inspected the saved schema before matching records. It
failed closed on missing or duplicate identities, non-finite values, shape
mismatches, status differences, parameter-ID differences, and unavailable
records. It matched the 78 individual and 6 joint records by
`(bin_index, kind, pair)`, then compared bounds, redshift, each 2x2 Fisher
matrix, `sigma_ap`, `sigma_at`, `correlation`, and the 12x12 combined Fisher
matrix. It also cross-checked the NPZ row identities (`bin_index`, kind, pair
indices, bounds, redshift, and availability) against the saved JSON records.
The comparison used exactly `rtol=1e-8`, `atol=0.0`.

Result: **passed**, with all 84 included records available and matched. The
maximum relative difference over all compared quantities was
`9.460748895689632e-11`; the maximum absolute difference was
`4.203866410534829e-07`. The complete per-record JSON, stdout, stderr, and
timing are `accuracy-comparison.json` and `accuracy-comparison.*`.

The independent audit found 505 comparison metrics. It verified that the
native NPZ Fisher rows equal the expanded raw Fisher arrays exactly and that
the saved settings contain 84 available records plus 12 excluded records.
The worst relative metric is the bin-0 `lya(qso)` auto individual Fisher
matrix (`9.460748895689632e-11`); the worst absolute metric is the bin-3 joint
Fisher matrix (`4.203866410534829e-07`).

The final-result hashes are retained in `sha256-final.txt`.

The immutable saved reference inputs were
`.validation/desi2-examples/accuracy/settings.json` (SHA-256
`d52dff6a6ba29b2406f824dc6378a1095d461a8dbd43ee41a000051251504be7`) and
`results.npz` (SHA-256
`760cdc1645f13446130c66d1d9ba4e36add469c8e985ac802a4946be4bad705c`). The
new native raw/result hashes are recorded in `sha256-final.txt` and can be
regenerated from the exact commands above.

## Scope and limitations

No production source, scientific defaults, tolerances, or saved historical
inputs were changed for this validation. No Slurm allocation was used. The
comparison establishes numerical reproduction of the retained accuracy
artifact under the stated environment; it does not establish scientific
acceptance, accuracy against an external survey requirement, or convergence
of the underlying forecast model. The live result is ready for independent
review only.
