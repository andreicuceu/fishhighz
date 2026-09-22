# Independent review: Step 5 live validation

## Review status

**Pass for numerical reproduction.** The three explicitly authorized login-node
calculations completed successfully. This verdict is implementation evidence; it
does not constitute scientific acceptance of the forecast or adoption of a
cosmological or weighting profile.

All calculations ran serially on `login07` with one thread per numerical
library and without Slurm. No production source, scientific default, tolerance,
or immutable reference input was changed by the validation.

## Evidence checked

The final native/reference CAMB comparison uses the historical reference
construction with `z_ref=2.3`, sorted redshift centres, and growth quantities
matched by explicit redshift. The native and reference Planck18 configuration
files are byte-identical. The maximum relative differences are:

| quantity | maximum relative difference |
| --- | ---: |
| H | 0.0 |
| D_M | 0.0 |
| sigma8 | 6.246761000721432e-11 |
| growth rate f | 3.334302704998877e-09 |

The combined CAMB comparison required 176.74 s wall time. Two earlier harness
outputs are preserved and marked superseded: one reversed the historical
damping/template roles and one associated unsorted reference output by array
position. Neither contributes to the accepted result. Their chronology and
hashes are recorded in `camb-attempts.json`.

The public call `Forecast("desi2_accuracy.ini").run()` required 131.76 s wall
time and produced 78 individual results, 6 joint results, and 12 excluded
records. I independently verified that the individual, joint, and combined
Fisher arrays in `native-result/results.npz` equal the expanded arrays in
`native-raw.json` exactly.

The strict comparison used `rtol=1e-8` and `atol=0.0`. All 84 included records
were available and matched by `(bin_index, kind, pair)`. Across 505 comparison
metrics, the maximum relative difference was `9.460748895689632e-11`, and the
maximum absolute difference was `4.203866410534829e-07`. The worst relative
metric was the bin-0 `lya(qso)` auto individual Fisher matrix; the worst
absolute metric was the bin-3 joint Fisher matrix.

## Conclusion

The live calculation reproduces the retained DESI-2 accuracy artifact at a
maximum relative difference approximately 106 times smaller than the requested
`1e-8` tolerance. The evidence is internally consistent and ready for the
researcher's review. No claim beyond this bounded numerical reproduction is
made.
