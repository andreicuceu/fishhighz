# Step 02 implementation handoff

Plan: Step 02, revision 1 (2026-09-11). This handoff is ready for user review;
it does not record user acceptance or authorize Step 03.

## Outcome and scope

Implemented a maintained standard-library capture/check controller and a separate
scientific worker. The controller requires explicit lyaforecast checkout and
Python paths, rejects existing destinations, resolves inputs through the
reference utilities, snapshots source/input evidence, runs every case
sequentially in a fresh isolated process, and validates typed JSON without NumPy
or the original checkout. The worker alone imports lyaforecast, NumPy, SciPy, and
CAMB and calls `NewForecast(effective_config_path).new_run_forecast()` unchanged.

No FishHighz science, runtime dependency, package API, sibling source, planning
document, or authoritative input was changed. `AGENTS.md` and
`IMPLEMENTATION_STEP.md` were already modified by the user and were preserved.

Maintained implementation changes at handoff:

| File | State | SHA-256 |
| --- | --- | --- |
| `scripts/lyaforecast_baseline.py` | added | `d85f5a7e412425556bdaa780c494493d3c3c085e7f60455f0ee8055d7ec6437f` |
| `scripts/lyaforecast_baseline_worker.py` | added | `63f3655e5785795aeb2a6fc0c97233d58c2d9be138de0e7ea85dad07abc698e0` |
| `tests/test_lyaforecast_baseline.py` | added | `dfb8605765601fc7d7ac815d7448572edf288541ec7112bf92f8d718deef61be` |
| `tests/fixtures/stub_lyaforecast_worker.py` | added | `9c5bf6a66d082786c7f3f866fe95d866809c5de11746082eefd04fb50ba23f65` |
| `README.md` | modified | `548f7ed6d6827ddb929b1ebf7ac1fbdcabf594d2f6d4c7d1bb1d9b69557d7707` |

Accepted scaffold files remained byte-identical: `scripts/check.sh`
`fc4a18f...9cc3d`, `pyproject.toml` `1c3fe1f...26c6`,
`fishhighz/__init__.py` `32084cf9...39c4`, and `tests/test_import.py`
`49298621...158`.

FishHighz was at `130d144` on `main`. Final `git status --short` contained the
two pre-existing modified governance files, modified `README.md`, and the four
new script/test paths above plus this report. `git diff --check` passed. No
commit, push, publish, environment installation, or Slurm action was performed.

## Reference and environment

The acceptance bundle is
`.validation/baseline/20260911T184946Z-7a6a795d/` (21 MiB), manifest SHA-256
`ee2186cea7a796f7aa95e8f655d9c04269e267db37a78b4879e81637ea57af1f`.
Its schema is `fishhighz.lyaforecast-baseline` version 1 and state is `complete`.

Reference checkout: `../lyaforecast`, clean branch `cleanup`, revision
`5abe8bcf8d12cc31d1f5ecd89c87e739c6a14b81`. The worker imported
`lyaforecast.forecast_new` from that checkout. The explicit interpreter was
`../lyaforecast/.validation/dev-env/bin/python`; the recorded executable retains
that virtual-environment path. It reported Python 3.13.0, lyaforecast 0.1.0,
NumPy 2.3.5, SciPy 1.15.3, and CAMB 2.0.1. All workers recorded the three thread
variables as `1` and disabled bytecode writes.

Twenty-six Python/package-metadata files were snapshotted. Six unique resolved
inputs comprising 28 files were snapshotted: Planck18.ini, three density tables,
and every one of the 12 entries in each DESI-2-LBG and DESI-2-QSO SNR directory.
All source/input before, snapshot, and after hashes matched.

## Forecast evidence

Each status below has exit code 0. Initialization and forecast timings come from
inside the worker; wall timing includes process and serialization overhead.

| Primary INI | Selected/all pairs | Init (s) | Forecast (s) | Wall (s) |
| --- | ---: | ---: | ---: | ---: |
| `lbg_lae_3x2pt.ini` | 3/3 | 124.901 | 0.279 | 125.861 |
| `lya_lbg_lae_3x2pt.ini` | 3/6 | 125.265 | 1.321 | 127.526 |
| `lya_lbg_lae_6x2pt.ini` | 6/6 | 121.061 | 1.398 | 123.520 |
| `lya_qso_2x2pt.ini` | 2/3 | 127.777 | 0.997 | 129.590 |
| `lya_qso_lbg_lae_4x2pt.ini` | 4/15 | 128.683 | 4.207 | 133.830 |
| `lya_qso_lbg_lae_8x2pt.ini` | 8/15 | 121.698 | 3.796 | 126.554 |
| `lya_qso_lbg_lae_15x2pt.ini` | 15/15 | 127.910 | 5.262 | 134.056 |
| separate 15x2pt repeat | 15/15 | 123.620 | 4.137 | 128.760 |

The checker verified exactly seven primaries, all artifact hashes, path-only INI
changes, complete input provenance, original grid/settings agreement, ordered
redshifts and enclosing edges, array shapes, finite values, positive selected and
combined uncertainties, bounded correlations, and zero-filled unselected pair
arrays. Pickle and typed-JSON round trips preserved complete dictionary order,
keys, shapes, and dtypes.

The independent repeat compared 302 stored numbers at `rtol=1e-10`,
`atol=1e-12`. Keys/order/shapes/dtypes were identical; maximum absolute and
relative differences were both exactly 0.0. The comparison records zero versus
zero as zero relative difference and nonzero versus a zero primary denominator
as infinity.

## Commands and acceptance checks

Capture, from the FishHighz root with its existing development environment:

```bash
.venv/bin/python scripts/lyaforecast_baseline.py capture \
  --reference-checkout ../lyaforecast \
  --python ../lyaforecast/.validation/dev-env/bin/python
```

Result: exit 0; seven primary cases and the repeat captured and internally
validated.

Direct checker, invoked from `/tmp`:

```bash
/global/homes/a/acuceu/desi_acuceu/vega_dev/lib/fishhighz/.venv/bin/python \
  /global/homes/a/acuceu/desi_acuceu/vega_dev/lib/fishhighz/scripts/lyaforecast_baseline.py \
  check /global/homes/a/acuceu/desi_acuceu/vega_dev/lib/fishhighz/.validation/baseline/20260911T184946Z-7a6a795d
```

Result: exit 0, `PASS: 7 primary cases and separate ...15x2pt.ini repeat`.

Relocation check, without the original reference source:

```bash
cp -a .validation/baseline/20260911T184946Z-7a6a795d \
  /tmp/fishhighz-step02-relocated-7a6a795d
cd /tmp
/global/homes/a/acuceu/desi_acuceu/vega_dev/lib/fishhighz/.venv/bin/python \
  /global/homes/a/acuceu/desi_acuceu/vega_dev/lib/fishhighz/scripts/lyaforecast_baseline.py \
  check /tmp/fishhighz-step02-relocated-7a6a795d
```

Result: exit 0 with the same seven-case/repeat PASS result.

Final development check:

```bash
PATH="$PWD/.venv/bin:$PATH" ./scripts/check.sh
```

Result: Python 3.13.15; 13 pytest tests passed in 8.44 seconds; Ruff lint passed;
Ruff format verification passed (13 files already formatted). The 12 new tests
use only synthetic fixtures and stub subprocesses and cover valid/relocated
bundles, selected and zero-filled pairs, missing/duplicate cases, worker failure,
missing/corrupt artifacts, malformed/non-finite results, unauthorized science
changes, destination refusal, consecutive fresh captures, and repeat numerical
or structural failure.

## Failure history, limitations, and unresolved findings

An audit of the first ignored capture
`.validation/baseline/20260911T182955Z-afed8d5e/` found that the controller had
resolved the requested virtual-environment Python symlink to its base Conda
executable. That bundle accurately recorded NumPy 2.4.6, but it is not the
acceptance artifact. The controller now preserves the lexical interpreter path,
a regression assertion checks worker argv, and the complete acceptance capture
was rerun with NumPy 2.3.5. The earlier bundle remains as diagnostic evidence and
was neither reused nor overwritten.

Native pickles are provided for exact local reference use and should be loaded
only as trusted generated evidence; the checker deliberately reads typed JSON.
This step validates reproducible legacy lyaforecast BAO outputs only. FishHighz
scientific numerical implementation has not begun and is not claimed validated.
No unresolved Step 02 acceptance finding remains for implementation review.
