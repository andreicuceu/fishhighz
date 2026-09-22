# Step 02 revision 3 implementation handoff

Plan: Step 02, revision 3 (2026-09-11). This implementation is ready for user
review. It does not record acceptance or authorize Step 03.

## Outcome and scope

Implemented the default quick reference workflow and closed all four findings in
the [revision-1 review](step-02-review-r1.md). The original implementation report
is retained unchanged as [step-02-r1.md](step-02-r1.md).

The controller now defaults to `--suite quick`, requires a validated compatible
full baseline, resolves and runs only `lya_qso_lbg_lae_15x2pt.ini`, and embeds
the baseline manifest, primary result, compatibility evidence, and recomputed
comparison. `--suite full` remains an explicit opt-in with seven primaries plus
one independent repeat. New manifests use schema version 2 and record exact
suite coverage. The checker retains explicit schema-v1 full-bundle support and
implements `--require-suite` without launching a worker.

No FishHighz numerical code, runtime dependency, package API, sibling source,
scientific setting, or authoritative input was changed. The user-modified
`AGENTS.md` and `IMPLEMENTATION_STEP.md` were preserved. No commit, push,
publish, environment installation, Slurm action, or new real full capture was
performed.

Maintained implementation state:

| File | State | SHA-256 |
| --- | --- | --- |
| `scripts/lyaforecast_baseline.py` | added/revised | `f4a6c8c8ffd59f3a579698eb7f21bb96f7be2cbafc19199d250119fe9d37ee60` |
| `scripts/lyaforecast_baseline_worker.py` | added/unchanged from r1 | `63f3655e5785795aeb2a6fc0c97233d58c2d9be138de0e7ea85dad07abc698e0` |
| `tests/test_lyaforecast_baseline.py` | added/revised | `a26fcdf2d076c079460be138a1ef040f5b76bc305639ed58f620f50c6b02216f` |
| `tests/fixtures/stub_lyaforecast_worker.py` | added/revised | `1b79fe86ff1f467ccf1284e858cabcde765c5be0860b786c68e43588ff34c2f9` |
| `README.md` | modified | `50ef78d1f978c4ef9acb0f3a9f6f947983d449252f67e76f3312a3b48e5ac24e` |
| `reviews/step-02-r1.md` | retained r1 report | `3fc8e74ffbf21e24c0a89bb83d2666dfa387eb24ce1e4d2d196248c236749e7a` |

FishHighz remained at `130d144321c5cc7abea2f961b1f0f6d32b1ed1f5`
on `main`. The sibling reference remained clean at
`5abe8bcf8d12cc31d1f5ecd89c87e739c6a14b81` on `cleanup`.
`git diff --check` passed.

## Review findings addressed

| Finding | Resolution and regression coverage |
| --- | --- |
| Tool snapshots/shared evidence | The checker requires the schema-specific shared inventory, validates both nested tool snapshot paths and hashes, validates the nested Git-diff record, and rejects missing, corrupt, aliased, unexpected, or escaping paths. |
| Repeat worker evidence | One common validator now checks primary, repeat, and quick roles: manifest/status identity and success, response equality, round-trip state, thread state, request paths, required/distinct artifacts, configuration evidence, and results. Full primary/repeat directories and artifacts must be disjoint. |
| Coordinates and identities | Every k/mu coordinate, redshift center/edge, metadata array, and fiducial redshift is checked against the original INI with tight `5e-14` tolerances. Tracers, forest backgrounds, all pairs, and selected pairs are derived from the configuration. |
| Configuration/inventory preservation | Authoritative INIs are copied before resolution and used throughout. Schema-v2 captures compare pre/post INI hashes and exact authoritative inventory, bounded source inventory plus contents, and consumed directory inventories plus contents. |

The portable matrix refreshes semantic artifact hashes before adversarial checks
where appropriate. It covers suite action counts, default selection, compatibility
and relocation, numerical/structural mismatch, missing/invalid/quick baseline
rejection, source/config/input/environment mismatch before forecasting,
location-only changes, malformed coverage, nested artifacts, failed repeat state,
unverified round trips, artifact aliases, interior k/mu changes, redshift/bin
linkage, tracer/pair identities, and source/input/config content additions and
deletions during capture. Synthetic full orchestration runs eight stub `run`
actions; synthetic quick orchestration runs one `run` and one resolution action.

## Portable and offline acceptance checks

Final ordinary validation:

```bash
PATH="$PWD/.venv/bin:$PATH" ./scripts/check.sh
```

Result: exit 0 with Python 3.13.15; 52 tests passed in 22.34 seconds;
Ruff lint passed; Ruff format check reported 16 files already formatted.

The historical accepted full bundle was checked with the final controller:

```bash
.venv/bin/python scripts/lyaforecast_baseline.py check \
  .validation/baseline/20260911T184946Z-7a6a795d --require-suite full
```

Result: exit 0; schema v1 interpreted as full, all seven primaries and the
separate 15x2pt repeat passed. The checker reports, without inventing evidence,
that v1 records content preservation but lacks source, input, and configuration
inventory-preservation fields. Its manifest remains byte-identical at
`ee2186cea7a796f7aa95e8f655d9c04269e267db37a78b4879e81637ea57af1f`.

## Real quick reference evidence

Acceptance bundle:
`.validation/baseline/20260911T205607Z-81e76c36/` (20 MiB), schema 2,
suite `quick`, state `complete`; manifest SHA-256
`3ddb33c829203575941a0c7eb42eb4ce125d9e5d30fcad39ab6497545ca3efe0`.
Its coverage is exactly one quick-role primary, no repeat, and one scientific
run.

The bundle's capture-time tool snapshots are controller
`7da060d4441501914ab8df150ba58ca752cf2b6b00b7282a765050cbcb6ecd54`
and worker
`63f3655e5785795aeb2a6fc0c97233d58c2d9be138de0e7ea85dad07abc698e0`.
After capture, checker-only consistency checks were strengthened, producing the
maintained controller hash reported above; the final controller passed portable
tests and rechecked both the original and relocated quick bundle. Capture inputs,
execution, and result comparison were unchanged.

Command:

```bash
/usr/bin/time -f 'TOTAL_WALL_SECONDS=%e' \
  .venv/bin/python scripts/lyaforecast_baseline.py capture \
  --suite quick \
  --baseline .validation/baseline/20260911T184946Z-7a6a795d \
  --reference-checkout ../lyaforecast \
  --python ../lyaforecast/.validation/dev-env/bin/python
```

Result: exit 0 in 143.27 seconds total. The case wall time was 138.581 seconds;
worker initialization 132.121 seconds, forecast 5.602 seconds, and total worker
time 137.731 seconds. Serialization round trips passed.

The environment recorded Python 3.13.0, lyaforecast 0.1.0, NumPy 2.3.5,
SciPy 1.15.3, and CAMB 2.0.1. `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and
`MKL_NUM_THREADS` were each `1`; bytecode writes were disabled. Compatibility
passed for the exact recorded configuration, source, input, and environment
fields while excluding location-only paths. All 26 bounded source records,
six resolved inputs containing 28 files, and the one scoped authoritative INI
were unchanged with identical inventories and no additions or deletions.

The embedded historical manifest digest is
`ee2186cea7a796f7aa95e8f655d9c04269e267db37a78b4879e81637ea57af1f`.
The copied primary result digest is
`f2616fc928b34da3a97541ea2c7e133b3f5f50d03e928d46cae872146386046c`;
the fresh result has the same digest. The stored comparison artifact digest is
`b0e4e90300144ccabf3811020276087a72ffcf1c111cd6c3bcc7c56991d0873d`.
It compared 302 numbers at `rtol=1e-10`, `atol=1e-12`, with identical
keys/order/shapes/dtypes and zero maximum absolute and relative differences.

Relocation check:

```bash
cp -a .validation/baseline/20260911T205607Z-81e76c36 \
  /tmp/fishhighz-step02-r3-relocated-81e76c36
cd /tmp
/global/homes/a/acuceu/desi_acuceu/vega_dev/lib/fishhighz/.venv/bin/python \
  /global/homes/a/acuceu/desi_acuceu/vega_dev/lib/fishhighz/scripts/lyaforecast_baseline.py \
  check /tmp/fishhighz-step02-r3-relocated-81e76c36
```

Result: exit 0 using only the relocated quick bundle's embedded comparison
evidence. This later check does not claim to revalidate uncopied artifacts from
the external full bundle; that full bundle was validated at capture time.

## Retained diagnostic failure and limitations

The first quick attempt is retained at
`.validation/baseline/20260911T205206Z-cb9a056b/` with state `incomplete`.
Its one scientific worker completed, but the controller reported source inventory
failure after comparing identical membership in different path sort orders
(`lyaforecast/` versus `lyaforecast.egg-info/`). The evidence itself records
identical pre/post inventories, no additions/deletions, and unchanged hashes.
After changing the checker to exact order-independent membership and adding an
`.egg-info` regression fixture, that retained evidence passes every substantive
check under the internal incomplete-state diagnostic path. It was not overwritten
or promoted; the acceptance capture used a new destination. The failed command
took 139.00 seconds total and motivated the clean quick retry above.

No new real full schema-v2 suite was run, per the user-selected quick validation
policy. Full-mode orchestration and revised preservation are covered synthetically;
the real full scientific evidence remains the historical schema-v1 bundle. Quick
mode exercises only the full 15x2pt configuration and does not replace galaxy-only
or pair-subset full coverage. Native pickles remain trusted local evidence and are
not loaded by the offline checker. This step validates the lyaforecast reference
workflow; FishHighz scientific numerics have not begun and are not claimed valid.
