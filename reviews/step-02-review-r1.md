# Step 02 revision 1 review

Outcome: changes requested. The captured forecasts and supplied report check out,
but the evidence checker has gaps. The user's request for a faster development
validation path is incorporated in Step 02 revision 2. Step 03 is not authorized.

## Verified evidence

Reviewed both maintained scripts, synthetic tests/worker, README, current plan,
and `reviews/step-02.md`. All five full implementation-file hashes and the
acceptance manifest hash in that report match the current files. FishHighz HEAD
is `130d144`; implementation files remain uncommitted. The reference checkout
remains clean.

The acceptance bundle is `.validation/baseline/20260911T184946Z-7a6a795d/`.
The reviewer ran:

```bash
PATH="$PWD/.venv/bin:$PATH" ./scripts/check.sh
.venv/bin/python scripts/lyaforecast_baseline.py check \
  .validation/baseline/20260911T184946Z-7a6a795d
```

Results: 13 tests passed in 10.63 seconds, Ruff lint passed, Ruff format passed
(14 files), and the checker accepted all seven primaries plus the separate repeat.
Independent inspection of all eight stored results confirmed full k/mu arrays
and redshift centers/edges agree with their configurations and metadata. Both
stored tool snapshots match their recorded hashes. The stored repeat compares
302 numbers with zero maximum absolute/relative difference.

No real forecast was rerun for this review. Existing full capture evidence,
independent artifact checks, portable tests, and focused failure probes were
sufficient to assess this revision. Implementation source was not changed.

## Findings requiring revision

1. **P2 — tool snapshot integrity is not checked.**
   `scripts/lyaforecast_baseline.py:1060` verifies shared metadata, source, and
   input snapshots, but never verifies the tool snapshot entries in `tool.json`.
   Appending a comment to a copied controller snapshot, without changing any
   stored hash, still produces a successful bundle check. A missing snapshot is
   similarly not covered. Check every required tool snapshot and its recorded
   digest, using the existing bundle-relative path containment rules.

2. **P2 — repeat validation skips worker failure evidence.**
   `scripts/lyaforecast_baseline.py:1187` does not apply the primary-case status,
   exit-code, and round-trip checks to the repeat. A copied bundle with repeat
   exit code 17 and a stored `failed` status is accepted after updating that
   status artifact's checksum. Reuse a common case validator for primary,
   repeat, and future quick records; reject conflicting status/response records
   and require successful exits and verified round trips in every role.

3. **P2 — grid/redshift checks do not verify the recorded calculation fully.**
   `scripts/lyaforecast_baseline.py:865` checks redshift ordering/enclosure but
   does not compare centers/edges with metadata and configuration. At line 916,
   only grid lengths/endpoints are checked. On copied evidence the checker
   accepts an interior k value of -1000, and separately accepts a result redshift
   shifted by 0.001 while its metadata/configuration remain unchanged. Validate
   every coordinate and expected shape, plus redshift and tracer/pair identity,
   against the actual reference conventions. Add direct failure coverage.

4. **P2 — preservation checks miss configuration and inventory changes.**
   `scripts/lyaforecast_baseline.py:318` rehashes only previously listed source
   files; authoritative example INIs are not part of the before/after records.
   `scripts/lyaforecast_baseline.py:415` likewise rehashes existing input files
   without checking additions to the original directories. Synthetic probes
   adding a Python module or SNR file, or changing an example INI, still record
   unchanged evidence. Snapshot/hash authoritative INIs before preparation, run
   from those captured configurations, and compare bounded before/after source,
   INI, and consumed-directory inventories as well as file contents.

The mutation probes are under `.validation/step02-review-7qx1vtwy/`, including
`probe.py` and `results.json`. Preservation probes and their results are under
`.validation/step02-preservation-probe-cjhb_wgv/`. All mutations used reviewer
copies or synthetic directories. These findings concern missing validation;
they do not show that the original stored forecasts are incorrect.

## Validation-speed decision

The eight recorded worker wall times sum to 1029.70 seconds (17.16 minutes).
The full primary 15x2pt case took 134.06 seconds, with 127.91 seconds in
initialization and 5.26 seconds in the forecast itself. Ordinary tests already
avoid real forecasts and take roughly ten seconds.

| Option | Benefit | Cost/limitation |
| --- | --- | --- |
| Full-resolution 15x2pt once, compared with saved full baseline | About two minutes based on this capture; retains all fields and scientific settings | Does not cover galaxy-only and pair-subset workflows |
| Two concurrent single-thread workers | Idealized full-suite time around nine minutes if initialization scales independently | Same total work; more memory/CPU use and concurrent failure handling; speedup unmeasured |
| Reduced k/mu/redshift grids | Less forecast work | Initialization dominates here; changes the scientific reference and loses coverage |
| Share/cache initialized cosmology across cases | Could remove repeated initialization | Changes the independent reference execution and introduces cache correctness work |

Recommendation: keep fast ordinary tests, add explicit quick 15x2pt validation,
and reserve full captures for relevant review checkpoints. Keep execution serial
in this revision. Do not add cosmology caches or alter grids. A full capture is
needed once on the final revised capture implementation; fixes during development
use unit tests and, when relevant, the quick mode. Future FishHighz steps should
reuse the immutable reference bundle and test FishHighz outputs directly; merely
rerunning unchanged lyaforecast does not validate new FishHighz code.
