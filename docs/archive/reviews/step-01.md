# Step 01 handoff — revision 4

Implemented 2026-09-11 against `IMPLEMENTATION_STEP.md`, revision 4.
Status: ready for user review, not accepted. No subsequent step was started.

## Changes addressing review

**U1 (direct layout):** moved the unchanged minimal module from
`src/fishhighz/__init__.py` to `fishhighz/__init__.py`. Setuptools now discovers
only `fishhighz` and `fishhighz.*` from the repository root with namespace
discovery disabled. The egg-info ignore rule now applies at the root.
The editable installation was refreshed and its new origin checked with isolated
Python from `/tmp`. Version, runtime dependencies, smoke test, and runner behavior
are unchanged.

Before removing the obsolete root `src/`, its contents were inspected: the module,
generated egg-info, and a Python cache were the only files. The module was moved;
the remaining generated tree was preserved at
`.validation/revision-3-helper/obsolete-src-generated/`. No root `src/` remains.

**R1 (stale validation):** the helper allocates a new `.validation/r4-*` directory
for every invocation. Each contains a clean source staging tree copied from an
explicit list of maintained files, fresh artifact output directories, and two
new installation environments. Selection requires exactly one matching artifact.
Both environments are checked for absence of FishHighz before installation;
`--no-index --no-deps` installs the exact wheel. Each isolated probe verifies
quiet import, site-packages origin, matching metadata (including raw METADATA
bytes), and every installed package file against the selected wheel payload.
An unexpected installed package file also fails the probe. Bytecode caches are
excluded from that payload comparison.

The source archive is separately extracted under `/tmp` and rebuilt using declared
isolated build requirements. Its direct layout and executable runner are checked.
Only `fishhighz/` and its dist-info belong to either wheel. Helper/probe snapshots,
commands, logs, source snapshots, artifact paths, and hashes are retained per run.

Changed files: `pyproject.toml`, `.gitignore`, `README.md`,
`.validation/validate.py`, `.validation/probe.py`, and `reviews/step-01.md`.
Moved source: `src/fishhighz/__init__.py` to `fishhighz/__init__.py`.
Added archive: `reviews/step-01-r3.md` (exact copy of the previous handoff).
`MANIFEST.in`, `scripts/check.sh`, and `tests/test_import.py` remain unchanged.
The old helpers were archived under `.validation/revision-3-helper/` before
replacement. All original revision-3 final logs and distribution artifacts remain
at their original paths; all three old artifact hashes were rechecked successfully.
The reviewer's `reviews/step-01-review-r3.md` is unchanged.

## Setup and exact workflow

Used the existing isolated `.venv` created for revision 3 with:

```bash
/global/common/software/nersc/pe/conda-envs/26.8.0/python-3.13/nersc-python/bin/python -m venv .venv
```

For this revision, from `lib/fishhighz`, ran the following twice consecutively
with unchanged maintained source and helpers, with no intervening cleanup:

```bash
.venv/bin/python .validation/validate.py
.venv/bin/python .validation/validate.py
```

Each invocation reinstalls editable development support using
`<dev-python> -m pip install -e '.[dev]'`, then invokes the executable
`scripts/check.sh` from the repository root and an unrelated `/tmp` directory.
The runner executes `python -m pytest`, `python -m ruff check .`, and
`python -m ruff format --check .`. Every subprocess has PYTHONPATH removed,
`.venv/bin` first on PATH, and OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, and
MKL_NUM_THREADS set to 1. Checks remained lightweight on the login node.

The per-run helper and logs preserve full arguments and working directories for:

```text
<dev-python> -m build --outdir <run>/dist <run>/source
<dev-python> -m venv <run>/wheel-env
<wheel-python> -m pip install --no-index --no-deps <exact-wheel>
<wheel-python> -I <run>/probe.py <exact-wheel>
<dev-python> -m build --wheel --outdir <run>/rebuilt <temporary-extracted-source>
<dev-python> -m venv <run>/rebuilt-wheel-env
<rebuilt-python> -m pip install --no-index --no-deps <exact-rebuilt-wheel>
<rebuilt-python> -I <run>/probe.py <exact-rebuilt-wheel>
```

Each installation environment disables system site packages, starts with only pip,
and ends with only pip and FishHighz. Imports run from temporary directories outside
the checkout with isolated Python. Temporary extraction/import directories are
removed after use; logs and installed environments remain. Build staging avoids
all previous build/egg-info contents. Original `dist/` artifacts are revision-3
evidence; the new artifacts are exclusively in the run directories below.

Python 3.13.15 (conda-forge, GCC 14.4.0), pip 26.2.1, pytest 9.1.1,
Ruff 0.16.7, build 1.6.1, isolated setuptools 84.0.0. Other dev tools:
iniconfig 2.3.0, packaging 26.3, pluggy 1.6.0, Pygments 2.21.0,
pyproject_hooks 1.2.0. Each run records `versions.log` and `python-version.log`;
isolated backend versions appear in `build.log` and `sdist-build.log`.

## Acceptance evidence

Run A: `.validation/r4-5vcpk0qj/` (local-only path: `../.validation/r4-5vcpk0qj/`).
Run B: `.validation/r4-x26gxng2/` (local-only path: `../.validation/r4-x26gxng2/`).
Both completed all 18 logged commands with exit status 0. Each has `status.json`.
All four wheel installations report successful installation, with none skipped
as already installed. The two run roots and all four installation environment
paths are distinct; the helper asserts each environment path does not exist
before creation.

| Assessment | Result in both runs | Evidence in each run |
| --- | --- | --- |
| Layout/discovery | Direct package in both distributions; no src entries or unwanted wheel packages | `sdist-contents.log`, `wheel-contents.log`, `rebuilt-wheel-contents.log` |
| Editable development | Reinstallation succeeds; root and outside runner calls each pass one pytest test plus both Ruff checks | `editable-install.log`, `checks.log`, `checks-outside.log` |
| Fresh-process import | Quiet-import test passes; editable origin is the current `fishhighz/__init__.py` | check logs, `editable-origin.log` |
| Exact built wheel | Fresh install and isolated origin/metadata/payload checks pass | `wheel-{venv,empty,install,import,packages}.log` |
| Source archive rebuild | Isolated build succeeds from independent extraction; executable check runner included; rebuilt wheel passes same install checks | `sdist-build.log`, `sdist-contents.log`, `rebuilt-wheel-{venv,empty,install,import,packages}.log` |
| Repeat reliability | Two consecutive executions, unchanged source/helpers, separate outputs/environments, no installation skips | both `status.json` files, install logs and per-run manifests |
| Scope/evidence | Earlier report archived; reviewer/governance hashes unchanged; old artifact hashes still match | `reviews/step-01-r3.md`, `.validation/revision-3-helper/protected-sha256.txt`, original `.validation/artifact-sha256.txt` |

The probes report `fishhighz`, version `0.1.0.dev0`, Python `>=3.11`, and the
module under each new environment's `lib/python3.13/site-packages/fishhighz/`.
Package import itself is quiet; probes deliberately print verification evidence.

## SHA-256 manifests

The following manifests identify maintained source/configuration/tests and the
actual helper/probe snapshots used. Each was verified after both runs; both
snapshots match the current helpers. Reports are not hashed into themselves.

### r4-5vcpk0qj

```text
b68abc33791d9979d1266075d1c7adcaa1b59b266d0ce49f1518cc3243e2c3c6  .gitignore
8d3a7bfcd5bd654becff5cf1ae2ddfac34f030ca142b7fa0bd5f8c2585c96272  MANIFEST.in
f36acc5441f5f194d264a97d013d58dab56f7d4821ce532815d3caed43c355a6  README.md
1c3fe1f5196f5154fadae20e2005fc5c4fabd9b1a2719622e95c36b3de5e26c6  pyproject.toml
fc4a18f0a503279533f23556925e4bb6b031f7b259fc6c7c4cb14bfa9bf9cc3d  scripts/check.sh
32084cf9099038a201e89ed594484bfcf80ea76673c557e7a4949e7a5faa39c4  fishhighz/__init__.py
4929862120311eb39cedd3838ff7858decda84c597abd0f9d408d53cdd12b158  tests/test_import.py
b734542948a87e5032f81138409bbc603fda1813bb00e7ba8364101a9d90f017  .validation/r4-5vcpk0qj/validate.py
90de186a2d4d8d4e6bb27efc82061ee0936e345770797897bdf059641e5cdcc9  .validation/r4-5vcpk0qj/probe.py
```

Artifacts:

```text
30192312978c66e2a9672b6cfb23619f3fe8770935b839feb7fa68aa978fdaf7  .validation/r4-5vcpk0qj/dist/fishhighz-0.1.0.dev0-py3-none-any.whl
e2aa55c13e4dfb97370b6c5d9295ae22241b3c46c78b25383f19a542d4b4a065  .validation/r4-5vcpk0qj/dist/fishhighz-0.1.0.dev0.tar.gz
eac0fef4bad8a26c1590c676d001dfea282382d4fd7a448d735a177c1b06f8db  .validation/r4-5vcpk0qj/rebuilt/fishhighz-0.1.0.dev0-py3-none-any.whl
```

### r4-x26gxng2

```text
b68abc33791d9979d1266075d1c7adcaa1b59b266d0ce49f1518cc3243e2c3c6  .gitignore
8d3a7bfcd5bd654becff5cf1ae2ddfac34f030ca142b7fa0bd5f8c2585c96272  MANIFEST.in
f36acc5441f5f194d264a97d013d58dab56f7d4821ce532815d3caed43c355a6  README.md
1c3fe1f5196f5154fadae20e2005fc5c4fabd9b1a2719622e95c36b3de5e26c6  pyproject.toml
fc4a18f0a503279533f23556925e4bb6b031f7b259fc6c7c4cb14bfa9bf9cc3d  scripts/check.sh
32084cf9099038a201e89ed594484bfcf80ea76673c557e7a4949e7a5faa39c4  fishhighz/__init__.py
4929862120311eb39cedd3838ff7858decda84c597abd0f9d408d53cdd12b158  tests/test_import.py
b734542948a87e5032f81138409bbc603fda1813bb00e7ba8364101a9d90f017  .validation/r4-x26gxng2/validate.py
90de186a2d4d8d4e6bb27efc82061ee0936e345770797897bdf059641e5cdcc9  .validation/r4-x26gxng2/probe.py
```

Artifacts:

```text
f0180efd43bf58ffe329dc1196fe3df9c9823333e27d834ca7a0ef1030bc1a56  .validation/r4-x26gxng2/dist/fishhighz-0.1.0.dev0-py3-none-any.whl
50bafdc3d67bc297979ddb9f70e95af0cfd9decef9b0b6ba110b36e3c07dfc3f  .validation/r4-x26gxng2/dist/fishhighz-0.1.0.dev0.tar.gz
f424daf110cb947308d27ea650c42d8c68b22dffa9bec516d8d1f72370945f8f  .validation/r4-x26gxng2/rebuilt/fishhighz-0.1.0.dev0-py3-none-any.whl
```

## Limitations and remaining decisions

No required assessment is incomplete and no scientific work was introduced.
An initial sandboxed attempt failed during editable installation because PyPI DNS
resolution was unavailable (exit 1), retained in
`.validation/r4-ob5iwzyc/editable-install.log` and `status.json`. Both complete
runs used approved network access for declared build dependencies. No shared
software installation was modified.

Validation exercised Python 3.13.15 only; no multi-version matrix or scientific
baseline was required. Clean staging intentionally lists this scaffold's current
inputs; future package resources must be added to that list when introduced.
Separate builds need not have byte-identical wheel hashes.

The design, roadmap, AGENTS.md, current-step plan, and reviewer report were left
unchanged and their saved hashes rechecked. No siblings, scientific APIs, runtime
numerical dependencies, invented metadata, Git operations, publication, or Slurm
actions were introduced. Implementation is ready for the user's review; acceptance
and any progression remain user decisions.
