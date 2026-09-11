# Step 01 handoff — revision 3

Implemented 2026-09-11 against `IMPLEMENTATION_STEP.md`, revision 3.
Status: ready for user review; acceptance and progression remain with the user.

## Delivered behavior and files

FishHighz is an installable src-layout scaffold with a quiet import, no runtime
requirements, and no scientific API. The version has one definition in
`pyproject.toml`. Setuptools builds both distribution formats. The executable
check runner uses the active Python, changes to the package root, sets all three
thread limits to one, and stops on any failed command without modifying files.

Added files: `.gitignore`, `MANIFEST.in`, `README.md`, `pyproject.toml`,
`scripts/check.sh` (executable), `src/fishhighz/__init__.py`,
`tests/test_import.py`, and this report, `reviews/step-01.md`.
`MANIFEST.in` includes the documented runner in source distributions.
Generated environments, logs, validation helpers, and hashes are in `.validation/`
and `.venv/`; distribution/build metadata are in `dist/`, `build/`, and
`src/fishhighz.egg-info/`. These are ignored, not implementation source.

No design, roadmap, AGENTS.md, current-step plan, sibling package, scientific
reference, or parent ignore file was edited. No Git repository was initialized,
and nothing was committed, published, or submitted to Slurm. Metadata contains
no invented author, license, or repository information. Ruff uses E4/E7/E9/F/I,
Python 3.11 syntax targets, 88 columns, and excludes generated validation evidence.

## Environment and reproducible commands

Lightweight checks ran on the login node. The shell's default Python was 3.6.15.
`module load python` located the interpreter below but emitted shell-initialization
`ERROR:: command not found` messages; subsequent direct invocation worked.
No shared environment was modified. The workspace home path resolves to the CFS
path shown in the logs; all writes stayed within the authorized package directory.

From `lib/fishhighz`:

```bash
/global/common/software/nersc/pe/conda-envs/26.8.0/python-3.13/nersc-python/bin/python -m venv .venv
.venv/bin/python -m pip install -e '.[dev]' > .validation/install.log 2>&1
# Retry with approved network access after sandbox DNS failure:
.venv/bin/python -m pip install -e '.[dev]' > .validation/install-retry.log 2>&1
.venv/bin/python .validation/validate.py > .validation/validation.log 2>&1
```

The retained [validation driver](../.validation/validate.py) records every exact
command and working directory in its per-command logs. It clears `PYTHONPATH`,
prepends `.venv/bin` to PATH, and sets `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`,
and `MKL_NUM_THREADS` to 1. It runs `scripts/check.sh` from the package root and
from a temporary directory under `/tmp`, then `python -m build` with build isolation.
For each wheel it runs `python -m venv <environment>`, installs the exact absolute
wheel path with `python -m pip install --no-deps <wheel>`, and runs the retained
[import probe](../.validation/probe.py) using `python -I -c` from a fresh `/tmp`
directory. No system site packages or development extras enter those environments.

For the source distribution, the driver extracts the archive into a separate
`TemporaryDirectory` under `/tmp` using the tarfile data filter, then runs:

```text
<dev-python> -m build --wheel --outdir <package>/.validation/rebuilt <extracted-source>
```

This uses declared isolated build requirements and no adjacent checkout files.
Temporary extraction/import directories are removed after assessment; the command
logs, archive member lists, rebuilt wheel, and installed environments remain.
After adding the runner to the source manifest, the full driver was repeated.
The first successful install environments were moved to
`.validation/initial-wheel-env` and `.validation/initial-rebuilt-wheel-env` before
creating fresh final environments. Those archived environments are not final evidence.

Versions: Python 3.13.15 (conda-forge, GCC 14.4.0), pip 26.2.1, pytest 9.1.1,
Ruff 0.16.7, build 1.6.1, isolated setuptools 84.0.0. Other dev dependencies:
iniconfig 2.3.0, packaging 26.3, pluggy 1.6.0, Pygments 2.21.0, pyproject_hooks 1.2.0.
See [versions](../.validation/versions.log),
[Python details](../.validation/python-version.log), and build logs for provenance.

## Acceptance evidence

These are implementation results for review, not user acceptance decisions.

| Assessment | Observed result | Evidence under `.validation/` |
| --- | --- | --- |
| Isolated development install | Retry exited 0; editable FishHighz and declared dev tools installed | `install-retry.log`, `versions.log` |
| Source workflow | Exit 0: 1 pytest test passed; Ruff lint and format checks passed | `checks.log` |
| Runner from unrelated cwd | Exit 0 with the same three successful checks | `checks-outside.log` |
| Fresh-process import | Test passes from pytest temporary directory, with isolated Python and cleared PYTHONPATH; empty stdout/stderr; neither neighboring package in sys.modules | `checks.log`, `tests/test_import.py` |
| Distribution build | Exit 0; wheel and source archive built | `build.log` |
| Built wheel independence | All venv/install/import/list commands exit 0; only FishHighz and pip installed | `wheel-{venv,install,import,packages}.log` |
| Source archive completeness | Separate extraction and isolated rebuild exit 0; source, metadata, test, README, and executable runner included | `sdist-contents.log`, `sdist-build.log` |
| Rebuilt wheel independence | All venv/install/import/list commands exit 0; only FishHighz and pip installed | `rebuilt-wheel-{venv,install,import,packages}.log` |
| Contents and metadata | Both wheels contain `fishhighz/__init__.py` and four dist-info files; no runtime Requires-Dist; only dev-conditional requirements | `wheel-contents.log`, `rebuilt-wheel-contents.log` |

Both installed probes report name `fishhighz`, version `0.1.0.dev0`, Python
requirement `>=3.11`, and a module inside the respective environment's
`lib/python3.13/site-packages/fishhighz/__init__.py`, rather than `src/`.
Import itself is quiet; the probe deliberately prints metadata/path evidence.
[Final driver results](../.validation/validation.log) record all checks exiting 0.

## Hash manifests (SHA-256)

Implemented files, excluding this report to avoid hashing it into itself:

```text
e985bf814e8fcba1230853f13b53d6e1c1e07b49dbb40d32d54716df9955b057  .gitignore
8d3a7bfcd5bd654becff5cf1ae2ddfac34f030ca142b7fa0bd5f8c2585c96272  MANIFEST.in
92459c1fed5542d739b86ff0bcb50e2d922a8a42fa743044fdb8c33550dde1be  README.md
8db9df8215d5938f331c4212604da03bbe51f9a35e65cb2d299d06639c55fd28  pyproject.toml
fc4a18f0a503279533f23556925e4bb6b031f7b259fc6c7c4cb14bfa9bf9cc3d  scripts/check.sh
32084cf9099038a201e89ed594484bfcf80ea76673c557e7a4949e7a5faa39c4  src/fishhighz/__init__.py
4929862120311eb39cedd3838ff7858decda84c597abd0f9d408d53cdd12b158  tests/test_import.py
```

Built artifacts:

```text
24b344bfb5de6a978d11458abe450883149db45c1597df63ef2396f332def957  dist/fishhighz-0.1.0.dev0-py3-none-any.whl
8f425e1bd461720c70e445b74953f679e8804c9f01bb9de4eed34ed2dc626084  dist/fishhighz-0.1.0.dev0.tar.gz
8943c55c6cc5a777e3045d52bee023848689de8855f5ee269057ec2db5160d9d  .validation/rebuilt/fishhighz-0.1.0.dev0-py3-none-any.whl
```

Manifests also exist at `.validation/source-sha256.txt` and
`.validation/artifact-sha256.txt`. Wheel hashes identify artifacts, not an assertion
that separate builds are byte-for-byte reproducible.

## Limitations and review decisions

No scope departures or incomplete acceptance checks remain. Initial sandboxed
installation and build attempts exited 1 because PyPI DNS resolution failed;
errors are retained in `install.log` and `build-sandbox-failure.log`. Approved
network retries succeeded. Only Python 3.13.15 was exercised; the declared 3.11
floor is reflected in Ruff configuration but a multi-version test matrix was
not part of this step. No scientific baseline or scientific validation was run.

Implementation is ready for the user's review. Step 01 is not marked accepted,
and no subsequent step has been drafted or implemented.
