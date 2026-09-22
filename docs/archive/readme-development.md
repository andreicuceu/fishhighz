# Earlier development and reference-capture instructions

This is the historical README record. Commands and local paths describe the
original development environment, not required current contributor procedures.

## Development

Use Python 3.11 or newer. From this directory, create an isolated environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m pip install -e '.[dev,templates]'
./scripts/check.sh
```

Source lives directly in `fishhighz/`. Package discovery includes only `fishhighz`
and its package descendants, with implicit namespaces disabled. Reinstall the
editable package after changing its layout. Imports from the checkout root alone
do not verify installation; the smoke test imports in an isolated subprocess
outside the checkout.

The runner uses the active environment and works when invoked by absolute or
relative path from another directory. It stops on a failed check and does not
install tools or change formatting. Its individual commands are:

```bash
python -m pytest
python -m ruff check .
python -m ruff format --check .
```

Ruff targets Python 3.11 with an 88-character line length. Enabled rules are E4
(imports), E7 (statement errors), E9 (runtime/syntax errors), F (Pyflakes), and I
(import sorting). Pytest, Ruff, and the build frontend are declared in the `dev`
extra. Setuptools >=68 is the isolated build requirement. The version is defined
only in `pyproject.toml`; installed versions can be read with
`importlib.metadata.version("fishhighz")`.

Ruff excludes generated `.validation/` evidence even before Git is initialized.

## Distributions

```bash
python -m build
```

This creates a wheel and source distribution in `dist/`. The retained local
acceptance helper builds from a clean source snapshot into a new output directory
and creates two fresh installation environments on every invocation:

```bash
python .validation/validate.py
```

The helper is local review evidence, not a distributed package tool. It records
exact wheel paths and hashes, rejects ambiguous artifact selection, and compares
installed package bytes with each selected wheel. It repeats those checks for a
wheel rebuilt from a separately extracted source archive. Run directories and
helper/probe snapshots remain under `.validation/r4-*`; no cleanup is needed
between runs. Isolated builds may download the declared build requirements.

The [Step 01 report](reviews/step-01.md) records the exact isolated interpreter,
commands, artifact hashes, and independent installation checks used for verification.
Workspace contributors can consult the [design](planning/FISHHIGHZ_DESIGN.md) and
[roadmap](planning/FISHHIGHZ_IMPLEMENTATION_PLAN.md); neither is needed to install,
import, or test the package.

## Lyaforecast reference capture

Step 02's maintained tool defaults to a quick scientific check: one unchanged,
full-resolution 15x2pt forecast compared with a previously accepted full bundle.
Supply the baseline, reference checkout, and scientific Python explicitly. The
destination is created exclusively and must not already exist:

```bash
python scripts/lyaforecast_baseline.py capture \
  --suite quick \
  --baseline .validation/baseline/<existing-full-bundle> \
  --reference-checkout ../lyaforecast \
  --python ../lyaforecast/.validation/dev-env/bin/python
```

`--suite quick` is the default, but is shown for clarity. It requires a compatible
full baseline and never falls back to a full run. A complete capture remains
available only as explicit opt-in; it runs all seven authoritative DESI-2 cases
and an independent 15x2pt repeat:

```bash
python scripts/lyaforecast_baseline.py capture \
  --suite full \
  --reference-checkout ../lyaforecast \
  --python ../lyaforecast/.validation/dev-env/bin/python
```

The default destination is a unique, ignored
`.validation/baseline/<UTC-timestamp>-<suffix>/` directory. Each case runs
serially in a fresh isolated subprocess with single-thread settings. Bundles
contain original/effective INIs, configuration/input/source snapshots and
inventories, Git and environment provenance, complete pickle and typed-JSON
results, logs, timings, and checksums. Quick bundles copy the full reference
result, compatibility fields, and numerical comparison so offline checks survive
relocation. Existing destinations are never resumed or overwritten.

Validate a bundle without rerunning forecasts, NumPy, CAMB, or the original
checkout:

```bash
python scripts/lyaforecast_baseline.py check \
  .validation/baseline/<UTC-timestamp>-<suffix>
python scripts/lyaforecast_baseline.py check \
  .validation/baseline/<UTC-timestamp>-<suffix> --require-suite full
```

The checker supports the accepted schema-v1 full evidence and identifies
inventory guarantees unavailable in that legacy schema. Real reference forecasts
are explicit validation work and are not run by `scripts/check.sh`; ordinary
pytest coverage uses synthetic inputs and stub workers only. Run the real full
suite only when it is specifically requested.

