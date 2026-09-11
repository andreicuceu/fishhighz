# FishHighz

FishHighz is an early development scaffold for high-redshift Fourier-space
Gaussian covariances and Fisher forecasts. No forecasting API or scientific
calculation is implemented yet. The planned starting backend is NumPy with
dedicated Numba-friendly kernels; this scaffold has no runtime dependencies.

## Development

Use Python 3.11 or newer. From this directory, create an isolated environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m pip install -e '.[dev]'
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
Workspace contributors can consult the [design](../../FISHHIGHZ_DESIGN.md) and
[roadmap](../../FISHHIGHZ_IMPLEMENTATION_PLAN.md); neither is needed to install,
import, or test the package.
