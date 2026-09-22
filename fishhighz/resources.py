"""Access to data bundled with FishHighz.

The package data are accessed through :mod:`importlib.resources` so that the
same API works from a source checkout, an installed wheel, and an installed
zip-style distribution.  Callers that need a filesystem path must keep the
``bundled_path`` context manager open while using the returned path.
"""

from __future__ import annotations

from contextlib import ExitStack, contextmanager
from importlib import resources
from pathlib import Path, PurePosixPath
from typing import Iterator

_DATA_ROOT = resources.files("fishhighz").joinpath("data")


def _relative_name(name: str | Path) -> str:
    """Validate and normalize a package-relative resource name."""
    if not isinstance(name, (str, Path)):
        raise TypeError("resource name must be str or pathlib.Path")
    text = str(name).replace("\\", "/")
    path = PurePosixPath(text)
    if (
        not text
        or path.is_absolute()
        or any(part in ("", ".", "..") for part in path.parts)
    ):
        raise ValueError(f"resource name must be a relative data path: {name!r}")
    return path.as_posix()


def bundled_resource(name: str | Path):
    """Return a bundled data file as an ``importlib.resources`` traversable.

    Parameters
    ----------
    name : str or pathlib.Path
        Path relative to the package ``data/`` directory, for example
        ``"camb_configs/Planck18.ini"``.

    Raises
    ------
    FileNotFoundError
        If the requested package data file is not present.
    """
    resource = _DATA_ROOT.joinpath(_relative_name(name))
    if not resource.is_file():
        raise FileNotFoundError(f"bundled FishHighz resource not found: {name!r}")
    return resource


@contextmanager
def bundled_path(name: str | Path) -> Iterator[Path]:
    """Yield a temporary filesystem path for one bundled data file.

    ``importlib.resources.as_file`` may materialize a temporary file when the
    package is imported from a wheel or another non-filesystem distribution;
    therefore the yielded path is valid only inside this context manager.
    """
    with resources.as_file(bundled_resource(name)) as path:
        yield path


@contextmanager
def bundled_paths(directory: str | Path) -> Iterator[tuple[Path, ...]]:
    """Materialize each file in a bundled directory for Python 3.11 safety.

    ``importlib.resources.as_file`` did not support directory traversables on
    all supported Python 3.11 releases.  Enumerating the traversable and
    materializing files individually works for source trees, wheels and zip
    imports; yielded paths remain valid until this context exits.
    """

    relative = _relative_name(directory)
    resource = _DATA_ROOT.joinpath(relative)
    if not resource.is_dir():
        raise FileNotFoundError(
            f"bundled FishHighz resource directory not found: {directory!r}"
        )
    names = sorted(child.name for child in resource.iterdir() if child.is_file())
    if not names:
        raise FileNotFoundError(
            f"bundled FishHighz resource directory is empty: {directory!r}"
        )
    with ExitStack() as stack:
        paths = tuple(
            stack.enter_context(bundled_path(f"{relative}/{name}")) for name in names
        )
        yield paths


def resolve_input_path(path: str | Path, survey_ini_path: str | Path) -> Path:
    """Resolve a user input path relative to its survey INI, not the CWD.

    Absolute user paths are retained (apart from normal path resolution).
    Relative paths are anchored at the directory containing
    ``survey_ini_path``.  This helper deliberately performs no existence
    check, allowing the caller to issue a context-specific diagnostic.
    """
    candidate = Path(path).expanduser()
    ini_path = Path(survey_ini_path).expanduser().resolve(strict=False)
    if not candidate.is_absolute():
        candidate = ini_path.parent / candidate
    return candidate.resolve(strict=False)


__all__ = [
    "bundled_path",
    "bundled_paths",
    "bundled_resource",
    "resolve_input_path",
]
