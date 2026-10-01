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
    """Validate and normalize a package-relative resource name.

    Parameters
    ----------
    name : str or pathlib.Path
        File or directory name relative to package data.

    Returns
    -------
    relative : str
        Normalized POSIX resource name.

    Raises
    ------
    TypeError
        If name is neither str nor pathlib.Path.
    ValueError
        If name is empty, absolute or contains traversal components.
    """
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
    """Return a bundled data file as an importlib resource.

    Parameters
    ----------
    name : str or pathlib.Path
        Path relative to package data, such as camb_configs/Planck18.ini.

    Returns
    -------
    resource : importlib.resources.abc.Traversable
        Bundled file supporting resource access.

    Raises
    ------
    FileNotFoundError
        If the package data file is absent.
    """
    resource = _DATA_ROOT.joinpath(_relative_name(name))
    if not resource.is_file():
        raise FileNotFoundError(f"bundled FishHighz resource not found: {name!r}")
    return resource


@contextmanager
def bundled_path(name: str | Path) -> Iterator[Path]:
    """Yield a temporary filesystem path for one bundled data file.

    Parameters
    ----------
    name : str or pathlib.Path
        File name relative to package data.

    Notes
    -----
    ``importlib.resources.as_file`` may materialize a temporary file when the
    package is imported from a wheel or another non-filesystem distribution;
    therefore the yielded path is valid only inside this context manager.

    Yields
    ------
    path : pathlib.Path
        Materialized file path, valid only within the context.
    """
    with resources.as_file(bundled_resource(name)) as path:
        yield path


@contextmanager
def bundled_paths(directory: str | Path) -> Iterator[tuple[Path, ...]]:
    """Materialize each file in a bundled directory for Python 3.11 safety.

    Parameters
    ----------
    directory : str or pathlib.Path
        Directory name relative to package data.

    Raises
    ------
    FileNotFoundError
        If the directory is absent or contains no files.

    Notes
    -----
    ``importlib.resources.as_file`` did not support directory traversables on
    all supported Python 3.11 releases.  Enumerating the traversable and
    materializing files individually works for source trees, wheels and zip
    imports; yielded paths remain valid until this context exits.

    Yields
    ------
    paths : tuple of pathlib.Path
        Materialized files in sorted filename order, valid within the context.
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

    Parameters
    ----------
    path : str or pathlib.Path
        Input path, possibly relative to the INI.
    survey_ini_path : str or pathlib.Path
        Survey INI path defining the relative-path base.

    Returns
    -------
    resolved : pathlib.Path
        Absolute expanded input path; existence is not checked.

    Notes
    -----
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
