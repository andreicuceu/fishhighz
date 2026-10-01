"""Seven explicit reference inventories, without a general INI translator."""

import configparser
import hashlib
import json
from pathlib import Path

import numpy as np

from ..fields import ObservedField, PairSelection
from ._recipes import RECIPES

CASE_IDS = tuple(RECIPES)


def recipe(case):
    """Return an owned copy of one known original recipe, rejecting unknown IDs.

    Parameters
    ----------
    case : str
        Identifier of one of the seven original DESI-2 validation
        configurations.

    Returns
    -------
    recipe : dict
        Independent nested configuration mapping for the original case.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    if case not in RECIPES:
        raise ValueError(f"unknown DESI-2 case {case!r}")
    return json.loads(json.dumps(RECIPES[case]))


def selection(case):
    """Preserve original field indices and reference upper-triangle pair order.

    Parameters
    ----------
    case : str
        Identifier of one of the seven original DESI-2 validation
        configurations.

    Returns
    -------
    selection : PairSelection
        Original fields and selected spectra, preserving upstream ordering.
    """
    configuration = recipe(case)
    fields = []
    for key, tracer_settings in configuration.items():
        if not key.startswith("tracer "):
            continue
        forest = tracer_settings["tracer_type"] == "continuous"
        name = (
            f"lya({tracer_settings['background_tracer']})"
            if forest
            else tracer_settings["tracer"]
        )
        fields.append(
            ObservedField(
                name,
                "forest" if forest else "galaxy",
                tracer_settings["tracer"],
                background=tracer_settings.get("background_tracer"),
            )
        )
    wanted = configuration["control"]["correlations"].split()
    pairs = [
        (a.id, b.id)
        for i, a in enumerate(fields)
        for b in fields[i:]
        if "all" in wanted or f"{a.id}_{b.id}" in wanted or f"{b.id}_{a.id}" in wanted
    ]
    return PairSelection(fields, pairs)


def bins(case):
    """Original equal edges from explicit limits/count, retaining centre labels.

    Parameters
    ----------
    case : str
        Identifier of one of the seven original DESI-2 validation
        configurations.

    Returns
    -------
    bounds : list of tuple of float
        Lower and upper dimensionless redshift edges in bin order.
    """
    survey_settings = recipe(case)["survey"]
    redshift_edges = np.linspace(
        float(survey_settings["z bin min"]),
        float(survey_settings["z bin max"]),
        int(survey_settings["num z bins"]) + 1,
    )
    return [
        (float(lo), float(hi))
        for lo, hi in zip(redshift_edges[:-1], redshift_edges[1:])
    ]


def verify_inventory(directory):
    """Read-only exact seven-INI verification, including unknown effective options.

    Parameters
    ----------
    directory : str or pathlib.Path
        Directory containing exactly the seven reference INI files.

    Returns
    -------
    inventory : dict
        Paths, SHA-256 hashes, configurations, selected pairs and bin bounds
        keyed by case.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    paths = sorted(Path(directory).glob("*.ini"))
    if [p.stem for p in paths] != list(CASE_IDS):
        raise ValueError("original INI inventory differs from seven explicit recipes")
    result = {}
    for p in paths:
        configuration = configparser.ConfigParser()
        configuration.read(p)
        actual = {s: dict(configuration[s]) for s in configuration.sections()}
        if actual != recipe(p.stem):
            raise ValueError(f"{p.stem}: unknown or changed recipe options")
        result[p.stem] = dict(
            path=str(p.resolve()),
            sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
            original=actual,
            selected_pairs=selection(p.stem).selected_pairs.tolist(),
            bins=bins(p.stem),
        )
    return result
