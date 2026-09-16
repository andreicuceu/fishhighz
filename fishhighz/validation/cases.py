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
    """Return an owned copy of one known original recipe, rejecting unknown IDs."""
    if case not in RECIPES:
        raise ValueError(f"unknown DESI-2 case {case!r}")
    return json.loads(json.dumps(RECIPES[case]))


def selection(case):
    """Preserve original field indices and reference upper-triangle pair order."""
    r = recipe(case)
    fields = []
    for key, t in r.items():
        if not key.startswith("tracer "):
            continue
        forest = t["tracer_type"] == "continuous"
        name = f"lya({t['background_tracer']})" if forest else t["tracer"]
        fields.append(
            ObservedField(
                name,
                "forest" if forest else "galaxy",
                t["tracer"],
                background=t.get("background_tracer"),
            )
        )
    wanted = r["control"]["correlations"].split()
    pairs = [
        (a.id, b.id)
        for i, a in enumerate(fields)
        for b in fields[i:]
        if "all" in wanted or f"{a.id}_{b.id}" in wanted or f"{b.id}_{a.id}" in wanted
    ]
    return PairSelection(fields, pairs)


def bins(case):
    """Original equal edges from explicit limits/count, retaining centre labels."""
    s = recipe(case)["survey"]
    e = np.linspace(
        float(s["z bin min"]), float(s["z bin max"]), int(s["num z bins"]) + 1
    )
    return [(float(lo), float(hi)) for lo, hi in zip(e[:-1], e[1:])]


def verify_inventory(directory):
    """Read-only exact seven-INI verification, including unknown effective options."""
    paths = sorted(Path(directory).glob("*.ini"))
    if [p.stem for p in paths] != list(CASE_IDS):
        raise ValueError("original INI inventory differs from seven explicit recipes")
    result = {}
    for p in paths:
        c = configparser.ConfigParser()
        c.read(p)
        actual = {s: dict(c[s]) for s in c.sections()}
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
