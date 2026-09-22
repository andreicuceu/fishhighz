"""Explicit recipe identity and scientifically selected forecast observables."""

from ..accuracy import ADAPTIVE, REFERENCE, REVISION, STOPPING
from ..fields import PairSelection
from .cases import selection

__all__ = [
    "ADAPTIVE",
    "REFERENCE",
    "REVISION",
    "STOPPING",
    "forecast_selection",
    "identity",
]

PROFILES = ("full-compatibility", "fixed-compatibility", "accuracy")


def forecast_selection(case, index):
    """Exclude weak bin-1 tracers before any noise or Fisher preparation."""
    original = selection(case)
    if index != 0:
        return original
    excluded = {"lbg", "lae", "lya(lbg)"}
    pairs = [
        p
        for p in original.selected_pairs
        if all(original.fields[i].id.lower() not in excluded for i in p)
    ]
    return PairSelection(original.fields, pairs)


def identity(profile, method=None):
    """Identity used for new records and cache comparisons."""
    profile = "full-compatibility" if profile == "compatibility" else profile
    if profile not in PROFILES:
        raise ValueError("unknown forecast profile")
    return dict(
        recipe_revision=REVISION,
        profile=profile,
        method=method
        or ("legacy" if profile == "full-compatibility" else "early_lyaforecast"),
        reference=REFERENCE,
        stopping=STOPPING,
        selection="bin1-qso-only; bins2-6-all",
    )
