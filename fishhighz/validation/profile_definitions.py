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
    """Exclude weak bin-1 tracers before any noise or Fisher preparation.

    Parameters
    ----------
    case : str
        Identifier of one of the seven original DESI-2 validation
        configurations.
    index : int
        Zero-based redshift-bin index.

    Returns
    -------
    selection : PairSelection
        Retained observables with original field indices and required covariance
        pairs.
    """
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
    """Identity used for new records and cache comparisons.

    Parameters
    ----------
    profile : str
        Scientific profile name; compatibility is normalized to full-
        compatibility.
    method : str
        Forest-weight prescription: legacy, inverse_variance, early_lyaforecast
        or mcdonald, as applicable. Default is ``None``.

    Returns
    -------
    identity : dict
        Revision, profile, weight method, reference coordinates, stopping rules
        and redshift selection.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
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
