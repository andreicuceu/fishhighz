"""Explicit recipe identity and scientifically selected forecast observables."""

from ..fields import PairSelection
from .cases import selection

REVISION = "early-lyaforecast-2026-09-18"
PROFILES = ("full-compatibility", "fixed-compatibility", "accuracy")
ADAPTIVE = ("early_lyaforecast", "mcdonald")
STOPPING = dict(rtol=1e-4, min_updates=3, stable_steps=3, max_updates=96)
REFERENCE = dict(
    k_t_deg=2.4,
    k_p_velocity=0.00035,
    convention="fiducial_auto_p3d_and_p1d_times_response_squared",
)


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
