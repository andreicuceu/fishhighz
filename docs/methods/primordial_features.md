# Direct cosmology and primordial features

This sibling preserves the original compressed-AP forecast mode. New direct
cosmology calculations use physical AP factors DM/DMfid and Hfid/H, a single
volume Jacobian, and fixed reference h_fid units. They must not simultaneously
fit independent compressed AP/growth parameters. There is no additional
sound-horizon factor when the trial linear spectrum is regenerated.

`CambSnapshotProvider` now prepares immutable linear/smooth spectra, growth,
distances and sound horizon for all requested redshifts in one actual transfer
solve. The five shared coordinates are omega_b,omega_c,h,n_s,ln(10^10 A_s).
Neutrino, radiation, helium and dark-energy content are held at the explicit
base model. Scalar total-matter growth and EdS EFT kernels are approximations;
neutrino mass is not a supported varied coordinate. Cache keys include the
full source, solver, smoothing, base cosmology and requested grids.

CAMB can silently clamp a multi-redshift spectrum below its transfer minimum.
This provider rejects that input by default. `low_k_policy="matched_power_law"`
explicitly matches the value and slope just above the true transfer boundary;
it is a documented soft-integral continuation, not additional computed modes.
High-k extrapolation is never silently accepted.

`CosmologicalEFTAdapter` refreshes the Gaussian one-loop spectrum and damping
with each cosmology. Shared observed forests can use one physical coefficient
set. Known survey responses/noise remain external. Forecast-specific residual
stochastic/FoG terms require an explicit additive callback. One-loop primordial
features are deliberately rejected by this Gaussian adapter: the companion
campaign constructs their validated null response with feature-specific IR.
There is no relabeling of an old BAO forecast as a primordial-feature result.

Focused checks are `tests/test_physical_ap.py`, `tests/test_camb_snapshot.py`
and `tests/test_cosmological_eft.py`; real-transfer campaign checks are separate
from the deterministic solver fixtures in unit tests. Companion EFT branch:
codex/lya-galaxy-primordial-features, based on
f5f8427b49814dd0e910a7a5aa60ed53e74cff09.
