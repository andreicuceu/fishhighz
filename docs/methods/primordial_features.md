# Direct cosmology and primordial features

This sibling preserves the original compressed-AP forecast mode. New direct
cosmology calculations use physical AP factors DM/DMfid and Hfid/H, a single
volume Jacobian, and fixed reference h_fid units. They must not simultaneously
fit independent compressed AP/growth parameters. There is no additional
sound-horizon factor when the trial linear spectrum is regenerated.

Initial support is a tested physical-coordinate helper. The transfer-aware
cosmology provider and feature-aware one-loop forecast are pending their
scientific gates. No output from the existing BAO-only mode is relabeled as
a primordial-feature forecast. Companion EFT branch:
codex/lya-galaxy-primordial-features, based on
f5f8427b49814dd0e910a7a5aa60ed53e74cff09.
