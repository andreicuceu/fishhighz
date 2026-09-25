# Fixed-growth DESI-2 BAO package implementation

The opt-in `bao_marginalized` native mode forecasts bin-local
`(alpha_iso, phi)` wiggle dilation with identity smooth dilation. Its mapping is
`ap=alpha_iso*phi**(-2/3)` and `at=alpha_iso*phi**(1/3)`; the pre-existing
`alpha_phi` convention remains distinct. The galaxy Kaiser rate is a fixed
fiducial `f(z)` shared by QSO, LBG and LAE. Active nuisance parameters are one
bias for each selected galaxy tracer, one Lyα bias and one Lyα beta shared by
both forest fields. The full group has two targets and five nuisances; bin 1
of the standard selection has two targets and three nuisances. No prior is
applied. Each bin retains independent parameter identities.

The mode reuses the fixed fiducial survey preparation, interval Fisher
accumulation, covariance closure and derivative engine. Template support is
checked over the actual smooth and wiggle scaling bases and derivative stencils.
`Forecast.run(individuals=False)` and CLI `--joint-only` omit individual Fisher
calculations. Empty selected bins retain an excluded joint record. Saving uses
`fishhighz-bao-marginalized-result` version 1 and the named matrix layout of
the existing full-shape output. `target_ids` are ordered
`alpha_iso_N,phi_N`, with zero-based original bin index `N`. The resolved
settings include the numerical fixed galaxy growth rate by bin.

Validation on 2026-09-24: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1 .venv/bin/python -m pytest -q
tests/test_bao_marginalized.py tests/test_full_shape.py
tests/test_public_forecast.py tests/test_survey_config.py` passed 74 tests.
Ruff check and format-check passed for all eight modified or new Python files.
The synthetic checks cover dilation mapping, active registries, empty bins,
dense marginalization, singular status, joint-only serialization and fixed
growth binding. They do not establish numerical convergence or physical
adequacy at the requested `k_max=0.5 h/Mpc`; those require the planned
interactive-node calculations and independent review.
