# S5.1 and S5.2 documentation handoff

The package documentation now identifies the tested S2--S4 accuracy recipe as
the recommended research baseline while preserving the two compatibility
profiles as reproduction and weighting-isolation controls. No numerical recipe,
profile identity, public default, historical evidence, or forecast result changed
in this documentation subtask. The coordinated S5.2 implementation separately
updates public docstrings and three input-error messages without changing the
scientific calculation.

## Documentation changes

- `RESEARCH_BASELINE.md` is a standalone scientific guide. It specifies signal,
  response, physical FWHM conversion, mixed-pair damping, full wiggle AP
  derivatives, early-lyaforecast weights and units, adaptive failure behavior,
  forest noise, Gaussian covariance, selection, result metadata and the exact
  three-profile identities.
- The guide states that the validated calculation fits exactly `ap` and `at` per
  bin with no unspecified nuisance marginalization or priors. It distinguishes
  this recipe from the general parameter/prior facilities of the public API.
- `README.md` now reflects implemented survey preparation, default P1D and
  optional compiled Fisher contraction. It points to the selected baseline and
  to `examples/research_bao_forecast.py:run()`. The old Step-12 fixed-reference
  account is explicitly marked historical rather than silently rewritten.
- Strict nonnegative public inputs are separated from signed/floored
  compatibility adapters. Adaptive iteration convergence is separated from the
  S2 finite-grid qualification. Profile names and revision identity are retained.

## Scientific evidence and limits

The guide reuses the reviewed S2 maxima of 0.0011411% individual and 0.0012685%
joint BAO-error refinement change, the S3 accuracy/fixed joint changes, and the
S4 endpoint closure and attribution. It adds no calculation. The selected
baseline adopts the tested assumptions; it does not establish that the
mixed-pair reconstruction physics is physically preferred. D1/D2 and future
bin-1 selection remain identified without adding a new prescription.

## Validation

Focused checks passed for local Markdown links, referenced public names,
required baseline/profile content and whitespace. The coordinated compact
example runs successfully from the package environment and its importable
`run()` returns a converged early-weight, rank-two joint AP result. Ruff check
and format check pass for the example. `git diff --check` passes for the S5
documentation, example and manifest changes. The coordinator is performing the
installed wheel/sdist validation. No scientific forecast, Slurm action, commit,
or push is part of this handoff.

```bash
.venv/bin/python examples/research_bao_forecast.py
.venv/bin/python -c '<import example; assert rank, method, status and AP shape>'
.venv/bin/python -c '<check local Markdown targets>'
.venv/bin/python -m ruff check examples/research_bao_forecast.py
.venv/bin/python -m ruff format --check examples/research_bao_forecast.py
git diff --check -- README.md RESEARCH_BASELINE.md \
  reviews/s5-documentation.md MANIFEST.in examples/research_bao_forecast.py
```

All commands passed. The direct script reported three selected spectra, a
rank-two joint result, converged early-lyaforecast weights after 10 updates, and
the expected fixed three-spectrum covariance preparation.

Ready for independent scientific, software, and documentation review.
