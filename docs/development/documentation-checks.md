# Documentation implementation checks

The documentation consolidation on 2026-09-22 preserves the adopted scientific
prescription, numerical code, examples, inputs, and existing tests. It does not
constitute a new scientific qualification of the forecasts.

## Migration and preservation

The [inventory](migration.md) accounts for 163 moved research/development
records and two derived README archive records. The README is now 116 lines;
its scientific tutorials remain in the methods pages and its historical
account remains in the archives.

Before/after SHA-256 comparisons found package source, bundled data, examples,
and tests byte-identical. Nine diagnostic scripts changed only live document
paths. Saved evidence, recorded hashes, historical commands, and snapshot
filenames were not regenerated. Operational agent instructions and bundled
provenance remain in place.

## Bounded checks

These checks passed without real-survey calculations:

```bash
python -m sphinx -n -W --keep-going -b html docs docs/_build/html
python -m pytest -q tests/test_survey_config.py tests/test_public_forecast.py tests/test_resources.py tests/test_desi2_accuracy_example.py
```

The focused selection passed **51 tests**. Both native and annotated INIs were
parsed without preparing a forecast. The synthetic covariance and Fisher Python
blocks were executed directly from the methods pages. All 113 explicitly
selected autodoc objects resolved. Ruff lint and formatting checks passed for
the Sphinx configuration and nine scripts with relocated document paths.

The integrated HTML build passed with warnings treated as errors. The generated
HTML audit checked **18,545 local links and anchors** across 202 HTML pages with
zero errors, including this check record. Representative HTML contained
mathematical elements, tables, highlighted INI examples, and API signatures.

## Distributions

Wheel and source-distribution builds passed. Documentation sources and included
examples ship in the source archive; research/development archives and generated
HTML are excluded from the runtime wheel. A strict Sphinx build from the
extracted source archive passed outside the checkout. Autodoc imported installed
FishHighz with neither Vega nor lyaforecast installed.

The temporary environment used Python 3.13, Sphinx 8.2.3, MyST-Parser 4.0.1, and
sphinx-rtd-theme 3.1.0. Logs, original-file snapshots, and the HTML audit are
local-only evidence under `/tmp/fishhighz-docs-20260922/`, not distributed assets.

## Independent review

Sol (high) independently reviewed the INI/Python/result documentation against
source, sampled the active research records against preserved originals,
reviewed API coverage and packaging, and inspected validation evidence.
**Pass: no outstanding scientifically or operationally consequential findings.**

The review requested clarification of the 163 moved files versus 165 manifest
entries and correction of duplicate/counting errors in archive navigation.
Both were corrected; all 162 archived destinations occur exactly once.
Scientific qualifications and the distinction between individual and joint
constraints were retained. The review did not run new forecasts or confer
scientific acceptance on historical results.

## Limits

Read the Docs and CI configurations are prepared; no hosted build or publication
was requested. No real-survey forecast, full scientific validation suite, Slurm
action, commit, or push was performed. Numerical results and their scientific
acceptance status remain those of the preserved research records.
