# Results and saved output

`Forecast.run()` returns a `SurveyResult`. Keep the distinction between
individual-spectrum and joint constraints when interpreting BAO precision:

| Attribute | Content |
| --- | --- |
| `individual` / `individual_spectra` | One constraint per selected spectrum, calculated with its own covariance. |
| `joint` / `joint_constraints` | One constraint per bin from the full selected-spectrum covariance. |
| `excluded` | Explicit records for pairs omitted by the bin selection. |
| `constraints` | Individual records followed by joint records; excludes omitted pairs. |
| `fisher_results` | Fisher objects in `constraints` order. |
| `combined` | Combined joint Fisher result across independent bins, in the full named parameter basis. |
| `convergence` | Recorded forest-weight iteration diagnostics by bin and field. |
| `resolved_settings` | Expanded configuration, input identities, parameter order and execution settings. |
| `config`, `prepared` | Parsed configuration and retained prepared state. |

Each `SpectrumConstraint` gives `bin_index` (zero-based), `bin_id`, redshift
`bounds`, `z_eval`, `kind`, `pair`, `parameter_ids`, `status`, `fisher` (also
`result`), `sigma_ap`, `sigma_at` and `correlation`. Joint records have
`pair=None`. INI bin sections are one-based; parameter IDs and result indices
are zero-based (`ap_0`, `at_0` for the first bin).

`status="available"` means the two-parameter Fisher information has full rank;
`available` is the corresponding Boolean property. Rank-deficient constraints
have `status="unavailable"`, retain their Fisher object, and report errors and
correlation as `None`. Excluded records have `status="excluded"` and no Fisher
object or uncertainties. Do not interpret either status as a zero uncertainty.

The two BAO errors are marginalized within each record's AP pair. Individual
constraints are not blocks of the inverse joint covariance. Correlations between
spectra remain in the joint calculation even though the native sampling-noise
model is diagonal. Across bins, `combined` retains separate AP parameters; it
does not impose a shared dilation parameter.

```python
for item in result.joint:
    print(item.bin_index, item.status, item.sigma_ap, item.sigma_at)

output = result.save("accuracy-desi2")
```

`save()` creates a new directory, including missing parents, and returns its
resolved path. Existing destinations raise an error. Two files are written:

- `settings.json`: result schema `fishhighz-forecast-result`, version 1;
  expanded settings and input identities/checksums captured during construction
  and preparation; convergence records; constraint identities, statuses,
  uncertainties and counts. Non-finite JSON numbers are rejected.
- `results.npz`: `individual_fisher`, `joint_fisher` and `combined_fisher` arrays
  of data Fisher information. Individual and joint arrays follow their record
  order; excluded records have no matrix. Pair matrices use their recorded
  `parameter_ids`; the combined matrix uses saved `parameter_order`.

```python
import json
import numpy as np

with open("accuracy-desi2/settings.json") as stream:
    metadata = json.load(stream)
with np.load("accuracy-desi2/results.npz", allow_pickle=False) as arrays:
    combined = arrays["combined_fisher"]
```

This is an inspectable saved result, not a serialization of the full live
`SurveyResult`: no public inverse `load()` reconstructs prepared objects. Input
hashes identify the inputs used; the output directory does not contain copies
of every source table. Retain external inputs separately for reproducibility.

For direct Python analyses, `FisherResult` keeps data and prior information
separate and exposes named fixing, marginalization and rank diagnostics. Bounds
and derivative steps are not priors. A converged forest iteration does not
establish quadrature convergence or physical accuracy; consult the
[research prescription](../research/RESEARCH_BASELINE.md) for the qualified
calculation and remaining limitations.
