# Documentation migration

This page records the inventory-driven Markdown migration completed on 2026-09-22. The fixed inventory contains 163 moved-file records: 160 archived records and three research records. `migration.json` now contains 165 records in total after adding two derived README archive records.

The 163 source files were checked against both their recorded hashes and the preserved originals before moving. Scientific prose, equations, historical commands, numerical values and original paths were retained. After the move, only navigational links were repaired, with two minimal MyST rendering adjustments: an H1 title was added to the cumulative forest-weight report, and the archived design's Mermaid source was labelled and rendered as text. No figures, plots or numerical evidence were copied; `.validation/` remains outside the documentation tree.

The coordinator extracted the former README sections into two additional archive records. They are included at the end of `migration.json` as section-qualified extensions; their hashes cover the final extracted Markdown records after link portability repair and are not immutable source hashes:

| Former README material | Documentation destination |
| --- | --- |
| Introduction and native workflow | [`getting-started.md`](../getting-started.md) and the user guides |
| Development, distributions and reference capture | [`development/contributing.md`](contributing.md) and [`archive/readme-development.md`](../archive/readme-development.md) |
| Preparing synthetic arrays | [`methods/arrays.md`](../methods/arrays.md) |
| Covariance | [`methods/covariance.md`](../methods/covariance.md) |
| Fisher and compiled execution | [`methods/fisher.md`](../methods/fisher.md) |
| External models | [`methods/external.md`](../methods/external.md) |
| Templates | [`methods/templates.md`](../methods/templates.md) |
| Kaiser model | [`methods/kaiser.md`](../methods/kaiser.md) |
| Geometry and P1D | [`methods/geometry.md`](../methods/geometry.md) |
| Weights and noise | [`methods/weights.md`](../methods/weights.md) |
| Survey and raw inputs | [`methods/survey.md`](../methods/survey.md) |
| Historical Step 12 and Step 13 | [`archive/readme-validation.md`](../archive/readme-validation.md) |

Relative links between migrated Markdown records now point to their archive or research destinations. Links to external literature and versioned upstream sources are retained. Links to `.validation/` outputs, absolute workstation paths, source files and other evidence outside the Sphinx documentation tree are rendered as explicit `local-only path` text with their original path, so the HTML build does not expose broken local hyperlinks.

The archive indexes are [original implementation](../archive/original-implementation.md), [weighting diagnostics](../archive/weighting-diagnostics.md), [compatibility](../archive/compatibility.md), [S1--S5](../archive/stages-s1-s5.md) and [standalone](../archive/standalone.md). They preserve reachability without treating archived plans as current instructions.
