# Step 01 review of implementation revision 4

Reviewed 2026-09-11 against `IMPLEMENTATION_STEP.md` revision 4 and the current
[implementation report](step-01.md). The user reported no further comments.

## Conclusion

No outstanding findings. The requested direct package layout and the previous
validation reliability issue are resolved. The reviewer recommends closing
Step 01 and is ready to prepare Step 02 when the user requests it. This review
does not dispatch another agent or introduce a next-step assignment.

## Resolution of prior feedback

- **U1, direct layout:** source is in `fishhighz/__init__.py`; the root `src/`
  directory is gone. Setuptools discovery explicitly includes only `fishhighz`
  and its descendants, with implicit namespaces disabled. Ignore rules and
  README describe the new layout, and the refreshed editable installation
  resolves to the correct source path from outside the checkout.
- **R1, stale validation:** each helper invocation allocates a unique run
  directory, stages maintained inputs without old build metadata, rejects
  ambiguous artifact selection, and creates new environments. Probes check
  installed metadata and package bytes against the exact wheel, including
  unexpected installed files. The two reported completed runs use distinct
  outputs/environments and successfully install all four wheels without skips.

## Evidence checked independently

- Verified all 24 SHA-256 entries in the revision-4 report, including both sets
  of source/helper snapshots and six artifacts. Repeated source entries match
  the current maintained files; both helper snapshots match the current helpers.
- Inspected the validation helper and probe, package metadata/discovery, README,
  ignore/manifest files, check runner, import test, and implementation report.
- Checked both reported `status.json` files: all 18 logged commands succeeded
  in each run. Inspected successful installation evidence and archive contents.
- Reran the executable check runner from an unrelated temporary directory with
  the development environment on PATH and PYTHONPATH cleared: one pytest test
  passed, Ruff lint passed, and Ruff formatting verification passed.
- Independently checked the editable module origin using isolated Python from
  outside the checkout.
- Created four new temporary environments and installed both direct wheels and
  both source-distribution-rebuilt wheels using `--no-index --no-deps`. All four
  isolated import, metadata, and full package-payload probes passed.
- Checked that wheel contents include only the intended package and dist-info,
  that source archives use the direct layout, and that their check runner is
  executable. Archived source and wheel module bytes match current source.

Fresh reviewer logs are retained in
`../.validation/review-r4-noe6srg2/` (local-only path: `../.validation/review-r4-noe6srg2/`).
The temporary independent installation environments were removed after review;
the implementer's retained environments and earlier evidence were left intact.

## Limits and scope

Independent execution used Python 3.13.15. The review reused the hashed
distribution artifacts rather than rerunning their network-dependent builds;
the two original isolated build/rebuild runs are documented in retained evidence.
No Python-version matrix or scientific validation was required by Step 01.

No implementation source, tests, helper scripts, or artifact files were changed
during review. The current detailed step remains revision 4. Numerical work and
DESI-2 scientific baseline capture have not started as part of this review.
