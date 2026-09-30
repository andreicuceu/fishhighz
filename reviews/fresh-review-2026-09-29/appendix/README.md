# Appendix: fishhighz vs lyaforecast review (2026-09-29)

**Contents of this directory:**
- `runs_summary.csv` — one row per E-set case: 108 rows (107 completed plus the failed `fh_t_w_mcdonald`) and 92 columns. Columns:
  - `case`, `code`, `group`, `desc`, `wall_s`, `peak_rss_mb`;
  - combined σ∥, σ⊥ and their Δ against the own-code baseline;
  - per-bin joint σ;
  - the key spectra per bin and combined: lya(qso)², lya(lbg)², qso², lya(qso)×qso and lbg².

  The 12 E3 corrective cases are not in this file. Their tables are in `../work/E3_corrections.md`, and their raw output is in the run directories listed below.
- `pytest_ruff.log` — pytest summary (1576 passed, 7 failed; see summary §3.1, fh-2) and the Ruff result (clean).
- `figures/`, produced from the E set; they do not include the E3 cases:
  - `fig1_perbin_joint_fh_vs_lf.png` — per-bin joint σ for fh, lf, fixed-compat and the full chain, with residuals against lf.
  - `fig2_chain_waterfall.png` — the cumulative convention chain fh → lf, combined, with the Q = 1 and SE07 endpoints.
  - `fig3_delta_sigma_all_variations.png` — Δσ/σ (joint, combined) for all fh variations, coloured by group (symlog).
  - `fig4_forest_auto_weighting.png` — lya(qso) and lya(lbg) auto σ per bin for the early, inverse-variance, McDonald, legacy-3, ME07-mode and lf weightings.
  - `fig5_fh_over_lf_ratio_matched.png` — σ_fh/σ_lf over the matched variations.

**Supporting material in `../work/`:**
- `A1_ME07.md`, `A2_MW11.md`, `A3_FR14.md` — paper equation digests.
- `map_fishhighz.md`, `map_lyaforecast.md` — unverified seed maps, corrected in B.
- `B_inventory.md` — the stage-by-stage inventory.
- `C1_weights.md`, `C2_signal_fisher.md`, `C3_survey_noise.md` — scientific assessments.
- `D1_numerical_core.md`, `D2_orchestration_inputs.md` — implementation audits.
- `E_results.md` — all E-set tables and checks.
- `F1_verification.md` — adversarial verification. It is authoritative where it conflicts with earlier files.
- `E3_corrections.md` — the 12 corrective runs:
  - k⊥ Nyquist with the per-bin sightline density, and the MW11 n̄_eff variant;
  - range-averaged, truncated LBG-forest density;
  - b·D bias runs and DESI DR1 b_F/β_F;
  - per-row cell-width normalisation;
  - the grid-converged lyaforecast reference.

**Harness and raw runs** (`/pscratch/sd/a/acuceu/fishhighz-review/`; `$SCRATCH`, subject to purge):
- `harness/` — the drivers `run_fishhighz*.py` and `run_lyaforecast*.py`, and `run_case.sh`.
- Case definitions: `cases.py`, `cases2.py`, `cases3.py`.
- Toggles: `toggles.py`, `toggles2.py`, `toggles3.py`, documented in `TOGGLES.md`.
- Job lists: `jobs*.txt`.
- Analysis scripts: `analyse*.py`.
- Environment and schema: `README.md`.
- `runs/<case>/` — `result.json` (per-bin joint and spectra σ, Fisher matrices, forest-weight records, timing, toggle and INI-edit logs), `effective.ini`, and logs in `runs/logs/`.
- `analysis2/` — the intermediate report fragments behind E.
- `papers/` — PDFs of ME07, MW11 and FR14.

