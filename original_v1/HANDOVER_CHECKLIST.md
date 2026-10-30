# Handover checklist — "Same scan, different pipeline"

## Included and verified present in this package

| Item requested | Included? | Location |
|---|---|---|
| Complete manuscript (docx, unaltered) | Yes | `Manuscript/Same_scan_different_pipeline_manuscript.docx` |
| Manuscript as PDF | Yes | `Manuscript/Same_scan_different_pipeline_manuscript.pdf` |
| Manuscript as plain text (all sections + tables, for machine parsing) | Yes | `Manuscript/manuscript_full_text.txt` |
| Exact experimental methodology, formulas | Yes | `Reproducibility/METHODOLOGY.md` |
| All 20 transition results, unrounded | Yes | `Experimental_Results/transitions_main.csv` |
| All 18 non-control + 2 control results | Yes | same file, `control` column distinguishes them |
| Coverage per source-target pair | Yes | `transitions_main.csv` columns `cov_src/raw/std/sitestd/wcp/wcp_reg` |
| Nominal coverage level | Yes | documented as α=0.10 (90%) in METHODOLOGY.md, fixed throughout |
| Interval widths | Yes | `width_src`, `width_src_site`, `recal_main.csv` `mean_width` |
| Infinite-interval percentages | Yes | `wcp_frac_inf`, `wcp_reg_frac_inf` in `transitions_main.csv` |
| Sample sizes per comparison | Yes | `N` column in every transitions/lgbm/mixed/loso file; `matched_N_*.csv` |
| P-values (uncorrected) | Yes | `mcnemar_p_raw`, `mcnemar_p_std` in `transitions_main.csv` |
| Adjusted p-values | **Partially** | Bonferroni comparison (α/18) is stated as a count in the manuscript text and in `CLAIMS_VS_EVIDENCE.md`, but a full per-transition adjusted-p column was **not** written to CSV — recomputable in one line from `mcnemar_p_std` |
| Confidence intervals | Yes | `cov_*_lo/hi`, `d_*_lo/hi` (bootstrap 95% CI) in `transitions_main.csv`; Wilson CI in `strata_main.csv` |
| Effect sizes | Partial | reported as coverage differences (`d_raw`, `d_std`, `d_sitestd`) and MAE differences; no standardized effect size (e.g. Cohen's d) was computed — **not included, not invented** |
| Subgroup analysis (age, diagnosis) | Yes | `strata_main.csv` |
| Leave-one-site-out results | Yes | `loso_by_site.csv`, `loso_pooled.csv` |
| QC sensitivity results | Yes | `transitions_qc0.csv`, `transitions_qc2.csv` |
| Mixed-pipeline training results | Yes | `mixed_pipeline.csv` |
| Failed/excluded/inconclusive experiments | Yes, disclosed | weighted CP's near-total failure documented in data and text; strata with n<20 silently absent (documented in `Experimental_Results/README.md`); no experiment was fully "failed" in the sense of being discarded |
| All Python scripts used to generate results | Yes | `Source_Code/` (abide_data.py, swap_engine.py, run_all.py, make_figs.py, build_docx.py) |
| Jupyter notebooks | **None exist** — all analysis was plain `.py` scripts, not notebooks |
| requirements.txt | Yes | `Source_Code/requirements.txt` |
| Execution instructions | Yes | `Source_Code/README.md` |
| Random seeds | Yes, exact values | `Source_Code/README.md` and `METHODOLOGY.md` §4 |
| Dataset download instructions | Yes | `Dataset_Documentation/README.md` |
| GitHub repo URL + commit hash for the analysis code | **Not available — honestly stated, not invented.** No git repo was ever initialized for this code. |
| Exact commit of the `compare-surf-tools` data repo used | **Not available — not recorded at the time.** Repo name/URL given; commit hash not logged. |
| Figures (all 4) | Yes | `Figures_and_Tables/fig1..4_*.png` |
| Summary/table CSVs used to build figures | Yes | `Figures_and_Tables/table_transitions.csv`, `table_recal_summary.csv` |
| Literature audit: all ~41 searches, exact queries, dates | Yes | `Literature_Review/LITERATURE_AUDIT.md` |
| Six closest competing papers, full detail | Yes | same file |
| ReForm 2026 and Samanta 2026 specifics | Yes, with evidence tags | same file |
| Preregistration / kill-criterion record | **No formal preregistration exists — stated explicitly, not fabricated.** Closest artifact (a dated research-agent report) is described and its limits as a "preregistration" are spelled out. |
| Claims-vs-evidence table | Yes | `Literature_Review/CLAIMS_VS_EVIDENCE.md` |
| Train/calibration/test subject identifiers | **Not separately exported** — per-subject outcome files include subject IDs, but fold/seed role assignment was not saved to disk (reconstructible by re-running the code with the same seeds; documented in `Dataset_Documentation/README.md`) |
| Reproducibility check (second independent run) | Yes | `Reproducibility/rerun_diff.txt`, `Reproducibility/full_run_log.txt` |

## Known issues to resolve before submission (carried over honestly, not smoothed over)

1. **Reference list gaps.** Several manuscript references ([2], [4]–[9], [11]) still need final author/title/DOI verification; partially resolved across the two research passes (see `Literature_Review/LITERATURE_AUDIT.md`) but not fully reconciled into the manuscript's reference list as of this handover.
2. **Denominator inconsistency** (20 vs 18 transitions) in a couple of ridge/LightGBM comparison sentences — see `CLAIMS_VS_EVIDENCE.md` note.
3. **ReForm journal full text never obtained** — both the MIT Press page and the PMC copy were blocked (bot detection / reCAPTCHA) on multiple independent attempts. The calibration-size rule (40/80/120 vs a conflicting ≥120 snippet) was never reconciled.
4. **No second dataset.** Generalizability beyond ABIDE I is not demonstrated; OASIS-3 was investigated as a fallback and found to have unresolved access/approval terms and no ready multi-pipeline derivative set.
5. **Forward-citation chasing failed entirely** (Semantic Scholar 429 on every attempt, both passes) — scoop risk from very recent/unindexed work cannot be fully ruled out.
6. **No formal preregistration.** State this plainly if asked by reviewers.

## What was explicitly NOT manufactured

No files, results, citations, preregistration records, commit hashes, or experimental evidence were invented to fill gaps in this package. Every "not available" line above reflects something that genuinely was not produced or recorded during the project, stated as such rather than worked around.
