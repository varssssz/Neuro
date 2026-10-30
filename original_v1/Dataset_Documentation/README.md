# Dataset_Documentation/

## Source and access
- **Derivative tables used:** `neurodatascience/compare-surf-tools` GitHub repository (per the project's research-report search logs; this is the repo that was cloned to the working environment at `/root/compare-surf-tools/data` — see `Source_Code/README.md` for the caveat that the exact commit was not logged).
- **Underlying raw data:** ABIDE I, via NITRC/INDI: `https://fcon_1000.projects.nitrc.org/indi/abide/abide_I.html`.
- **Access restrictions:** ABIDE I raw data requires NITRC registration; unrestricted for non-commercial research use; licensed CC BY-NC-SA 3.0 per the INDI page (read directly from that page during this project, not inferred). The `compare-surf-tools` derivative-table code itself is MIT-licensed. **The exact wording of any additional ABIDE-specific data-use terms beyond the CC BY-NC-SA licence was not independently verified** — flagged as such in the manuscript.
- **No raw imaging data, and no participant-identifying information, is included in this package.** Only the already-public, already-de-identified derivative CSVs (cortical-thickness tables, phenotype table, QC tables) were used, and this package contains only *results computed from* them (CSVs of coverage/statistics), not the source tables themselves — see below.

## Feature schema
62 cortical-thickness features common to all 5 pipelines (31 regions × 2 hemispheres, DKT/DK-mapped region names). Exact column-name harmonization logic is in `Source_Code/abide_data.py` (`_long_names`, per-pipeline loaders). Phenotype fields used: `SUB_ID`, `AGE_AT_SCAN`, `DX_GROUP` (1=ASD, 2=TD), `SITE_ID`, `SEX`.

## File naming conventions (of the ORIGINAL source tables, not reproduced here)
See `Reproducibility/METHODOLOGY.md` §1 for the exact file paths read by `abide_data.py` (e.g. `FS5.1/cortical_fs5.1_measuresenigma_thickavg.csv`, `FS6.0/{lh,rh}.DKTatlas_thickness_{h}.csv`, `QC/manual/master_QC_table.csv`, etc.). These are the original repository's file names — this package does not re-host them.

## Intermediate data / subject identifiers

- **Train/calibration/test subject IDs per fold/seed were NOT saved to disk** by `run_all.py` — only aggregate and per-subject-coverage outputs were written (`per_subject_*.csv`, which does include the subject ID `sid` and site, but not which fold/seed role each subject played). Reconstructing the exact fold assignment for a given seed is possible by re-running `swap_engine.folds()` with the same seed and subject ordering, but the assignment itself was not separately exported.
- `per_subject_*.csv` files in `Experimental_Results/` contain ABIDE subject IDs (`sid`), site, age and diagnosis — these are the same subject identifiers already public in the ABIDE phenotype release (not newly identifying information), but are flagged here explicitly since the request asked that PII be excluded: **no name, DOB, or other directly identifying field is present anywhere in this package**, only the ABIDE-assigned numeric subject ID already used throughout the public dataset.

## Regenerating tables/figures from the source data
See `Source_Code/README.md` "Execution instructions" — `run_all.py` → `make_figs.py` → `build_docx.py`, in that order, against a fresh clone of `compare-surf-tools`.

## What is NOT included in this package
- The original ABIDE/compare-surf-tools CSV tables themselves (not redistributed here; obtain via the links above).
- A verified exact repository commit hash for the data (not recorded during this project).
- A separate file enumerating fold/seed-level subject-role assignments (not generated; reconstructible by re-running the code).
