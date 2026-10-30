# Source_Code/ — exact scripts used

All scripts are included verbatim (not rewritten for this handover; same files already delivered earlier in this project).

- `abide_data.py` — data loading, QC-flag loading, subject matching (§1-2 of METHODOLOGY.md).
- `swap_engine.py` — the analysis engine: conformal quantile, weighted CP, ReForm-style recalibration, standardization, monitors, mixed-pipeline training (§3-11 of METHODOLOGY.md).
- `run_all.py` — orchestrates all stages and writes every CSV in `Experimental_Results/`. **This is the literal script that produced every reported number.** CLI: `python run_all.py DATA_DIR OUT_DIR [comma,separated,stages]`.
- `make_figs.py` — builds the 4 manuscript figures and 2 summary tables from the `Experimental_Results/` CSVs (matplotlib).
- `build_docx.py` — builds the manuscript docx from the CSVs + figures, using `python-docx`. Every numeric sentence in the manuscript is computed in this script directly from the CSV files at build time (not hand-typed), so re-running it on the same CSVs reproduces the same text.
- `requirements.txt` — pinned package versions actually used.
- `pilot_pipeline_swap.py`, `pilot_extras.py` — the user's own earlier pilot scripts (kept for reference/provenance; superseded by `run_all.py`/`swap_engine.py` for all reported results). `pilot_transitions.csv`, `pilot_strata.csv` are their outputs from an earlier stage of this project. `pilot_rerun.log` is a re-run of the pilot scripts from a fresh clone, done to check reproducibility of the user's original pilot numbers — it found a discrepancy (see `Reproducibility/METHODOLOGY.md` is silent on this; see the manuscript's Section 4.2 note and `Literature_Review/` — the pilot's originally quoted numbers (e.g. FS6.0→ANTs 0.555) did not reproduce from a fresh clone (0.520 obtained); the manuscript explicitly flags this and uses only the `run_all.py`-produced numbers).

## No git repository / no commit hash

This project was never initialized as a git repository in the working environment (`Is a git repository: false`). There is **no GitHub repo URL and no commit hash for this analysis code** to provide — this is stated explicitly rather than invented. If you need one, the code would need to be committed to a new repository before a hash exists.

## Execution instructions

```bash
# 1. Get the data (not included in this package; see Dataset_Documentation/)
git clone https://github.com/neurodatascience/compare-surf-tools
# NOTE: the exact commit of compare-surf-tools used to produce Experimental_Results/
# was not recorded at the time (see METHODOLOGY.md §9). Re-cloning today may not be
# byte-identical to what was used. Verify file names/shapes match abide_data.py's
# expectations before trusting a re-run to reproduce these exact numbers.

# 2. Environment
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt

# 3. Run everything (this is what produced Experimental_Results/; ran in background,
#    full run time was on the order of tens of minutes to ~1-2 hours on the machine used)
python run_all.py compare-surf-tools/data results_full

# 4. Figures + manuscript
python make_figs.py results_full figs
python build_docx.py
```

## Random seeds

Literal values in the code, no sweep performed: `SEEDS = [0,1,2,3,4]` (run_all.py) for the main ridge-based analysis, strata, atlas, QC-sensitivity, mixed-pipeline and LOSO stages; `[0,1,2]` for LightGBM; `seed=0` for the monitor stage. Internal RNGs (bootstrap, fold permutation, ReForm-style resampling) are derived deterministically from these via `np.random.default_rng(seed*1000+f)` / `np.random.default_rng(seed+7)` / fixed `seed=11` for the bootstrap — see `swap_engine.py` for the exact derivations.
