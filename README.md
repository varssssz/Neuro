# Same scan, different pipeline — analysis code

Code and aggregate results for *"Same scan, different pipeline: coverage of frozen conformal brain-age intervals when structural MRI is re-processed"* (manuscript in submission).

A ridge brain-age model with a 90% split-conformal interval is trained on one structural MRI pipeline (FreeSurfer 5.1/5.3/6.0, ANTs, CIVET 2.1; ABIDE I cortical thickness) and applied, frozen, to the same subjects processed by another pipeline. The code evaluates coverage, label-free adaptations, labelled recalibration, a drift monitor, and age-conditional coverage.

## Data (not included)
ABIDE I derivatives come from [compare-surf-tools](https://github.com/neurodatascience/compare-surf-tools), commit `debcc78`. ABIDE data are CC BY-NC-SA 3.0. Subject-level outputs are regenerated locally and are excluded from this repository (`.gitignore`).

## Setup
```bash
git clone https://github.com/neurodatascience/compare-surf-tools
cd compare-surf-tools
git checkout debcc78
cd ..
git clone [this repository URL] neurovision
cd neurovision
python -m venv .venv
# Windows: .venv\Scripts\activate    macOS/Linux: source .venv/bin/activate
pip install -r original_v1/Source_Code/requirements.txt pytest
python -m pytest code_v2/tests        # expect 44 passed
```

## Reproduce all results (≈15 min on 2 cores)
Run from the repository root; `D` = `../compare-surf-tools/data`.
```bash
python code_v2/run_recalibration_v2.py ../compare-surf-tools/data results_v2/recalibration 2
python code_v2/summarise_recalibration.py results_v2/recalibration
python code_v2/run_frozen_v2.py ../compare-surf-tools/data results_v2/frozen 2
python code_v2/run_secondary_textbook.py original_v1/Source_Code ../compare-surf-tools/data results_v2/secondary_textbook
python code_v2/phase3_inference.py results_v2 original_v1/Experimental_Results audit/phase1/out
python code_v2/phase4_std_monitor.py results_v2 original_v1 ../compare-surf-tools/data
python code_v2/phase5_age.py results_v2
python code_v2/make_figs_v2.py results_v2 manuscript_v2/figures
python code_v2/build_numbers.py results_v2 original_v1/Experimental_Results manuscript_v2/numbers.json
python code_v2/render_manuscript.py results_v2 manuscript_v2   # needs pandoc
```
Optional: `python code_v2/run_lgbm_textbook.py original_v1/Source_Code ../compare-surf-tools/data results_v2/lgbm_textbook` (needs lightgbm).

`manuscript_v2/numbers.json` records every number reported in the manuscript with its source file and column.

## Layout
| Path | Contents |
|---|---|
| `code_v2/` | Corrected analysis engine, analysis scripts, tests |
| `results_v2/` | Aggregate results (recalibration, inference, standardisation/monitor, age, sensitivity) |
| `original_v1/` | Initial analysis code and aggregate outputs (superseded; kept for transparency) |
| `audit/phase1/` | Baseline reproduction and split-integrity checks |
| `provenance/` | Environment, data provenance and result manifests |
| `*.md` reports | Audit and correction reports (recalibration, inference, standardisation, monitoring, age, references, changelog) |

## Licence
Code: MIT (see `LICENSE`). Data: not included; see the ABIDE and compare-surf-tools licences.
