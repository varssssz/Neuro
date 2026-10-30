# AUDIT_BASELINE.md — Phase 1: frozen and validated baseline (v1)

Status labels used throughout:
**VERIFIED-RERUN**: re-executed from raw data in this audit and matched the stored output.
**VERIFIED-CODE**: established by reading or instrumenting the v1 code, without running new experiments.
**NOT RE-EXECUTED**: stored output accepted as-is for now.
**OPEN**: unresolved; assigned to a later phase.

Nothing under `original_v1/` was modified. It is write-protected, and every file's hash is in `provenance/ORIGINAL_V1_SHA256.txt`.

---

## 1. What was frozen

| Item | Value |
|---|---|
| Git | Repository initialised **during this audit** (2026-10-30). No earlier git history exists, and none is implied. Baseline commit tagged `v1-baseline`. |
| Original package | `original_v1/` (unzipped, read-only) and `original_v1_package.zip` (byte copy of the upload) |
| Data | `neurodatascience/compare-surf-tools` @ `debcc78a346d0640295245ee331e4da47dda9e46` (commit date 2023-01-30, clean tree). Kept **outside** the repo because ABIDE derivatives are not redistributed. SHA-256 of all 13 input files: `provenance/DATA_PROVENANCE.md`. |
| Data-commit status | The commit used for the original run was **not logged**. Re-running the v1 code on this commit reproduces the stored v1 outputs exactly (§3), so `debcc78` is the best-supported identification, not a logged fact. |
| Environment | `provenance/ENVIRONMENT.txt`. numpy 2.5.3, scipy 1.18.1, scikit-learn 1.9.1, pandas 3.0.5, matplotlib 3.11.2 and python-docx 1.2.0 all equal the pins in `requirements.txt`. **LightGBM 4.7.0 could not be installed** (PyPI blocked from this session), so `lgbm.csv` cannot be re-executed here. |

## 2. Inventory

- **Code:** `abide_data.py` (load, QC, matching), `swap_engine.py` (models, conformal prediction, adaptations, monitors, mixed training), `run_all.py` (9 stages that write every CSV), `make_figs.py`, `build_docx.py`. There are also pilot scripts, which no reported number uses.
- **Results:** 23 CSV/TXT files in `Experimental_Results/`, 4 figures, 2 summary tables.
- **Manuscript:** docx, pdf, txt.
- **File-to-code map:** `provenance/RESULTS_PROVENANCE_MANIFEST.csv`. This maps each result file to its stage call, code path, seeds, QC level, units and Phase 1 status.
- **Number-to-file map:** `provenance/HEADLINE_NUMBER_PROVENANCE.csv`. This maps each manuscript number to its file, column, aggregation and denominator.

## 3. Reproduction of the baseline (VERIFIED-RERUN)

| Output | Check | Result |
|---|---|---|
| `transitions_main.csv`, all 43 numeric columns, 20 transitions | Re-executed v1 `run_transition` + unmodified `summarize()` (recalibration draws skipped) | max \|diff\| = 4×10⁻¹⁵ |
| `recal_main.csv` (reform_raw, all k) | FS60→ANTS, CIVET→FS51 | max \|diff\| = 1×10⁻¹⁶ |
| `monitors_b50.csv` (MMD raw/std detection) | FS60→ANTS, CIVET→FS51 | identical |
| `matched_N_summary.csv` | 26 pairs × 3 QC levels | identical |
| `per_subject_main.csv` | subject IDs, order, age, site, diagnosis vs phenotype, all 20 transitions | identical |
| `mixed_pipeline.csv` N | 669 | identical |

Not re-executed: atlas, QC0/QC2 and LOSO values; LightGBM; strata. Strata derive deterministically from `per_subject_main.csv`, and manuscript Table 2 was recomputed from that file and matches.

## 4. Transitions and denominators (VERIFIED-CODE)

- **Main arm:** 20 directed transitions = 5 pipelines × 4 targets (`PAIRS` in `run_all.py`). **None was excluded or dropped.**
- **Labelled "control":** `CONTROL = {(FS60,FS53),(FS53,FS60)}`. The 2 control rows + 18 = 20. "18 non-control" is therefore a reporting choice, not an exclusion.
- **Atlas arm:** 6 additional pairs involving FS60DK, in a separate file and not part of the 20.
- **Inconsistent denominators found in the manuscript:**
  - The **LightGBM sentence** (0.716/0.854) and its ridge comparator (0.531/0.859) use **all 20**. Over 18: LightGBM 0.695/0.848, ridge 0.489/0.853.
  - The **LOSO sentence** (0.867/0.504/0.818/0.910/0.919; worst site 0.672) uses **all 20**. Over 18: 0.865/0.462/0.810/0.908/0.917; worst site 0.666.
  - Everything else (headline, McNemar counts, weighted-CP infinite share 0.896/0.794, recalibration, Table 2, MMD 0.139) uses **18**.

## 5. Data-flow table (VERIFIED-CODE)

Notation: in each cross-fitting fold, **tr** = proper-training (≈60%), **ca** = source calibration (≈20%), **te** = test (≈20%). Xs and Xt are the source and target features of the *same* subjects.

| Step / condition | Fit or estimated on | Labels used | Applied to | Notes |
|---|---|---|---|---|
| Subject matching | intersection of both pipelines' subjects, with phenotype; QC exclusions | — | — | `matched_ids`; no imputation (zero rows with NaN in any pipeline) |
| Source scaler `sc_s` | Xs[tr] | no | Xs, Xt everywhere ("raw" = target features in source scaling) | |
| Ridge model + λ | Xs[tr], y[tr] | tr only | — | RidgeCV efficient LOO inside tr; fixed 30-value grid; **ca is not used for fitting or for λ** |
| Conformal q | residuals on Xs[ca] | ca | all test conditions | one rank more conservative than the textbook rank (§8) |
| **src / raw** coverage | — | te scoring only | Xs[te] / Xt[te] | same te subjects for every condition in a fold |
| **std** (per-pipeline) | `sc_t` on **Xt[tr]** | no | Xt[te] | target statistics come from the *same subjects* as source training: population-matched, an optimistic case |
| **sitestd** (per-site) | site means/SDs from X[tr] by site; sites with <15 tr subjects use global tr stats | refit on tr, q on ca | Zt[te] | **refits the model**; 2–8 sites per fold fall back to global stats, covering 12.3% of test subjects |
| **WCP** | logistic density ratio: Xs[ca] vs Xt[tr] | ca residuals | weights evaluated at Xt[te] | test features enter only through each test point's own weight (standard WCP) |
| **Recalibration** (reform_raw / std_then_reform) | k subjects drawn from **pool = tr ∪ ca** using Xt[pool] | y[pool] | te | pool is 74.9% training subjects; pool ∩ te = ∅ |
| **Monitors** | reference = 200 tr subjects (Xs); threshold = 95th percentile over 200 batches from ca (Xs); false-alarm check on te (Xs) | no | alternative batches = Xt[te] (raw or std) | uses the **same te subjects** as coverage; 40 overlapping batches of 50 per fold. KS reference uses *in-sample* training predictions (threshold built the same way, so consistent). Details in Phase 4. |
| **LOSO** | tr/ca = 75/25 random split of the other sites; te = held-out site | tr, ca | te | **site statistics use ref = tr ∪ te**, i.e. the held-out site's own unlabelled test features (transductive) |
| **Mixed-pipeline** | pooled Xd[p][tr] across pipelines; one pooled scaler; q on pooled ca | tr, ca | Xd[t][te] | N = 669 (all-5 intersection); "mix_all" includes the target pipeline's training rows; pooled q gives per-target coverage that is empirical only |

## 6. Split-integrity table (VERIFIED-CODE; reconstructed from v1 code + seeds)

Exported: `audit/phase1/out/subject_roles_main_and_atlas_QC1.csv.gz`. It holds every subject's role (train/cal/test) in each of the 5 folds, for 26 pairs × 5 seeds. Per-fold checks: `split_integrity_by_fold.csv`, `loso_integrity.csv`.

| Check | Main + atlas (26 pairs × 5 seeds × 5 folds = 650 folds) | LOSO (20 pairs × 6 sites × 5 seeds) |
|---|---|---|
| train ∩ cal | 0 | 0 |
| train ∩ test | 0 | 0 |
| cal ∩ test | 0 | 0 |
| each subject tested exactly once per seed | yes (asserted) | n/a (whole site is test) |
| fold sizes | train 414–580, cal 139–194, test 138–194 | test = site n (44–167) |
| recalibration pool ∩ test | 0 | — |
| recalibration pool ∩ train | **74.9% of the pool** | — |
| test features used for normalisation statistics | no | **yes** (site_z ref includes te) |

## 7. Leakage vs. transductive vs. design choice

| Finding | Classification | Consequence |
|---|---|---|
| Test **labels** used anywhere except scoring | **None found** | — |
| Recalibration pool contains training subjects (`pool = np.r_[tr, ca]`) | **Protocol violation (role overlap), not test leakage.** Training-subject residuals are partly in-sample. | Phase 2 rerun. Pre-audit check: effect ≤0.4 points on 3 transitions; to be confirmed. |
| §3.8 of the manuscript ("the same subject never appears in train, calibration and test") | **True for train/cal/test; false for the recalibration role** | text fix |
| LOSO per-site standardisation uses the held-out site's test features | **Transductive** (label-free, uses the test batch) | disclose as transductive; contradicts §3.8 ("site statistics are computed from training subjects") |
| Per-pipeline std uses Xt of the *same subjects* as source training | **Design choice, optimistic** (population-matched) | Phase 4 single sensitivity analysis |
| Per-site std refits the model | **Not frozen-model adaptation** | reframe (Phase 4) |
| WCP weights at test points; monitors on test batches | **Inherent to the method** (label-free, per-point or per-batch) | describe accurately |
| Site effects shared across folds (non-LOSO) | **Design choice**, already disclosed in §3.8 | keep |

## 8. Other baseline checks

- **Feature schema:** 62 common features (31 per hemisphere), identical names and order in every pipeline (`feature_schema_62.csv`). Loader accounting (`loader_accounting.csv`): no rows lacked a subject ID, no duplicate IDs, no NaN rows. FS5.x and FS60-DK load 68 regions (DK atlas); FS60 loads 63 (DKT); ANTs and CIVET load 62.
- **Feature alignment** (`feature_alignment.csv`):
  - Within FreeSurfer, every region's best-correlated counterpart is its namesake (same-region r = 0.84–0.97).
  - For ANTs pairs, 24–32 of 62 regions are better correlated with a *neighbouring or contralateral* region than with the namesake (mean same-region r = 0.42–0.46, minimum 0.09). CIVET pairs show 9–13 such regions.
  - There is **no systematic pattern indicating a mapping bug** (no consistent hemisphere swap or offset). This is consistent with weak regional agreement plus a global thickness factor.
  - **Implication:** ANTs/CIVET "same-name" features are not equivalent measurements. This is a caveat for the paper, not a code defect.
- **QC:** QC1 exclusion sets are applied as described: manual fail 22, high-dimensional flag 105, auto low-dimensional 55–71 per pipeline. Pipeline names in the QC table match the code's mapping, so automatic flags are genuinely applied. Two caveats:
  - The high-dimensional flag is not pipeline-specific; a subject is excluded if flagged at all.
  - FS60-DK reuses the FS60 flags.
- **Conformal quantile:** `cq()` uses `np.quantile(r, ceil((n+1)(1-α))/n, method="higher")`. numpy's "higher" method positions at (n−1)·level, so this returns the **(⌈(n+1)(1−α)⌉+1)-th** smallest score, one rank above the manuscript's equation, except when capped at the maximum. The intervals are valid but slightly conservative. The same construction in the recalibration code explains why coverage tracks the discrete finite-sample levels. Full treatment and tests: Phase 2.
- **Same held-out subjects across conditions:**
  - src, raw, std, sitestd, WCP and recalibration all score the **same te subjects** within a fold and seed. Paired comparisons are therefore legitimate.
  - LOSO and mixed-pipeline use different subject sets (by design).
  - Different transitions use different but heavily overlapping subject sets (N 692–967 from about 1,100 subjects), so **transitions are not independent units**. This is for Phase 3.

## 9. Open items handed to later phases

| # | Item | Phase |
|---|---|---|
| 1 | Clean four-way split for recalibration; separate pool effect from quantile-rank effect | 2 |
| 2 | Control relabelling (FS53↔FS60 vs FS60DK↔FS60 r = 0.974) | 3 |
| 3 | Dependence across seeds, transitions and shared subjects in inference; seed-0-only McNemar | 3 |
| 4 | Per-pipeline std with disjoint subjects for target statistics; transductive/refit taxonomy | 4 |
| 5 | Monitor validation (overlapping batches, same test subjects, in-sample KS reference) | 4 |
| 6 | Age-conditional baseline miscoverage (≥25: about 48 subjects per transition, 0.277 at source) | 5 |
| 7 | LightGBM not re-executable in this environment | report as a limitation, or rerun where PyPI is reachable |

## 10. Reproduce this audit

```bash
cd neurovision
D=/path/to/compare-surf-tools/data   # @ debcc78
python -I audit/phase1/split_integrity.py      original_v1/Source_Code original_v1/Experimental_Results $D audit/phase1/out
python -I audit/phase1/reproduce_main.py       original_v1/Source_Code original_v1/Experimental_Results $D audit/phase1/out
python -I audit/phase1/reproduce_recal_monitor.py original_v1/Source_Code original_v1/Experimental_Results $D
python -I audit/phase1/build_manifest.py       original_v1 provenance
```
