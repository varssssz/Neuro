# Exact Experimental Methodology — "Same scan, different pipeline"

Source of truth: `Source_Code/abide_data.py`, `Source_Code/swap_engine.py`, `Source_Code/run_all.py` (verbatim, in this package). Every statement below is read directly off that code, not reconstructed from memory.

## 1. Dataset

- **Name / source:** ABIDE I derivative tables distributed with the `compare-surf-tools` repository.
  - Primary GitHub: `neurodatascience/compare-surf-tools` (per earlier research-report search log; this was the clone used — exact remote URL and commit hash were **not recorded** at clone time, see §9 "What is NOT available").
  - Underlying raw data: ABIDE I, via NITRC/INDI (`fcon_1000.projects.nitrc.org/indi/abide/abide_I.html`).
- **Version:** not pinned to a specific repo tag/commit (see §9). Files read, verbatim from `abide_data.py`:
  - `FS5.1/cortical_fs5.1_measuresenigma_thickavg.csv`
  - `FS5.3/ABIDE_fs5.3_thickness.csv`
  - `FS6.0/{lh,rh}.DKTatlas_thickness_{h}.csv` and `FS6.0/{lh,rh}.aparc.thickness.table.csv` (DK, for the atlas-confound arm)
  - `ANTs/ABIDE_ants_thickness_data.csv`
  - `CIVET/ABIDE_civet2.1_thickness.csv` + `CIVET/DKT/DKT_parcel_map_FS_CIVET.csv` (region-name mapping)
  - `ABIDE_Phenotype.csv` (columns used: `SUB_ID`, `AGE_AT_SCAN`, `DX_GROUP`, `SITE_ID`, `SEX`)
  - `QC/manual/master_QC_table.csv` (`QC_maarten`, `QC_gleb`), `QC/auto/software/software_outlier_data.csv` (`outlier_ind`), `QC/auto/high-dim/high_dim_outlier_data.csv` (`OD_high_dim`)
- **Pipelines / software versions:** FreeSurfer 5.1, FreeSurfer 5.3, FreeSurfer 6.0 (DKT atlas; also DK atlas extracted separately as "FS60DK" to isolate the atlas confound), ANTs, CIVET 2.1. No reprocessing was performed — all tables are pre-computed outputs from the repository.
- **Prediction target:** `AGE_AT_SCAN` (chronological age, continuous).
- **Features:** cortical thickness, 62 regions common to all five pipelines (31 regions × 2 hemispheres), after name harmonization in `abide_data.py` (`_long_names`, per-pipeline column parsing). No volume, area, curvature or subcortical features were used.

### Inclusion/exclusion criteria (exact, from `abide_data.py`)
- A subject is included in a directed pair (source, target) if: present (non-NaN, non-duplicated index) in both pipeline tables, present in the phenotype table with non-missing age/diagnosis/site.
- **QC1 (primary):** additionally excluded if flagged by ANY of: `QC_maarten != 0` (manual lenient fail), `outlier_ind==1` in either pipeline of the pair (auto low-dimensional outlier), `OD_high_dim==1` (auto high-dimensional outlier). A subject absent from a QC table is **not** excluded by that table.
- **QC2 (stringent):** QC1 conditions AND `QC_gleb` in `{"1","-+1"}` (explicit manual stringent pass); unrated subjects are dropped.
- **QC0:** no exclusions at all.
- These three levels were defined in `abide_data.py` before any coverage result was inspected (per the docstring) and all three are reported (`transitions_main.csv`=QC1, `transitions_qc0.csv`, `transitions_qc2.csv`).

### Missing data
No imputation anywhere. A subject missing a feature table, a phenotype field, or failing any applicable QC filter is simply excluded from that specific directed pair's analysis (pairwise, not listwise across all 20 pairs).

## 2. Subject matching across pipelines

`abide_data.matched_ids(P, ph, qc, s, t, level)`: intersection of subject IDs present in source table, target table and phenotype table, minus QC-excluded IDs, intersected with the stringent-pass set if level=2. IDs are normalized by regex-extracting the last ≥5-digit run from whatever ID format each table uses (`sid()` function), so ID formats differing across tables (e.g. zero-padded vs not, site-prefixed vs not) are reconciled this way — **this normalization was not independently audited against a ground-truth ID crosswalk.**

## 3. Model

- **Primary:** `sklearn.linear_model.RidgeCV`, `alphas=np.logspace(-2, 4, 30)` (30 candidate penalties, 0.01 to 10,000), internal CV picks the penalty. Features standardized with `StandardScaler` fit on the fold's training subjects only.
- **Secondary (robustness check):** LightGBM (`lightgbm.LGBMRegressor`), `n_estimators=300, learning_rate=0.03, num_leaves=15, min_child_samples=20, subsample=0.8, subsample_freq=1, colsample_bytree=0.8, n_jobs=4`. Run for seeds {0,1,2} only (not all 5), on QC1 data, all 20 transitions.

## 4. Splits, cross-fitting, seeds

- 5-fold cross-fitting: `sklearn.model_selection.KFold(5, shuffle=True, random_state=seed)` splits all matched subjects into 5 folds; in each fold, the 80% not held out as test is further split 75/25 into proper-train / calibration via `np.random.default_rng(seed*1000+f).permutation`. So overall ≈60% train / 20% calibration / 20% test, and **every subject is a test subject in exactly one fold per seed.**
- **Seeds:** `SEEDS = [0,1,2,3,4]` for the ridge-based main analysis (`run_all.py`), strata, atlas, QC-sensitivity, mixed-pipeline training and LOSO stages. LightGBM used seeds `[0,1,2]`. Monitor stage used `seed=0` only. These are the literal values in the committed scripts; no other seeds were tried or discarded.
- No separate held-out "final" test set beyond this 5-fold, 5-seed scheme; all reported numbers are averages over the 5×5 (or 5×3) repeats as implemented in `summarize()`/`stage_*` functions.

## 5. Conformal prediction — exact formula and code

Absolute-residual split-conformal score, with the finite-sample-corrected quantile level:

```
q = Quantile( {|y_i − f(x_i)|}_{i in calibration} ;  ⌈(n+1)(1−α)⌉ / n ),   α = 0.10
Ĉ(x) = [ f(x) − q,  f(x) + q ]
```

Code (`swap_engine.cq`):
```python
def cq(res, alpha=ALPHA):
    n = len(res)
    lev = min(1.0, np.ceil((n + 1) * (1 - alpha)) / n)
    return np.quantile(res, lev, method="higher")
```
`ALPHA = 0.10` (nominal 90% coverage) throughout. `method="higher"` in `np.quantile` implements the "higher" interpolation, i.e. the finite-sample-valid empirical quantile for split conformal.

**Coverage** = fraction of test subjects with `|y_true − prediction| <= q`, where `prediction` and `q` are from the frozen source-trained model/calibration, applied to whichever pipeline's features are being scored (source or target).

**Interval width** = `2q` (years), per fold, averaged across folds/seeds.

### Conditions evaluated on the same test subjects (all paired)
| Code name | Meaning |
|---|---|
| `src` | source-pipeline features, frozen model, frozen q (reference) |
| `raw` | target-pipeline features, frozen model, frozen q (the swap, no fix) |
| `std` | target features z-scored with mean/SD of **unlabelled** target-pipeline training features of that fold; model and q still frozen |
| `src_site` / `sitestd` | per-site z-scoring (site stats from fold's training subjects; sites with <15 training subjects fall back to pipeline-global stats via `MIN_SITE=15`); model **re-fit** and q **re-computed** on site-standardized source features |
| `wcp` / `wcp_reg` | weighted conformal prediction (density-ratio reweighting), two settings (see §6) |
| ReForm-style (separate table) | k-labelled-subject recalibration (see §7) |

### Weighted conformal prediction (exact)
Density ratio from a logistic-regression classifier distinguishing source-calibration features from target-training features (both in the frozen source scaling):
```
w(x) = p̂(target|x) / (1 − p̂(target|x))  ×  n_cal / n_train
```
Weighted quantile (`swap_engine.wq`): cumulative-weight search for the (1−α)-weighted-mass point among calibration residuals; returns `inf` when the target weight exceeds what the calibration set's weighted mass can reach — **this `inf` is the literal mechanism by which "infinite interval" is defined and counted** (`wcp_frac_inf` = fraction of test subjects for which this search returns `inf`).
- `wcp`: `LogisticRegression(C=1.0)`.
- `wcp_reg` (sensitivity variant): `LogisticRegression(C=0.01)` (stronger L2) with predicted probabilities clipped to `[0.05, 0.95]` instead of `[0.001, 0.999]`.

### ReForm-style recalibration (exact — and explicitly an approximation, not the published algorithm)
For k ∈ {10, 20, 40, 80, 160} labelled target-pipeline subjects, drawn without replacement from the fold's pooled train+calibration subjects (100 draws per fold per seed, `rrep=100`):
```
e_lo = (f(x) − q) − y        # how far the true age fell below the lower bound
e_hi = y − (f(x) + q)        # how far above the upper bound
level = ⌈(k+1)(1 − α/2)⌉ / k      # one-sided, α/2 per side
d_lo = Quantile(e_lo over the k draw; level, method="higher")
d_hi = Quantile(e_hi over the k draw; level, method="higher")
Ĉ_recal(x) = [ f(x) − q − d_lo,  f(x) + q + d_hi ]
```
Two variants reported: `reform_raw` (recalibrating the raw-feature prediction) and `std_then_reform` (standardize features first with the `std` scaler, then apply the same recalibration). **This is the authors' own reading of the published ReForm method's description (separate lower/upper empirical-quantile corrections); the authors' own code was not obtained or run, so this is disclosed in the manuscript as "ReForm-style," not ReForm.**

## 6. Negative control

`CONTROL = {("FS60","FS53"), ("FS53","FS60")}` — FreeSurfer 6.0 ↔ FreeSurfer 5.3, chosen as a near-identity pair a priori (closest software versions among the five), flagged in every output table's `control` column, and reported separately in the manuscript rather than pooled with the other 18 transitions.

## 7. Atlas-confound isolation

FreeSurfer 6.0 Desikan-Killiany (DK) thickness was extracted separately from the same FS6.0 run that produced the DKT table (`load_fs6_dk`, reading `{h}.aparc.thickness.table.csv`), giving a "FS60DK" pseudo-pipeline. `ATLAS_PAIRS = [(FS60DK,FS60), (FS60,FS60DK), (FS60DK,FS53), (FS53,FS60DK), (FS60DK,FS51), (FS51,FS60DK)]` isolates: same software + different atlas (FS60DK↔FS60) vs same atlas + different software version (FS60DK↔FS5.x). Run with `stage_main(1, False, ATLAS_PAIRS, "atlas")` — i.e. without the wcp/reform/site extras (`extras=False`), ridge only.

## 8. Statistics

- **Paired bootstrap:** subject-level, seed-averaged coverage indicator, `B=2000` resamples, 95% CI from the 2.5th/97.5th percentiles (`swap_engine.boot_mean`), fixed `seed=11` for the bootstrap RNG.
- **McNemar:** exact binomial test (`scipy.stats.binomtest`) on discordant pairs, using **seed 0 only** (`cov["src"][0], cov[k][0]` — the `[0]` index selects the first seed's per-subject coverage array).
- **Wilson intervals:** standard formula (`swap_engine.wilson`, z=1.96), applied in `stage_strata` to the seed-mean per-subject coverage proportion for each (transition, stratum) cell, using the stratum's raw subject count as n — **this is an approximation** (it treats the seed-averaged coverage as if it were a single-seed binomial proportion with the full-stratum n; no multi-seed correction was applied to the Wilson interval itself). Strata with fewer than 20 subjects are dropped (`if n < 20: continue`).
- **Multiple comparisons:** no correction is applied inside the code; the manuscript reports both raw McNemar p-values and a Bonferroni comparison (α/18) computed afterward from the 18 non-control p-values, not inside `run_all.py`.
- **Seed variability:** `cov_{src,raw,std}_seed_sd` = standard deviation across the 5 seeds' per-transition mean coverage (a sensitivity diagnostic, not a formal CI).
- **Fold-level noise floor:** `frac_fold_below_088_{src,raw,std}` = fraction of individual fold-level coverage draws (25 draws = 5 folds × 5 seeds) below 0.88, reported including for the **source** pipeline itself, to show the sampling-noise floor given test folds of ≈170 subjects.

## 9. Monitoring (label-free shift detection; descriptive only)

Unbiased RBF-kernel MMD² (`swap_engine._rbf_mmd`; median-heuristic bandwidth from a 150-point subsample) and a two-sample KS test on predicted ages, computed per fold on batches of size `bsz` (50 primary; 25 and 100 as sensitivity, see `stage_monitor`). Thresholds set at the 95th percentile of the statistic's distribution on source-pipeline calibration-partition batches (`mmd_cal`/`ks_cal`); false-alarm rate checked on source-pipeline held-out test batches (`mmd_nul`/`ks_nul`); detection rate = fraction of alternative (target-pipeline) batches exceeding that threshold (`mmd_alt`/`ks_alt`). `seed=0` only; `nb=40` batches per condition per fold; `nref=200` reference subjects per fold. The transition (not the subject or batch) is the unit of analysis for any cross-transition comparison (20 points), and the manuscript treats this as descriptive, not inferential.

## 10. Mixed-pipeline training (Korbmacher-style)

Trained on the intersection of subjects present in **all five** pipelines simultaneously (N=669, see `stage_mixed`). "all" pools all five pipelines' features (including the target) with **one StandardScaler fit on the pooled mixture** (not per-pipeline); "loo" pools the four non-target pipelines only. Ridge, same `ALPHAS` grid, same 5-fold/5-seed scheme.

## 11. Leave-one-site-out (LOSO)

Run on the six largest sites by subject count in the all-pipeline-QC1 intersection (`stage_loso`, `top=6`). For each held-out site, train/calibrate on the other sites (75/25 split of the remaining subjects, 5 seeds), test on the held-out site's subjects; `src_site`/`sitestd` conditions use the held-out site's own (unlabelled) features for its site statistics (`ref = np.r_[tr, te]`).

## 12. Reproducibility actually performed

A second, independent full run of `run_all.py` on the same data directory (`results_full/`) was diffed against the first (`results/`): **every primary numeric table matched to 0.0 maximum absolute difference** for `transitions_main`, `strata_main`, `recal_main`, `lgbm`, `loso_by_site`, `loso_pooled`, `mixed_pipeline`, `monitors_b50`, `matched_N_summary`. (Per-subject and atlas/QC-sensitivity tables differed only in **column count**, not values, between runs — see `Reproducibility/rerun_diff.txt` for the exact diff output; this was not investigated further before this handover and should be checked by the auditor.) This confirms the pipeline is deterministic given the fixed seeds, not that the seeds themselves were chosen without hindsight (they were fixed before any run, per the script's literal content, but this cannot be externally verified from the code alone).

## What was NOT done / NOT available

- No git repository was initialized for this analysis code; there is **no commit hash**, because the environment is not a git repo (`Is a git repository: false`). Code is packaged here as plain files.
- The exact clone command/URL/commit of `compare-surf-tools` used at data-prep time was **not logged**; the data directory path used was `/root/compare-surf-tools/data` in this environment, and the repository identity should be re-confirmed by the auditor against `neurodatascience/compare-surf-tools` before relying on an exact file-for-file match.
- No participant-level identifiers are included in this package (see `Dataset_Documentation/`).
- Random seeds are the literal integers in the code (0–4 for the main analysis); no seed-selection or seed-sweep beyond what is listed was performed or discarded.
