# STANDARDISATION_SENSITIVITY.md — Phase 4a

Code: `code_v2/swap_engine_v2.run_frozen` (conditions `std`, `std_disjoint`, `std_same_sub`) and `code_v2/phase4_std_monitor.py`.
Results: `results_v2/standardisation_monitor/std_sensitivity_*.csv`.
Textbook rank; 18 non-control transitions; joint site-cluster bootstrap CIs.

## 1. Who contributes to each standardisation estimate (VERIFIED-CODE)

| Condition | Statistics estimated from | Same people as model training? | Uses test-fold features? | Labels? | Model refit? |
|---|---|---|---|---|---|
| Per-pipeline std (v1, `std`) | target features of the fold's **training** subjects (n ≈ 500) | **yes**: the same individuals, re-processed | no | no | no |
| `std_disjoint` (new sensitivity) | target features of the **source-calibration** subjects (n ≈ 167) | **no** | no | no (labels of these subjects are used only for source q, not here) | no |
| `std_same_sub` (new, size-matched) | random training subset of size n(cal) | yes | no | no | no |
| Per-site std (`sitestd`) | per-site means/SDs of each pipeline from training subjects; sites with <15 → global | yes | no | no | **yes** (ridge refit + new q on site-standardised source) |
| LOSO per-site std | training subjects **plus the held-out site's own test subjects** | partly | **yes (transductive)** | no | yes |

## 2. Focused sensitivity analysis: does population matching inflate the std result?

The question: v1 estimated the target standardisation from the *same people* the model was trained on. That is a best case, because the target sample's age and site mix is identical. Does coverage hold when the statistics come from different people?

| | Per-pipeline std (v1 definition) | Size-matched, same people | **Disjoint people** |
|---|---|---|---|
| Mean coverage | 0.845 | 0.844 | **0.843** |
| Median coverage | 0.854 | 0.853 | **0.854** |
| Minimum | 0.742 | 0.741 | **0.743** |
| Transitions < 0.88 | 14 | 14 | **14** |
| Mean MAE (years) | 3.27 | 3.28 | **3.28** |
| Interval width | unchanged (frozen q, 11.8 y) | same | same |

Paired differences, mean over 18 [joint site CI]:
- Disjoint − v1 std: **−0.0015 [−0.0025, −0.0005]**, per-transition range −0.008 to +0.008.
- Size-matched − v1 std: −0.0006 [−0.0022, +0.0006]. The sample-size effect is nil.
- Disjoint − size-matched: −0.0009 [−0.0019, +0.0004]. The "different people" effect is not distinguishable from zero.

**Conclusion:**
- The v1 per-pipeline standardisation result is **not an artefact of using the same subjects**.
- Estimating target statistics from about 170 different, unlabelled subjects gives the same coverage to within about 0.2 points.
- The residual cross-family loss (§3 of INFERENCE_AUDIT) remains.

**Caveat:** "disjoint" subjects still come from the same ABIDE sites and age mix. This sensitivity analysis removes *individual* overlap, not distributional matching. A target sample with a different age range or site mix could behave differently. That is not testable within ABIDE without changing the design.

## 3. Taxonomy to use in the manuscript

| Class | Arms in this study | Needs at deployment |
|---|---|---|
| **A. Frozen model, inductive label-free feature adaptation** | per-pipeline standardisation (v1 and disjoint variants) | unlabelled target-pipeline features from a separate sample |
| **B. Transductive label-free adaptation** (uses the features of the very batch being predicted) | LOSO per-site standardisation (held-out site's own features); weighted CP (density ratio evaluated at test points; fitted on training-subject target features) | the test batch itself |
| **C. Adaptation requiring refitting** (label-free at the target, labelled source refit) | per-site standardisation with refit; mixed-pipeline training ("mix_loo": target pipeline unseen) | source labels + retraining; site labels |
| **D. Adaptation requiring labelled target data** | ReForm-style recalibration (k = 20–160); mixed training including the target pipeline ("mix_all") | k labelled target subjects, or a labelled target-pipeline training set |

**Wording rule:** per-site standardisation is reported as "site-wise standardisation with model refitting (class C)". It should never be called "standardising a frozen predictor". Its comparison baseline is its own refitted source condition (0.903), not the frozen source (0.904). It pays a 23% width cost: 14.5 vs 11.8 years under the textbook rank (v1 rank: 14.9 vs 12.1).
