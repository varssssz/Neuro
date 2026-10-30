# RECALIBRATION_CORRECTION.md — Phase 2

Code: `code_v2/swap_engine_v2.py` (`run_recal`, `folds_v2`, `conformal_quantile`), `code_v2/run_recalibration_v2.py`, `code_v2/summarise_recalibration.py`.
Tests: `code_v2/tests/test_v2.py` (44 passing).
Results: `results_v2/recalibration/`.
All aggregates are over the 18 non-control transitions unless stated otherwise.

## 1. What the v1 code did (VERIFIED-CODE, then VERIFIED-RERUN)

| Question | Finding |
|---|---|
| How were the k labelled target subjects chosen? | Uniformly at random (100 draws per fold and seed) from `pool = np.r_[tr, ca]`, i.e. training plus source-calibration subjects, using their **target-pipeline** features and true ages. |
| Did training residuals contaminate the pool? | Yes. 74.9% of the pool were model-training subjects. Their target-pipeline residuals are partly in-sample, because the model was fitted on the same subjects' source features. |
| Did source-calibration subjects influence fitting or λ? | No. RidgeCV selects λ by leave-one-out inside `tr` only, and `ca` is used only for q. |
| Test-label access? | None. The pool never intersects the test fold (verified for all 650 folds). |
| Quantile formula implemented | `np.quantile(r, min(1, ⌈(k+1)(1−α/2)⌉/k), method="higher")`. This returns the (r+1)-th order statistic, where r = ⌈(k+1)(1−α/2)⌉, **capped at the maximum**. The textbook rule is the r-th order statistic, or **+∞ when r > k**. |
| Why does coverage follow theory so exactly, for every transition? | Recalibration on target-pipeline residuals is itself a split-conformal procedure. Under exchangeability its expected coverage is fixed by the rank alone: 1 − 2(k+1−r)/(k+1), whatever the shift between pipelines. The v1 values (0.813, 0.900, 0.948, 0.922, 0.909) equal the theoretical values of the v1 rank rule (0.818, 0.905, 0.951, 0.926, 0.913) up to the small contamination bias. The non-monotone curve and the "conservative k = 40" are rank artefacts, not empirical phenomena. |
| k = 10 | ⌈11 × 0.95⌉ = 11 > 10, so the textbook threshold is **infinite**: no finite interval can carry the guarantee. v1 silently used the maximum residual and reported 0.813 coverage. That number is below nominal and carries no guarantee. |

The design D0 (the v1 code path re-executed) reproduces `original_v1/Experimental_Results/recal_main.csv` exactly: max |diff| = 0 for coverage and width.

## 2. Corrected design

Four disjoint partitions inside each of v1's outer folds. The test folds are **identical to v1**, so the same subjects are tested.

| Partition | Share of all subjects | Mean n (18 transitions) | Used for |
|---|---|---|---|
| Model training | 40% | 335 | ridge fit + λ (source features) |
| Source calibration | 16% | 134 | frozen q (source features) |
| Target recalibration pool | 24% | 201 (min 166) | k labelled target subjects drawn from here |
| Final test | 20% | 167 | scoring only |

The model stays frozen during recalibration. The conformal ranks are textbook ranks for both q and the recalibration corrections. Each (transition, k) cell has 2,500 draws (5 seeds × 5 folds × 100).

Repartitioning costs little accuracy:
- source MAE 2.774 vs 2.741 years;
- source coverage 0.902 vs 0.904 under the same rule;
- source width 11.96 vs 11.82 years.

## 3. Original versus corrected (`results_v2/recalibration/original_vs_corrected.csv`)

Effects are attributed step by step:
1. **D0 → D1:** quantile rule only.
2. **D1 → D3:** remove training subjects from the pool. Pool = source-calibration subjects, so the training size is unchanged.
3. **D3 → D5:** move to the full four-way split.

| k | v1 coverage | v1 width | **Corrected coverage** (range over 18) | Theory | **Corrected width** mean / median (y) | Effect of rule | Effect of pool | Effect of repartition |
|---|---|---|---|---|---|---|---|---|
| 10 | 0.813 | 13.6 | **1.000 (∞ interval, all draws)** | 1.000 | ∞ | +0.187 (→ ∞) | 0 | 0 |
| 20 | 0.900 | 16.6 | **0.904** (0.899–0.909) | 0.905 | 16.8 / 15.8 | 0.000 | +0.005 | −0.001 |
| 40 | 0.948 | 19.2 | **0.903** (0.896–0.907) | 0.902 | 15.8 / 15.2 | −0.050 | +0.004 | +0.001 |
| 80 | 0.922 | 16.4 | **0.902** (0.895–0.906) | 0.901 | 15.3 / 14.8 | −0.025 | +0.004 | +0.001 |
| 160 | 0.909 | 15.4 | **0.901** (0.893–0.906) | 0.901 | 15.1 / 14.7 | −0.013 | +0.006 | n/a* |

\*The pool = calibration-only bridge (D3) can draw k = 160 in only 10 of 18 transitions, so the repartition step is not separable at k = 160.

**Readings:**
- **Almost all of the v1 behaviour is the quantile rule.** The conservative 0.948 at k = 40 and the 19.2-year width disappear: textbook coverage is 0.903 at 15.8 years.
- **Training-subject contamination was real but small.** Removing it raises coverage by 0.004–0.006, bringing it onto the theoretical value. The direction matches expectation: in-sample residuals are slightly too small, so the v1-pool intervals slightly under-covered relative to their own theory.
- **Repartitioning has a negligible effect** (≤0.0014).
- **Standardise-then-recalibrate (hybrid, D5):** coverage is the same (0.904/0.903/0.902/0.901), with narrower intervals of 15.5/14.6/14.2/14.0 years at k = 20/40/80/160.
- **Draw-to-draw variability:** the share of individual draws with coverage below 0.88 is 0.24–0.30. That is about the binomial noise floor for test folds of about 167 subjects, so it is not evidence of invalidity.

### Width cost by transition (corrected, k = 40, `reform_raw`)

The ratio to the source width (about 12 years) has median 1.25 and range 0.99–1.82.
- The largest costs are in FreeSurfer → ANTs (1.51–1.82) and → CIVET (1.43–1.53).
- With standardisation first, the median ratio is 1.23 and the maximum 1.53. The hybrid mainly helps the ANTs-target transitions (e.g. FS60→ANTS 1.82 → 1.29).

## 4. How the manuscript should present recalibration

- **Not a finding:** "labelled recalibration restores coverage". Under exchangeability this holds by construction, and the coverage values above are what theory predicts.
- **The informative quantity is the deployment trade-off:**
  - at least 20 labelled target subjects are needed for any finite 90% interval with separate α/2 corrections;
  - with 20–160 subjects the intervals are typically 1.2–1.3× (up to 1.8×) the source width;
  - label-free standardisation before recalibration reduces that cost.
- **Delete or replace these v1 statements:**
  - "conservative coverage at k = 40, which we attribute to rounding…";
  - "about 40 labelled subjects gave conservative coverage at roughly 1.6 times the source width";
  - "its 40-subject rule is borne out".

  Replace them with the corrected table and its explanation. ReForm's own guidance (about 40 for 90% intervals) is about estimating quantiles stably, not about the minimum k for a finite interval. Do not present our results as confirming or refuting it.
- **Methods text:** describe the four-way split and the textbook rank, and state that test subjects are never in the pool.
