# AGE_CONDITIONAL_VALIDITY.md — Phase 5

Code: `code_v2/phase5_age.py`. Results: `results_v2/age/`.
Textbook rank; 18 non-control transitions. CIs come from a joint site-cluster bootstrap over transitions. No age-adaptive method was run (see §5).

## 1. Who is in each age bin

| Bin | Unique subjects | Sites | Age range | Per-transition n (min–median–max) |
|---|---|---|---|---|
| < 12 | 294 | 13 | 6.5–12.0 | 164–225–288 |
| 12–15 | 264 | 16 | 12.0–15.0 | 188–236–264 |
| 15–18 | 176 | 15 | 15.0–17.9 | 134–167–176 |
| 18–25 | 189 | 11 | 18.0–25.0 | 154–177–188 |
| **≥ 25** | **54** | **9** | **25.0–35.0** | **44–49–54** |

- The "≥ 25" bin is **not an older-adult bin**. It contains 54 young adults aged 25–35, and 28 of them come from two sites (USM 14, MAX_MUN 14).
- The same 54 people recur in every transition, so the 18 per-transition estimates are **not independent replicates**.

## 2. Manuscript-ready subgroup table (mean over 18 transitions [joint site-cluster 95% CI])

| Age bin | Source | Raw | Per-pipeline std | Std − source | Site-std (refit) | Mean signed error at source (y) |
|---|---|---|---|---|---|---|
| < 12 | 0.907 [0.892, 0.918] | 0.424 | 0.813 [0.773, 0.842] | **−0.094 [−0.131, −0.064]** | 0.899 | +2.2 |
| 12–15 | 0.979 [0.965, 0.989] | 0.524 | 0.952 [0.936, 0.970] | −0.026 [−0.044, −0.009] | 0.995 | +1.1 |
| 15–18 | 0.988 [0.981, 0.994] | 0.546 | 0.970 [0.956, 0.983] | −0.018 [−0.027, −0.009] | 0.992 | +0.1 |
| 18–25 | 0.908 [0.872, 0.935] | 0.500 | 0.822 [0.779, 0.861] | **−0.085 [−0.113, −0.060]** | 0.867 | −2.4 |
| ≥ 25 (n ≈ 49) | **0.253 [0.160, 0.374]** | 0.282 | 0.166 [0.106, 0.246] | −0.087 [−0.137, −0.038] | 0.120 | **−7.5** |

- **Seed variation is small.** The SD across the 5 seeds of pooled bin coverage is ≤ 0.016 for every bin and condition (`age_bin_seed_variation.csv`).
- **Variation across transitions** is reported as min/max in `age_bin_coverage.csv`. For example, source coverage in the ≥ 25 bin ranges from 0.155 to 0.396 across transitions.

## 3. Interpretation

1. **The baseline is miscalibrated conditionally, not marginally.**
   - Marginal source coverage (0.904) averages over-coverage at 12–18 (0.98–0.99) with severe under-coverage of the young-adult tail (0.25).
   - The cause is regression to the mean in a cohort concentrated in late childhood and adolescence: ridge under-predicts the age of 25–35-year-olds by 7.5 years on average, and over-predicts the youngest by 2.2.
   - This is a property of the model and cohort, present before any pipeline change. It is consistent with the general fact that split-conformal guarantees are marginal only.
2. **The pipeline effect is not explained by the adult tail.**
   - Adults are 5.8% of subjects and contribute only −0.005 (median −0.004) of the −0.059 marginal std loss.
   - Restricted to subjects under 25, where source coverage is 0.944, the std − source loss is **−0.058 [−0.072, −0.044]**, essentially unchanged.
   - 8 of 18 transitions still lose ≥ 5 points under 25, against 9 with all ages.
   - The raw-transfer loss is also undiminished under 25 (median −0.243 vs −0.224).
3. **The pipeline effect is age-dependent in the opposite direction from what the baseline would suggest.**
   - Losses after standardisation are largest in the bins with no slack: < 12 (−0.094) and 18–25 (−0.085).
   - They are smallest where the source interval is over-wide (12–18: −0.018 to −0.026).
   - A marginal-coverage summary therefore understates the deterioration for children under 12 and for young adults.
4. **Adults (≥ 25) after transfer.**
   - Raw transfer changes adult coverage by +0.030 [−0.072, +0.128], i.e. not detectably. In some transitions the global offset moves predictions upward, which by chance helps these under-predicted subjects.
   - Standardisation lowers adult coverage further (−0.087 [−0.137, −0.038]).
   - With about 49 correlated subjects from 9 sites, these estimates are imprecise, and they should not be generalised to older adults.

## 4. Statements for the manuscript

- **Report** the subgroup table above with n per bin and site-cluster CIs. State that the oldest bin spans only 25–35 years (n = 54 unique subjects, 9 sites).
- **Replace** "Marginal coverage hid severe conditional failure: coverage for ages 25 and above was below 0.30 even at the source" with a precise, two-part statement:
  - (i) the source model under-covers young adults (0.25, n ≈ 49), a baseline limitation independent of pipeline;
  - (ii) the pipeline-induced loss persists, at the same size, when this tail is excluded.
- **Do not** claim the study shows conditional coverage behaviour for older adults, or for any group other than ABIDE's 6–35-year-olds.

## 5. Why no age-adaptive conformal arm was run

The existing analysis is interpretable without one: the decomposition in §3.2 already separates baseline age miscalibration from the pipeline effect. A Mondrian or normalised-residual arm would need a separate calibration quantile for a bin with about 10 calibration subjects per fold (about 49 adults × 20%). That is too few for a stable 90% conformal quantile; by the rank rule in Phase 2, at least 9 are needed for a finite threshold. The adult-bin results would therefore be dominated by noise. This is listed as future work, not as a gap the current evidence must fill.
