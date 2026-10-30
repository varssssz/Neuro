# NOVELTY_MATRIX.md — Phase 6

Sources: REFERENCES_VERIFIED.md. ✔ = done in that work; ✘ = not done; ◐ = partly.

| Dimension | Bhagwat 2021 | Korbmacher 2024 | Delzant 2025 | Sanz-Robinson 2026 (PP) | Li 2022 (PP) | Ernsting 2023 (PP) | Samanta 2026 (PP) | ReForm (Chen, PP) | Olsson 2022 | **This study** |
|---|---|---|---|---|---|---|---|---|---|---|
| Same scans processed by several pipelines | ✔ (same ABIDE, same 5 tools) | ✔ (FS 5.3 / 7.4.1) | ✔ | ✔ | ✔ | ✘ | ✘ | ✘ (FS 6.0.1 only) | ✘ (scanner) | ✔ |
| Age / brain-age prediction | ✔ (within pipeline) | ✔ | ✔ (29 traits) | ✔ | ◐ (classification) | ✔ | ✔ | ✘ (reference intervals) | ✘ | ✔ |
| Model trained on one pipeline, applied to another | ✘ | ✔ | ✘ | ◐ | ✔ | ✘ | ✘ | ✘ | ✔ (scanner/lab) | ✔ (all 20 directed pairs) |
| Prediction intervals / conformal | ✘ | ✘ | ✘ | ✘ | ✘ | ✔ | ✔ | ✔ | ✔ | ✔ |
| Coverage measured under a processing change | ✘ | ✘ | ✘ | ✘ | ✘ | ✘ | ✘ | ✘ (study/scanner shift) | ◐ (scanner) | **✔** |
| Label-free adaptation evaluated for coverage | ✘ | ✘ | ✘ | ✘ | ✘ | ✘ | ✘ | ◐ (harmonisation baselines) | ✘ | ✔ |
| Labelled recalibration vs k, with width | ✘ | ✘ | ✘ | ✘ | ✘ | ✘ | ✘ | ✔ | ✘ | ✔ (applied, not new) |
| Drift monitoring vs coverage loss | ✘ | ✘ | ✘ | ✘ | ✘ | ✘ | ✘ | ✘ | ◐ (KS on predictions) | ✔ (descriptive) |
| Age-conditional coverage | ✘ | ✘ | ✘ | ✘ | ✘ | ◐ | ✔ | ✘ (listed as future work) | ✘ | ✔ |

## What is new, what is combination, what is not ours

| Category | Item |
|---|---|
| **Genuinely new (to the extent the searches can show)** | A per-transition measurement of split-conformal coverage when a frozen brain-age model and its calibration are moved between processing pipelines, holding scans fixed, for all 20 directed transitions among five pipelines. The **software-family result**: label-free standardisation reduces loss to <1 point within FreeSurfer, but leaves 6.5–8.3 points across families. The descriptive observation that a standard MMD/KS feature-drift monitor rarely alarmed after standardisation (batches of 50) while coverage loss remained. |
| **New combination of existing methods** | Split conformal (Lei 2018) + same-scan multi-pipeline data (Bhagwat 2021) + label-free standardisation + weighted CP (Tibshirani 2019) + ReForm-style recalibration (Chen) + MMD (Gretton 2012), evaluated in one paired design |
| **Replication / extension** | Pipeline choice changes brain-age predictions (Korbmacher 2024; Bhagwat 2021; Delzant 2025): replicated here in MAE terms. Conformal age-bin under-coverage (Samanta 2026): replicated at the source pipeline in a different cohort. Conformal validity fails under acquisition/processing change (Olsson 2022): extended to neuroimaging pipelines. |
| **Expected by theory, not a finding** | Coverage restored by labelled recalibration (exchangeability). Marginal-but-not-conditional validity. Weighted CP producing infinite intervals when the density ratio is nearly separable. |
| **Engineering only** | Pipeline code, reproducibility package, tests |
| **Not supported by our experiments** | Any claim about older adults, dementia, other datasets, other feature types (volumes, area), deep-learning models, other conformal scores, other drift detectors, or clinical use |

## Novelty statement to use

> "To our knowledge, based on searches up to 30 October 2026, no published study or preprint has measured the coverage of conformal brain-age intervals when the same scans are processed by different structural MRI pipelines. Closely related work has measured pipeline effects on brain-age point predictions [Korbmacher 2024; Bhagwat 2021; Delzant 2025], conformal brain-age intervals within a single pipeline [Ernsting 2023; Samanta 2026], and conformal recalibration of MRI reference intervals across studies [Chen ReForm]. The absence of search results does not establish absence of such work."

Do **not** claim novelty for the general observation that conformal coverage fails under distribution shift (Tibshirani 2019; Barber 2023; Olsson 2022).

**Scoop risk:** the Chen group (ReForm, Samanta 2026) works on both conformal brain-age and pipeline/recalibration questions. A coverage-under-pipeline-change preprint from that group is plausible. Submit a preprint early, e.g. to bioRxiv, once the manuscript is final.
