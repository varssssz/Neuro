# INFERENCE_AUDIT.md — Phase 3: controls and statistical inference

Code: `code_v2/phase3_inference.py`.
Inputs: `results_v2/frozen/per_subject_seed.csv.gz`, from `code_v2/run_frozen_v2.py`: v1 design, per-seed and per-subject outputs, under both the v1 rank and the textbook rank.
Outputs: `results_v2/inference/`.

**Reproduction:** under the v1 rank, the new run reproduces `per_subject_main.csv` exactly (16,718 rows; max |diff| = 0 for src/raw/std/sitestd).
**Primary results now use the textbook rank**, which matches the equation in the paper. v1-rank values are kept alongside.

## 1. Problems with the v1 inference

| v1 practice | Problem | v2 replacement |
|---|---|---|
| Exact McNemar on **seed 0 only** | Discards seeds 1–4. The result depends on the seed: of the 18 non-control std-vs-source tests, only 11 are p < 0.05 on all five seeds, and 4 that are significant on seed 0 are not significant on every seed. McNemar also assumes independent subjects. | Cluster-robust test on seed-averaged paired indicators |
| Subject-level bootstrap | Ignores site clustering. Subjects from one acquisition site share scanner and protocol effects. | Site-cluster bootstrap (B = 2,000); subject bootstrap kept for comparison |
| Bonferroni count "11 of 18" stated in the text only | No adjusted p in any file | Holm-adjusted (and Bonferroni) p-values per transition, per pre-defined family |
| Mean over transitions as the headline | Treats transitions as independent, although the same subjects recur in every transition; CIVET collapses dominate the mean | Mean **and median** with a joint site bootstrap (the same resampled sites applied to all transitions); min/max; counts |
| "Near-identity control" justified as the "closest software versions" | Not supported (§3) | Two controls with different meanings, reported separately |

## 2. Inference design (assumptions stated)

- **Unit:** subject. Each subject's covered indicator is averaged over the 5 cross-fitting seeds, which integrates out split randomness.
- **Paired contrast** per transition: d = coverage(condition) − coverage(reference), on the same test subjects.
- **Dependence:** subjects are clustered by acquisition site (18–20 sites per transition), and clusters are assumed independent.
  - Fold-level dependence (subjects in one fold share a fitted model and q) is not modelled separately; averaging over 5 seeds × 5 folds dilutes it.
  - Inference is therefore **conditional on the fitting procedure** and is not a guarantee about new datasets.
- **Per-transition test:** CR1 cluster-robust t-test of mean d = 0, df = G − 1. 95% CI by site-cluster bootstrap.
- **Families for multiplicity:** one family per contrast, each of the 18 non-control transitions, with Holm correction (Bonferroni also given). The families, all confirmatory and fixed by the original design, are:
  - (a) raw vs source;
  - (b) per-pipeline std vs source;
  - (c) per-site std vs its own refitted source.

  Controls, the atlas arm, LOSO, mixed-pipeline training, LightGBM, monitors and strata are **exploratory**: no multiplicity correction, CIs only.

## 3. Controls (`results_v2/inference/controls_table.csv`, textbook rank)

| Pair | Differs in | Same-region r (mean / min) | Offset (mm) | SD ratio | Raw − source [site-cluster 95% CI] | Std − source [CI] |
|---|---|---|---|---|---|---|
| FS60-DK → FS60 | atlas only (same FS 6.0 run) | 0.974 / 0.780 | 0.000 | 0.98 | +0.003 [−0.004, +0.008] | +0.001 [−0.005, +0.005] |
| FS60 → FS60-DK | atlas only | 0.974 / 0.780 | 0.000 | 1.02 | −0.003 [−0.010, +0.006] | −0.002 [−0.006, +0.004] |
| FS5.3 → FS6.0 | version **and** atlas (DK→DKT) | 0.884 / 0.673 | +0.031 | 0.82 | −0.004 [−0.029, +0.012] | +0.001 [−0.023, +0.018] |
| FS6.0 → FS5.3 | version and atlas | 0.884 / 0.673 | −0.031 | 1.22 | −0.009 [−0.025, +0.014] | −0.004 [−0.019, +0.015] |
| FS5.1 → FS5.3 *(comparator)* | version only (both DK) | 0.867 / 0.760 | −0.124 | 0.84 | −0.104 [−0.152, −0.061] | −0.022 [−0.049, −0.002] |
| FS5.3 → FS5.1 *(comparator)* | version only | 0.867 / 0.760 | +0.124 | 1.18 | −0.039 [−0.065, −0.005] | +0.011 [−0.005, +0.031] |

Conclusions:
- **FS60-DK↔FS60 is the appropriate negative control.** It is the same software run on the same scans, with only the parcellation changed, and it shows no detectable loss (CIs within about ±0.01).
- **FS5.3↔FS6.0 is not a "near-identity" pair by construction.** It differs in software version *and* atlas, and its feature agreement (r = 0.884) is about the same as FS5.1↔5.3 (0.867), which loses up to 10 points.
  - Its benign behaviour is an *empirical observation*: the shift is small in offset (0.03 mm) and evidently harmless here.
  - It should be reported as an "observed low-impact transition", not as a negative control chosen a priori.
  - The 18 vs 20 split can be kept as a reporting convention, but it must be described this way.
- The atlas arm also shows that FS60-DK → FS5.3 (same atlas, version 6.0 vs 5.3) loses nothing (−0.006 [−0.013, +0.003]). FS60-DK → FS5.1 loses −0.112 raw, which standardisation repairs (−0.010). This replicates the v1 conclusion that the FS5.1 effect is a version effect, not an atlas effect.

## 4. Aggregate results, 18 non-control transitions (`aggregate_18_noncontrol.csv`)

| Quantity | v1 rank (as published) | **Textbook rank (v2)** |
|---|---|---|
| Source coverage, mean (range) | 0.910 (0.901–0.919) | **0.904** (0.895–0.915) |
| Raw coverage, mean / median (range) | 0.489 / 0.690 (0.001–0.876) | **0.482 / 0.678** (0.001–0.871) |
| Std coverage, mean / median (range) | 0.853 / 0.862 (0.758–0.926) | **0.845 / 0.854** (0.742–0.921) |
| Std transitions < 0.88 | 12 / 18 | **14 / 18** |
| Std − source, mean [joint site CI] | −0.057 [−0.068, −0.044] | **−0.059 [−0.072, −0.046]** |
| Std − source, median [joint site CI] | −0.048 [−0.068, −0.035] | **−0.050 [−0.072, −0.036]** |
| Std transitions with loss ≥ 5 points | 9 | **9** |
| Std: Holm-significant (cluster-robust) | 11 | **10** |
| Raw − source, median [CI] | −0.218 [−0.288, −0.173] | **−0.224 [−0.292, −0.178]** |
| Raw: Holm-significant | 18 | **18** |
| Site-std coverage, mean (range) | 0.897 (0.875–0.923) | **0.890** (0.867–0.919) |
| Site-std − refitted source, mean [CI] | −0.012 [−0.018, −0.007] | **−0.013 [−0.019, −0.007]**; 2 Holm-significant |

**Which results change:** absolute coverages fall by about 0.006–0.008. This is purely the rank correction: v1 intervals were one rank wide. The std count below 0.88 rises from 12 to 14, because two transitions sat just above the threshold in v1 and fall just below it: FS51→FS53 (0.884 → 0.879) and CIVET→FS60 (0.881 → 0.871). *Thresholds such as 0.88 are fragile and should not carry the argument; report the paired differences and their CIs.*

**Which results do not change:** the size and sign of every paired difference, the Holm-significant set (10 vs 11), and every qualitative conclusion.

## 5. Software-family sensitivity (`category_sensitivity.csv`, textbook; joint site CI)

| Category (transitions) | Raw − source | Std − source |
|---|---|---|
| Atlas only, FS60-DK↔FS60 (2) | −0.000 [−0.003, +0.003] | −0.000 [−0.003, +0.002] |
| FS version 6.0 vs 5.3, same DK atlas (2) | −0.006 [−0.013, +0.003] | −0.002 [−0.007, +0.004] |
| FS 5.3↔6.0, designated control (2) | −0.007 [−0.015, +0.002] | −0.001 [−0.009, +0.006] |
| FS version 6.0 vs 5.1, same DK atlas (2) | −0.112 [−0.160, −0.067] | −0.010 [−0.016, −0.003] |
| Within FreeSurfer, involving 5.1 (4) | −0.095 [−0.123, −0.065] | −0.007 [−0.011, −0.004] |
| FreeSurfer↔ANTs (6) | −0.194 [−0.245, −0.151] | **−0.069 [−0.092, −0.050]** |
| FreeSurfer↔CIVET (6) | −0.847 [−0.904, −0.780] | **−0.083 [−0.110, −0.055]** |
| ANTs↔CIVET (2) | −0.487 [−0.506, −0.466] | **−0.065 [−0.084, −0.044]** |

**The clearest single result of the study.** Within FreeSurfer (version changes, atlas changes), per-pipeline standardisation reduces the loss to under 1 point. Across software families (ANTs, CIVET), a 6.5–8.3 point loss remains after standardisation. The cross-family pairs are also the ones whose same-name regions agree poorly (mean r 0.44–0.50; Phase 1).

Direction matters too. The std loss for CIVET → FreeSurfer targets is −0.02 to −0.03, but for FreeSurfer → CIVET it is −0.12 to −0.17. So "pipeline difference" is not symmetric, and per-transition reporting is needed.

## 6. Per-transition table

`results_v2/inference/per_transition_table_textbook.csv` gives, for all 20 transitions:
- N and number of sites;
- coverage per condition;
- each contrast with its site-cluster CI and Holm/Bonferroni-adjusted p;
- seed SD;
- the seed-0 McNemar p and the number of seeds with p < 0.05.

The per-transition seed SD of standardised coverage is 0.002–0.024.
