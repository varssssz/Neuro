# DRIFT_MONITOR_VALIDATION.md — Phase 4b

Code: the original v1 `swap_engine.run_monitor` (unmodified), re-executed by `code_v2/phase4_std_monitor.py`.
Results: `results_v2/standardisation_monitor/monitor_diagnostic_by_transition.csv`, `monitor_summary.csv`, `fig_monitor_alarm_vs_coverage_loss.png`.
Coverage losses come from Phase 3 (textbook rank).

## 1. Verification

Re-running the v1 monitor code for all 20 transitions at batch size 50 reproduces the stored `monitors_b50.csv` **exactly** (max |diff| = 0 across MMD/KS × raw/std detection rates).

## 2. What the monitor actually is (VERIFIED-CODE)

| Element | Implementation | Independent of the evaluation data? |
|---|---|---|
| Purpose | Detect **feature-distribution shift** in a batch of target subjects. It was not designed or calibrated as a coverage-loss estimator. | — |
| Statistics | (i) unbiased RBF-kernel MMD² between the batch and a reference set, in 62-feature space; median-heuristic bandwidth from 150 reference subjects. (ii) two-sample KS statistic between the batch's predicted ages and the reference subjects' predicted ages. | — |
| Reference set | 200 random **training** subjects, source features, source scaling. The KS reference uses their in-sample fitted ages. | yes (training subjects only) |
| Null distribution / threshold | 95th percentile of the statistic over 40 batches per fold (200 in total) of **source-calibration** subjects' source features | **yes**: calibration subjects never enter the test assessment |
| False-alarm check | batches of **test** subjects' *source* features | — |
| Alternative (what is being monitored) | batches of **the same test subjects** whose coverage is scored, using target-pipeline features (raw, or standardised with target statistics from training subjects) | the same subjects as the coverage assessment, which is appropriate: the monitor watches the batch whose intervals would be affected. No labels are used. |
| Batch size / sampling | 50 subjects drawn **without replacement within a batch** from a test fold of 138–194; 40 batches per fold, which **overlap heavily**. Alarm rates are therefore proportions of correlated batches, not independent trials. | — |

**Observed false-alarm rate** on source test batches (nominal 0.05):
- mean 0.067 (MMD raw), 0.071 (MMD std), 0.072 (KS std);
- per-transition range 0.01–0.16.

The thresholds come from only 200 correlated batches, so per-transition alarm rates below about 0.15 cannot be distinguished from false alarms.

## 3. Results (non-control transitions, batch 50)

- **Raw features:**
  - The MMD alarms on **every** batch of every transition (rate 1.00), including the two harmless FS5.3↔6.0 transitions (rate 1.00).
  - The KS test on predicted ages is more selective for those (0.00 and 0.045).
- **After per-pipeline standardisation:**
  - MMD alarm rate mean 0.139, median 0.103, maximum 0.32. KS mean 0.092.
  - Coverage loss remains: **15 of 18** transitions have a site-cluster CI for (std − source) entirely below 0, and 9 lose 5 points or more.
  - **None of those 15** reaches a majority-alarm rate (≥ 0.5) on either MMD or KS. Several sit at the false-alarm level: FS60→CIVET loses 11.7 points with an MMD alarm rate of 0.08, and ANTS→CIVET loses 8.5 points with 0.085.
  - Across the 18 transitions, the MMD alarm rate correlates weakly with the size of the loss (Spearman ρ = −0.49, p = 0.04; 18 correlated points). This is descriptive only.
- **Batch size:** the low sensitivity is about **power**, not complete blindness.
  - MMD std alarm rates are 0.07 / 0.14 / 0.22 (maximum 0.18 / 0.32 / 0.72) at batch sizes 25 / 50 / 100 (stored v1 files).
  - Batches of 100 are drawn from test folds of about 167, so they overlap even more.

## 4. Interpretation the manuscript may make (survives verification)

- In this design, standardisation removes the first and second moments of each feature. That is what a full-space kernel MMD with a median-heuristic bandwidth mainly responds to at n = 50.
- What remains is a change in the *joint* structure of the features: for ANTs/CIVET, same-name regions correlate weakly with their FreeSurfer counterparts (r ≈ 0.44–0.50; Phase 1). That still mis-calibrates a frozen linear model, but it produces only small MMD shifts.
- **Supported claim:** "After standardisation, a standard label-free feature-drift monitor (MMD, KS on predicted ages) rarely alarmed in batches of 50 subjects, although 15 of 18 transitions still lost coverage. Alarm rates rose with batch size."
- **Not supported:** "monitors cannot detect harmful shift", or any claim about other detectors, other bandwidths, or larger samples. Only one detector configuration was tested, and the analysis is descriptive at the level of 18 correlated transitions.

No additional detector was implemented, because the existing analysis is valid as specified.

## 5. Transition-level diagnostic table

The full table is `monitor_diagnostic_by_transition.csv`. For each transition it gives:
- coverage loss and its CI;
- thresholds and false-alarm rates;
- MMD/KS alarm rates for the raw and standardised conditions;
- harmful/missed flags.

The figure `fig_monitor_alarm_vs_coverage_loss.png` plots alarm rate against coverage change for raw (left) and standardised (right) features, with the mean false-alarm level marked.
