# Complete Novelty and Literature Audit — raw material, not conclusions

This reproduces the actual search logs and paper-level findings from the two prior research-agent passes run in this project (task IDs `wf-0a70cf30-cb62-59b6-b73b-75e0b7647571`, dated 25 Oct 2026, and `wf-d7dbc55e-7ca1-5cd4-a7cb-86ed82d43b7d`, dated 25 Oct 2026, a follow-up/extension of the first). Nothing below has been re-verified in this handover pass; it is copied forward, with tags preserved. **Total searches run: 23 (pass 1) + 18 successful + 1 failed (pass 2, numbered #21–#39) + 6 subagent full-text fetches (pass 2, F1–F8, 2 failed) + ~10 subagent searches (S1–S10, pass 1 and 2 combined) ≈ 41 search/fetch actions**, matching the "approximately 41 searches" figure. PubMed, Scopus and Web of Science were **never directly queried** — both passes ran web-search equivalents instead and say so explicitly; this is disclosed, not hidden.

---

## Pass 1 search log (23 entries, verbatim query text where recorded)

| # | Query (abridged, as logged) | Tool | Results used |
|---|---|---|---|
| 1 | Bhagwat preprocessing pipelines cortical surface GigaScience | web search | PMC7821710, GigaDB 100832, Nat Commun 2024 commentary |
| 2 | Delzant Colliot Couvy-Duchesne pipelines HBM | web search | Wiley hbm.70372, PMC12572822, UQ profile |
| 3 | Samanta conformal brain age bioRxiv 2026 | web search | Pitt profile, GitHub repo, BayesAge conformal (adjacent) |
| 4 | Sanz-Robinson segmentation variability brain age | web search | bioRxiv 2026.07.23.740400, Epics lab page |
| 5 | "Same brain, different prediction" EEG | web search | arXiv 2605.07212 abs/HTML |
| 6 | Korbmacher FreeSurfer version-shuffling | web search | ScienceDirect S2666956024000205, PMC12172796 |
| 7 | Normative modelling deviation FreeSurfer/FastSurfer centiles | web search | Sci Rep 2025 reference curves; FastSurfer paediatric database; CentileBrain |
| 8 | Harmful shift / label-free performance estimation medical imaging | web search | arXiv 2506.17442, 2409.17800, 2412.12910, npj DM 2024, Nat Commun 2023 |
| 9 | Numerical variability MCA FreeSurfer brain age Parkinson | web search | bioRxiv 2026.01.09.698203; arXiv 2509.05238; FS-versions PD preprint |
| 10 | (fetch) Sanz-Robinson landing page | fetch | Abstract |
| 11 | Samanta Shafiei Satterthwaite conformal coverage | web search | ResearchGate abstract (RBC, ages 6-22) |
| 12 | Brain age software packages test-retest Dörfel | web search | HBM hbm.26502; motion-reliability paper |
| 13 | Normative transfer new site FreeSurfer version PCNtoolkit | web search | arXiv 2509.07237; PCNtoolkit protocol update |
| 14 | Conformal covariate shift neuroimaging brain age | web search | Tibshirani et al. PDF; Ernsting supplement |
| 15 | compare-surf-tools ABIDE FS ANTs CIVET | web search | companat and neurodatascience GitHub; RIO BrainHack report |
| 16 | DL segmentation seed variability FastSurfer SynthSeg | web search | arXiv 2509.05238; arXiv 2504.15931 |
| 17 | Sulcal shape pipeline reproducibility FreeSurfer versions | web search | calcSulc (Madan); HBM FS 7.1/6.0/5.3 compatibility |
| 18 | Brain charts centiles FreeSurfer vs SynthSeg | web search | medRxiv 2023 clinical-controls charts; Brain MoNoCle |
| 19 | Drift detection MRI software update monitoring | web search | Nat Commun 2023 (simulated software updates); MMC+ |
| 20 | Conformal normative centiles neuroimaging | web search | **FAILED: search budget exhausted** |
| S1-S4 | Subagent: Samanta full text; "conformal FreeSurfer version"; "conformal coverage preprocessing pipeline MRI"; "normative calibration software version" | subagent | Chen et al. 2025 ReForm (full text); nothing direct on pipeline-shift conformal |

**Pass 1 self-reported shortfall:** 19 direct + ~4 subagent ≈ 23 searches, against a 40-search target. PubMed/Scopus/WoS Boolean strings were written but **not executed natively** in those databases.

## Pass 2 search log (continues numbering from #20; 18 succeeded, 1 failed, plus 6 full-text fetches)

| # | Query / action | Tool | Outcome |
|---|---|---|---|
| F1 | arxiv.org/pdf/2302.05304 (Ernsting) | fetch | **Success — full text + Supplementary Material** |
| F2 | nature.com/articles/s41467-022-34945-8 (Olsson) | fetch | Success main text; SI not fetchable (anchor-only link) |
| F3 | ncbi.nlm.nih.gov/.../PMC12572822 (Delzant) | fetch | Success, truncated at §4.4; SI not read |
| F4 | pmc.ncbi.nlm.nih.gov/.../PMC12172796 (Korbmacher) | fetch | Success, truncated before limitations paragraph |
| F5 | arxiv.org/html/2605.07212v1 (Hou) | fetch | Success, appendix partially truncated |
| F6 | arxiv.org/pdf/2412.12910v1 (Amoukou) | fetch | Success through §5.2; appendices not read |
| 21 | conformal prediction brain age FreeSurfer version coverage | search | Korbmacher, HCP conformal chapter, Ernsting; no pipeline-swap conformal paper |
| 22 | ReForm reference interval calibration MRI preprocessing pipeline | search | ReForm bioRxiv v1/v2, PMC12697551, PubMed 41394588/42444711, journal version |
| F7 | direct.mit.edu ReForm journal page | fetch | **FAILED (bot detection)** |
| 23 | Satterthwaite Chen conformal brain age bioRxiv 2026 | search | Samanta preprint page; no newer preprint from the group |
| F8 | pmc.ncbi.nlm.nih.gov/.../PMC12697551 (ReForm) | fetch | **FAILED (reCAPTCHA)** — also re-attempted and failed again later in this project |
| 24 | github andy1764 ReForm reference interval conformal | search | ReForm README (40/80/120 calibration-size rule), journal snippets |
| 25 | Antonopoulos VBM pipelines NeuroImage 2023 120292 | search | manuscript ref [3] verified |
| 26 | Nevin Lees Groth mismatched cleaning pipelines PeerJ | search | manuscript ref [12] verified |
| 27 | Lei G'Sell Rinaldo Tibshirani Wasserman JASA 113 523 | search | manuscript ref [7] verified, DOI found |
| 28 | Di Martino ABIDE Mol Psychiatry 2014 | search | manuscript ref [13] verified |
| 29 | Klein Tourville 101 labeled brain images | search | manuscript ref [14] verified |
| 30 | weighted conformal covariate shift neuroimaging site scanner brain age | search | methods papers only; no MRI-pipeline weighted-CP application found |
| 31 | Nipoppy derivatives FreeSurfer versions PPMI QPN | search | Nipoppy bioRxiv 2026; FS 6.0.1/7.3.2 configs |
| 32 | traveling subjects structural MRI FreeSurfer reproducibility older adults | search | OASIS-1, SRPBS traveling subjects, FastSurfer Traveling Human Phantom |
| 33 | OASIS-3 FreeSurfer 5.3 and 7 derivatives dementia | search | OASIS-3 FS5.3 derivatives (n=1038); DL+DiReCT vs FS vs ANTs on OASIS-3 |
| 34 | compare-surf-tools ABIDE csv | search | both GitHub repos, file names, GigaDB MIT licence |
| 35 | conformal coverage diagnosis subgroup autism brain age (G9 falsification) | search | AD subgroup under-coverage preprint (arXiv 2608.04254); nothing on ASD × pipeline |
| 36 | QC exclusion choice changes brain age FreeSurfer (G10 falsification) | search | Qoala-T; pediatric QC comparison; none on coverage specifically |
| 37 | conformal coverage preprocessing pipeline shift MRI segmentation (G1/G5 falsification) | search | COMPASS; "Conformalized LLMs under Configuration Shift"; pitfalls-of-CP paper; no MRI-pipeline conformal study |
| 38 | normative model deviation FreeSurfer version pretrained transfer | search | Alyas practical guide; no version-swap centile study |
| 39 | label-free performance estimation regression shift medical imaging | search | **FAILED: search budget exhausted** |
| S5-S10 | Subagent: ReForm/Samanta citing papers; Chen/Satterthwaite preprints Apr-Oct 2026; andy1764 repos; weighted CP for MRI regression; FastSurfer/SynthSeg centile coverage | subagent | No scoop found. Samanta bioRxiv page fetch failed (URL not permitted). Semantic Scholar citation API failed (HTTP 429), repeatedly. |

**Pass 2 self-reported shortfall:** 18 of 20 additional searches requested succeeded (+1 failed); forward-citation chasing (Google Scholar/Semantic Scholar "cited by") **failed for every one of the 7 seed papers** due to a 429 rate-limit — this was never resolved, so scoop risk from forward citations is genuinely unverified, not just "not found."

---

## The six closest competing papers — full detail

### 1. Chen et al. (ReForm)
- **Title:** "Calibration of MRI-based reference intervals to new samples."
- **Authors (as verified):** Chen AA, Seidlitz J, Gardner M, Bethlehem RAI, Dorfschmidt L, Kafadar E, Lifespan Brain Chart Consortium, Benitez A, Jensen JH, Vandekar S, Satterthwaite TD, Alexander-Bloch AF.
- **Venue/year/DOI:** *Imaging Neuroscience* 4:IMAG.a.1290, 2026, doi `10.1162/IMAG.a.1290` (peer-reviewed; PubMed record 10 Jul 2026). Preprint: bioRxiv `10.1101/2025.11.25.690453` (v1, v2).
- **Evidence tag:** preprint [FULL-TEXT, read directly in this project — see the project's earlier full-text fetch of the bioRxiv PDF]; journal version [ABSTRACT-ONLY] — the Imaging Neuroscience / MIT Press and PMC copies both failed to open (bot block, reCAPTCHA respectively), confirmed on **two separate attempts in this project.**
- **What it investigated:** recalibrates normative GAMLSS/quantile-regression reference intervals (CQR-style, separate lower/upper constants from calibration residuals) to a new sample. Experiments varied **study and field strength**, not processing pipeline.
- **Quote (≤15 words, verified from the full-text fetch):** "All data analyzed in this work used FreeSurfer version 6.0.1."
- **What it did NOT investigate:** no pipeline swap; no FreeSurfer-version comparison; no frozen-model cross-pipeline transfer.
- **Calibration-size rule:** README states ~40/80/120 calibration subjects for 90%/95%/99% intervals; the journal-page snippet (not full text) suggested ≥120 needed and that 160 gave no significant difference for 90%/95% — **this conflict between README and journal snippet was never resolved** because the journal full text could not be opened.
- **Difference from this study:** ReForm's guarantee holds whenever calibration and new observation are exchangeable; pipeline-matched calibration restores coverage *by construction* under its own theory. This study measures what happens *without* that pipeline-matched calibration, and tests the "ReForm-style" fix (an approximation — the real method's exact algorithm beyond what the preprint describes was not independently confirmed) as one baseline among several, under a pipeline swap specifically.

### 2. Samanta et al. 2026
- **Title/authors:** not independently re-verified beyond "Samanta U, Shafiei G, Karim H, Satterthwaite TD, Chen A. 'Quantifying Uncertainty in Brain Age Predictions via Conformal Prediction.'" — this title is as given by the research agent, **tagged [UNVERIFIED]** in the manuscript's reference list because it was not re-checked against the actual bioRxiv page in the final verification pass (a direct fetch of the Samanta bioRxiv page failed — "URL not permitted" — in both the subagent attempt and a prior session attempt).
- **Venue:** bioRxiv, doi `10.64898/2026.09.04.749442`, preprint, not peer-reviewed.
- **What it investigated (per the project's earlier full-text read, not re-verified here):** single dataset (RBC), single pipeline stack (FreeSurfer 6.0.1 + sMRIPrep), split/CV+/CQR conformal methods, ages 6-22. Marginal coverage ~91%, age-bin coverage as low as 79.7% (18-20) and 48.3% (20-22).
- **What it did NOT investigate:** no pipeline swap of any kind — single pipeline throughout.
- **Difference from this study:** confirms (within one pipeline) that marginal coverage can hide severe age-conditional failure — this study finds the same pattern (age ≥25 coverage 0.277 at the **source** pipeline, before any pipeline swap) and extends it across pipelines.

### 3. Korbmacher et al. 2024
- **Title:** "FreeSurfer version-shuffling can enhance brain age predictions."
- **Venue/DOI:** *NeuroImage: Reports* 4(3):100214, 2024 (published 16 Jul 2024), doi `10.1016/j.ynirp.2024.100214`. [FULL-TEXT main read via PMC12172796, truncated before the limitations paragraph; supplement not read.]
- **What it investigated:** UK Biobank, N=4395, ages 49.33-82.84, **same T1 scans** processed with FS 5.3.0 and FS 7.4.1 (Desikan-Killiany). Cross-version brain-age differences up to 1.42 years; version-shuffled training improved cross-version prediction; cross-version predictions changed the sign of a reported sex difference.
- **What it did NOT investigate:** no prediction intervals or coverage at all — point predictions (MAE, R²) only.
- **Difference from this study:** this is the closest "same scans, multiple pipelines, frozen/transferred model" precedent, but entirely without conformal intervals — this study adds the interval/coverage dimension that Korbmacher's design structurally lacks.

### 4. Sanz-Robinson et al. 2026
- **Venue:** bioRxiv `10.64898/2026.07.23.740400`, preprint.
- **What it investigated (per earlier full-text read):** 4 segmentation pipelines, trains per-pipeline or concatenated models with nested group CV; no frozen-model transfer, no intervals; inter-pipeline disagreement depends on age and sex; uses NeuroCI/Nipoppy tooling.
- **What it did NOT investigate:** no conformal/interval framework at all.

### 5. Ernsting et al. 2023
- **Title:** "From Group-Differences to Single-Subject Probability: Conformal Prediction-based Uncertainty Estimation for Brain-Age Modeling."
- **Venue:** arXiv:2302.05304, 2023, preprint. [FULL-TEXT incl. Supplementary Material, read.]
- **What it investigated:** NAKO training (N=10,691), coverage (PICP) evaluated on **independent cohorts** (BiDirect, MACS, IXI, ADNI) — all processed with **one pipeline** (CAT12 build 1450 / SPM12 v7487, DARTEL, 8mm FWHM). Calibration n=1000.
- **What it did NOT investigate:** cross-pipeline anything — every dataset used the identical processing stack; the shift tested is cohort/population shift, not software shift.

### 6. Olsson et al. 2022 (non-MRI analogue, included because it is the closest design match found anywhere)
- **Title:** "Estimating diagnostic uncertainty in artificial intelligence assisted pathology using conformal prediction."
- **Venue/DOI:** *Nat Commun* 13:7761, 2022, doi `10.1038/s41467-022-34945-8`. [FULL-TEXT main read; Supplementary Information not fetchable — anchor-only link, not a retrievable URL.]
- **What it investigated:** the **same 449 pathology slides**, digitized on two different scanners (Aperio and Hamamatsu); Mondrian inductive conformal prediction trained on one scanner's images, tested on the other's.
- **Verified quote (≤15 words):** "The prediction regions were non-valid when evaluated on the new scanner (Hamamatsu)."
- **Why it matters:** this is a same-sample, device/software-change, frozen-model, conformal-coverage design — structurally identical to this study's design, just in pathology/scanner-software rather than MRI/processing-pipeline. It is the single strongest piece of **indirect** evidence for what this study's results in fact show (coverage collapse under a frozen conformal model moved to a new acquisition/processing device), and it is the paper a reviewer is most likely to cite as "this has essentially been shown before, just not in MRI."

---

## Other papers that could potentially invalidate the novelty claim (full list, not just the headline six)

- **Li et al. 2022, arXiv:2208.12909** — "Pipeline-invariant representation learning for neuroimaging." Same 2000 UK Biobank T1 scans, each preprocessed by **3 pipelines** (C-PAC:Default, C-PAC:fMRIPrep, UKB FSL-SPM), frozen-model transfer across pipelines for a 10-class age×sex classification task. **No conformal intervals, no coverage — classification accuracy only.** Closest same-scan/multi-pipeline/frozen-model precedent found in either pass; does not use conformal prediction, so it does not invalidate the coverage-specific novelty claim, but it must be cited.
- **"Conformalized Large Language Models under Configuration Shift," arXiv:2608.01460 (2026)** — non-MRI (LLMs), reports that "configuration shift consistently erodes CP validity" across 9 LLMs/4 datasets/4 nonconformity scores. [ABSTRACT-ONLY] A non-MRI precedent for the general phenomenon, not MRI-specific.
- **COMPASS (Cheung, Veeraraghavan, Balakrishnan), arXiv:2509.22240** — weighted conformal prediction for medical **segmentation metrics** under covariate shift, not a pipeline swap. [ABSTRACT-ONLY]
- **Mehrtens, Bucher, Brinker, arXiv:2506.18162** — "Pitfalls of Conformal Predictions for Medical Image Classification," reports CP "unreliable under distributional shifts," restored by target-domain recalibration. [ABSTRACT-ONLY] Generic, not pipeline-specific.
- **"When Is a Conformal Guarantee Fair? Auditing Silent Subgroup Under-Coverage in Alzheimer's Disease Longitudinal Prediction," arXiv:2608.04254 (2026)** — subgroup under-coverage in AD prediction, **no pipeline dimension.** [ABSTRACT-ONLY/snippet; authors not resolved]
- No Chen/Satterthwaite-group preprint on pipeline or FreeSurfer-version shift was found between April and October 2026 in either pass, despite targeted searches of bioRxiv and the `andy1764` GitHub profile. **This is "not found," not "confirmed absent"** — forward-citation chasing failed entirely (Semantic Scholar 429 in every attempt), so a very recent or not-yet-indexed scoop cannot be ruled out.

---

## Preregistration / kill criterion — honest statement

**No formal, independently timestamped preregistration document exists for this study**, in the sense of a document committed or published *before* any results were seen, separate from the analysis code and research-report text. What does exist:
1. A QC-definition docstring in `abide_data.py` stating the three QC levels "were decided before looking at any coverage result" — this is a claim inside the code comments, not an externally verifiable timestamped artifact.
2. A kill criterion stated in the second research-agent report (dated 25 Oct 2026, task `wf-d7dbc55e-...`, Task G): *"not worth writing as a positive-result paper if (i) every non-control transition has coverage ≥0.87 with CI including 0, (ii) no age bin n≥30 falls below 0.80, (iii) a competing preprint appears before submission."* This report's internal date (25 Oct 2026) predates the final clean full re-run of `run_all.py` used for the manuscript (also 25 Oct 2026, later the same session) — so sequence is plausible but **there is no cryptographic timestamp, repository commit, or third-party registry (e.g., OSF) confirming the kill criterion was fixed before the full results were computed.** The research report itself is the closest thing to a dated protocol, and it is included in full above; it is not a formal preregistration in the OSF/AsPredicted sense.

**State this plainly in the manuscript and to any journal:** no formal preregistration was filed. If reviewers or the target journal require one, none can be retroactively produced honestly — this should be disclosed, not worked around.
