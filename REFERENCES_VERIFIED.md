# REFERENCES_VERIFIED.md — Phase 6

Verification was done on 2026-10-30 against primary or publisher-indexed records (journal page, PMC/PubMed record, arXiv/bioRxiv record, or university repository).

Status:
- **V-PR**: verified, peer-reviewed.
- **V-PP**: verified, preprint, not peer-reviewed.
- **NV**: could not be verified; do not cite as is.

PMC full texts were blocked by a CAPTCHA, so bioRxiv full texts were used where content checks were needed.

## Core references (keep)

| # | Reference (as it should appear) | Status | What was checked |
|---|---|---|---|
| 1 | Chen AA, Seidlitz J, Gardner M, Bethlehem RAI, Dorfschmidt L, Kafadar E, Lifespan Brain Chart Consortium, Benitez A, Jensen JH, Vandekar S, Satterthwaite TD, Alexander-Bloch AF. *Calibration of MRI-based reference intervals to new samples.* **Imaging Neuroscience 2026, doi:10.1162/IMAG.a.1290** (online 10 Jul 2026; preprint doi:10.1101/2025.11.25.690453) | **V-PR** (updated 31 Oct 2026: journal version confirmed via MIT Press and PMC records; author list from the PubMed preprint record, check against the journal page) | The full text confirms several things. Calibration-set guidance is about 40 / 80 / at least 120 subjects for 90 / 95 / 99% intervals (this resolves the earlier "≥120" conflict: ≥120 applies to 99%). Separate lower/upper constants are used. Data were processed with FreeSurfer 6.0.1. Conditional (age/sex) coverage is listed as a limitation/future work. Notes: the full author list was not extracted, so copy it from the bioRxiv page. The CQR section uses rank ⌈(1−α/2)(n+1)⌉, while the ReForm section writes (1−α)(1+1/n); check the published version. **The manuscript's "Imaging Neuroscience, doi:10.1162/IMAG.a.1290" could not be confirmed:** the DOI resolver and the MIT Press pages were not reachable from this session. Cite the preprint unless the authors confirm the journal version themselves. |
| 2 | Samanta U, Shafiei G, Karim H, Satterthwaite TD, Chen A. *Quantifying Uncertainty in Brain Age Predictions via Conformal Prediction.* bioRxiv, posted 11 Sep 2026. doi:10.64898/2026.09.04.749442 | **V-PP** | Data: Reproducible Brain Charts (ages 6–22, n = 3,847), FreeSurfer 6.0.1 / sMRIPrep / C-PAC. Methods: CV+ with RidgeCV, compared with split conformal and CQR. It reports below-nominal coverage in some age bins and **no cross-pipeline evaluation**. It shares a senior author (A. Chen) with ReForm, so the scoop risk from this group is real. |
| 3 | Korbmacher M, Westlye LT, Maximov II. *FreeSurfer version-shuffling can enhance brain age predictions.* NeuroImage: Reports 2024;4(3):100214. doi:10.1016/j.ynirp.2024.100214 | **V-PR** | UK Biobank; the same scans processed with FS 5.3 and 7.4.1; cross-version brain-age differences of up to 1.42 y; version-shuffled training. Point predictions only, with no intervals. |
| 4 | Sanz-Robinson J, Glatard T, Poline J-B. *Leveraging Segmentation Variability to Improve Brain Age Prediction.* bioRxiv 2026. doi:10.64898/2026.07.23.740400 | **V-PP** (authors confirmed 31 Oct 2026) | Title, DOI and abstract confirmed: four segmentation pipelines on five datasets; multi-pipeline aggregation improves prediction; differences relate to age and sex. No intervals. |
| 5 | Ernsting J, Winter NR, Leenings R, et al. *From Group-Differences to Single-Subject Probability: Conformal Prediction-based Uncertainty Estimation for Brain-Age Modeling.* arXiv:2302.05304 (2023) | **V-PP** | Title and authors confirmed. Conformal brain-age framework; N = 16,903. |
| 6 | Delzant E, Colliot O, Couvy-Duchesne B. *Choice of Processing Pipelines for T1-Weighted Brain MRI Impacts Association and Prediction Analyses.* Hum Brain Mapp 2025;46(16):e70372. doi:10.1002/hbm.70372 | **V-PR** | UK Biobank, 39,655 scans; FSL, CAT12 and FreeSurfer representations; prediction accuracy differs by pipeline. No intervals. |
| 7 | Olsson H, Kartasalo K, Mulliqi N, et al. *Estimating diagnostic uncertainty in artificial intelligence assisted pathology using conformal prediction.* Nat Commun 2022;13:7761. doi:10.1038/s41467-022-34945-8 | **V-PR** | Conformal predictions became non-valid on an external scanner. This is the closest cross-domain precedent for "processing or acquisition change breaks conformal validity". |
| 10 | Bhagwat N, Barry A, Dickie EW, et al. *Understanding the impact of preprocessing pipelines on neuroimaging cortical surface analyses.* GigaScience 2021;10(1):giaa155. doi:10.1093/gigascience/giaa155 | **V-PR** | This is the source of the compare-surf-tools data. It uses the same ABIDE scans and the same five pipelines. Its ML component was a random-forest age prediction **within** each pipeline (RMSE 5.7–7.2 y). It did **not** test cross-pipeline transfer, and it reported **no** intervals or conformal prediction. |
| 11 | Tibshirani RJ, Barber RF, Candès EJ, Ramdas A. *Conformal prediction under covariate shift.* NeurIPS 32 (2019), pp. 2530–2540 | **V-PR** | NeurIPS proceedings page. |
| new | Lei J, G'Sell M, Rinaldo A, Tibshirani RJ, Wasserman L. *Distribution-free predictive inference for regression.* JASA 2018;113(523):1094–1111. doi:10.1080/01621459.2017.1307116 | **V-PR** | Split-conformal reference for the rank formula |
| new | Barber RF, Candès EJ, Ramdas A, Tibshirani RJ. *Conformal prediction beyond exchangeability.* Ann Stat 2023;51(2):816–845. doi:10.1214/23-AOS2276 | **V-PR** | Coverage gap under non-exchangeability; this should replace the unverified "harmful shift" citations |
| new | Gretton A, Borgwardt KM, Rasch MJ, Schölkopf B, Smola A. *A kernel two-sample test.* JMLR 2012;13:723–773 | **V-PR** | MMD |
| new | Di Martino A, Yan CG, Li Q, et al. *The autism brain imaging data exchange…* Mol Psychiatry 2014;19:659–667. doi:10.1038/mp.2013.78 | **V-PR** (from a reference list in an indexed HBM record) | ABIDE I |
| new | Li X, Fedorov A, Mathur M, Abrol A, Kiar G, Plis S, Calhoun V. *Pipeline-Invariant Representation Learning for Neuroimaging.* arXiv:2208.12909 (ML4H 2022 extended abstract) | **V-PP** | Pipeline-invariant deep representations (UK Biobank). Not intervals. |
| new (optional) | Haddad E, et al. *Multisite test–retest reliability and compatibility of brain metrics derived from FreeSurfer versions 7.1, 6.0, and 5.3.* Hum Brain Mapp 2023;44(4):1515–1532. doi:10.1002/hbm.26147 | **V-PR** | Between-version incompatibility of regional thickness |

## Remove or replace

| v1 ref | Problem | Action |
|---|---|---|
| [8] Hou et al., arXiv:2605.07212 (conformal EEG) | **NV**: the arXiv record could not be retrieved, and searches found no matching paper | Remove; it is not needed for the argument |
| [9] Amoukou et al., arXiv:2412.12910 ("conformal under harmful shift") | **NV**: the record could not be retrieved | Replace with Barber et al. 2023 (verified) |

Software citations (scikit-learn, NumPy, SciPy, pandas, LightGBM) must still be added by the authors from each project's citation page. The session did not verify them.
