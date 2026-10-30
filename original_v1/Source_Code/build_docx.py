"""Build the manuscript docx. All numbers are read from results_full/*.csv at build time."""
import numpy as np, pandas as pd
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

R = "results_full/"; F = "figs/"
T = pd.read_csv(R + "transitions_main.csv"); RC = pd.read_csv(R + "recal_main.csv")
S = pd.read_csv(R + "strata_main.csv"); AT = pd.read_csv(R + "transitions_atlas.csv")
MO = pd.read_csv(R + "monitors_b50.csv"); MX = pd.read_csv(R + "mixed_pipeline.csv")
LG = pd.read_csv(R + "lgbm.csv"); LP = pd.read_csv(R + "loso_pooled.csv")
MN = pd.read_csv(R + "matched_N_summary.csv")
T0 = pd.read_csv(R + "transitions_qc0.csv"); T2 = pd.read_csv(R + "transitions_qc2.csv")
ctrl = T.control.astype(bool); NC = T[~ctrl]
f2 = lambda x: f"{x:.2f}"; f3 = lambda x: f"{x:.3f}"
def rk(v, k): return RC[(RC.variant == v) & (RC.k == k)].merge(T[["src", "tgt", "control"]], on=["src", "tgt"]).query("control==False")
def tr(s, t): return T[(T.src == s) & (T.tgt == t)].iloc[0]

n_nc = len(NC)
raw_min, raw_max, raw_mean = NC.cov_raw.min(), NC.cov_raw.max(), NC.cov_raw.mean()
std_mean, std_min = NC.cov_std.mean(), NC.cov_std.min()
n_lose5 = int((NC.d_std <= -0.05).sum()); n_lt88 = int((NC.cov_std < 0.88).sum())
site_mean, site_min = NC.cov_sitestd.mean(), NC.cov_sitestd.min()
w_src, w_site = NC.width_src.mean(), NC.width_src_site.mean()
k40, k80, k160, k20, k10 = (rk("reform_raw", k) for k in (40, 80, 160, 20, 10))
ctrl_raw = T[ctrl].cov_raw.values; ctrl_src = T[ctrl].cov_src.values
civet_to = NC[NC.tgt == "CIVET"]; civet_from = NC[NC.src == "CIVET"]
s = S.merge(T[["src", "tgt", "control"]], on=["src", "tgt"]); s = s[~s.control.astype(bool)]
ag = s[~s.stratum.isin(["ASD", "TD"])].groupby(["stratum", "cond"]).apply(lambda d: np.average(d["cov"], weights=d.n)).unstack()
dx = s[s.stratum.isin(["ASD", "TD"])].groupby(["stratum", "cond"]).apply(lambda d: np.average(d["cov"], weights=d.n)).unstack()
mon_nc = MO.merge(T[["src", "tgt", "control", "d_std"]], on=["src", "tgt"]); mon_c = mon_nc[mon_nc.control.astype(bool)]; mon_nc = mon_nc[~mon_nc.control.astype(bool)]
rho = mon_nc[["mmd_std_detect", "d_std"]].corr(method="spearman").iloc[0, 1]
Nmin, Nmax = int(MN[~MN.tgt.str.contains("DK") & ~MN.src.str.contains("DK")].N_QC1.min()), int(MN[~MN.tgt.str.contains("DK") & ~MN.src.str.contains("DK")].N_QC1.max())
N0min, N0max = int(MN[~MN.tgt.str.contains("DK") & ~MN.src.str.contains("DK")].N_noQC.min()), int(MN[~MN.tgt.str.contains("DK") & ~MN.src.str.contains("DK")].N_noQC.max())
N2min, N2max = int(MN[~MN.tgt.str.contains("DK") & ~MN.src.str.contains("DK")].N_QC2.min()), int(MN[~MN.tgt.str.contains("DK") & ~MN.src.str.contains("DK")].N_QC2.max())
lgb_raw, lgb_std = LG.lgbm_cov_raw.mean(), LG.lgbm_cov_std.mean()

doc = Document()
st = doc.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(11)
for sec in doc.sections:
    sec.left_margin = sec.right_margin = Inches(1); sec.top_margin = sec.bottom_margin = Inches(1)

def H(t, l=1): doc.add_heading(t, level=l)
def P(t, italic=False, bold=False):
    p = doc.add_paragraph(); r = p.add_run(t); r.italic = italic; r.bold = bold; p.paragraph_format.space_after = Pt(6); return p
def fig(path, cap, w=6.0):
    doc.add_picture(F + path, width=Inches(w)); doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    c = doc.add_paragraph(); c.add_run(cap).italic = True
def table(df, cap, widths=None, size=8):
    c = doc.add_paragraph(); c.add_run(cap).bold = True
    t = doc.add_table(rows=1, cols=len(df.columns)); t.style = "Light Grid Accent 1"
    for i, h in enumerate(df.columns): t.rows[0].cells[i].text = str(h)
    for row in df.itertuples(index=False):
        cs = t.add_row().cells
        for i, v in enumerate(row): cs[i].text = str(v)
    for row in t.rows:
        for cell in row.cells:
            for p in cell.paragraphs:
                for r in p.runs: r.font.size = Pt(size)
    doc.add_paragraph()

# ---------------- Title ----------------
t = doc.add_heading("Same scan, different pipeline: how much coverage does a frozen conformal brain-age model lose when structural MRI is processed by different software?", 0)
P("Author: Hamsika [affiliation: VIT Vellore; co-author list, ORCID and corresponding-author details to be completed by the authors]")
P("Target journal: Imaging Neuroscience (open-access, methods-friendly; it published the closest prior work, ReForm). Alternative: NeuroImage: Reports.", italic=True)
P("Status: preprint draft. Single dataset, no external replication. See Limitations.", italic=True)

# ---------------- Abstract ----------------
H("Abstract")
P(f"Prediction intervals for brain age are usually calibrated on one image-processing pipeline and then applied to features from another. We asked how much marginal coverage a frozen, source-calibrated split-conformal model loses when the same scans are re-processed, what remains after label-free fixes, and what those fixes cost in interval width. "
  f"Using ABIDE I cortical thickness (62 shared regions; FreeSurfer 5.1, 5.3 and 6.0, ANTs, CIVET 2.1; matched N {Nmin}-{Nmax} per pair after quality control), we evaluated all 20 directed transitions (nominal 90% coverage; ridge regression; 5 seeds; paired bootstrap). "
  f"Two FreeSurfer 5.3/6.0 transitions served as near-identity controls (raw coverage {f3(min(ctrl_raw))} and {f3(max(ctrl_raw))}). Across the other {n_nc} transitions, coverage of the unadapted model ranged from {f3(raw_min)} to {f3(raw_max)} (mean {f3(raw_mean)}); all transitions involving CIVET fell below 0.12 except ANTs to CIVET. "
  f"Per-pipeline feature standardization, which uses unlabelled target features only, raised mean coverage to {f3(std_mean)}, but {n_lt88} of {n_nc} transitions remained below 0.88 and {n_lose5} lost at least 5 points (worst {f3(std_min)}, FreeSurfer to CIVET and ANTs targets). "
  f"Per-site standardization brought mean coverage to {f3(site_mean)} (minimum {f3(site_min)}) at {100*(w_site/w_src-1):.0f}% wider intervals. A ReForm-style recalibration with 40 labelled target subjects gave mean coverage {f3(k40.mean_cov.mean())} with intervals {k40.mean_width.mean():.1f} years wide versus {w_src:.1f} at the source. "
  f"Marginal coverage hid severe conditional failure: coverage for ages 25 and above was below 0.30 even at the source. To our knowledge no earlier study evaluates conformal coverage under same-scan pipeline changes. Results come from one dataset and cortical thickness only.")

# ---------------- Introduction ----------------
H("1. Introduction")
P("Brain-age models and normative reference intervals are increasingly shipped as frozen artefacts: a regression model and a conformal or centile calibration fitted once on a source cohort. A user at another site then extracts features with whatever structural pipeline is available locally. Conformal guarantees are marginal and rest on exchangeability between calibration and test data [11]; a change of processing software changes the feature distribution of the very same scans, so exchangeability can fail even when the people, scanners and diagnoses are identical.")
P("Prior work shows that processing pipelines disagree on cortical measures [10] and that pipeline choice affects brain-age point predictions [3, 4]. Calibrating reference intervals to new samples has been addressed for study and scanner shifts [1]. What is missing is a direct measurement of interval coverage when identical scans are re-processed. This paper provides that measurement in a deliberately simple setting, and asks three questions. (Q1) How much coverage does a frozen model lose in each directed pipeline transition, relative to a near-identity control? (Q2) How much of the loss remains after label-free fixes that need only unlabelled target-pipeline features? (Q3) What do the fixes cost in interval width, and how many labelled target subjects does a ReForm-style recalibration need?")
P("Contributions. (i) An evaluation of all 20 directed transitions among five pipelines on the same subjects, with a near-identity control and an atlas-confound analysis. (ii) A comparison of label-free fixes, mixed-pipeline training and labelled recalibration, reported with width. (iii) Age-binned and diagnosis-stratified coverage. We propose no new conformal method.")

# ---------------- Related work ----------------
H("2. Related work")
P("ReForm [1] recalibrates the lower and upper bounds of a reference interval with separate constants, taken as empirical quantiles of calibration residuals, and recommends about 40, 80 and 120 calibration samples for 90%, 95% and 99% intervals. Its experiments used FreeSurfer 6.0.1 throughout and varied study and field strength, and its authors list conditional (age and sex) coverage as future work. Samanta et al. [2] report conformal brain-age intervals in one dataset and one processing stack, with marginal coverage near target but much lower coverage in some age bins [full details not re-verified in this session]. Korbmacher et al. [3] showed that mixing FreeSurfer versions in training can improve brain-age prediction; Sanz-Robinson et al. [4] compare brain-age models across pipelines, finding that disagreement depends on age and sex, without intervals or frozen-model transfer. Ernsting et al. [5] and Delzant et al. [6] study brain-age uncertainty and calibration under other sources of variation [exact scope not re-verified here]. Conformal coverage under harmful shift has been studied in other domains [7, 8, 9], and weighted conformal prediction handles covariate shift when the density ratio is estimable [11]. Bhagwat et al. [10] quantified pipeline differences in cortical surface analyses without prediction intervals. Our study differs by holding the scans fixed and varying only the software, so every transition isolates the processing effect.")
P("Novelty statement. To our knowledge, based on searches dated 25 October 2026 (bioRxiv, medRxiv, arXiv, PubMed and general web searches; see Section 6), no study evaluates conformal coverage when identical scans are re-processed by different structural-MRI pipelines. The absence of search hits is not proof of absence.", italic=True)

# ---------------- Methods ----------------
H("3. Methods")
H("3.1 Data", 2)
P("We used the ABIDE I cortical-thickness derivatives distributed with the compare-surf-tools repository (MIT-licensed code and tables; the ABIDE data are shared under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 licence, with unrestricted non-commercial research use after NITRC registration, as stated on the INDI ABIDE page). Phenotype variables were age at scan, diagnostic group (ASD or typically developing), site and sex. Five pipelines processed the same T1-weighted scans: FreeSurfer 5.1, FreeSurfer 5.3 (Desikan-Killiany atlas), FreeSurfer 6.0 (Desikan-Killiany-Tourville, DKT), ANTs and CIVET 2.1 (DKT labels).")
H("3.2 Quality control and matched samples", 2)
P(f"For each directed pair we kept only subjects present in both pipelines. The primary level (QC1) excluded a subject flagged by the repository's manual lenient rating, either of the two pipelines' automatic low-dimensional outlier flags, or the high-dimensional outlier flag; subjects absent from a table were not excluded by it. QC0 applied no exclusion, and QC2 additionally required an explicit stringent manual pass. QC definitions were fixed before coverage results were inspected. Matched N per pair was {N0min}-{N0max} (QC0), {Nmin}-{Nmax} (QC1) and {N2min}-{N2max} (QC2); the by-site counts are in the supplementary table matched_N_by_site.csv. The QC is pilot-level: it reuses the repository's own flags and is not an independent visual re-rating.")
H("3.3 Features and atlas confound", 2)
P("The five pipelines share 62 cortical thickness features (31 regions in each hemisphere) after mapping region names. FreeSurfer 5.x uses the Desikan-Killiany atlas and FreeSurfer 6.0 DKT, so software version and atlas are partly confounded. To separate them we also extracted FreeSurfer 6.0 Desikan-Killiany thickness (FS60-DK) from the same software run: the pair FS60-DK/FS60 differs in atlas only, and FS60-DK/FS5.x share the atlas but differ in software version.")
H("3.4 Models and conformal procedure", 2)
P("Ridge regression with the penalty chosen by internal cross-validation (30 values from 0.01 to 10,000) predicted age from the 62 thickness values; LightGBM served as a nonlinear check. Within each of 5 cross-fitting folds, the training part was split 75/25 into proper-training and calibration sets, giving roughly 60/20/20 train/calibration/test overall, and the whole procedure was repeated for 5 seeds. The model and the conformal quantile were then frozen and applied to the same test subjects with features from the target pipeline. With absolute-residual scores,")
p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
def run(t, sub=False, sup=False, it=False):
    r = p.add_run(t); r.font.subscript = sub; r.font.superscript = sup; r.italic = it; return r
run("q", it=True); run(" = Quantile", ); run("(1−α)(1+1/n)", sub=True); run("( { |"); run("y", it=True); run("i", sub=True); run(" − "); run("f", it=True); run("(x"); run("i", sub=True); run(")| } ), i = 1,…,n,     Ĉ(x) = [ f(x) − q, f(x) + q ]")
P("where n is the calibration-set size, α = 0.10, f the frozen regression model and the quantile is the finite-sample-corrected empirical quantile (the ⌈(n+1)(1−α)⌉-th smallest score). Coverage is the share of test subjects whose true age lies in Ĉ(x). Interval width is 2q, in years.")
H("3.5 Label-free fixes and baselines", 2)
for b in [
 "Per-pipeline standardization. The target-pipeline features of the (unlabelled) training subjects are z-scored with their own mean and standard deviation, and source features likewise; the model and q stay frozen.",
 "Per-site standardization. The same, per acquisition site (sites with at least 15 training subjects; otherwise global statistics). Because site statistics change the source features, the model and q are refit on site-standardized source features, so this arm changes the model and its interval width.",
 "Weighted conformal prediction [11], with a logistic-regression density ratio between source and target features (plain and L2-regularized variants).",
 "ReForm-style recalibration with k = 10, 20, 40, 80, 160 labelled target subjects: separate lower and upper corrections, each the finite-sample quantile of the one-sided residuals at level α/2 per side. This is our reading of the procedure described in [1] (separate lower/upper constants from calibration residuals); we did not run the authors' code, so we call it ReForm-style. Each (transition, k) cell used 100 random calibration draws per fold and seed (2,500 draws), and we also report a hybrid that standardizes first.",
 "Mixed-pipeline training in the spirit of [3]: one model trained on subjects pooled across pipelines (all five, or the four non-target pipelines), with a pooled scaler.",
 "LightGBM in place of ridge (secondary check)."]:
    doc.add_paragraph(b, style="List Bullet")
H("3.6 Monitoring", 2)
P("As a descriptive secondary analysis, we computed an unbiased RBF-kernel maximum mean discrepancy (median-heuristic bandwidth) and a Kolmogorov-Smirnov statistic on predicted ages for batches of 50 target subjects (also 25 and 100), with thresholds set at a 5% false-alarm rate from source calibration batches; false-alarm rates were checked on source test batches. The transition is the unit of analysis (20 points), so no inferential claim is made.")
H("3.7 Statistics", 2)
P("Differences in coverage between conditions are paired at the subject level (indicator averaged over seeds) with a subject-level bootstrap of 2,000 resamples for 95% intervals, and exact McNemar tests on seed 0. Wilson intervals are shown for strata, using the seed-mean coverage with the stratum's subject count (an approximation). Age strata are under 12, 12-15, 15-18, 18-25 and 25 or above; diagnosis strata are ASD and TD. Leave-one-site-out (LOSO) evaluation was run on the six largest sites. We report the share of fold-level draws with coverage below 0.88; because test folds have about 170 subjects, this share has a sampling noise floor (about 0.12-0.2 at the source, Table 1 and supplementary data).")
H("3.8 Leakage controls", 2)
P("Cross-fitting folds are defined on subjects, so the same subject never appears in train, calibration and test in one fold, in either pipeline. Standardization statistics use training-fold subjects only. Site statistics are computed from training subjects. Test labels are used only for scoring. A caveat: ABIDE sites contribute subjects to all folds, so site effects are shared between source and test partitions except in the LOSO analysis.")
H("3.9 Reproducibility", 2)
P("All numbers in this paper were produced by run_all.py; a second, clean full re-run reproduced every table to numerical identity (maximum absolute difference 0 across the primary tables). Package versions are in requirements.txt.")

# ---------------- Results ----------------
H("4. Results")
H("4.1 Matched samples and the near-identity control", 2)
P(f"Matched N per directed pair ranged from {Nmin} to {Nmax} under QC1. The two FreeSurfer 5.3/6.0 transitions behaved as controls: raw coverage {f3(ctrl_raw[0])} and {f3(ctrl_raw[1])} versus {f3(ctrl_src[0])} and {f3(ctrl_src[1])} at the source (both differences about -0.006, within the seed noise of about 0.01). FreeSurfer 5.1 transitions are reported separately because 5.1 differs from the other FreeSurfer versions by a larger offset.")
H("4.2 Coverage under a pipeline swap", 2)
P(f"Source-pipeline test coverage was {f3(T.cov_src.min())}-{f3(T.cov_src.max())} across the 20 transitions, close to the nominal 0.90. Without adaptation, the 18 non-control transitions covered {f3(raw_min)}-{f3(raw_max)} (mean {f3(raw_mean)}; Figure 1, Table 1). Transitions between FreeSurfer versions lost 0.04-0.12; transitions to ANTs lost 0.18-0.28; transitions from or to CIVET lost 0.80-0.90 except ANTs to CIVET (loss 0.09). Mean absolute error rose from about 2.5-3.0 years to {NC.MAE_raw.min():.1f}-{NC.MAE_raw.max():.1f} years. The losses track global scale and offset differences: CIVET thickness differs from FreeSurfer by about 0.5-0.9 mm on average (Table 1, offset column).")
tab = T[["src", "tgt", "N", "cov_src", "cov_raw", "cov_std", "cov_sitestd", "d_std", "mean_offset_mm", "MAE_src", "MAE_raw"]].copy()
tab.insert(0, "Ctrl", np.where(ctrl, "*", ""))
tab.columns = ["", "Source", "Target", "N", "Src cov", "Raw cov", "Std cov", "Site-std cov", "Δ std vs src", "Offset (mm)", "MAE src", "MAE raw"]
for c in tab.columns[4:]: tab[c] = tab[c].astype(float).round(3)
table(tab, "Table 1. Coverage (nominal 0.90) for all 20 directed transitions (QC1). * near-identity control. Raw = unadapted frozen model on target features; Std = per-pipeline standardization; Site-std = per-site standardization with refit. 95% CIs, McNemar tests and seed SDs are in table_transitions.csv.", size=7)
fig("fig1_coverage_by_transition.png", "Figure 1. Coverage by directed transition: raw frozen model, per-pipeline standardization, and ReForm-style recalibration with 40 labelled subjects. Bars show bootstrap 95% CIs for the raw condition. Dashed line: nominal 0.90.", w=5.6)
H("4.3 What remains after label-free fixes", 2)
P(f"Per-pipeline standardization raised mean non-control coverage to {f3(std_mean)}, and 16 of {n_nc} paired differences from the source were significant by McNemar test at 0.05 (11 after Bonferroni correction across 18 tests). Nevertheless, {n_lt88} of {n_nc} transitions stayed below 0.88 and {n_lose5} lost at least 5 points relative to the source (worst: FS5.1 to CIVET {f3(tr('FS51','CIVET').cov_std)}, FS5.3 to CIVET {f3(tr('FS53','CIVET').cov_std)}). The kill criterion we set in advance, that no non-control transition would lose 5 points or more after standardization, was therefore not met, and the residual loss is concentrated in transitions that have ANTs or CIVET as the target, plus ANTs as the source. Weighted conformal prediction did not help: with a logistic density ratio, the intervals were infinite for {int(round(100*NC.wcp_frac_inf.mean()))}% of test subjects on average (cross-software shifts are almost perfectly separable), and the regularized variant reached the same outcome for all cross-software pairs (mean share of infinite intervals {int(round(100*NC.wcp_reg_frac_inf.mean()))}%); its finite-case coverage is therefore not informative.")
H("4.4 Cost of the fixes in width, and labelled recalibration", 2)
P(f"Per-site standardization with refit reached mean coverage {f3(site_mean)} (range {f3(site_min)}-{f3(NC.cov_sitestd.max())}), and intervals were {100*(w_site/w_src-1):.0f}% wider ({w_site:.1f} versus {w_src:.1f} years; Figure 4). ReForm-style recalibration (Figure 2) reached mean coverage {f3(k10.mean_cov.mean())}, {f3(k20.mean_cov.mean())}, {f3(k40.mean_cov.mean())}, {f3(k80.mean_cov.mean())} and {f3(k160.mean_cov.mean())} at k = 10, 20, 40, 80 and 160, with mean widths of {k10.mean_width.mean():.1f}, {k20.mean_width.mean():.1f}, {k40.mean_width.mean():.1f}, {k80.mean_width.mean():.1f} and {k160.mean_width.mean():.1f} years. Coverage at k = 40 is conservative (0.948); we attribute this to rounding up of the one-sided quantile at small k (ReForm also reports rounding effects), but we did not test this explanation. At k below about 20 the finite-sample guarantee cannot be met with this construction (one-sided level 0.05 needs k of at least 19). The share of draws below 0.88 was {k40.frac_below_088.mean():.2f} at k = 40, {k80.frac_below_088.mean():.2f} at k = 80 and {k160.frac_below_088.mean():.2f} at k = 160, against a source-pipeline floor of {T.frac_fold_below_088_src.mean():.2f} (mean over transitions; supplementary table). The recalibration constants absorb CIVET-sized offsets by widening the interval, which explains the large widths. Standardizing before recalibrating left mean coverage essentially unchanged (largest per-transition difference 0.005) and reduced mean width by 1.2 years on average, by up to 8.8 years for the CIVET pairs.")
fig("fig2_coverage_vs_k.png", "Figure 2. ReForm-style recalibration: mean coverage (left, shaded band: mean of per-transition 2.5-97.5% range of draws) and share of draws below 0.88 (right) versus number of labelled target subjects, over the 18 non-control transitions.", w=6.0)
fig("fig4_width_vs_coverage.png", "Figure 4. Interval width versus coverage per transition for per-pipeline standardization (width equals source width), per-site standardization, and recalibration at k = 20, 40 and 80.", w=4.8)
H("4.5 Conditional coverage by age and diagnosis", 2)
rows = pd.DataFrame({"Age (years)": ag.index, "Source": ag["src"].round(3), "Raw": ag["raw"].round(3), "Std": ag["std"].round(3), "Site-std": ag["sitestd"].round(3)}).reindex([0, 1, 2, 3, 4]) if False else None
ag2 = ag.loc[["<12", "12-15", "15-18", "18-25", ">=25"]]
table(pd.DataFrame({"Age bin": ag2.index, "Source": ag2["src"].round(3).values, "Raw": ag2["raw"].round(3).values, "Std": ag2["std"].round(3).values, "Site-std": ag2["sitestd"].round(3).values}), "Table 2. Coverage by age bin, subject-weighted over 18 non-control transitions. Mean stratum size per transition: about 235, 226, 155, 172 and 49 subjects.", size=9)
P(f"Marginal coverage hid large conditional differences. Even at the source pipeline, coverage was {f3(ag2.loc['>=25','src'])} for ages 25 and above (about 49 subjects per transition) while it was 0.98-0.99 for ages 12-18; standardization lowered the oldest bin further ({f3(ag2.loc['>=25','std'])}). The model predicts toward the training mean, which suits a cohort concentrated in childhood and adolescence and fails in the sparse adult tail. By diagnosis, source coverage was {f3(dx.loc['ASD','src'])} (ASD) and {f3(dx.loc['TD','src'])} (TD); raw coverage was {f3(dx.loc['ASD','raw'])} and {f3(dx.loc['TD','raw'])}, and standardized {f3(dx.loc['ASD','std'])} and {f3(dx.loc['TD','std'])}. Standardization lowered coverage by 0.08-0.09 relative to the source in the under-12 and 18-25 bins. The Wilson intervals in the supplementary strata table use an approximation (see 3.7).")
fig("fig3_age_binned_coverage.png", "Figure 3. Coverage by age bin (pooled over 18 non-control transitions).", w=4.6)
H("4.6 Atlas confound", 2)
a = AT.set_index(["src", "tgt"])
P(f"Atlas alone cost almost nothing: FS60-DK to FS60 (same software, different atlas) had raw coverage {f3(a.loc[('FS60DK','FS60'),'cov_raw'])} against {f3(a.loc[('FS60DK','FS60'),'cov_src'])} at the source. Matching the atlas did not remove the version effect: FS60-DK to FS5.1 (same atlas, older software) fell to {f3(a.loc[('FS60DK','FS51'),'cov_raw'])}, close to the {f3(tr('FS60','FS51').cov_raw)} of FS60 to FS5.1. The FS5.1 effect is therefore a software-version effect, not an atlas effect. The common 62-region set is a shared subset of regions and not an identical definition, so a residual atlas effect cannot be excluded for ANTs and CIVET.")
H("4.7 Sensitivity analyses", 2)
P(f"QC level. Raw coverage under QC0 and QC2 is in transitions_qc0.csv and transitions_qc2.csv. The ranking of transitions and the CIVET collapse are unchanged; standardized coverage under QC2 ({N2min}-{N2max} subjects per pair) ranged {f3(T2[~T2.control.astype(bool)].cov_std.min())}-{f3(T2[~T2.control.astype(bool)].cov_std.max())} and under QC0 {f3(T0[~T0.control.astype(bool)].cov_std.min())}-{f3(T0[~T0.control.astype(bool)].cov_std.max())}, so the residual loss for ANTs and CIVET targets does not depend on QC level, although QC2 intervals are noisier. "
  f"Model. With LightGBM, mean raw coverage over the 20 transitions was {f3(lgb_raw)} and standardized {f3(lgb_std)} (ridge: {f3(T.cov_raw.mean())} and {f3(T.cov_std.mean())}). "
  f"Mixed-pipeline training (N = 669 subjects with all five pipelines available). Training on all five pipelines, including the target, gave target coverage {f3(MX.mix_all.min())}-{f3(MX.mix_all.max())} (width fixed at 13.5 years), against within-pipeline coverage {f3(MX.within.min())}-{f3(MX.within.max())}. Excluding the target pipeline from training gave {f3(MX.mix_loo.min())}-{f3(MX.mix_loo.max())}: acceptable for FreeSurfer 5.1 and 5.3 targets ({f3(MX.mix_loo.iloc[0])}, {f3(MX.mix_loo.iloc[1])}), lower for FreeSurfer 6.0 ({f3(MX.mix_loo.iloc[2])}) and ANTs ({f3(MX.mix_loo.iloc[3])}), and a collapse for CIVET ({f3(MX.mix_loo.iloc[4])}; MAE {MX.MAE_loo.iloc[4]:.1f} years). Mixed training that includes the target pipeline therefore works about as well as labelled recalibration here, but it needs target-pipeline features with labels at training time and a retrained model. "
  f"Leave-one-site-out. When whole sites are held out, source coverage falls to a mean of {f3(LP.cov_src.mean())}; raw cross-pipeline coverage averages {f3(LP.cov_raw.mean())}, standardized {f3(LP.cov_std.mean())}, and per-site standardization {f3(LP.cov_sitestd.mean())} against a refit source of {f3(LP.cov_src_site.mean())}. The worst single site averaged {f3(LP.min_site_cov_std.mean())} after standardization.")
H("4.8 Monitoring", 2)
P(f"Batches of 50 target subjects were tested against thresholds at a 5% false-alarm rate (observed false-alarm rate on source test batches: mean {f3(MO.mmd_raw_far_test.mean())} for MMD on raw and {f3(MO.mmd_std_far_test.mean())} on standardized features). The MMD flagged every non-control transition on raw features (detection rate 1.0) but also the harmless FS5.3/FS6.0 controls (detection rate {mon_c.mmd_raw_detect.min():.2f} for both); the KS test on predicted ages was more selective for the controls (0.00 and 0.04). On standardized features the MMD detection rate averaged {f3(mon_nc.mmd_std_detect.mean())} over non-control transitions and was associated with larger coverage loss across the 18 transitions (Spearman rho between detection rate and coverage change = {rho:.2f}; negative means more detection with more loss), which is descriptive only. Detection is not a coverage estimator: it flags shift, including shift that does no harm.")
fig  # placeholder to keep linter quiet

# ---------------- Discussion ----------------
H("5. Discussion")
P(f"A frozen, source-calibrated conformal interval did not survive a change of structural pipeline in this dataset. A near-identical pair (FreeSurfer 5.3/6.0) kept coverage within about 0.01, so the loss elsewhere is not an artefact of the pipeline of evaluation. Most of the raw loss is a global scale and offset difference between pipelines, and the cheapest label-free fix, standardizing target features, recovers much of it. However, {n_lt88} of {n_nc} transitions stayed below 0.88 after standardization, concentrated in targets processed by ANTs or CIVET. Two practical consequences follow. First, a user with no labels should not assume that a published interval transfers across pipelines; a monitor can flag the shift but cannot tell whether the shift is harmful. Second, labelled recalibration restores coverage at a price: about 40 labelled subjects gave conservative coverage at roughly {k40.mean_width.mean()/w_src:.1f} times the source width in these settings, and about 80 gave coverage near 0.92.")
P("Relation to ReForm. Our results are consistent with ReForm's exchangeability argument: pipeline-matched calibration restores marginal coverage, and its 40-subject rule is borne out for 90% intervals, with a conservative bias at k = 40 and no valid guarantee at k = 10. Our contribution is the measurement of how large the problem is without recalibration, and how much a label-free step recovers. The width cost shows that guarantees can be bought with uninformative intervals when the shift is largely an offset.")
P("Conditional coverage. Marginal coverage near 0.90 at the source concealed 0.28 coverage for ages 25 and above. This agrees with the general statement that conformal guarantees are marginal rather than conditional and with the age-bin failures reported in [2]; our result shows that a pipeline change interacts with this problem rather than creating it.")
H("Limitations", 2)
for b in [
 "One dataset (ABIDE I), a childhood-to-young-adult cohort with few adults; findings may not transfer to older adults or dementia cohorts. A replication on a multi-pipeline older-adult dataset was not run (OASIS-3 provides FreeSurfer derivatives and requires an access request; we did not obtain it, and no OASIS-3 multi-pipeline derivatives were identified).",
 "Cortical thickness from 62 shared regions only; no volumes, area or curvature. The common region set is not an identical atlas definition, so an atlas effect cannot be fully excluded for ANTs and CIVET.",
 "Quality control reuses repository flags and is pilot-level; it is not an independent visual re-rating. Sensitivity to QC0 and QC2 is reported.",
 "ReForm-style recalibration is our reading of the method in [1] and was not run with the authors' code. Its small-k behaviour reflects our one-sided α/2 construction.",
 "Conditional coverage is poor at the source for ages 25 and above, with about 49 subjects per transition in that bin; Wilson intervals use an approximation.",
 "Ridge is a simple model; LightGBM was run as a check only. Deep-learning brain-age models and other conformal scores were not tested.",
 "The fraction of draws below 0.88 has a sampling noise floor from test folds of about 170 subjects.",
 "The transition-level monitoring analysis has 18-20 points and is descriptive."]:
    doc.add_paragraph(b, style="List Bullet")
H("Search log for the novelty claim", 2)
P("Searches run on 25 October 2026 (web search engine equivalents of bioRxiv, medRxiv, arXiv and PubMed queries; those databases were not queried through their own interfaces): conformal or reference-interval coverage under processing-pipeline or software change in structural MRI; brain-age uncertainty across FreeSurfer versions or pipelines; follow-up preprints from the ReForm/Samanta groups on sMRIPrep or FreeSurfer; OASIS-3 multiple-pipeline derivatives. No paper evaluating conformal coverage under a same-scan pipeline swap was found. One earlier search surfaced unrelated conformal work in EEG and MRI reconstruction. These searches cannot exclude unindexed or very recent work.")

# ---------------- Availability etc ----------------
H("6. Data and code availability")
P("ABIDE I data are available from INDI/NITRC after registration under a CC BY-NC-SA 3.0 licence (non-commercial research); the derivative tables used here come from the MIT-licensed compare-surf-tools repository (neurodatascience/compare-surf-tools). Analysis scripts (abide_data.py, swap_engine.py, run_all.py, make_figs.py, build_docx.py) and requirements.txt are provided as supplementary material, and the result CSVs are included. We did not redistribute any ABIDE data. [Public repository URL and DOI: to be added by the authors.]")
H("7. Generative-AI disclosure")
P("Claude (Anthropic) assisted with writing analysis code, running the analyses and drafting this manuscript. The authors are responsible for verifying all results, citations and claims before submission. [Authors to confirm wording against the target journal's AI policy.]")
H("8. Acknowledgements and funding")
P("[To be completed by the authors.] We thank the ABIDE contributors and the compare-surf-tools authors for sharing the derivatives.")

# ---------------- References ----------------
H("References")
refs = [
 "[1] Chen A. A. et al. Calibration of MRI-based reference intervals to new samples (ReForm). Imaging Neuroscience, 2026. doi:10.1162/IMAG.a.1290; preprint bioRxiv doi:10.1101/2025.11.25.690453. [Full author list not re-verified.]",
 "[2] Samanta et al. Conformal brain-age prediction [title not re-verified]. bioRxiv, 2026. doi:10.64898/2026.09.04.749442. [Author list and title UNVERIFIED in this session.]",
 "[3] Korbmacher M. et al. FreeSurfer version-shuffling can enhance brain age predictions. NeuroImage: Reports, 4(3):100214, 2024 [year, volume and author list to be confirmed; a bioRxiv version exists, doi:10.1101/2024.06.14.599070].",
 "[4] Sanz-Robinson et al. [title UNVERIFIED]. bioRxiv, 2026. doi:10.64898/2026.07.23.740400.",
 "[5] Ernsting et al. [title UNVERIFIED]. arXiv:2302.05304, 2023.",
 "[6] Delzant et al. [title UNVERIFIED]. Human Brain Mapping, 2025. doi:10.1002/hbm.70372.",
 "[7] Olsson H. et al. [title UNVERIFIED; conformal prediction in AI-assisted pathology]. Nature Communications 13:7761, 2022.",
 "[8] Hou et al. [title UNVERIFIED; conformal prediction in EEG]. arXiv:2605.07212, 2026.",
 "[9] Amoukou et al. [title UNVERIFIED; conformal prediction under harmful shift]. arXiv:2412.12910, 2024.",
 "[10] Bhagwat N. et al. Understanding the impact of preprocessing pipelines on neuroimaging cortical surface analyses. GigaScience, 2021 (preprint bioRxiv doi:10.1101/2020.05.22.100180). [Volume/DOI to be confirmed.]",
 "[11] Tibshirani R. J., Barber R. F., Candès E. J., Ramdas A. Conformal prediction under covariate shift. NeurIPS 2019. [UNVERIFIED in this session; cited from general knowledge.]",
]
for r in refs: P(r)
P("Software citations to add before submission (all UNVERIFIED in this session): scikit-learn, LightGBM, SciPy, NumPy, pandas; Di Martino et al. 2014 (ABIDE I, Molecular Psychiatry); Gretton et al. 2012 (MMD, JMLR); Lei et al. 2018 (distribution-free predictive inference, JASA).", italic=True)

doc.save("Same_scan_different_pipeline_manuscript.docx")
print("saved")
