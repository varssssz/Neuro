"""Build RESULTS_CHANGELOG.md: v1 manuscript claims vs v2 values (from numbers.json, with source).
Usage: python -I code_v2/build_changelog.py manuscript_v2/numbers.json RESULTS_CHANGELOG.md"""
import sys, json
N = json.load(open(sys.argv[1])); OUT = sys.argv[2]
v = lambda *ks: " / ".join(N[k]["value"] for k in ks)
src = lambda *ks: "<br>".join(sorted({N[k]["source"].split(":")[0] for k in ks}))
SCRIPT = {"inference": "code_v2/phase3_inference.py", "recalibration": "code_v2/run_recalibration_v2.py + summarise_recalibration.py",
          "frozen": "code_v2/run_frozen_v2.py", "age": "code_v2/phase5_age.py", "standardisation_monitor": "code_v2/phase4_std_monitor.py",
          "secondary_textbook": "code_v2/run_secondary_textbook.py", "Experimental_Results": "original_v1 run_all.py (unchanged)"}
def script(*ks):
    s = set()
    for k in ks:
        f = N[k]["source"]
        s |= {sc for key, sc in SCRIPT.items() if key in f}
    return "<br>".join(sorted(s)) or "—"
rows = [
 ("Source coverage, mean (18)", "0.910", ["cov_src_mean"], "rank correction (v1 used one rank above the stated equation)"),
 ("Raw coverage, mean; range", "0.489; 0.001–0.876", ["cov_raw_mean", "cov_raw_min", "cov_raw_max"], "rank correction; median now also reported"),
 ("Raw coverage, median", "not reported", ["cov_raw_median"], "added"),
 ("Std coverage, mean", "0.853", ["cov_std_mean"], "rank correction"),
 ("Std transitions < 0.88", "12 of 18", ["cov_std_n_below_088"], "two transitions crossed the 0.88 threshold; threshold de-emphasised"),
 ("Std worst transition", "0.757 (FS5.1→CIVET)", ["fs51_civet_cov_std"], "rank correction"),
 ("Std transitions losing ≥ 5 points", "9", ["d_std_vs_src_n_le_minus005"], "unchanged"),
 ("Std vs source significance", "16/18 McNemar p<0.05 (seed 0); 11 after Bonferroni", ["d_std_vs_src_n_holm_sig"], "method replaced: site-cluster CR1 test, Holm"),
 ("Std − source, median [CI]", "not reported", ["d_std_vs_src_median", "d_std_vs_src_median_lo", "d_std_vs_src_median_hi"], "added (joint site bootstrap)"),
 ("Per-site std (refit) coverage mean (range)", "0.897 (0.875–0.923)", ["cov_sitestd_mean", "cov_sitestd_min", "cov_sitestd_max"], "rank correction"),
 ("Per-site std width vs source (y)", "14.9 vs 12.1 (+23%)", ["width_site", "width_src", "width_site_pct"], "rank correction"),
 ("Recalibration coverage k=10/20/40/80/160", "0.813/0.900/0.948/0.922/0.909", ["rc_cov_10", "rc_cov_20", "rc_cov_40", "rc_cov_80", "rc_cov_160"], "four-way split + textbook rank; k=10 now infinite interval (coverage 1)"),
 ("Recalibration width k=20/40/80/160 (y)", "16.6/19.2/16.4/15.4 (k=10: 13.6)", ["rc_w_20", "rc_w_40", "rc_w_80", "rc_w_160"], "k=10 width now ∞"),
 ("Std-then-recalibrate width k=20/40/80/160", "≈1.2 y narrower on average", ["rch_w_20", "rch_w_40", "rch_w_80", "rch_w_160"], "recomputed"),
 ("Source width in recalibration design (y)", "12.1", ["rc_width_src"], "repartitioned training (40%) + rank"),
 ("Age bins source coverage <12/12–15/15–18/18–25/≥25", "0.916/0.981/0.990/0.915/0.277 (subject-weighted)", ["age_lt12_cov_src", "age_b1215_cov_src", "age_b1518_cov_src", "age_b1825_cov_src", "age_ge25_cov_src"], "rank correction + aggregation now mean over transitions with site CIs"),
 ("Age ≥25 std coverage", "0.183", ["age_ge25_cov_std"], "as above"),
 ("Diagnosis strata (ASD/TD)", "0.923/0.899 source", ["cov_src_mean"], "removed from main text (not needed for conclusions); v1 values remain in original_v1 strata_main.csv"),
 ("MMD std detection (non-control, b=50)", "0.139", ["mon_mmd_std_detect_mean_nc"], "re-executed; unchanged (no conformal rank involved)"),
 ("Spearman detection vs coverage change", "−0.49", ["mon_spearman_mmd_std_detect_vs_d_std"], "coverage change now textbook rank; unchanged"),
 ("LightGBM raw/std mean", "0.716/0.854 (over 20)", ["v1_lgbm_cov_raw", "v1_lgbm_cov_std"], "denominator corrected to 18; still v1 rank (not re-executable)"),
 ("Ridge comparator in LightGBM sentence", "0.531/0.859 (over 20)", ["cov_raw_mean", "cov_std_mean"], "denominator corrected to 18; sentence removed"),
 ("LOSO src/raw/std/site-std/refit", "0.867/0.504/0.818/0.910/0.919 (over 20)", ["loso_cov_src", "loso_cov_raw", "loso_cov_std", "loso_cov_sitestd", "loso_cov_src_site"], "denominator 18 + textbook rank; site-std labelled transductive"),
 ("LOSO worst site after std", "0.672", ["loso_min_site_cov_std"], "as above"),
 ("Mixed training, all pipelines (range)", "0.880–0.914", ["mix_all_min", "mix_all_max"], "textbook rank"),
 ("Mixed training, target left out (FS5.1/5.3/6.0/ANTs/CIVET)", "0.906/0.894/0.853/0.784/0.216", ["mix_loo_FS51", "mix_loo_FS53", "mix_loo_FS60", "mix_loo_ANTS", "mix_loo_CIVET"], "textbook rank"),
 ("QC0 std range", "0.771–0.923", ["qc0_std_min", "qc0_std_max"], "textbook rank"),
 ("QC2 std range", "0.719–0.957", ["qc2_std_min", "qc2_std_max"], "textbook rank (larger change: QC2 calibration sets ≈45, where one rank matters more)"),
 ("Weighted CP infinite-interval share", "90% / 79%", ["wcp_frac_inf", "wcp_reg_frac_inf"], "unchanged (weighted quantile not affected)"),
 ("MAE source / raw range (y)", "≈2.5–3.0 / 3.1–18.5", ["MAE_src_mean", "MAE_raw_min", "MAE_raw_max"], "unchanged (point predictions)"),
 ("Near-identity control", "FS5.3↔6.0 'near-identity control'; raw 0.899/0.908", ["fs53_fs60_cov_raw", "fs60_fs53_cov_raw"], "relabelled 'low-impact pair'; new negative control FS6.0-DK↔FS6.0"),
 ("Atlas-only control change", "0.914 vs 0.914", ["fs60dk_fs60_d_raw_vs_src", "fs60dk_fs60_d_raw_vs_src_clu_lo", "fs60dk_fs60_d_raw_vs_src_clu_hi"], "now with site-cluster CI"),
 ("Kill criterion", "'set in advance': ≥5-point loss", ["cov_std_mean"], "removed; no preregistration stated"),
]
L = ["# RESULTS_CHANGELOG.md", "",
     "Every number that changed between the original manuscript (v1, `original_v1/Manuscript/manuscript_full_text.txt`) and the revised manuscript (v2, `manuscript_v2/manuscript_v2.md`), "
     "with the numbers.json key(s), the generating script and the result file. The value of every key is traceable in `manuscript_v2/numbers.json` (field `source`). "
     "Generated by `code_v2/build_changelog.py`.", "",
     "Main causes of change: (1) conformal rank corrected to the stated equation (≈ −0.006 on frozen-model coverages); (2) recalibration re-run in a disjoint four-way split; "
     "(3) 18-vs-20 denominators unified to 18; (4) inference replaced (site-clustered, Holm); (5) controls relabelled.", "",
     "| Quantity | v1 manuscript | v2 manuscript | Reason | numbers.json key(s) | Script | Result file |", "|---|---|---|---|---|---|---|"]
for name, v1v, ks, why in rows:
    v2v = v(*ks) if "removed" not in why or name.startswith("Kill") is False else "—"
    if name in ("Diagnosis strata (ASD/TD)", "Kill criterion"): v2v = "—"
    L.append(f"| {name} | {v1v} | {v2v} | {why} | {', '.join(ks)} | {script(*ks)} | {src(*ks)} |")
L += ["", "## All numbers used in the v2 manuscript", "", "`manuscript_v2/numbers_used.json` lists the 160 keys substituted into the text; each key's value and source are in `manuscript_v2/numbers.json`."]
open(OUT, "w").write("\n".join(L) + "\n")
print(len(rows), "rows")
