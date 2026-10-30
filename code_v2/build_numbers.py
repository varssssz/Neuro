"""Collect every number cited in manuscript_v2 into numbers.json, with its source.
Usage: python -I code_v2/build_numbers.py RESULTS_V2 ORIGINAL_RESULTS OUT_JSON"""
import sys, json
import numpy as np, pandas as pd

RV2, RV1, OUT = sys.argv[1:4]
N = {}
CONTROL = {("FS60", "FS53"), ("FS53", "FS60")}
nc = lambda d, a="src", b="tgt": d[[(x, y) not in CONTROL for x, y in zip(d[a], d[b])]]


def put(key, val, src, fmt="{:.3f}"):
    v = float(val) if not isinstance(val, str) else val
    N[key] = dict(value=(fmt.format(v) if not isinstance(v, str) else v), raw=None if isinstance(v, str) else v, source=src)


# --- samples
mn = pd.read_csv(f"{RV1}/matched_N_summary.csv"); mm = mn[[(s, t) in [(a, b) for a in "FS51 FS53 FS60 ANTS CIVET".split() for b in "FS51 FS53 FS60 ANTS CIVET".split() if a != b] for s, t in zip(mn.src, mn.tgt)]]
for c, k in (("N_QC1", "n_qc1"), ("N_noQC", "n_qc0"), ("N_QC2", "n_qc2")):
    put(f"{k}_min", mm[c].min(), f"{RV1}/matched_N_summary.csv:{c}:min(20 main)", "{:.0f}"); put(f"{k}_max", mm[c].max(), f"{RV1}/matched_N_summary.csv:{c}:max(20 main)", "{:.0f}")
# --- aggregates
A = pd.read_csv(f"{RV2}/inference/aggregate_18_noncontrol.csv"); a = A[A.rule == "textbook"].set_index("quantity"); a1 = A[A.rule == "v1"].set_index("quantity")
src = f"{RV2}/inference/aggregate_18_noncontrol.csv"
for q in ("cov_src", "cov_raw", "cov_std", "cov_src_site", "cov_sitestd"):
    for s in ("mean", "median", "min", "max"):
        put(f"{q}_{s}", a.loc[q, s], f"{src}:{q}:{s}")
        put(f"v1_{q}_{s}", a1.loc[q, s], f"{src}[rule=v1]:{q}:{s}")
    put(f"{q}_n_below_088", a.loc[q, "n_below_088"], f"{src}:{q}:n_below_088", "{:.0f}")
    put(f"v1_{q}_n_below_088", a1.loc[q, "n_below_088"], f"{src}[rule=v1]:{q}:n_below_088", "{:.0f}")
for q in ("d_raw_vs_src", "d_std_vs_src", "d_sitestd_vs_srcsite"):
    for s in ("mean", "mean_lo", "mean_hi", "median", "median_lo", "median_hi", "min", "max"):
        put(f"{q}_{s}", a.loc[q, s], f"{src}:{q}:{s}")
    for s in ("n_le_minus005", "n_holm_sig", "n_clusterCI_excl0"):
        put(f"{q}_{s}", a.loc[q, s], f"{src}:{q}:{s}", "{:.0f}")
put("v1_d_std_n_holm_sig", a1.loc["d_std_vs_src", "n_holm_sig"], f"{src}[rule=v1]", "{:.0f}")
# --- per-transition extremes / controls
T = pd.read_csv(f"{RV2}/inference/per_transition_inference.csv"); t = T[(T.rule == "textbook")].set_index(["src", "tgt"])
ps = f"{RV2}/inference/per_transition_inference.csv[rule=textbook]"
for (s_, t_), k in ((("FS51", "CIVET"), "fs51_civet"), (("FS53", "CIVET"), "fs53_civet"), (("FS60", "CIVET"), "fs60_civet"), (("ANTS", "CIVET"), "ants_civet"),
                    (("FS53", "FS60"), "fs53_fs60"), (("FS60", "FS53"), "fs60_fs53"), (("FS60DK", "FS60"), "fs60dk_fs60"), (("FS60", "FS60DK"), "fs60_fs60dk"),
                    (("FS51", "FS53"), "fs51_fs53"), (("FS60", "ANTS"), "fs60_ants"), (("CIVET", "FS51"), "civet_fs51")):
    r = t.loc[(s_, t_)]
    for c in ("cov_src", "cov_raw", "cov_std", "d_raw_vs_src", "d_raw_vs_src_clu_lo", "d_raw_vs_src_clu_hi", "d_std_vs_src", "d_std_vs_src_clu_lo", "d_std_vs_src_clu_hi", "mean_same_region_r", "mean_offset_mm"):
        if c in r and pd.notna(r[c]):
            put(f"{k}_{c}", r[c], f"{ps}:{s_}>{t_}:{c}")
# --- categories
Cg = pd.read_csv(f"{RV2}/inference/category_sensitivity.csv")
cmap = {"atlas only (same FS6.0 run, DK vs DKT)": "cat_atlas", "FS version, same DK atlas (6.0 vs 5.3)": "cat_v53", "FS version, same DK atlas (6.0 vs 5.1)": "cat_v51",
        "within FreeSurfer: 5.3<->6.0 (designated control)": "cat_ctrl", "within FreeSurfer: involving 5.1": "cat_fs51", "FreeSurfer<->ANTs": "cat_fsants",
        "FreeSurfer<->CIVET": "cat_fscivet", "ANTs<->CIVET": "cat_antscivet"}
for _, r in Cg.iterrows():
    k = cmap[r.category]; c = "raw" if r.contrast == "raw_vs_src" else "std"
    for s in ("mean_d", "lo", "hi", "min_d", "max_d"):
        put(f"{k}_{c}_{s}", r[s], f"{RV2}/inference/category_sensitivity.csv:{r.category}:{r.contrast}:{s}")
    put(f"{k}_n", r.n_transitions, "category_sensitivity.csv:n_transitions", "{:.0f}"); put(f"{k}_r", r.mean_feature_r, "category_sensitivity.csv:mean_feature_r", "{:.2f}")
# --- widths
W = pd.read_csv(f"{RV2}/frozen/widths_by_transition.csv"); w = nc(W[(W.arm == "main") & (W.rule == "textbook")])
put("width_src", w.width_src.mean(), f"{RV2}/frozen/widths_by_transition.csv:width_src:mean(18,textbook)", "{:.1f}")
put("width_site", w.width_site.mean(), f"{RV2}/frozen/widths_by_transition.csv:width_site:mean(18,textbook)", "{:.1f}")
put("width_site_pct", 100 * (w.width_site.mean() / w.width_src.mean() - 1), "derived width_site/width_src", "{:.0f}")
# --- MAE (v1 per-transition table; identical point predictions under either rank)
tm = nc(pd.read_csv(f"{RV1}/transitions_main.csv"))
for c in ("MAE_src", "MAE_raw", "MAE_std", "MAE_sitestd"):
    put(f"{c}_mean", tm[c].mean(), f"{RV1}/transitions_main.csv:{c}:mean(18)", "{:.2f}"); put(f"{c}_min", tm[c].min(), f"{RV1}/transitions_main.csv:{c}:min(18)", "{:.1f}"); put(f"{c}_max", tm[c].max(), f"{RV1}/transitions_main.csv:{c}:max(18)", "{:.1f}")
put("wcp_frac_inf", tm.wcp_frac_inf.mean(), f"{RV1}/transitions_main.csv:wcp_frac_inf:mean(18)", "{:.2f}")
put("wcp_reg_frac_inf", tm.wcp_reg_frac_inf.mean(), f"{RV1}/transitions_main.csv:wcp_reg_frac_inf:mean(18)", "{:.2f}")
# --- recalibration
S = pd.read_csv(f"{RV2}/recalibration/summary_by_design_k_18noncontrol.csv")
Bm = pd.read_csv(f"{RV2}/recalibration/base_model_by_design_18noncontrol.csv").set_index("design")
rs = f"{RV2}/recalibration/summary_by_design_k_18noncontrol.csv"
for k in (10, 20, 40, 80, 160):
    for dn, pre in (("D5_v2split_tb", "rc"), ("D0_v1_reproduction", "v1rc")):
        for var, vp in (("reform_raw", ""), ("std_then_reform", "h")):
            r = S[(S.design == dn) & (S.variant == var) & (S.k == k)].iloc[0]
            put(f"{pre}{vp}_cov_{k}", r.mean_cov, f"{rs}:{dn}:{var}:k={k}:mean_cov")
            if np.isfinite(r.mean_width):
                put(f"{pre}{vp}_w_{k}", r.mean_width, f"{rs}:{dn}:{var}:k={k}:mean_width", "{:.1f}")
            put(f"{pre}{vp}_min_{k}", r.min_cov, f"{rs}:{dn}:{var}:k={k}:min_cov"); put(f"{pre}{vp}_max_{k}", r.max_cov, f"{rs}:{dn}:{var}:k={k}:max_cov")
            put(f"{pre}{vp}_theory_{k}", r.expected_cov_theory, f"{rs}:{dn}:{var}:k={k}:expected_cov_theory")
oc = pd.read_csv(f"{RV2}/recalibration/original_vs_corrected.csv").set_index("k")
for k in (20, 40, 80, 160):
    put(f"rc_pool_eff_{k}", oc.loc[k, "effect_pool_cleaning"], f"{RV2}/recalibration/original_vs_corrected.csv:effect_pool_cleaning:k={k}", "{:+.3f}")
    put(f"rc_rule_eff_{k}", oc.loc[k, "effect_quantile_rule"], f"{RV2}/recalibration/original_vs_corrected.csv:effect_quantile_rule:k={k}", "{:+.3f}")
put("rc_width_src", Bm.loc["D5_v2split_tb", "width_src"], f"{RV2}/recalibration/base_model_by_design_18noncontrol.csv:D5:width_src", "{:.1f}")
put("rc_mae_src", Bm.loc["D5_v2split_tb", "MAE_src"], "base_model_by_design_18noncontrol.csv:D5:MAE_src", "{:.2f}")
put("rc_n_train", Bm.loc["D5_v2split_tb", "n_train_mean"], "base_model_by_design_18noncontrol.csv:D5:n_train_mean", "{:.0f}")
put("rc_n_cal", Bm.loc["D5_v2split_tb", "n_cal_mean"], "base_model_by_design_18noncontrol.csv:D5:n_cal_mean", "{:.0f}")
put("rc_n_pool", Bm.loc["D5_v2split_tb", "n_pool_mean"], "base_model_by_design_18noncontrol.csv:D5:n_pool_mean", "{:.0f}")
put("rc_n_pool_min", Bm.loc["D5_v2split_tb", "n_pool_min"], "base_model_by_design_18noncontrol.csv:D5:n_pool_min", "{:.0f}")
put("rc_n_test", Bm.loc["D5_v2split_tb", "n_test_mean"], "base_model_by_design_18noncontrol.csv:D5:n_test_mean", "{:.0f}")
R = pd.read_csv(f"{RV2}/recalibration/recal_by_design_transition_k.csv"); r40 = nc(R[(R.design == "D5_v2split_tb") & (R.k == 40)])
B = pd.read_csv(f"{RV2}/recalibration/base_model_by_design.csv"); b5 = B[B.design == "D5_v2split_tb"].set_index(["src", "tgt"]).width_src
for var, vp in (("reform_raw", ""), ("std_then_reform", "h")):
    g = r40[r40.variant == var].set_index(["src", "tgt"]); ratio = g.mean_width_finite / b5.loc[g.index]
    put(f"rc{vp}_ratio40_median", ratio.median(), "recal_by_design_transition_k.csv / base_model_by_design.csv (D5, k=40)", "{:.2f}")
    put(f"rc{vp}_ratio40_min", ratio.min(), "same", "{:.2f}"); put(f"rc{vp}_ratio40_max", ratio.max(), "same", "{:.2f}")
# --- standardisation sensitivity
SS = pd.read_csv(f"{RV2}/standardisation_monitor/std_sensitivity_aggregate.csv", index_col=0).value
for k in ("cov_std_disjoint_mean", "d_disjoint_minus_std_mean", "d_disjoint_minus_std_lo", "d_disjoint_minus_std_hi", "cov_std_disjoint_n_below_088"):
    put(f"ss_{k}", SS[k], f"{RV2}/standardisation_monitor/std_sensitivity_aggregate.csv:{k}", "{:.0f}" if "n_below" in k else ("{:+.4f}" if k.startswith("d_") else "{:.3f}"))
# --- monitors
MS = pd.read_csv(f"{RV2}/standardisation_monitor/monitor_summary.csv", index_col=0).value
for k, f in (("mmd_std_detect_mean_nc", "{:.2f}"), ("mmd_std_detect_median_nc", "{:.2f}"), ("ks_std_detect_mean_nc", "{:.2f}"), ("far_mmd_raw_mean", "{:.3f}"),
             ("far_mmd_std_mean", "{:.3f}"), ("n_harmful_ci", "{:.0f}"), ("n_harmful_ci_missed_by_mmd_majority", "{:.0f}"), ("spearman_mmd_std_detect_vs_d_std", "{:.2f}"), ("spearman_p", "{:.2f}")):
    put(f"mon_{k}", float(MS[k]), f"{RV2}/standardisation_monitor/monitor_summary.csv:{k}", f)
M = pd.read_csv(f"{RV2}/standardisation_monitor/monitor_diagnostic_by_transition.csv")
put("mon_far_max", M[[c for c in M if c.endswith("far_test")]].max().max(), "monitor_diagnostic_by_transition.csv:max far_test", "{:.2f}")
put("mon_std_detect_max", nc(M).mmd_std_detect.max(), "monitor_diagnostic_by_transition.csv:max mmd_std_detect(18)", "{:.2f}")
put("mon_ctrl_ks_raw_max", M[M.control].ks_raw_detect.max(), "monitor_diagnostic_by_transition.csv:control ks_raw_detect max", "{:.3f}")
mfa = M.set_index(["src", "tgt"])
put("mon_fs60_civet_detect", mfa.loc[("FS60", "CIVET"), "mmd_std_detect"], "monitor_diagnostic_by_transition.csv:FS60>CIVET:mmd_std_detect", "{:.2f}")
for b in (25, 100):
    mb = nc(pd.read_csv(f"{RV1}/monitors_b{b}.csv")); put(f"mon_b{b}_std", mb.mmd_std_detect.mean(), f"{RV1}/monitors_b{b}.csv:mmd_std_detect:mean(18)", "{:.2f}")
    put(f"mon_b{b}_std_max", mb.mmd_std_detect.max(), f"{RV1}/monitors_b{b}.csv:mmd_std_detect:max(18)", "{:.2f}")
# --- age
AG = pd.read_csv(f"{RV2}/age/age_bin_coverage.csv"); DE = pd.read_csv(f"{RV2}/age/bin_description.csv").set_index("bin")
bk = {"<12": "lt12", "12-15": "b1215", "15-18": "b1518", "18-25": "b1825", ">=25": "ge25"}
for _, r in AG.iterrows():
    k = bk[r.bin]; f = "{:+.1f}" if "error" in r.quantity else "{:.3f}"
    for s, col in (("", "mean_over_transitions"), ("_lo", "lo"), ("_hi", "hi"), ("_min", "min_transition"), ("_max", "max_transition")):
        put(f"age_{k}_{r.quantity}{s}", r[col], f"{RV2}/age/age_bin_coverage.csv:{r.bin}:{r.quantity}:{col}", f)
for b_, k in bk.items():
    put(f"age_{k}_n", DE.loc[b_, "n_unique_subjects"], f"{RV2}/age/bin_description.csv:{b_}:n_unique_subjects", "{:.0f}")
    put(f"age_{k}_sites", DE.loc[b_, "n_sites"], f"{RV2}/age/bin_description.csv:{b_}:n_sites", "{:.0f}")
    put(f"age_{k}_npt", DE.loc[b_, "n_per_transition_median"], f"{RV2}/age/bin_description.csv:{b_}:n_per_transition_median", "{:.0f}")
put("age_ge25_agemax", DE.loc[">=25", "age_max"], "bin_description.csv:>=25:age_max", "{:.0f}")
DS = pd.read_csv(f"{RV2}/age/decomposition_summary.csv", index_col=0)
for c, f in (("share_ge25", "{:.1%}"), ("contribution_ge25", "{:+.3f}"), ("d_std_lt25", "{:+.3f}"), ("cov_src_lt25", "{:.3f}"), ("d_raw_lt25", "{:+.3f}")):
    put(f"dec_{c}_mean", DS.loc["mean", c], f"{RV2}/age/decomposition_summary.csv:mean:{c}", f); put(f"dec_{c}_median", DS.loc["median", c], f"{RV2}/age/decomposition_summary.csv:median:{c}", f)
dbt = pd.read_csv(f"{RV2}/age/decomposition_by_transition.csv"); put("dec_n_lt25_ge5", (dbt.d_std_lt25 <= -0.05).sum(), "decomposition_by_transition.csv:count d_std_lt25<=-0.05", "{:.0f}")
# --- secondary (textbook rank re-execution)
SB = f"{RV2}/secondary_textbook"
for f_ in ("qc0", "qc2"):
    q = nc(pd.read_csv(f"{SB}/transitions_{f_}.csv"))
    put(f"{f_}_std_min", q.cov_std.min(), f"{SB}/transitions_{f_}.csv:cov_std:min(18)"); put(f"{f_}_std_max", q.cov_std.max(), f"{SB}/transitions_{f_}.csv:cov_std:max(18)")
    put(f"{f_}_d_std_median", q.d_std.median(), f"{SB}/transitions_{f_}.csv:d_std:median(18)", "{:+.3f}"); put(f"{f_}_src_mean", q.cov_src.mean(), f"{SB}/transitions_{f_}.csv:cov_src:mean(18)")
mx = pd.read_csv(f"{SB}/mixed_pipeline.csv").set_index("target")
put("mix_N", mx.N.iloc[0], f"{SB}/mixed_pipeline.csv:N", "{:.0f}"); put("mix_all_min", mx.mix_all.min(), f"{SB}/mixed_pipeline.csv:mix_all:min"); put("mix_all_max", mx.mix_all.max(), f"{SB}/mixed_pipeline.csv:mix_all:max")
put("mix_within_min", mx.within.min(), f"{SB}/mixed_pipeline.csv:within:min"); put("mix_within_max", mx.within.max(), f"{SB}/mixed_pipeline.csv:within:max")
put("mix_width_all", mx.width_all.mean(), f"{SB}/mixed_pipeline.csv:width_all", "{:.1f}")
for tg in ("FS51", "FS53", "FS60", "ANTS", "CIVET"):
    put(f"mix_loo_{tg}", mx.loc[tg, "mix_loo"], f"{SB}/mixed_pipeline.csv:{tg}:mix_loo")
put("mix_mae_loo_CIVET", mx.loc["CIVET", "MAE_loo"], f"{SB}/mixed_pipeline.csv:CIVET:MAE_loo", "{:.1f}")
lp = nc(pd.read_csv(f"{SB}/loso_pooled.csv"))
for c in ("cov_src", "cov_raw", "cov_std", "cov_src_site", "cov_sitestd", "min_site_cov_std"):
    put(f"loso_{c}", lp[c].mean(), f"{SB}/loso_pooled.csv:{c}:mean(18)")
lg = nc(pd.read_csv(f"{RV1}/lgbm.csv"))
for c in ("lgbm_cov_src", "lgbm_cov_raw", "lgbm_cov_std"):
    put(f"v1_{c}", lg[c].mean(), f"{RV1}/lgbm.csv:{c}:mean(18) [v1 rank; not re-executed]")
put("v1_lgbm_std_n_below_088", (lg.lgbm_cov_std < 0.88).sum(), f"{RV1}/lgbm.csv:count lgbm_cov_std<0.88 (18)", "{:.0f}")
json.dump(N, open(OUT, "w"), indent=1)
print(len(N), "numbers written")
