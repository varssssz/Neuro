"""Phase 2 summary tables. Usage: python -I code_v2/summarise_recalibration.py results_v2/recalibration"""
import sys, numpy as np, pandas as pd
D = sys.argv[1]
R = pd.read_csv(f"{D}/recal_by_design_transition_k.csv"); B = pd.read_csv(f"{D}/base_model_by_design.csv")
nc = R[~R.control]
rows = []
for (dn, var, k), g in nc.groupby(["design", "variant", "k"]):
    g = g[g.n_draws > 0]
    rows.append(dict(design=dn, variant=var, k=k, n_transitions=len(g), mean_cov=g.mean_cov.mean(), median_cov=g.mean_cov.median(),
                     min_cov=g.mean_cov.min(), max_cov=g.mean_cov.max(), expected_cov_theory=g.expected_cov_theory.iloc[0] if len(g) else np.nan,
                     mean_frac_draws_below_088=g.frac_below_088.mean(), frac_inf_width=g.frac_inf_width.mean(),
                     mean_width=g.mean_width_finite.mean(), median_width=g.median_width.median()))
S = pd.DataFrame(rows); S.to_csv(f"{D}/summary_by_design_k_18noncontrol.csv", index=False)
b = B[~B.control].groupby("design")[["cov_src", "cov_raw", "MAE_src", "MAE_raw", "width_src", "n_train_mean", "n_cal_mean", "n_pool_mean", "n_pool_min", "n_test_mean", "pool_test_overlap_max"]].mean()
b.to_csv(f"{D}/base_model_by_design_18noncontrol.csv")
# original-vs-corrected headline table (reform_raw)
s = S[S.variant == "reform_raw"].set_index(["design", "k"])
comp = []
for k in (10, 20, 40, 80, 160):
    v1, q, pool, rep, fin = (s.loc[(d, k)] for d in ("D0_v1_reproduction", "D1_v1split_trca_tb", "D3_v1split_ca_tb", "D5_v2split_tb", "D5_v2split_tb"))
    comp.append(dict(k=k, v1_cov=v1.mean_cov, v1_width=v1.mean_width, corrected_cov=fin.mean_cov, corrected_cov_min=fin.min_cov,
                     corrected_cov_max=fin.max_cov, corrected_theory=fin.expected_cov_theory, corrected_width_mean=fin.mean_width,
                     corrected_width_median=fin.median_width, corrected_frac_inf=fin.frac_inf_width,
                     effect_quantile_rule=q.mean_cov - v1.mean_cov, effect_pool_cleaning=pool.mean_cov - q.mean_cov,
                     effect_repartition=rep.mean_cov - pool.mean_cov if pool.n_transitions == 18 else np.nan,
                     width_effect_quantile_rule=q.mean_width - v1.mean_width))
C = pd.DataFrame(comp); C.to_csv(f"{D}/original_vs_corrected.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(C.round(4).to_string()); print(b.round(3).to_string())
print(S[(S.design == "D5_v2split_tb")].round(3).to_string())
print(S[S.design.isin(["D2_v1split_ca_v1", "D3_v1split_ca_tb"]) & (S.k == 160)][["design", "n_transitions"]])
