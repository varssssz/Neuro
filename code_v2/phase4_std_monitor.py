"""
Phase 4: standardisation sensitivity + drift-monitor validation.
Usage: python -I code_v2/phase4_std_monitor.py RESULTS_V2 ORIGINAL_V1_DIR DATA_DIR
- Standardisation: uses results_v2/frozen (textbook rank). Conditions:
    std          target stats from Xt of the model-training subjects (v1 definition; n ~ 500; SAME people as training)
    std_disjoint target stats from Xt of the source-calibration subjects (n ~ 167; DIFFERENT people from training)
    std_same_sub target stats from a random subset of training subjects of the same size as calibration (same people, n ~ 167)
  => std_same_sub vs std isolates sample size; std_disjoint vs std_same_sub isolates "same people" (population matching).
- Monitors: re-executes ORIGINAL v1 run_monitor (unmodified) for all 20 transitions at batch 50 and compares with stored output;
  then relates alarm rate to the coverage loss under per-pipeline standardisation.
"""
import sys, os
import numpy as np
import pandas as pd
from scipy import stats

RV2, V1, DATA = sys.argv[1:4]
OUT = f"{RV2}/standardisation_monitor"; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, f"{V1}/Source_Code")
import abide_data as AD
import swap_engine as E1          # original v1 engine, unmodified

CONTROL = {("FS60", "FS53"), ("FS53", "FS60")}
d = pd.read_csv(f"{RV2}/frozen/per_subject_seed.csv.gz")
d = d[(d.arm == "main") & (d.rule == "textbook")]
rng = np.random.default_rng(99)
sites_all = np.array(sorted(d.site.unique())); pos = {s: i for i, s in enumerate(sites_all)}
BM = rng.multinomial(len(sites_all), np.ones(len(sites_all)) / len(sites_all), size=2000)

rows, boots = [], {}
for (s, t), g in d.groupby(["src", "tgt"], sort=False):
    sm = g.groupby("sid").agg(site=("site", "first"), **{c: (c, "mean") for c in ["cov_src", "cov_std", "cov_std_disjoint", "cov_std_same_sub"]},
                              **{f"ae_{k}": (f"err_{k}", lambda x: np.abs(x).mean()) for k in ["src", "std", "std_disjoint", "std_same_sub"]})
    r = dict(src=s, tgt=t, control=(s, t) in CONTROL, N=len(sm))
    for k in ("src", "std", "std_disjoint", "std_same_sub"):
        r[f"cov_{k}"] = sm[f"cov_{k}"].mean(); r[f"MAE_{k}"] = sm[f"ae_{k}"].mean()
    r["width"] = g.width.mean()
    for nm, (a, b) in {"disjoint_minus_std": ("cov_std_disjoint", "cov_std"), "samesub_minus_std": ("cov_std_same_sub", "cov_std"),
                       "disjoint_minus_samesub": ("cov_std_disjoint", "cov_std_same_sub"), "disjoint_minus_src": ("cov_std_disjoint", "cov_src")}.items():
        x = (sm[a] - sm[b]).values
        gs = pd.DataFrame(dict(x=x, s=sm.site.values)).groupby("s").x.agg(["sum", "count"])
        idx = np.array([pos[q] for q in gs.index]); w = BM[:, idx]
        bt = (w * gs["sum"].values).sum(1) / (w * gs["count"].values).sum(1)
        r[f"d_{nm}"] = x.mean(); r[f"d_{nm}_lo"], r[f"d_{nm}_hi"] = np.nanpercentile(bt, [2.5, 97.5])
        boots[(s, t, nm)] = bt
    rows.append(r)
S = pd.DataFrame(rows); S.to_csv(f"{OUT}/std_sensitivity_by_transition.csv", index=False)
nc = S[~S.control]
agg = dict(n=len(nc))
for k in ("src", "std", "std_disjoint", "std_same_sub"):
    agg[f"cov_{k}_mean"] = nc[f"cov_{k}"].mean(); agg[f"cov_{k}_median"] = nc[f"cov_{k}"].median()
    agg[f"cov_{k}_min"] = nc[f"cov_{k}"].min(); agg[f"cov_{k}_n_below_088"] = int((nc[f"cov_{k}"] < 0.88).sum()); agg[f"MAE_{k}_mean"] = nc[f"MAE_{k}"].mean()
for nm in ("disjoint_minus_std", "samesub_minus_std", "disjoint_minus_samesub", "disjoint_minus_src"):
    E_ = np.vstack([boots[(s, t, nm)] for s, t in zip(nc.src, nc.tgt)]).mean(0)
    agg[f"d_{nm}_mean"] = nc[f"d_{nm}"].mean(); agg[f"d_{nm}_lo"], agg[f"d_{nm}_hi"] = np.percentile(E_, [2.5, 97.5])
    agg[f"d_{nm}_min"], agg[f"d_{nm}_max"] = nc[f"d_{nm}"].min(), nc[f"d_{nm}"].max()
pd.Series(agg).to_csv(f"{OUT}/std_sensitivity_aggregate.csv", header=["value"])

# ------------------------------------------------------------------ monitors: re-execute v1 for all 20 at b=50
P, ph, common, qc = AD.build(DATA)
stored = pd.read_csv(f"{V1}/Experimental_Results/monitors_b50.csv").set_index(["src", "tgt"])
inf = pd.read_csv(f"{RV2}/inference/per_transition_table_textbook.csv").set_index(["src", "tgt"])
mrows = []
for s, t in [(a, b) for a in AD.PIPES5 for b in AD.PIPES5 if a != b]:
    ids = AD.matched_ids(P, ph, qc, s, t, 1)
    Xs, Xt, y = P[s].loc[ids].values, P[t].loc[ids].values, ph.loc[ids, "AGE_AT_SCAN"].values
    m = E1.run_monitor(Xs, Xt, y, seed=0, bsz=50)
    r = dict(src=s, tgt=t, control=(s, t) in CONTROL, N=len(ids), test_fold_mean=len(ids) / 5, batches_per_condition=len(m["raw"]["mmd_alt"]))
    for cond in ("raw", "std"):
        for st in ("mmd", "ks"):
            thr = np.quantile(m[cond][f"{st}_cal"], 0.95)
            r[f"{st}_{cond}_threshold"] = thr
            r[f"{st}_{cond}_far_test"] = float((np.array(m[cond][f"{st}_nul"]) > thr).mean())
            r[f"{st}_{cond}_detect"] = float((np.array(m[cond][f"{st}_alt"]) > thr).mean())
            r[f"{st}_{cond}_detect_stored"] = stored.loc[(s, t), f"{st}_{cond}_detect"]
            r[f"{st}_{cond}_alt_median_stat_over_thr"] = float(np.median(m[cond][f"{st}_alt"]) / thr) if thr > 0 else np.nan
    for k in ("cov_src", "cov_raw", "cov_std", "d_raw_vs_src", "d_std_vs_src", "d_std_vs_src_clu_lo", "d_std_vs_src_clu_hi", "p_std_vs_src_holm"):
        r[k] = inf.loc[(s, t), k]
    r["harmful_after_std"] = bool(r["d_std_vs_src"] <= -0.05)
    r["harmful_after_std_ci"] = bool(r["d_std_vs_src_clu_hi"] < 0)
    r["mmd_std_alarm_majority"] = bool(r["mmd_std_detect"] >= 0.5)
    mrows.append(r); print(s, t, "monitor done", flush=True)
M = pd.DataFrame(mrows); M.to_csv(f"{OUT}/monitor_diagnostic_by_transition.csv", index=False)
mnc = M[~M.control]
summ = dict(max_abs_diff_detect_vs_stored=float(max((M[f"{st}_{c}_detect"] - M[f"{st}_{c}_detect_stored"]).abs().max() for st in ("mmd", "ks") for c in ("raw", "std"))),
            mmd_raw_detect_mean_nc=mnc.mmd_raw_detect.mean(), mmd_std_detect_mean_nc=mnc.mmd_std_detect.mean(), mmd_std_detect_median_nc=mnc.mmd_std_detect.median(),
            ks_std_detect_mean_nc=mnc.ks_std_detect.mean(), far_mmd_raw_mean=M.mmd_raw_far_test.mean(), far_mmd_std_mean=M.mmd_std_far_test.mean(),
            far_ks_std_mean=M.ks_std_far_test.mean(), far_max_any=float(M[[c for c in M if c.endswith("far_test")]].max().max()),
            n_harmful_after_std=int(mnc.harmful_after_std.sum()), n_harmful_ci=int(mnc.harmful_after_std_ci.sum()),
            n_harmful_ci_missed_by_mmd_majority=int((mnc.harmful_after_std_ci & ~mnc.mmd_std_alarm_majority).sum()),
            n_harmful_ci_missed_by_ks_majority=int((mnc.harmful_after_std_ci & ~(mnc.ks_std_detect >= 0.5)).sum()),
            spearman_mmd_std_detect_vs_d_std=stats.spearmanr(mnc.mmd_std_detect, mnc.d_std_vs_src)[0],
            spearman_p=stats.spearmanr(mnc.mmd_std_detect, mnc.d_std_vs_src)[1],
            controls_mmd_raw_detect=M[M.control].mmd_raw_detect.tolist(), controls_ks_raw_detect=M[M.control].ks_raw_detect.tolist())
pd.Series(summ).to_csv(f"{OUT}/monitor_summary.csv", header=["value"])

# ------------------------------------------------------------------ figure
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), sharey=False)
for a, cond, title in ((ax[0], "raw", "Raw target features"), (ax[1], "std", "After per-pipeline standardisation")):
    dd = M[f"d_{cond}_vs_src"]
    a.axhline(0.5, color="#999", lw=0.8, ls=":"); a.axvline(-0.05, color="#999", lw=0.8, ls=":")
    a.scatter(dd[~M.control], M.loc[~M.control, f"mmd_{cond}_detect"], s=36, color="#2b6cb0", label="non-control transition")
    a.scatter(dd[M.control], M.loc[M.control, f"mmd_{cond}_detect"], s=40, marker="s", facecolor="white", edgecolor="#c05621", label="FS5.3↔FS6.0")
    a.set_xlabel("Coverage change vs source pipeline"); a.set_title(title, fontsize=10); a.set_ylim(-0.04, 1.04)
    far = M[f"mmd_{cond}_far_test"].mean(); a.axhline(far, color="#c05621", lw=0.8, ls="--"); a.text(a.get_xlim()[0], far + 0.02, f" mean false-alarm rate {far:.2f}", fontsize=7, color="#c05621")
ax[0].set_ylabel("MMD alarm rate (batches of 50, 5% threshold)"); ax[0].legend(fontsize=8, loc="center right", frameon=False)
fig.tight_layout(); fig.savefig(f"{OUT}/fig_monitor_alarm_vs_coverage_loss.png", dpi=350); plt.close(fig)
pd.set_option("display.width", 250)
print(pd.Series(agg).round(4).to_string()); print(pd.Series(summ).to_string())
print(M[["src", "tgt", "d_std_vs_src", "d_std_vs_src_clu_hi", "mmd_std_detect", "ks_std_detect", "mmd_std_far_test", "harmful_after_std_ci"]].round(3).to_string())
