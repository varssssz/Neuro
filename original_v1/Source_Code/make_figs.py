"""Build figures and summary tables from run_all.py outputs. Usage: python make_figs.py RESULTS_DIR OUT_DIR"""
import sys, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R, O = sys.argv[1], sys.argv[2]
os.makedirs(O, exist_ok=True)
C = dict(src="#444444", raw="#D55E00", std="#0072B2", rec="#009E73", site="#CC79A7")
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})

T = pd.read_csv(f"{R}/transitions_main.csv")
RC = pd.read_csv(f"{R}/recal_main.csv")
S = pd.read_csv(f"{R}/strata_main.csv")
T["lab"] = T.src + "→" + T.tgt
T["pair_ctrl"] = T.control.astype(bool)

# ---------- Fig 1: coverage by transition ----------
order = T.sort_values("cov_raw").reset_index(drop=True)
rk = RC[(RC.variant == "reform_raw") & (RC.k == 40)].set_index(["src", "tgt"]).mean_cov
fig, ax = plt.subplots(figsize=(6.4, 6.2))
y = np.arange(len(order))
ax.axvline(0.90, color="k", lw=0.8, ls="--")
for i, r in order.iterrows():
    ax.plot([r.cov_raw_lo, r.cov_raw_hi], [i, i], color=C["raw"], lw=1, alpha=.5)
ax.scatter(order.cov_raw, y, color=C["raw"], label="Frozen model, raw target features", zorder=3, s=22)
ax.scatter(order.cov_std, y, color=C["std"], label="Per-pipeline standardization (label-free)", zorder=3, s=22, marker="s")
ax.scatter([rk[(r.src, r.tgt)] for r in order.itertuples()], y, color=C["rec"], label="ReForm-style, k = 40 labelled", zorder=3, s=26, marker="^")
ax.scatter(order.cov_src, y, color=C["src"], label="Source-pipeline test (reference)", zorder=2, s=14, marker="|")
ax.set_yticks(y)
ax.set_yticklabels([l + (" *" if c else "") for l, c in zip(order.lab, order.pair_ctrl)])
ax.set_xlabel("Empirical coverage (nominal 0.90)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5,-0.09), ncol=2, fontsize=7, frameon=False)
ax.set_title("Coverage after a same-scan pipeline swap (* near-identity control)", fontsize=9)
fig.tight_layout(); fig.savefig(f"{O}/fig1_coverage_by_transition.png"); plt.close(fig)

# ---------- Fig 2: coverage vs k ----------
nc = RC.merge(T[["src", "tgt", "control"]], on=["src", "tgt"])
nc = nc[~nc.control.astype(bool)]
g = nc.groupby(["variant", "k"]).agg(cov=("mean_cov", "mean"), lo=("cov_lo", "mean"), hi=("cov_hi", "mean"),
                                     fb=("frac_below_088", "mean"), w=("mean_width", "mean")).reset_index()
fig, ax = plt.subplots(1, 2, figsize=(7, 3.0))
for v, c, lab in [("reform_raw", C["rec"], "ReForm-style on raw features"), ("std_then_reform", C["std"], "Standardize, then ReForm-style")]:
    d = g[g.variant == v]
    ax[0].plot(d.k, d["cov"], "o-", color=c, label=lab)
    ax[0].fill_between(d.k, d.lo, d.hi, color=c, alpha=.12)
    ax[1].plot(d.k, d.fb, "o-", color=c)
ax[0].axhline(.9, color="k", ls="--", lw=.8)
ax[0].set_xscale("log"); ax[1].set_xscale("log")
ax[0].set_xticks([10, 20, 40, 80, 160]); ax[0].set_xticklabels([10, 20, 40, 80, 160])
ax[1].set_xticks([10, 20, 40, 80, 160]); ax[1].set_xticklabels([10, 20, 40, 80, 160])
ax[0].set_xlabel("Labelled target subjects k"); ax[1].set_xlabel("Labelled target subjects k")
ax[0].set_ylabel("Mean coverage (18 non-control transitions)"); ax[1].set_ylabel("Share of draws with coverage < 0.88")
ax[0].legend(fontsize=7.5, frameon=False, loc="lower right")
fig.tight_layout(); fig.savefig(f"{O}/fig2_coverage_vs_k.png"); plt.close(fig)

# ---------- Fig 3: age-binned coverage ----------
bins = ["<12", "12-15", "15-18", "18-25", ">=25"]
s = S.merge(T[["src", "tgt", "control"]], on=["src", "tgt"])
s = s[~s.control.astype(bool)]
fig, ax = plt.subplots(figsize=(5.2, 3.2))
off = {"src": -.27, "raw": -.09, "std": .09, "sitestd": .27}
nm = {"src": "Source test", "raw": "Raw target", "std": "Standardized", "sitestd": "Per-site standardized"}
col = {"src": C["src"], "raw": C["raw"], "std": C["std"], "sitestd": C["site"]}
present = [b for b in bins if b in set(s.stratum)] or sorted(set(s.stratum))
for cnd in ["src", "raw", "std", "sitestd"]:
    d = s[s.cond == cnd]
    vals = []
    for b in present:
        x = d[d.stratum == b]
        vals.append(np.average(x["cov"], weights=x.n) if len(x) else np.nan)
    ax.bar(np.arange(len(present)) + off[cnd], vals, width=.17, color=col[cnd], label=nm[cnd])
ax.axhline(.9, color="k", ls="--", lw=.8)
ax.set_xticks(range(len(present))); ax.set_xticklabels(present)
ax.set_xlabel("Age at scan (years)"); ax.set_ylabel("Coverage (pooled, 18 transitions)")
ax.legend(fontsize=7, frameon=False, ncol=2, loc="lower left")
fig.tight_layout(); fig.savefig(f"{O}/fig3_age_binned_coverage.png"); plt.close(fig)

# ---------- Fig 4: width vs coverage ----------
fig, ax = plt.subplots(figsize=(5.2, 3.8))
nT = T[~T.pair_ctrl]
ax.scatter(nT.width_src, nT.cov_std, color=C["std"], s=18, label="Standardized (width = source)")
ax.scatter(nT.width_src_site, nT.cov_sitestd, color=C["site"], s=18, marker="s", label="Per-site standardized")
for k, mk in zip([20, 40, 80], ["v", "^", "D"]):
    d = nc[(nc.variant == "reform_raw") & (nc.k == k)]
    ax.scatter(d.mean_width, d.mean_cov, color=C["rec"], s=16, marker=mk, alpha=.7, label=f"ReForm-style k={k}")
ax.axhline(.9, color="k", ls="--", lw=.8)
ax.set_xlabel("Mean interval width (years)"); ax.set_ylabel("Coverage")
ax.legend(fontsize=7, frameon=False, loc="lower right")
fig.tight_layout(); fig.savefig(f"{O}/fig4_width_vs_coverage.png"); plt.close(fig)

# ---------- Tables ----------
nc40 = nc[(nc.variant == "reform_raw")]
tab = T[["src", "tgt", "N", "cov_src", "cov_raw", "cov_std", "cov_sitestd", "d_raw", "d_raw_lo", "d_raw_hi",
         "d_std", "d_std_lo", "d_std_hi", "mcnemar_p_raw", "mcnemar_p_std", "width_src", "MAE_src", "MAE_raw", "MAE_std", "control"]]
tab.round(4).to_csv(f"{O}/table_transitions.csv", index=False)
g.round(4).to_csv(f"{O}/table_recal_summary.csv", index=False)
print(g.round(3).to_string())
grp = T.assign(kind=np.where(T.pair_ctrl, "control", np.where((T.src == "CIVET") | (T.tgt == "CIVET"), "CIVET", np.where(T.tgt == "ANTS", "to ANTs", np.where((T.src == "ANTS"), "from ANTs", "FS-family")))))
print(grp.groupby("kind")[["cov_src", "cov_raw", "cov_std", "cov_sitestd", "d_raw", "d_std"]].agg(["mean", "min"]).round(3).to_string())
