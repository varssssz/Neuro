"""Revised manuscript figures from results_v2 (textbook rank). Usage: python -I code_v2/make_figs_v2.py RESULTS_V2 OUT_DIR"""
import sys, os, shutil
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

RV2, O = sys.argv[1:3]; os.makedirs(O, exist_ok=True)
C = dict(src="#2a78d6", raw="#eb6834", std="#1baf7a", site="#eda100", ink="#2b2b2b", muted="#8a8a85", grid="#e4e4e0")
plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": C["muted"], "axes.labelcolor": C["ink"], "xtick.color": C["ink"], "ytick.color": C["ink"],
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 350, "font.family": "DejaVu Sans"})
NAME = {"FS51": "FS5.1", "FS53": "FS5.3", "FS60": "FS6.0", "ANTS": "ANTs", "CIVET": "CIVET"}

# ---------------- Figure 1: per-transition coverage
T = pd.read_csv(f"{RV2}/inference/per_transition_inference.csv")
t = T[(T.rule == "textbook") & (T.arm == "main")].copy()
order = ["within FreeSurfer: 5.3<->6.0 (designated control)", "within FreeSurfer: involving 5.1", "FreeSurfer<->ANTs", "FreeSurfer<->CIVET", "ANTs<->CIVET"]
lab = {order[0]: "FS5.3↔6.0 (low-impact pair)", order[1]: "Within FreeSurfer (5.1)", order[2]: "FreeSurfer↔ANTs", order[3]: "FreeSurfer↔CIVET", order[4]: "ANTs↔CIVET"}
t["cat_i"] = t.category.map({c: i for i, c in enumerate(order)})
t = t.sort_values(["cat_i", "cov_std"]).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(6.6, 6.2))
ypos, y, ticks, groups = [], 0, [], []
for i, r in t.iterrows():
    if i > 0 and r.cat_i != t.cat_i[i - 1]:
        y += 0.8
    ypos.append(y); y += 1
t["y"] = ypos
for c in order:
    g = t[t.category == c]; groups.append((g.y.mean(), lab[c]))
ax.axvline(0.9, color=C["muted"], lw=1, ls="--", zorder=0)
ax.text(0.9, -1.4, " nominal 0.90", color=C["muted"], fontsize=7.5, va="center")
for _, r in t.iterrows():
    ax.plot([r.cov_raw, r.cov_src], [r.y, r.y], color=C["grid"], lw=1, zorder=1)
ax.errorbar(t.cov_std, t.y, xerr=[t.cov_std - t.cov_std_clu_lo, t.cov_std_clu_hi - t.cov_std], fmt="none", ecolor=C["std"], elinewidth=1.4, capsize=0, zorder=2)
ax.scatter(t.cov_raw, t.y, marker="X", s=42, color=C["raw"], edgecolor="white", linewidth=0.6, zorder=3, label="Raw target features (frozen model)")
ax.scatter(t.cov_std, t.y, marker="D", s=34, color=C["std"], edgecolor="white", linewidth=0.6, zorder=4, label="Per-pipeline standardisation (95% site-cluster CI)")
ax.scatter(t.cov_src, t.y, marker="o", s=38, color=C["src"], edgecolor="white", linewidth=0.6, zorder=5, label="Source pipeline")
ax.set_yticks(t.y); ax.set_yticklabels([f"{NAME[a]} → {NAME[b]}" for a, b in zip(t.src, t.tgt)], fontsize=7.5)
ax.invert_yaxis(); ax.set_xlim(-0.02, 1.0); ax.set_xlabel("Empirical coverage of the 90% interval")
ax.grid(axis="x", color=C["grid"], lw=0.6); ax.set_axisbelow(True)
for yy, l in groups:
    ax.text(1.005, yy, l, transform=ax.get_yaxis_transform(), fontsize=7.5, color=C["ink"], va="center")
ax.legend(loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=1, frameon=False, fontsize=7.5)
fig.subplots_adjust(left=0.2, right=0.74, top=0.86, bottom=0.08)
fig.savefig(f"{O}/fig1_coverage_by_transition.png"); plt.close(fig)

# ---------------- Figure 2: recalibration (v1 vs corrected) and width
S = pd.read_csv(f"{RV2}/recalibration/summary_by_design_k_18noncontrol.csv")
B = pd.read_csv(f"{RV2}/recalibration/base_model_by_design_18noncontrol.csv").set_index("design")
ks = [10, 20, 40, 80, 160]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0))
a = axs[0]
for dn, lbl, col, ls, mk in (("D0_v1_reproduction", "v1 (training subjects in pool; rank +1)", C["muted"], "--", "o"),
                             ("D5_v2split_tb", "Corrected (disjoint pool; textbook rank)", C["src"], "-", "o")):
    s = S[(S.design == dn) & (S.variant == "reform_raw")].set_index("k").loc[ks]
    xs = [k for k in ks if np.isfinite(s.loc[k, "mean_width"])]
    a.plot(xs, s.loc[xs, "mean_cov"], ls=ls, color=col, lw=2, marker=mk, ms=6, mfc=col if ls == "-" else "white", label=lbl)
th = S[(S.design == "D5_v2split_tb") & (S.variant == "reform_raw")].set_index("k").loc[ks[1:], "expected_cov_theory"]
a.plot(ks[1:], th, ls=":", color=C["ink"], lw=1, label="Theoretical expected coverage")
a.axhline(0.9, color=C["grid"], lw=1, zorder=0)
a.plot([10], [0.964], marker="^", ms=7, color=C["src"], clip_on=False); a.text(10.9, 0.964, "k = 10, corrected: infinite interval (no finite 90% bound)", fontsize=6.8, color=C["ink"], va="center")
a.set_xscale("log"); a.minorticks_off(); a.set_xticks(ks); a.set_xticklabels(ks); a.set_xlabel("Labelled target subjects (k)"); a.set_ylabel("Mean coverage (18 transitions)")
a.set_ylim(0.79, 0.97); a.legend(fontsize=6.8, frameon=False, loc="lower right"); a.set_title("a  Coverage", loc="left", fontsize=9)
a = axs[1]
for var, lbl, col, mk in (("reform_raw", "Recalibrate raw target features", C["raw"], "o"), ("std_then_reform", "Standardise, then recalibrate", C["std"], "D")):
    s = S[(S.design == "D5_v2split_tb") & (S.variant == var)].set_index("k").loc[ks[1:]]
    a.plot(ks[1:], s.mean_width, color=col, lw=2, marker=mk, ms=6, label=lbl)
a.axhline(B.loc["D5_v2split_tb", "width_src"], color=C["src"], lw=1.2, ls="--")
a.text(165, B.loc["D5_v2split_tb", "width_src"] + 0.2, "source interval", fontsize=7, color=C["src"], ha="right")
a.set_xscale("log"); a.minorticks_off(); a.set_xticks(ks[1:]); a.set_xticklabels(ks[1:]); a.set_xlabel("Labelled target subjects (k)"); a.set_ylabel("Mean interval width (years)")
a.set_ylim(10, 18); a.legend(fontsize=6.8, frameon=False, loc="upper right"); a.set_title("b  Width (corrected design)", loc="left", fontsize=9)
fig.tight_layout(); fig.savefig(f"{O}/fig2_recalibration.png"); plt.close(fig)

# ---------------- Figure 3: age bins
A = pd.read_csv(f"{RV2}/age/age_bin_coverage.csv"); D = pd.read_csv(f"{RV2}/age/bin_description.csv").set_index("bin")
bins = ["<12", "12-15", "15-18", "18-25", ">=25"]
fig, ax = plt.subplots(figsize=(6.4, 3.0))
w = 0.18
for j, (q, lbl, col, mk) in enumerate((("cov_src", "Source pipeline", C["src"], "o"), ("cov_raw", "Raw transfer", C["raw"], "X"),
                                         ("cov_std", "Per-pipeline standardisation", C["std"], "D"), ("cov_sitestd", "Per-site std. (refit)", C["site"], "s"))):
    s = A[A.quantity == q].set_index("bin").loc[bins]; x = np.arange(5) + (j - 1.5) * w
    ax.errorbar(x, s.mean_over_transitions, yerr=[s.mean_over_transitions - s.lo, s.hi - s.mean_over_transitions], fmt=mk, color=col, ms=6, mec="white", mew=0.6, elinewidth=1.2, label=lbl)
ax.axhline(0.9, color=C["muted"], lw=1, ls="--", zorder=0)
ax.set_xticks(range(5)); ax.set_xticklabels([f"{b.replace('>=', '≥')}\n(n={int(D.loc[b, 'n_unique_subjects'])})" for b in bins])
ax.set_ylabel("Coverage (mean of 18 transitions)"); ax.set_xlabel("Age bin, years (unique subjects)"); ax.set_ylim(0, 1.02)
ax.grid(axis="y", color=C["grid"], lw=0.6); ax.set_axisbelow(True)
ax.legend(fontsize=7, frameon=False, ncol=2, loc="lower left")
fig.tight_layout(); fig.savefig(f"{O}/fig3_age_bins.png"); plt.close(fig)

shutil.copy(f"{RV2}/standardisation_monitor/fig_monitor_alarm_vs_coverage_loss.png", f"{O}/fig4_monitor_vs_coverage.png")
print("figures written")
