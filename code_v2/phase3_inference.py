"""
Phase 3: controls and inference. Usage: python -I code_v2/phase3_inference.py RESULTS_V2 ORIGINAL_RESULTS PHASE1_OUT
Inputs: results_v2/frozen/per_subject_seed.csv.gz (v1 + textbook rank), v1 transitions_main/atlas, phase-1 feature_alignment.
Writes results_v2/inference/.
Inference design (documented in INFERENCE_AUDIT.md):
  unit = subject; indicator averaged over 5 seeds (cross-fitting randomness integrated out);
  dependence: subjects clustered by acquisition site (site effects shared within site); clusters assumed independent;
  per-transition test: cluster-robust (CR1) t-test of the mean paired difference, df = G-1; CI by site-cluster bootstrap;
  multiplicity: Holm within each pre-defined family of 18 non-control transitions (one family per contrast);
  aggregate over transitions: joint site-cluster bootstrap (same resampled sites applied to every transition), which
  respects re-use of the same subjects across transitions.
"""
import sys, os, json
import numpy as np
import pandas as pd
from scipy import stats

RV2, RV1, P1 = sys.argv[1:4]
OUT = f"{RV2}/inference"; os.makedirs(OUT, exist_ok=True)
B = 2000
CONTROL = {("FS60", "FS53"), ("FS53", "FS60")}
d = pd.read_csv(f"{RV2}/frozen/per_subject_seed.csv.gz")

# ---------------------------------------------------------------- 0. reproduction check of v1 per-subject values
if os.path.exists(f"{RV1}/per_subject_main.csv"):   # subject-level v1 file is not in the public repository
    v1ps = pd.read_csv(f"{RV1}/per_subject_main.csv")
    m = d[(d.arm == "main") & (d.rule == "v1")].groupby(["src", "tgt", "sid"])[["cov_src", "cov_raw", "cov_std", "cov_sitestd"]].mean().reset_index()
    mm = m.merge(v1ps, on=["src", "tgt", "sid"], suffixes=("_v2", "_v1"))
    repro = {c: float((mm[f"{c}_v2"] - mm[f"{c}_v1"]).abs().max()) for c in ("cov_src", "cov_raw", "cov_std", "cov_sitestd")}
    repro["rows_matched"] = len(mm); repro["rows_v1"] = len(v1ps)
else:
    repro = {"skipped": "original_v1/Experimental_Results/per_subject_main.csv not present"}

# ---------------------------------------------------------------- helpers
def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); n = len(p); adj = np.empty(n)
    run = 0
    for i, j in enumerate(o):
        run = max(run, (n - i) * p[j]); adj[j] = min(1.0, run)
    return adj


def cr1_test(x, g):
    """Cluster-robust t-test that mean(x)=0; x per subject, g cluster labels."""
    x = np.asarray(x, float); n = len(x); mu = x.mean()
    cl = pd.Series(x - mu).groupby(np.asarray(g)).sum().values; G = len(cl)
    se = np.sqrt(G / (G - 1) * (cl ** 2).sum()) / n
    if se == 0:
        return mu, se, 1.0 if mu == 0 else 0.0, G
    t = mu / se
    return mu, se, float(2 * stats.t.sf(abs(t), G - 1)), G


def site_sums(x, g):
    s = pd.DataFrame(dict(x=x, g=g)).groupby("g").x.agg(["sum", "count"])
    return s.index.values, s["sum"].values, s["count"].values.astype(float)


rng = np.random.default_rng(2026)
ALL_SITES = np.array(sorted(d.site.unique()))
BOOT_M = rng.multinomial(len(ALL_SITES), np.ones(len(ALL_SITES)) / len(ALL_SITES), size=B)   # B x G site multiplicities (shared)
site_pos = {s: i for i, s in enumerate(ALL_SITES)}
rng_s = np.random.default_rng(7)


def cluster_boot_mean(x, g):
    sites, S, N = site_sums(x, g); idx = np.array([site_pos[s] for s in sites])
    w = BOOT_M[:, idx]; den = (w * N).sum(1)
    est = (w * S).sum(1) / np.where(den > 0, den, np.nan)
    return est


def subj_boot(x):
    x = np.asarray(x, float); i = rng_s.integers(0, len(x), (B, len(x)))
    return x[i].mean(1)


def mcnemar(a, b):
    n10 = int(((a == 1) & (b == 0)).sum()); n01 = int(((a == 0) & (b == 1)).sum())
    return stats.binomtest(n10, n10 + n01, 0.5).pvalue if n10 + n01 else 1.0


CONTRASTS = {"raw_vs_src": ("cov_raw", "cov_src"), "std_vs_src": ("cov_std", "cov_src"),
             "sitestd_vs_srcsite": ("cov_sitestd", "cov_src_site"), "sitestd_vs_src": ("cov_sitestd", "cov_src")}
fa = pd.read_csv(f"{P1}/feature_alignment.csv").set_index("pair")
tm = pd.read_csv(f"{RV1}/transitions_main.csv").set_index(["src", "tgt"])
ta = pd.read_csv(f"{RV1}/transitions_atlas.csv").set_index(["src", "tgt"])


def category(s, t):
    pr = {s, t}
    if "FS60DK" in pr:
        o = (pr - {"FS60DK"}).pop()
        return {"FS60": "atlas only (same FS6.0 run, DK vs DKT)", "FS53": "FS version, same DK atlas (6.0 vs 5.3)",
                "FS51": "FS version, same DK atlas (6.0 vs 5.1)"}[o]
    if pr <= {"FS51", "FS53", "FS60"}:
        return "within FreeSurfer: 5.3<->6.0 (designated control)" if pr == {"FS53", "FS60"} else "within FreeSurfer: involving 5.1"
    if pr == {"ANTS", "CIVET"}:
        return "ANTs<->CIVET"
    return "FreeSurfer<->ANTs" if "ANTS" in pr else "FreeSurfer<->CIVET"


def feat(s, t):
    key = f"{s}-{t}" if f"{s}-{t}" in fa.index else f"{t}-{s}"
    tab = ta if (s, t) in ta.index else tm
    r = tab.loc[(s, t)]
    return dict(mean_same_region_r=fa.loc[key, "mean_same_region_r"], min_same_region_r=fa.loc[key, "min_same_region_r"],
                regions_namesake_not_best=int(fa.loc[key, "n_regions_namesake_not_best_match"]),
                mean_offset_mm=r.mean_offset_mm, sd_ratio=r.sd_ratio)


rows, agg_inputs = [], {}
for rule in ("textbook", "v1"):
    dr = d[d.rule == rule]
    for (s, t, arm), g0 in dr.groupby(["src", "tgt", "arm"], sort=False):
        cols = [c for c in ("cov_src", "cov_raw", "cov_std", "cov_src_site", "cov_sitestd") if c in g0 and (arm == "main" or "site" not in c)]
        sm = g0.groupby("sid").agg({**{c: "mean" for c in cols}, "site": "first"})
        seed_cov = g0.groupby("seed")[cols].mean()
        r = dict(rule=rule, arm=arm, src=s, tgt=t, control=(s, t) in CONTROL, category=category(s, t), N=len(sm), n_sites=sm.site.nunique(), **feat(s, t))
        for c in cols:
            r[c] = sm[c].mean(); r[f"{c}_seed_sd"] = seed_cov[c].std(ddof=1)
            cb = cluster_boot_mean(sm[c].values, sm.site.values); r[f"{c}_clu_lo"], r[f"{c}_clu_hi"] = np.nanpercentile(cb, [2.5, 97.5])
        for nm, (a, b) in CONTRASTS.items():
            if a not in cols or b not in cols:
                continue
            x = (sm[a] - sm[b]).values
            mu, se, p, G = cr1_test(x, sm.site.values)
            sb = subj_boot(x); cb = cluster_boot_mean(x, sm.site.values)
            r[f"d_{nm}"] = mu; r[f"d_{nm}_cr1_se"] = se; r[f"p_{nm}_cr1"] = p
            r[f"d_{nm}_subj_lo"], r[f"d_{nm}_subj_hi"] = np.percentile(sb, [2.5, 97.5])
            r[f"d_{nm}_clu_lo"], r[f"d_{nm}_clu_hi"] = np.nanpercentile(cb, [2.5, 97.5])
            if nm in ("raw_vs_src", "std_vs_src"):
                ps = [mcnemar(g0[g0.seed == k].set_index("sid")[b].values, g0[g0.seed == k].set_index("sid")[a].values) for k in range(5)]
                r[f"mcnemar_{nm}_seed0"] = ps[0]; r[f"mcnemar_{nm}_n_seeds_p05"] = int(sum(p_ < 0.05 for p_ in ps))
                r[f"mcnemar_{nm}_p_min"], r[f"mcnemar_{nm}_p_max"] = min(ps), max(ps)
            agg_inputs[(rule, arm, s, t, nm)] = (x, sm.site.values)
        rows.append(r)
T = pd.DataFrame(rows)
# Holm within families (18 non-control main transitions, per rule, per contrast); Bonferroni shown for comparison
for rule in ("textbook", "v1"):
    sel = (T.rule == rule) & (T.arm == "main") & (~T.control)
    for nm in CONTRASTS:
        T.loc[sel, f"p_{nm}_holm"] = holm(T.loc[sel, f"p_{nm}_cr1"].values)
        T.loc[sel, f"p_{nm}_bonf"] = np.minimum(1, T.loc[sel, f"p_{nm}_cr1"].values * sel.sum())
T.to_csv(f"{OUT}/per_transition_inference.csv", index=False)

# ---------------------------------------------------------------- aggregates over the 18 (joint site bootstrap)
def joint_boot(keys):
    """mean and median over transitions of the per-transition mean difference, with jointly resampled sites."""
    ests = []
    for key in keys:
        x, g = agg_inputs[key]; ests.append(cluster_boot_mean(x, g))
    E_ = np.vstack(ests)                         # transitions x B
    point = np.array([agg_inputs[k][0].mean() for k in keys])
    return point, E_


agg = []
for rule in ("textbook", "v1"):
    base = T[(T.rule == rule) & (T.arm == "main") & (~T.control)]
    for c in ("cov_src", "cov_raw", "cov_std", "cov_src_site", "cov_sitestd"):
        agg.append(dict(rule=rule, quantity=c, n_transitions=len(base), mean=base[c].mean(), median=base[c].median(), min=base[c].min(), max=base[c].max(),
                        n_below_088=int((base[c] < 0.88).sum())))
    for nm in CONTRASTS:
        keys = [(rule, "main", s, t, nm) for s, t in zip(base.src, base.tgt)]
        point, E_ = joint_boot(keys)
        mean_b, med_b = np.nanmean(E_, 0), np.nanmedian(E_, 0)
        sub = base[f"d_{nm}"]
        agg.append(dict(rule=rule, quantity=f"d_{nm}", n_transitions=len(keys), mean=point.mean(), mean_lo=np.percentile(mean_b, 2.5), mean_hi=np.percentile(mean_b, 97.5),
                        median=np.median(point), median_lo=np.percentile(med_b, 2.5), median_hi=np.percentile(med_b, 97.5), min=point.min(), max=point.max(),
                        n_le_minus005=int((point <= -0.05).sum()), n_holm_sig=int((base[f"p_{nm}_holm"] < 0.05).sum()),
                        n_clusterCI_excl0=int(((base[f"d_{nm}_clu_hi"] < 0) | (base[f"d_{nm}_clu_lo"] > 0)).sum())))
A = pd.DataFrame(agg); A.to_csv(f"{OUT}/aggregate_18_noncontrol.csv", index=False)

# ---------------------------------------------------------------- category sensitivity (textbook rule)
cat = []
tt = T[(T.rule == "textbook")]
for cname, g in tt.groupby("category"):
    for nm in ("raw_vs_src", "std_vs_src"):
        keys = [("textbook", a, s, t_, nm) for a, s, t_ in zip(g.arm, g.src, g.tgt)]
        point, E_ = joint_boot(keys); mb = np.nanmean(E_, 0)
        cat.append(dict(category=cname, contrast=nm, n_transitions=len(keys), transitions=";".join(f"{s}>{t_}" for s, t_ in zip(g.src, g.tgt)),
                        mean_d=point.mean(), lo=np.percentile(mb, 2.5), hi=np.percentile(mb, 97.5), min_d=point.min(), max_d=point.max(),
                        mean_feature_r=g.mean_same_region_r.mean()))
C = pd.DataFrame(cat); C.to_csv(f"{OUT}/category_sensitivity.csv", index=False)

# ---------------------------------------------------------------- controls table
ctrl_pairs = [("FS53", "FS60"), ("FS60", "FS53"), ("FS60DK", "FS60"), ("FS60", "FS60DK"), ("FS51", "FS53"), ("FS53", "FS51")]
ct = tt.set_index(["src", "tgt"]).loc[ctrl_pairs].reset_index()
ct = ct[["src", "tgt", "category", "N", "n_sites", "mean_same_region_r", "min_same_region_r", "mean_offset_mm", "sd_ratio", "cov_src", "cov_raw", "cov_std",
         "d_raw_vs_src", "d_raw_vs_src_clu_lo", "d_raw_vs_src_clu_hi", "p_raw_vs_src_cr1", "d_std_vs_src", "d_std_vs_src_clu_lo", "d_std_vs_src_clu_hi"]]
ct.to_csv(f"{OUT}/controls_table.csv", index=False)
json.dump(repro, open(f"{OUT}/v1_reproduction_check.json", "w"), indent=2)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(json.dumps(repro)); print(A.round(4).to_string()); print(C.round(4).drop(columns="transitions").to_string()); print(ct.round(4).to_string())
