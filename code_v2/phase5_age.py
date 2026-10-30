"""
Phase 5: age-conditional coverage. Usage: python -I code_v2/phase5_age.py RESULTS_V2
Textbook rank; 18 non-control transitions. Writes results_v2/age/.
"""
import sys, os
import numpy as np
import pandas as pd

RV2 = sys.argv[1]; OUT = f"{RV2}/age"; os.makedirs(OUT, exist_ok=True)
CONTROL = {("FS60", "FS53"), ("FS53", "FS60")}
BINS = [("<12", 0, 12), ("12-15", 12, 15), ("15-18", 15, 18), ("18-25", 18, 25), (">=25", 25, 200)]
d = pd.read_csv(f"{RV2}/frozen/per_subject_seed.csv.gz")
d = d[(d.arm == "main") & (d.rule == "textbook")]
d = d[[(s, t) not in CONTROL for s, t in zip(d.src, d.tgt)]].copy()
d["bin"] = pd.cut(d.age, [b[1] for b in BINS] + [200], right=False, labels=[b[0] for b in BINS])
conds = ["src", "raw", "std", "src_site", "sitestd"]
sm = d.groupby(["src", "tgt", "sid"]).agg(site=("site", "first"), age=("age", "first"), bin=("bin", "first"),
                                         **{f"c_{c}": (f"cov_{c}", "mean") for c in conds}, err_src=("err_src", "mean")).reset_index()

# ---- dataset-level description of bins (unique subjects across the 18 transitions)
u = sm.drop_duplicates("sid")
desc = u.groupby("bin", observed=True).agg(n_unique_subjects=("sid", "size"), n_sites=("site", "nunique"), age_min=("age", "min"), age_max=("age", "max")).reset_index()
per_tr = sm.groupby(["src", "tgt", "bin"], observed=True).size().unstack()
desc["n_per_transition_min"] = per_tr.min().values; desc["n_per_transition_median"] = per_tr.median().values; desc["n_per_transition_max"] = per_tr.max().values
adult_sites = u[u.bin == ">=25"].site.value_counts()
desc.to_csv(f"{OUT}/bin_description.csv", index=False)

# ---- joint site bootstrap machinery
rng = np.random.default_rng(5)
sites = np.array(sorted(sm.site.unique())); pos = {s: i for i, s in enumerate(sites)}
BM = rng.multinomial(len(sites), np.ones(len(sites)) / len(sites), size=2000)


def boot_bin(df, col):
    """mean over transitions of within-transition bin mean of col; joint site bootstrap."""
    ests, pts = [], []
    for (s, t), g in df.groupby(["src", "tgt"]):
        gs = g.groupby("site")[col].agg(["sum", "count"]); idx = np.array([pos[q] for q in gs.index]); w = BM[:, idx]
        den = (w * gs["count"].values).sum(1)
        ests.append(np.where(den > 0, (w * gs["sum"].values).sum(1) / np.maximum(den, 1), np.nan)); pts.append(g[col].mean())
    E_ = np.nanmean(np.vstack(ests), 0)
    return float(np.mean(pts)), float(np.nanpercentile(E_, 2.5)), float(np.nanpercentile(E_, 97.5)), float(np.min(pts)), float(np.max(pts))


rows = []
for b, _, _ in BINS:
    g = sm[sm.bin == b].copy()
    for c in conds:
        m, lo, hi, mn, mx = boot_bin(g, f"c_{c}")
        rows.append(dict(bin=b, quantity=f"cov_{c}", mean_over_transitions=m, lo=lo, hi=hi, min_transition=mn, max_transition=mx))
    for nm, (a, r) in {"d_raw": ("raw", "src"), "d_std": ("std", "src"), "d_sitestd": ("sitestd", "src_site")}.items():
        g[nm] = g[f"c_{a}"] - g[f"c_{r}"]; m, lo, hi, mn, mx = boot_bin(g, nm)
        rows.append(dict(bin=b, quantity=nm, mean_over_transitions=m, lo=lo, hi=hi, min_transition=mn, max_transition=mx))
    m, lo, hi, mn, mx = boot_bin(g, "err_src")
    rows.append(dict(bin=b, quantity="mean_signed_error_src_years", mean_over_transitions=m, lo=lo, hi=hi, min_transition=mn, max_transition=mx))
A = pd.DataFrame(rows); A.to_csv(f"{OUT}/age_bin_coverage.csv", index=False)

# ---- seed variation (pooled over transitions, subject-weighted)
sv = d.groupby(["bin", "seed"], observed=True)[[f"cov_{c}" for c in ("src", "raw", "std")]].mean().groupby("bin", observed=True).agg(["mean", "std", "min", "max"])
sv.to_csv(f"{OUT}/age_bin_seed_variation.csv")

# ---- how much of the marginal loss is the adult tail?
dec = []
for (s, t), g in sm.groupby(["src", "tgt"]):
    tot = (g.c_std - g.c_src).mean(); young = g[g.age < 25]; old = g[g.age >= 25]
    w_old = len(old) / len(g)
    dec.append(dict(src=s, tgt=t, N=len(g), n_ge25=len(old), share_ge25=w_old, d_std_all=tot, d_std_lt25=(young.c_std - young.c_src).mean(),
                    d_std_ge25=(old.c_std - old.c_src).mean(), contribution_ge25=w_old * (old.c_std - old.c_src).mean(),
                    cov_src_all=g.c_src.mean(), cov_src_lt25=young.c_src.mean(), cov_std_lt25=young.c_std.mean(),
                    d_raw_all=(g.c_raw - g.c_src).mean(), d_raw_lt25=(young.c_raw - young.c_src).mean()))
D = pd.DataFrame(dec); D.to_csv(f"{OUT}/decomposition_by_transition.csv", index=False)
sumd = D[["share_ge25", "d_std_all", "d_std_lt25", "d_std_ge25", "contribution_ge25", "cov_src_all", "cov_src_lt25", "cov_std_lt25", "d_raw_all", "d_raw_lt25"]].agg(["mean", "median", "min", "max"])
sumd.to_csv(f"{OUT}/decomposition_summary.csv")
g = sm[sm.age < 25].copy(); g["dd"] = g["c_std"] - g["c_src"]
lt25 = boot_bin(g, "dd")
pd.set_option("display.width", 250)
print(desc.to_string()); print("adult sites:", adult_sites.to_dict())
print(A.round(3).to_string()); print(sv.round(3).to_string()); print(sumd.round(4).to_string()); print("d_std <25 mean [CI]:", np.round(lt25, 4))
print("n transitions with d_std_lt25 <= -0.05:", int((D.d_std_lt25 <= -0.05).sum()), " vs all ages:", int((D.d_std_all <= -0.05).sum()))
