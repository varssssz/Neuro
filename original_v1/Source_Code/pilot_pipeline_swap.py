"""
Pilot: frozen ridge + split conformal (90%) under a same-scan pipeline swap.
Data: ABIDE I derivative tables from neurodatascience/compare-surf-tools (cloned locally).

Usage:  python pilot_pipeline_swap.py /path/to/compare-surf-tools/data

Design notes (see report):
 * 5-fold cross-fitting: every participant is a test subject exactly once, so paired
   participant-level bootstrap CIs are valid. Within each fold the remaining 80% is split
   75/25 into train / calibration (i.e. 60/20/20 overall).
 * Model and conformal quantile are fit on the SOURCE pipeline and frozen; the same test
   participants are then scored with TARGET-pipeline features.
 * Common features = 31 DKT cortical regions x 2 hemispheres shared by all five pipelines.
   FS5.1/5.3 tables are Desikan-Killiany (ENIGMA); only the DKT-shared regions are used.
 * No QC exclusions and no site harmonisation in this pilot (add in the full analysis).
"""
import re
import sys
import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

DATA = sys.argv[1] if len(sys.argv) > 1 else "compare-surf-tools/data"
ALPHA, SEEDS, BOOT = 0.10, [0, 1, 2], 2000
RNG = np.random.default_rng(123)


def sid(x):
    m = re.findall(r"\d{5,}", str(x))
    return int(m[-1]) if m else np.nan


def long_names(df, pattern, hemi_fn):
    out = {}
    for c in df.columns:
        m = re.match(pattern, str(c))
        if m:
            h, r = hemi_fn(m)
            out[c] = f"{h}_{r.lower()}"
    return df[list(out)].rename(columns=out)


def load_fs5(path):                      # FS 5.1 / 5.3 (ENIGMA, DK atlas)
    d = pd.read_csv(path)
    X = long_names(d, r"^([LR])_(\w+)_thickavg$", lambda m: (m[1].lower(), m[2]))
    X.index = d["SubjID"].map(sid)
    return X


def load_fs6(dirpath):                   # FS 6.0 DKT, tab-separated, lh and rh files
    parts = []
    for h in ["lh", "rh"]:
        d = pd.read_csv(f"{dirpath}/DKTatlas_thickness_{h}.csv", sep="\t")
        idx = d.iloc[:, 0].map(sid)
        X = long_names(d, rf"^{h}_(\w+?)_thickness$", lambda m, h=h: (h[0], m[1]))
        X = X.drop(columns=[c for c in X.columns if "meanthickness" in c], errors="ignore")
        X.index = idx
        parts.append(X)
    return pd.concat(parts, axis=1)


def load_civet(path, mapfile):
    d = pd.read_csv(path)
    m = pd.read_csv(mapfile)
    inv = {r.CIVET: r.Freesurfer for r in m.itertuples() if r.Freesurfer != "unknown"}
    cols = {}
    for c in d.columns[1:]:
        name, h = c.rsplit(".", 1)
        if name in inv:
            cols[c] = f"{h.lower()}_{inv[name]}"
    X = d[list(cols)].rename(columns=cols)
    X.index = d["SubjID"].map(sid)
    return X


def load_ants(path):
    a = pd.read_csv(path, header=None, low_memory=False)
    names = a.iloc[2, 2:].astype(str).tolist()
    codes = a.iloc[1, 2:].astype(str).tolist()
    body = a.iloc[3:, :].reset_index(drop=True)
    ids = body.iloc[:, 1].map(sid)
    data = body.iloc[:, 2:].apply(pd.to_numeric, errors="coerce")
    cols = {}
    for j, (n, c) in enumerate(zip(names, codes)):
        if c.isdigit() and int(c) >= 1000 and (n.startswith("left ") or n.startswith("right ")):
            h = "l" if n.startswith("left") else "r"
            cols[j] = f"{h}_{n.split(' ', 1)[1].replace(' ', '')}"
    X = data.iloc[:, list(cols)]
    X.columns = list(cols.values())
    X.index = ids.values
    return X


def build():
    P = {
        "FS51": load_fs5(f"{DATA}/FS5.1/cortical_fs5.1_measuresenigma_thickavg.csv"),
        "FS53": load_fs5(f"{DATA}/FS5.3/ABIDE_fs5.3_thickness.csv"),
        "FS60": load_fs6(f"{DATA}/FS6.0"),
        "ANTS": load_ants(f"{DATA}/ANTs/ABIDE_ants_thickness_data.csv"),
        "CIVET": load_civet(f"{DATA}/CIVET/ABIDE_civet2.1_thickness.csv",
                            f"{DATA}/CIVET/DKT/DKT_parcel_map_FS_CIVET.csv"),
    }
    for k, v in P.items():
        v = v[~v.index.isna()]
        v = v[~v.index.duplicated()]
        P[k] = v.astype(float)
    common = sorted(set.intersection(*[set(v.columns) for v in P.values()]))
    ph = pd.read_csv(f"{DATA}/ABIDE_Phenotype.csv", low_memory=False)
    ph = ph.set_index("SUB_ID")[["AGE_AT_SCAN", "DX_GROUP", "SITE_ID"]].dropna()
    return {k: v[common] for k, v in P.items()}, ph, common


def cq(res, alpha=ALPHA):
    n = len(res)
    lev = min(1.0, np.ceil((n + 1) * (1 - alpha)) / n)
    return np.quantile(res, lev, method="higher")


def boot(x, B=BOOT):
    x = np.asarray(x, float)
    idx = RNG.integers(0, len(x), (B, len(x)))
    return np.percentile(x[idx].mean(1), [2.5, 97.5])


def run_pair(Xs, Xt, y, seed):
    ids = Xs.index.values
    n = len(ids)
    cs, ct, lo, hi, ts, tt, yy, sidx = [np.zeros(n) for _ in range(8)]
    pos = {s: i for i, s in enumerate(ids)}
    for tr_all, te in KFold(5, shuffle=True, random_state=seed).split(ids):
        rng = np.random.default_rng(seed)
        perm = rng.permutation(tr_all)
        k = int(0.75 * len(perm))
        tr, ca = perm[:k], perm[k:]
        model = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 4, 30)))
        model.fit(Xs.values[tr], y[tr])
        q = cq(np.abs(y[ca] - model.predict(Xs.values[ca])))
        ps, pt = model.predict(Xs.values[te]), model.predict(Xt.values[te])
        cs[te], ct[te] = np.abs(y[te] - ps) <= q, np.abs(y[te] - pt) <= q
        lo[te], hi[te] = pt - q, pt + q
        ts[te], tt[te], yy[te] = np.abs(y[te] - ps), np.abs(y[te] - pt), y[te]
    return cs, ct, lo, hi, ts, tt


def reform_style(y, lo, hi, ks=(10, 20, 40, 80, 160), R=500, alpha=ALPHA):
    rows = []
    n = len(y)
    for k in ks:
        if k > n - 30:
            continue
        covs = []
        for _ in range(R):
            p = RNG.permutation(n)
            cal, ev = p[:k], p[k:]
            lev = min(1.0, np.ceil((k + 1) * (1 - alpha / 2)) / k)
            dlo = np.quantile(lo[cal] - y[cal], lev, method="higher")
            dhi = np.quantile(y[cal] - hi[cal], lev, method="higher")
            covs.append(np.mean((y[ev] >= lo[ev] - dlo) & (y[ev] <= hi[ev] + dhi)))
        covs = np.array(covs)
        rows.append((k, covs.mean(), np.mean(covs < 0.88)))
    return rows


def main():
    P, ph, common = build()
    print(f"Common DKT features across 5 pipelines: {len(common)}")
    print({k: len(v) for k, v in P.items()}, "participants per table")
    names = list(P)
    out, bins_out = [], []
    for s in names:
        for t in names:
            if s == t:
                continue
            ids = sorted(set(P[s].dropna().index) & set(P[t].dropna().index) & set(ph.index))
            if len(ids) < 150:
                print(f"{s}->{t}: matched N={len(ids)} (skipped)")
                continue
            Xs, Xt = P[s].loc[ids], P[t].loc[ids]
            y = ph.loc[ids, "AGE_AT_SCAN"].values
            dx = ph.loc[ids, "DX_GROUP"].values
            dcs, dcts, w, mae_s, mae_t = [], [], [], [], []
            for sd in SEEDS:
                cs, ct, lo, hi, ts, tt = run_pair(Xs, Xt, y, sd)
                dcs.append(cs.mean()); dcts.append(ct.mean())
                mae_s.append(ts.mean()); mae_t.append(tt.mean())
                if sd == SEEDS[0]:
                    cs0, ct0, lo0, hi0 = cs, ct, lo, hi
            d = ct0 - cs0
            b, c = int(((cs0 == 1) & (ct0 == 0)).sum()), int(((cs0 == 0) & (ct0 == 1)).sum())
            p_mc = binomtest(b, b + c, 0.5).pvalue if b + c else 1.0
            lo_ci, hi_ci = boot(d)
            out.append(dict(src=s, tgt=t, N=len(ids), cov_src=np.mean(dcs), cov_tgt=np.mean(dcts),
                            dcov=np.mean(dcts) - np.mean(dcs), dcov_seed0=d.mean(),
                            ci_lo=lo_ci, ci_hi=hi_ci, mcnemar_p=p_mc,
                            width=float(np.mean(hi0 - lo0)), MAE_src=np.mean(mae_s), MAE_tgt=np.mean(mae_t),
                            mean_offset_mm=float((Xt.values - Xs.values).mean())))
            for lab, mask in [("<12", y < 12), ("12-15", (y >= 12) & (y < 15)), ("15-18", (y >= 15) & (y < 18)),
                              ("18-25", (y >= 18) & (y < 25)), (">=25", y >= 25),
                              ("ASD", dx == 1), ("TD", dx == 2)]:
                if mask.sum() >= 20:
                    bins_out.append(dict(src=s, tgt=t, stratum=lab, n=int(mask.sum()),
                                         cov_src=cs0[mask].mean(), cov_tgt=ct0[mask].mean()))
            if (s, t) in [("FS60", "CIVET"), ("FS60", "ANTS"), ("FS51", "FS53"), ("FS60", "FS53")]:
                print(f"\nReForm-style recalibration {s}->{t}: (k, mean coverage, P[cov<0.88])")
                for r in reform_style(y, lo0, hi0):
                    print(f"  k={r[0]:>3}  cov={r[1]:.3f}  P(<0.88)={r[2]:.2f}")
    res = pd.DataFrame(out).sort_values("dcov")
    res.to_csv("pilot_transitions.csv", index=False)
    pd.DataFrame(bins_out).to_csv("pilot_strata.csv", index=False)
    pd.set_option("display.width", 200)
    print("\n", res.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
