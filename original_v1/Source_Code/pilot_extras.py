"""Extras: (1) label-free per-pipeline standardisation baseline, (2) width after ReForm-style recalibration."""
import sys, numpy as np, pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
import pilot_pipeline_swap as M

DATA = sys.argv[1]
M.DATA = DATA
P, ph, common = M.build()
pairs = [("FS60","FS53"),("FS53","FS60"),("FS51","FS53"),("FS53","FS51"),("FS51","FS60"),
         ("FS60","ANTS"),("FS51","ANTS"),("ANTS","CIVET"),("FS60","CIVET"),("CIVET","FS60")]
rows = []
for s, t in pairs:
    ids = sorted(set(P[s].dropna().index) & set(P[t].dropna().index) & set(ph.index))
    Xs, Xt = P[s].loc[ids].values, P[t].loc[ids].values
    y = ph.loc[ids, "AGE_AT_SCAN"].values
    n = len(ids); raw = np.zeros(n); zc = np.zeros(n); src = np.zeros(n)
    maer = np.zeros(n); maez = np.zeros(n); lo = np.zeros(n); hi = np.zeros(n)
    for tr_all, te in KFold(5, shuffle=True, random_state=0).split(ids):
        rng = np.random.default_rng(0); perm = rng.permutation(tr_all)
        k = int(0.75*len(perm)); tr, ca = perm[:k], perm[k:]
        sc_s = StandardScaler().fit(Xs[tr]); sc_t = StandardScaler().fit(Xt[tr])   # t-scaler uses UNLABELLED target-pipeline train features
        rg = RidgeCV(alphas=np.logspace(-2,4,30)).fit(sc_s.transform(Xs[tr]), y[tr])
        q = M.cq(np.abs(y[ca]-rg.predict(sc_s.transform(Xs[ca]))))
        ps = rg.predict(sc_s.transform(Xs[te]))
        p_raw = rg.predict(sc_s.transform(Xt[te]))     # frozen, source scaling
        p_z   = rg.predict(sc_t.transform(Xt[te]))     # target scaling
        src[te] = np.abs(y[te]-ps) <= q
        raw[te] = np.abs(y[te]-p_raw) <= q; zc[te] = np.abs(y[te]-p_z) <= q
        maer[te] = np.abs(y[te]-p_raw); maez[te] = np.abs(y[te]-p_z)
        lo[te], hi[te] = p_raw-q, p_raw+q
    r = {"src":s,"tgt":t,"N":n,"cov_src":src.mean(),"cov_raw":raw.mean(),"cov_zscore_target":zc.mean(),
         "MAE_raw":maer.mean(),"MAE_zscore":maez.mean()}
    for k in (40, 80):
        ws=[]
        for _ in range(300):
            p = M.RNG.permutation(n); cal, ev = p[:k], p[k:]
            lev = min(1.0, np.ceil((k+1)*(1-0.05))/k)
            dlo = np.quantile(lo[cal]-y[cal], lev, method="higher"); dhi = np.quantile(y[cal]-hi[cal], lev, method="higher")
            ws.append((hi[ev]-lo[ev]+dlo+dhi).mean())
        r[f"width_recal_k{k}"] = np.mean(ws)
    r["width_src"] = float((hi-lo).mean())
    rows.append(r)
pd.set_option("display.width",220)
print(pd.DataFrame(rows).round(3).to_string(index=False))
