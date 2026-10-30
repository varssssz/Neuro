"""
Analysis engine: frozen ridge/LightGBM + split conformal (90%) under a same-scan pipeline swap.

Conformal procedure (absolute residual, split):
    q = Quantile( {|y_i - f(x_i)|}_{i in calibration} ; ceil((n+1)(1-alpha)) / n ),  C(x) = [f(x) - q, f(x) + q]
Cross-fitting: 5 folds; in each fold the remaining 80% is split 75/25 into train / calibration
(60/20/20 overall). Every subject is a test subject exactly once per seed.

Conditions evaluated on the same test subjects:
  src      source pipeline features, frozen model, frozen q                (reference)
  raw      target pipeline features, frozen model, frozen q                (the swap)
  std      raw + label-free per-pipeline standardisation: target features z-scored with the mean/SD of
           UNLABELLED target-pipeline features of the fold's training subjects; model and q still frozen
  sitestd  label-free per-site standardisation of each pipeline (site stats from the fold's training
           subjects; sites with <15 training subjects fall back to the pipeline-global stats);
           model re-fit on source site-standardised features, q re-computed on source calibration
  wcp      weighted conformal prediction (Tibshirani et al. 2019) with a logistic-regression density ratio between
           source calibration features and target training features (both in the frozen source scaling)
  reform   ReForm-style recalibration: k labelled target-pipeline subjects (drawn from the fold's train+calibration
           subjects) used to widen/shift the frozen interval's lower and upper bounds separately with
           alpha/2 per side (CQR-type correction)
"""
import numpy as np
import pandas as pd
from scipy.stats import binomtest, ks_2samp
from sklearn.linear_model import RidgeCV, LogisticRegression
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

ALPHA = 0.10
ALPHAS = np.logspace(-2, 4, 30)
KS = (10, 20, 40, 80, 160)
MIN_SITE = 15


def cq(res, alpha=ALPHA):
    n = len(res)
    lev = min(1.0, np.ceil((n + 1) * (1 - alpha)) / n)
    return np.quantile(res, lev, method="higher")


def wilson(p, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, c - h), min(1.0, c + h))


def boot_mean(x, B=2000, seed=11):
    x = np.asarray(x, float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), (B, len(x)))
    m = x[idx].mean(1)
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def mcnemar_p(a, b):
    """exact McNemar: a, b boolean covered indicators for the same subjects."""
    n10 = int(((a == 1) & (b == 0)).sum())
    n01 = int(((a == 0) & (b == 1)).sum())
    return binomtest(n10, n10 + n01, 0.5).pvalue if n10 + n01 else 1.0


def folds(n, seed):
    for f, (tr_all, te) in enumerate(KFold(5, shuffle=True, random_state=seed).split(np.arange(n))):
        rng = np.random.default_rng(seed * 1000 + f)
        perm = rng.permutation(tr_all)
        k = int(0.75 * len(perm))
        yield f, perm[:k], perm[k:], te


def site_z(X, sites, ref_idx):
    """z-score each row with the stats of its own site computed on ref_idx rows (unlabelled); fallback global."""
    g_mu, g_sd = X[ref_idx].mean(0), X[ref_idx].std(0) + 1e-9
    out = np.empty_like(X)
    for s in np.unique(sites):
        rows = np.where(sites == s)[0]
        ref = np.intersect1d(rows, ref_idx)
        if len(ref) >= MIN_SITE:
            mu, sd = X[ref].mean(0), X[ref].std(0) + 1e-9
        else:
            mu, sd = g_mu, g_sd
        out[rows] = (X[rows] - mu) / sd
    return out


def wq(r, w, wt, alpha=ALPHA):
    """weighted split-conformal quantile (one value per test weight in wt); inf when not reached."""
    o = np.argsort(r)
    rs, ws = r[o], w[o]
    C = np.cumsum(ws)
    W = ws.sum()
    idx = np.searchsorted(C, (1 - alpha) * (W + wt), side="left")
    q = np.where(idx >= len(rs), np.inf, rs[np.minimum(idx, len(rs) - 1)])
    return q


def run_transition(Xs, Xt, y, sites, seed, rrep=100, do_wcp=True, do_reform=True, do_site=True):
    """Return per-subject arrays (length n) and recalibration draws for one directed transition and one seed."""
    n = len(y)
    names = ["src", "raw", "std", "src_site", "sitestd", "wcp", "wcp_reg"]
    cov = {k: np.zeros(n) for k in names}
    qsite = np.zeros(n)
    wcp_reg_w = np.zeros(n)
    foldcov = {"src": [], "raw": [], "std": []}
    ae = {k: np.zeros(n) for k in ["src", "raw", "std", "sitestd"]}
    qv = np.zeros(n)
    wcp_w = np.zeros(n)
    rng = np.random.default_rng(seed + 7)
    rec = {k: {"cov": [], "width": []} for k in KS}
    rec2 = {k: {"cov": [], "width": []} for k in KS}
    for f, tr, ca, te in folds(n, seed):
        sc_s, sc_t = StandardScaler().fit(Xs[tr]), StandardScaler().fit(Xt[tr])
        rg = RidgeCV(alphas=ALPHAS).fit(sc_s.transform(Xs[tr]), y[tr])
        q = cq(np.abs(y[ca] - rg.predict(sc_s.transform(Xs[ca]))))
        ps = rg.predict(sc_s.transform(Xs[te]))
        pr = rg.predict(sc_s.transform(Xt[te]))
        pz = rg.predict(sc_t.transform(Xt[te]))
        for nm, p in (("src", ps), ("raw", pr), ("std", pz)):
            cov[nm][te] = np.abs(y[te] - p) <= q
            ae[nm][te] = np.abs(y[te] - p)
        qv[te] = q
        for nm in ("src", "raw", "std"):
            foldcov[nm].append(cov[nm][te].mean())
        if do_site:
            Zs, Zt = site_z(Xs, sites, tr), site_z(Xt, sites, tr)
            rg2 = RidgeCV(alphas=ALPHAS).fit(Zs[tr], y[tr])
            q2 = cq(np.abs(y[ca] - rg2.predict(Zs[ca])))
            p2s, p2t = rg2.predict(Zs[te]), rg2.predict(Zt[te])
            qsite[te] = q2
            cov["src_site"][te] = np.abs(y[te] - p2s) <= q2
            cov["sitestd"][te] = np.abs(y[te] - p2t) <= q2
            ae["sitestd"][te] = np.abs(y[te] - p2t)
        if do_wcp:
            Z = np.vstack([sc_s.transform(Xs[ca]), sc_s.transform(Xt[tr])])
            lab = np.r_[np.zeros(len(ca)), np.ones(len(tr))]
            clf = LogisticRegression(max_iter=2000, C=1.0).fit(Z, lab)

            def w(X):
                pp = np.clip(clf.predict_proba(sc_s.transform(X))[:, 1], 1e-3, 1 - 1e-3)
                return pp / (1 - pp) * len(ca) / len(tr)
            r_cal = np.abs(y[ca] - rg.predict(sc_s.transform(Xs[ca])))
            qw = wq(r_cal, w(Xs[ca]), w(Xt[te]))
            cov["wcp"][te] = np.abs(y[te] - pr) <= qw
            wcp_w[te] = 2 * qw
            # regularised variant (sensitivity): stronger L2 (C=0.01) and weights bounded by clipping p to [0.05, 0.95]
            clf2 = LogisticRegression(max_iter=2000, C=0.01).fit(Z, lab)

            def w2(X):
                pp = np.clip(clf2.predict_proba(sc_s.transform(X))[:, 1], 0.05, 0.95)
                return pp / (1 - pp) * len(ca) / len(tr)
            qw2 = wq(r_cal, w2(Xs[ca]), w2(Xt[te]))
            cov["wcp_reg"][te] = np.abs(y[te] - pr) <= qw2
            wcp_reg_w[te] = 2 * qw2
        if do_reform:
            pool = np.r_[tr, ca]
            pp = rg.predict(sc_s.transform(Xt[pool]))
            lo_p, hi_p = pp - q, pp + q
            e_lo, e_hi = lo_p - y[pool], y[pool] - hi_p
            pp2 = rg.predict(sc_t.transform(Xt[pool]))  # hybrid: label-free standardisation first, then recalibrate
            e_lo2, e_hi2 = (pp2 - q) - y[pool], y[pool] - (pp2 + q)
            for k in KS:
                lev = min(1.0, np.ceil((k + 1) * (1 - ALPHA / 2)) / k)
                for _ in range(rrep):
                    sel = rng.choice(len(pool), k, replace=False)
                    dlo = np.quantile(e_lo[sel], lev, method="higher")
                    dhi = np.quantile(e_hi[sel], lev, method="higher")
                    rec[k]["cov"].append(np.mean((y[te] >= pr - q - dlo) & (y[te] <= pr + q + dhi)))
                    rec[k]["width"].append(2 * q + dlo + dhi)
                    dlo2 = np.quantile(e_lo2[sel], lev, method="higher")
                    dhi2 = np.quantile(e_hi2[sel], lev, method="higher")
                    rec2[k]["cov"].append(np.mean((y[te] >= pz - q - dlo2) & (y[te] <= pz + q + dhi2)))
                    rec2[k]["width"].append(2 * q + dlo2 + dhi2)
    return cov, ae, qv, wcp_w, rec, dict(qsite=qsite, wcp_reg_w=wcp_reg_w, foldcov=foldcov, rec2=rec2)


def run_lgbm(Xs, Xt, y, seed):
    import lightgbm as lgb
    n = len(y)
    cov = {k: np.zeros(n) for k in ["src", "raw", "std"]}
    for f, tr, ca, te in folds(n, seed):
        sc_s, sc_t = StandardScaler().fit(Xs[tr]), StandardScaler().fit(Xt[tr])
        m = lgb.LGBMRegressor(n_estimators=300, learning_rate=0.03, num_leaves=15, min_child_samples=20, subsample=0.8,
                              subsample_freq=1, colsample_bytree=0.8, random_state=seed, n_jobs=4, verbose=-1)
        m.fit(sc_s.transform(Xs[tr]), y[tr])
        q = cq(np.abs(y[ca] - m.predict(sc_s.transform(Xs[ca]))))
        for nm, p in (("src", m.predict(sc_s.transform(Xs[te]))), ("raw", m.predict(sc_s.transform(Xt[te]))),
                      ("std", m.predict(sc_t.transform(Xt[te])))):
            cov[nm][te] = np.abs(y[te] - p) <= q
    return cov


# ------------------------------------------------------------------ monitors
def _rbf_mmd(X, Y, sigma, kyy):
    def K(A, B):
        d2 = (A * A).sum(1)[:, None] + (B * B).sum(1)[None, :] - 2 * A @ B.T
        return np.exp(-np.maximum(d2, 0) / (2 * sigma ** 2))
    kxx = K(X, X)
    m = len(X)
    return (kxx.sum() - np.trace(kxx)) / (m * (m - 1)) + kyy - 2 * K(X, Y).mean()


def run_monitor(Xs, Xt, y, seed=0, bsz=50, nb=40, nref=200):
    """Label-free monitors per fold; returns pooled null (calibration), null (test) and alternative batch stats."""
    n = len(y)
    rng = np.random.default_rng(seed + 99)
    out = {c: {"mmd_cal": [], "mmd_nul": [], "mmd_alt": [], "ks_cal": [], "ks_nul": [], "ks_alt": []} for c in ("raw", "std")}
    for f, tr, ca, te in folds(n, seed):
        sc_s, sc_t = StandardScaler().fit(Xs[tr]), StandardScaler().fit(Xt[tr])
        rg = RidgeCV(alphas=ALPHAS).fit(sc_s.transform(Xs[tr]), y[tr])
        refi = rng.choice(tr, min(nref, len(tr)), replace=False)
        ref = sc_s.transform(Xs[refi])
        sub = ref[rng.choice(len(ref), min(150, len(ref)), replace=False)]
        d = np.sqrt(np.maximum(((sub[:, None] - sub[None]) ** 2).sum(-1), 0))
        sigma = float(np.median(d[np.triu_indices(len(sub), 1)]))
        kyy_mat = np.exp(-np.maximum((ref * ref).sum(1)[:, None] + (ref * ref).sum(1)[None, :] - 2 * ref @ ref.T, 0) / (2 * sigma ** 2))
        nr = len(ref)
        kyy = (kyy_mat.sum() - np.trace(kyy_mat)) / (nr * (nr - 1))
        pref = rg.predict(sc_s.transform(Xs[tr]))
        sets = {"cal": ca, "nul": te}
        for cond, scal in (("raw", sc_s), ("std", sc_t)):
            for tag, pool in (("cal", ca), ("nul", te), ("alt", te)):
                for _ in range(nb):
                    b = rng.choice(pool, bsz, replace=False)
                    if tag == "alt":
                        F = scal.transform(Xt[b])
                        P = rg.predict(F)
                    else:
                        F = sc_s.transform(Xs[b])
                        P = rg.predict(F)
                    out[cond][f"mmd_{tag}"].append(_rbf_mmd(F, ref, sigma, kyy))
                    out[cond][f"ks_{tag}"].append(ks_2samp(P, pref).statistic)
    return out


# ------------------------------------------------------------------ mixed-pipeline training
def run_mixed(Xd, y, seed, pipes):
    """Xd: dict pipeline -> array (n x p) for the SAME subjects. Returns per-target per-subject coverage arrays."""
    n = len(y)
    res = {t: {"within": np.zeros(n), "xfer_raw": np.zeros((len(pipes) - 1, n)), "xfer_std": np.zeros((len(pipes) - 1, n)),
               "mix_all": np.zeros(n), "mix_loo": np.zeros(n), "w_all": np.zeros(n), "w_loo": np.zeros(n),
               "ae_all": np.zeros(n), "ae_loo": np.zeros(n)} for t in pipes}
    for f, tr, ca, te in folds(n, seed):
        sc = {p: StandardScaler().fit(Xd[p][tr]) for p in pipes}
        mod, qq = {}, {}
        for p in pipes:
            mod[p] = RidgeCV(alphas=ALPHAS).fit(sc[p].transform(Xd[p][tr]), y[tr])
            qq[p] = cq(np.abs(y[ca] - mod[p].predict(sc[p].transform(Xd[p][ca]))))
        for t in pipes:
            res[t]["within"][te] = np.abs(y[te] - mod[t].predict(sc[t].transform(Xd[t][te]))) <= qq[t]
            for j, s in enumerate([p for p in pipes if p != t]):
                res[t]["xfer_raw"][j, te] = np.abs(y[te] - mod[s].predict(sc[s].transform(Xd[t][te]))) <= qq[s]
                res[t]["xfer_std"][j, te] = np.abs(y[te] - mod[s].predict(sc[t].transform(Xd[t][te]))) <= qq[s]
            for mode, tp in (("all", pipes), ("loo", [p for p in pipes if p != t])):
                Xtr = np.vstack([Xd[p][tr] for p in tp]); ytr = np.concatenate([y[tr]] * len(tp))
                Xca = np.vstack([Xd[p][ca] for p in tp]); yca = np.concatenate([y[ca]] * len(tp))
                scp = StandardScaler().fit(Xtr)  # one scaler on the pooled mixture (NOT per pipeline): Korbmacher-style shuffling
                m = RidgeCV(alphas=ALPHAS).fit(scp.transform(Xtr), ytr)
                q = cq(np.abs(yca - m.predict(scp.transform(Xca))))
                pt = m.predict(scp.transform(Xd[t][te]))
                res[t][f"mix_{mode}"][te] = np.abs(y[te] - pt) <= q
                res[t][f"w_{mode}"][te] = 2 * q
                res[t][f"ae_{mode}"][te] = np.abs(y[te] - pt)
    return res
