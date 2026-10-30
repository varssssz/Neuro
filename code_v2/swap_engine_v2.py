"""
swap_engine_v2 — corrected analysis engine (v2). The v1 engine (original_v1/Source_Code/swap_engine.py) is kept unchanged.

Changes relative to v1, and nothing else:
  1. Conformal quantile uses the conventional finite-sample rank r = ceil((n+1)*level): the r-th smallest score,
     or +inf when r > n. (v1: np.quantile(..., r/n, method="higher") gives the (r+1)-th smallest, capped at the maximum.)
  2. Labelled target recalibration draws its k subjects from a pool that is disjoint from model training,
     source calibration and the final test set (four-way split). v1 drew from training + calibration subjects.
  3. Extra, clearly named conditions for the Phase 4 standardisation sensitivity analysis.
Outer test folds are identical to v1 (same KFold seeds), so every condition scores the same test subjects as v1.
"""
import numpy as np
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

ALPHA = 0.10
ALPHAS = np.logspace(-2, 4, 30)          # same ridge grid as v1
KS = (10, 20, 40, 80, 160)
MIN_SITE = 15
V2_FRACTIONS = (0.50, 0.20, 0.30)        # of the non-test 80%: train / source-cal / target-recal pool -> 40/16/24/20 overall


# ------------------------------------------------------------------ conformal quantiles
def conformal_rank(n, level):
    """Conventional split-conformal rank: ceil((n+1)*level). May exceed n (=> infinite threshold)."""
    return int(np.ceil((n + 1) * level - 1e-12))


def conformal_quantile(scores, level):
    """r-th smallest score with r = ceil((n+1)*level); +inf if r > n."""
    s = np.sort(np.asarray(scores, float))
    r = conformal_rank(len(s), level)
    return np.inf if r > len(s) else float(s[r - 1])


def v1_quantile(scores, level):
    """Exact replica of the v1 rule (np.quantile at min(1, r/n), method='higher')."""
    n = len(scores)
    lev = min(1.0, np.ceil((n + 1) * level) / n)
    return float(np.quantile(scores, lev, method="higher"))


def v1_rank_used(n, level):
    """Rank (1-based) actually selected by v1_quantile, for documentation/tests."""
    lev = min(1.0, np.ceil((n + 1) * level) / n)
    return int(np.ceil((n - 1) * lev - 1e-12)) + 1


def expected_recal_coverage(k, alpha=ALPHA, rule="textbook"):
    """Exact expected coverage of the two-sided, separate-alpha/2 correction under exchangeability (continuous scores).
    Lower and upper misses are disjoint events, each with probability (k+1-r)/(k+1)."""
    if rule == "textbook":
        r = conformal_rank(k, 1 - alpha / 2)
        return 1.0 if r > k else 1 - 2 * (k + 1 - r) / (k + 1)
    r = v1_rank_used(k, 1 - alpha / 2)
    return 1 - 2 * (k + 1 - r) / (k + 1)


QUANT = {"textbook": conformal_quantile, "v1": v1_quantile}


# ------------------------------------------------------------------ splits
def folds_v1(n, seed):
    """Identical to v1 swap_engine.folds(): 5-fold outer CV; non-test split 75/25 train/cal."""
    for f, (tr_all, te) in enumerate(KFold(5, shuffle=True, random_state=seed).split(np.arange(n))):
        perm = np.random.default_rng(seed * 1000 + f).permutation(tr_all)
        k = int(0.75 * len(perm))
        yield f, perm[:k], perm[k:], np.array([], int), te


def folds_v2(n, seed, fr=V2_FRACTIONS):
    """Same outer test folds and same permutation as v1; non-test split train/cal/recal (disjoint)."""
    for f, (tr_all, te) in enumerate(KFold(5, shuffle=True, random_state=seed).split(np.arange(n))):
        perm = np.random.default_rng(seed * 1000 + f).permutation(tr_all)
        a = int(round(fr[0] * len(perm))); b = a + int(round(fr[1] * len(perm)))
        yield f, perm[:a], perm[a:b], perm[b:], te


# ------------------------------------------------------------------ per-site z-scoring (as v1)
def site_z(X, sites, ref_idx):
    g_mu, g_sd = X[ref_idx].mean(0), X[ref_idx].std(0) + 1e-9
    out = np.empty_like(X)
    for s in np.unique(sites):
        rows = np.where(sites == s)[0]
        ref = np.intersect1d(rows, ref_idx)
        mu, sd = (X[ref].mean(0), X[ref].std(0) + 1e-9) if len(ref) >= MIN_SITE else (g_mu, g_sd)
        out[rows] = (X[rows] - mu) / sd
    return out


# ------------------------------------------------------------------ frozen-model conditions (main table, v1 design)
def run_frozen(Xs, Xt, y, sites, seed, rule="textbook", do_site=True):
    """v1 60/20/20 design. Returns dict of per-subject arrays (length n) for one seed.
    Conditions: src, raw, std (target stats on Xt[tr], as v1), std_disjoint (target stats on Xt[ca]: unlabelled,
    disjoint from model-training subjects), std_same_sub (target stats on a random subset of tr of size |ca|:
    same subjects as training, matched sample size), src_site, sitestd (refit)."""
    Q = QUANT[rule]
    n = len(y)
    keys = ["src", "raw", "std", "std_disjoint", "std_same_sub", "src_site", "sitestd"]
    cov = {k: np.zeros(n) for k in keys}
    err = {k: np.zeros(n) for k in keys}    # signed error pred - y
    q = np.zeros(n); qs = np.zeros(n)
    rng = np.random.default_rng(seed + 31)
    for f, tr, ca, _, te in folds_v1(n, seed):
        sc_s = StandardScaler().fit(Xs[tr])
        rg = RidgeCV(alphas=ALPHAS).fit(sc_s.transform(Xs[tr]), y[tr])
        qq = Q(np.abs(y[ca] - rg.predict(sc_s.transform(Xs[ca]))), 1 - ALPHA)
        sub = rng.choice(tr, len(ca), replace=False)
        preds = {"src": rg.predict(sc_s.transform(Xs[te])),
                 "raw": rg.predict(sc_s.transform(Xt[te])),
                 "std": rg.predict(StandardScaler().fit(Xt[tr]).transform(Xt[te])),
                 "std_disjoint": rg.predict(StandardScaler().fit(Xt[ca]).transform(Xt[te])),
                 "std_same_sub": rg.predict(StandardScaler().fit(Xt[sub]).transform(Xt[te]))}
        for k, p in preds.items():
            cov[k][te] = np.abs(y[te] - p) <= qq; err[k][te] = p - y[te]
        q[te] = qq
        if do_site:
            Zs, Zt = site_z(Xs, sites, tr), site_z(Xt, sites, tr)
            rg2 = RidgeCV(alphas=ALPHAS).fit(Zs[tr], y[tr])
            q2 = Q(np.abs(y[ca] - rg2.predict(Zs[ca])), 1 - ALPHA)
            for k, p in (("src_site", rg2.predict(Zs[te])), ("sitestd", rg2.predict(Zt[te]))):
                cov[k][te] = np.abs(y[te] - p) <= q2; err[k][te] = p - y[te]
            qs[te] = q2
    return dict(cov=cov, err=err, width=2 * q, width_site=2 * qs)


# ------------------------------------------------------------------ labelled target recalibration
def run_recal(Xs, Xt, y, seed, design, rrep=100, ks=KS):
    """design: dict(split='v1'|'v2', pool='tr+ca'|'ca'|'recal', rule='v1'|'textbook').
    Returns per-k lists of draw-level coverage and width, plus base-model diagnostics for this split."""
    Q = QUANT[design["rule"]]
    n = len(y)
    fold_fn = folds_v1 if design["split"] == "v1" else folds_v2
    rng = np.random.default_rng(seed + 7)          # same draw RNG seeding as v1
    out = {k: dict(cov=[], width=[], cov_h=[], width_h=[]) for k in ks}
    base = dict(cov_src=np.zeros(n), cov_raw=np.zeros(n), ae_src=np.zeros(n), ae_raw=np.zeros(n), width=np.zeros(n),
                n_train=[], n_cal=[], n_pool=[], n_test=[], overlap=[])
    lev = 1 - ALPHA / 2
    for f, tr, ca, rc, te in fold_fn(n, seed):
        pool = {"tr+ca": np.r_[tr, ca], "ca": ca, "recal": rc}[design["pool"]]
        base["overlap"].append(len(set(pool) & set(te)))
        sc_s, sc_t = StandardScaler().fit(Xs[tr]), StandardScaler().fit(Xt[tr])
        rg = RidgeCV(alphas=ALPHAS).fit(sc_s.transform(Xs[tr]), y[tr])
        q = Q(np.abs(y[ca] - rg.predict(sc_s.transform(Xs[ca]))), 1 - ALPHA)
        ps, pr, pz = rg.predict(sc_s.transform(Xs[te])), rg.predict(sc_s.transform(Xt[te])), rg.predict(sc_t.transform(Xt[te]))
        base["cov_src"][te] = np.abs(y[te] - ps) <= q; base["cov_raw"][te] = np.abs(y[te] - pr) <= q
        base["ae_src"][te] = np.abs(y[te] - ps); base["ae_raw"][te] = np.abs(y[te] - pr); base["width"][te] = 2 * q
        for nm, v in (("n_train", tr), ("n_cal", ca), ("n_pool", pool), ("n_test", te)):
            base[nm].append(len(v))
        pp = rg.predict(sc_s.transform(Xt[pool])); e_lo, e_hi = (pp - q) - y[pool], y[pool] - (pp + q)
        pp2 = rg.predict(sc_t.transform(Xt[pool])); e_lo2, e_hi2 = (pp2 - q) - y[pool], y[pool] - (pp2 + q)
        for k in ks:
            if k > len(pool):
                continue
            for _ in range(rrep):
                sel = rng.choice(len(pool), k, replace=False)
                dlo, dhi = Q(e_lo[sel], lev), Q(e_hi[sel], lev)
                out[k]["cov"].append(np.mean((y[te] >= pr - q - dlo) & (y[te] <= pr + q + dhi)))
                out[k]["width"].append(2 * q + dlo + dhi)
                dlo2, dhi2 = Q(e_lo2[sel], lev), Q(e_hi2[sel], lev)
                out[k]["cov_h"].append(np.mean((y[te] >= pz - q - dlo2) & (y[te] <= pz + q + dhi2)))
                out[k]["width_h"].append(2 * q + dlo2 + dhi2)
    return out, base
