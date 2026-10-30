"""
Run the full analysis. Usage:  python run_all.py /path/to/compare-surf-tools/data [outdir]
All numbers in the manuscript come from the CSV files this script writes to `outdir` (default: results/).
"""
import sys, os, time, itertools
import numpy as np
import pandas as pd
import abide_data as AD
import swap_engine as E

DATA = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else "results"
os.makedirs(OUT, exist_ok=True)
SEEDS = [0, 1, 2, 3, 4]
P, ph, common, qc = AD.build(DATA)
PIPES = AD.PIPES5
PAIRS = [(s, t) for s in PIPES for t in PIPES if s != t]
ATLAS_PAIRS = [("FS60DK", "FS60"), ("FS60", "FS60DK"), ("FS60DK", "FS53"), ("FS53", "FS60DK"), ("FS60DK", "FS51"), ("FS51", "FS60DK")]
BINS = [("<12", lambda a: a < 12), ("12-15", lambda a: (a >= 12) & (a < 15)), ("15-18", lambda a: (a >= 15) & (a < 18)),
        ("18-25", lambda a: (a >= 18) & (a < 25)), (">=25", lambda a: a >= 25)]
CONTROL = {("FS60", "FS53"), ("FS53", "FS60")}
log = open(f"{OUT}/run_log.txt", "w")


def say(*a):
    msg = " ".join(str(x) for x in a)
    print(msg, flush=True)
    log.write(msg + "\n"); log.flush()


def get(s, t, level):
    ids = AD.matched_ids(P, ph, qc, s, t, level)
    return ids, P[s].loc[ids].values, P[t].loc[ids].values, ph.loc[ids, "AGE_AT_SCAN"].values, ph.loc[ids, "SITE_ID"].values, ph.loc[ids, "DX_GROUP"].values


def feature_stats(Xs, Xt):
    cors = [np.corrcoef(Xs[:, j], Xt[:, j])[0, 1] for j in range(Xs.shape[1])]
    return dict(mean_offset_mm=float((Xt - Xs).mean()), sd_ratio=float(Xt.std(0).mean() / Xs.std(0).mean()), mean_feature_r=float(np.mean(cors)))


# ---------------------------------------------------------------- A. matched N by site, before / after QC
def stage_matched_n():
    rows = []
    for s, t in PAIRS + ATLAS_PAIRS:
        for lvl in (0, 1, 2):
            ids = AD.matched_ids(P, ph, qc, s, t, lvl)
            ss = ph.loc[ids, "SITE_ID"].value_counts()
            rows.append(dict(src=s, tgt=t, qc_level=lvl, site="ALL", n=len(ids)))
            for site, c in ss.items():
                rows.append(dict(src=s, tgt=t, qc_level=lvl, site=site, n=int(c)))
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/matched_N_by_site.csv", index=False)
    piv = df[df.site == "ALL"].pivot_table(index=["src", "tgt"], columns="qc_level", values="n").reset_index()
    piv.columns = ["src", "tgt", "N_noQC", "N_QC1", "N_QC2"]
    piv.to_csv(f"{OUT}/matched_N_summary.csv", index=False)
    say("matched N summary (QC0/QC1/QC2) min-max:", piv[["N_noQC", "N_QC1", "N_QC2"]].agg(["min", "max"]).to_dict())


# ---------------------------------------------------------------- B. main transitions
def summarize(s, t, level, cov, ae, qv, wcpw, extras=True, seed_sd=None, wcpw_reg=None, qsite=None, foldcov=None):
    """cov/ae/qv: dicts of arrays [seed, n]."""
    out = dict(src=s, tgt=t, qc_level=level)
    c = {k: v.mean(0) for k, v in cov.items()}
    n = len(c["src"])
    out["N"] = n
    for k in ("src", "raw", "std"):
        out[f"cov_{k}"] = c[k].mean(); out[f"cov_{k}_lo"], out[f"cov_{k}_hi"] = E.boot_mean(c[k])
    for k in ("raw", "std"):
        d = c[k] - c["src"]
        out[f"d_{k}"] = d.mean(); out[f"d_{k}_lo"], out[f"d_{k}_hi"] = E.boot_mean(d)
        out[f"mcnemar_p_{k}"] = E.mcnemar_p(cov["src"][0], cov[k][0])
    if extras:
        for k in ("src_site", "sitestd", "wcp", "wcp_reg"):
            out[f"cov_{k}"] = c[k].mean()
        out["width_src_site"] = float((2 * qsite).mean())
        out["wcp_reg_frac_inf"] = float(np.isinf(wcpw_reg).mean())
        fin2 = wcpw_reg[np.isfinite(wcpw_reg)]
        out["wcp_reg_median_width"] = float(np.median(fin2)) if len(fin2) else np.nan
        for kk in ("wcp", "wcp_reg"):
            w_ = wcpw if kk == "wcp" else wcpw_reg
            fm = np.isfinite(w_)
            out[f"{kk}_cov_if_finite"] = float(cov[kk].reshape(-1)[fm].mean()) if fm.any() else np.nan
        d = c["sitestd"] - c["src_site"]
        out["d_sitestd"] = d.mean(); out["d_sitestd_lo"], out["d_sitestd_hi"] = E.boot_mean(d)
        wf = wcpw
        out["wcp_frac_inf"] = float(np.isinf(wf).mean())
        fin = wf[np.isfinite(wf)]
        out["wcp_median_width"] = float(np.median(fin)) if len(fin) else np.nan
    out["width_src"] = float((2 * qv).mean())
    if seed_sd is not None:
        for k in ("src", "raw", "std"):
            out[f"cov_{k}_seed_sd"] = float(seed_sd[k])
    if foldcov is not None:
        for k in ("src", "raw", "std"):
            out[f"frac_fold_below_088_{k}"] = float((np.array(foldcov[k]) < 0.88).mean())
    for k in ae:
        out[f"MAE_{k}"] = ae[k].mean()
    return out, c


def stage_main(level=1, extras=True, pairs=None, tag="main", seeds=SEEDS):
    rows, subj_rows, rec_rows = [], [], []
    for s, t in (pairs or PAIRS):
        t0 = time.time()
        ids, Xs, Xt, y, sites, dx = get(s, t, level)
        accs = {k: [] for k in ["src", "raw", "std", "src_site", "sitestd", "wcp", "wcp_reg"]}
        aes = {k: [] for k in ["src", "raw", "std", "sitestd"]}
        qvs, wcps, wcprs, qss = [], [], [], []
        fcs = {"src": [], "raw": [], "std": []}
        recs = {k: {"cov": [], "width": []} for k in E.KS}
        recs2 = {k: {"cov": [], "width": []} for k in E.KS}
        for sd in seeds:
            cov, ae, qv, wcpw, rec, ex = E.run_transition(Xs, Xt, y, sites, sd, rrep=100 if extras else 0, do_wcp=extras, do_reform=extras, do_site=extras)
            for k in accs: accs[k].append(cov[k])
            for k in aes: aes[k].append(ae[k])
            qvs.append(qv); wcps.append(wcpw); wcprs.append(ex["wcp_reg_w"]); qss.append(ex["qsite"])
            for k in fcs: fcs[k] += ex["foldcov"][k]
            if extras:
                for k in E.KS:
                    recs[k]["cov"] += rec[k]["cov"]; recs[k]["width"] += rec[k]["width"]
                    recs2[k]["cov"] += ex["rec2"][k]["cov"]; recs2[k]["width"] += ex["rec2"][k]["width"]
        cov = {k: np.array(v) for k, v in accs.items()}
        ae = {k: np.array(v).mean(0) for k, v in aes.items()}
        sdv = {k: float(np.std(cov[k].mean(1))) for k in ("src", "raw", "std")}
        out, c = summarize(s, t, level, cov, {k: v for k, v in ae.items()}, np.mean(qvs, 0), np.concatenate(wcps), extras,
                           seed_sd=sdv, wcpw_reg=np.concatenate(wcprs), qsite=np.mean(qss, 0), foldcov=fcs)
        out.update(feature_stats(Xs, Xt))
        out["control"] = (s, t) in CONTROL
        rows.append(out)
        sub = pd.DataFrame(dict(src=s, tgt=t, sid=ids, site=sites, age=y, dx=dx, **{f"cov_{k}": v for k, v in c.items()}))
        subj_rows.append(sub)
        if extras:
            for k in E.KS:
                cc, ww = np.array(recs[k]["cov"]), np.array(recs[k]["width"])
                rec_rows.append(dict(src=s, tgt=t, variant="reform_raw", k=k, mean_cov=cc.mean(), cov_lo=np.percentile(cc, 2.5), cov_hi=np.percentile(cc, 97.5),
                                     frac_below_088=float((cc < 0.88).mean()), mean_width=ww.mean(), n_draws=len(cc)))
                c2, w2 = np.array(recs2[k]["cov"]), np.array(recs2[k]["width"])
                rec_rows.append(dict(src=s, tgt=t, variant="std_then_reform", k=k, mean_cov=c2.mean(), cov_lo=np.percentile(c2, 2.5), cov_hi=np.percentile(c2, 97.5),
                                     frac_below_088=float((c2 < 0.88).mean()), mean_width=w2.mean(), n_draws=len(c2)))
        say(f"[{tag}] {s}->{t} N={len(ids)} src={out['cov_src']:.3f} raw={out['cov_raw']:.3f} std={out['cov_std']:.3f} ({time.time()-t0:.1f}s)")
    pd.DataFrame(rows).to_csv(f"{OUT}/transitions_{tag}.csv", index=False)
    pd.concat(subj_rows).to_csv(f"{OUT}/per_subject_{tag}.csv", index=False)
    if rec_rows:
        pd.DataFrame(rec_rows).to_csv(f"{OUT}/recal_{tag}.csv", index=False)


# ---------------------------------------------------------------- strata (age bins, diagnosis) from per-subject file
def stage_strata(tag="main"):
    d = pd.read_csv(f"{OUT}/per_subject_{tag}.csv")
    rows = []
    for (s, t), g in d.groupby(["src", "tgt"]):
        strata = [(nm, f(g.age.values)) for nm, f in BINS] + [("ASD", (g.dx == 1).values), ("TD", (g.dx == 2).values)]
        for nm, m in strata:
            n = int(m.sum())
            if n < 20:
                continue
            for cond in ("src", "raw", "std", "sitestd"):
                if f"cov_{cond}" not in g: continue
                p = g.loc[m, f"cov_{cond}"].mean()
                lo, hi = E.wilson(p, n)
                rows.append(dict(src=s, tgt=t, stratum=nm, n=n, cond=cond, cov=p, lo=lo, hi=hi))
    pd.DataFrame(rows).to_csv(f"{OUT}/strata_{tag}.csv", index=False)


# ---------------------------------------------------------------- LightGBM
def stage_lgbm(level=1):
    rows = []
    for s, t in PAIRS:
        ids, Xs, Xt, y, sites, dx = get(s, t, level)
        acc = {k: [] for k in ("src", "raw", "std")}
        for sd in [0, 1, 2]:
            cv = E.run_lgbm(Xs, Xt, y, sd)
            for k in acc: acc[k].append(cv[k])
        c = {k: np.mean(v, 0) for k, v in acc.items()}
        r = dict(src=s, tgt=t, N=len(ids))
        for k in c: r[f"lgbm_cov_{k}"] = c[k].mean()
        for k in ("raw", "std"):
            d = c[k] - c["src"]; r[f"lgbm_d_{k}"] = d.mean(); r[f"lgbm_d_{k}_lo"], r[f"lgbm_d_{k}_hi"] = E.boot_mean(d)
        rows.append(r); say(f"[lgbm] {s}->{t} src={r['lgbm_cov_src']:.3f} raw={r['lgbm_cov_raw']:.3f} std={r['lgbm_cov_std']:.3f}")
    pd.DataFrame(rows).to_csv(f"{OUT}/lgbm.csv", index=False)


# ---------------------------------------------------------------- monitors
def stage_monitor(level=1, bsz=50):
    rows = []
    for s, t in PAIRS:
        ids, Xs, Xt, y, sites, dx = get(s, t, level)
        m = E.run_monitor(Xs, Xt, y, seed=0, bsz=bsz)
        r = dict(src=s, tgt=t, N=len(ids), batch=bsz)
        for cond in ("raw", "std"):
            for stat in ("mmd", "ks"):
                thr = np.quantile(m[cond][f"{stat}_cal"], 0.95)
                r[f"{stat}_{cond}_far_test"] = float((np.array(m[cond][f"{stat}_nul"]) > thr).mean())
                r[f"{stat}_{cond}_detect"] = float((np.array(m[cond][f"{stat}_alt"]) > thr).mean())
        rows.append(r); say(f"[monitor] {s}->{t} mmd raw det={r['mmd_raw_detect']:.2f} std det={r['mmd_std_detect']:.2f} ks raw={r['ks_raw_detect']:.2f} std={r['ks_std_detect']:.2f}")
    pd.DataFrame(rows).to_csv(f"{OUT}/monitors_b{bsz}.csv", index=False)


# ---------------------------------------------------------------- mixed pipeline training
def stage_mixed(level=1):
    ids = sorted(set.intersection(*[set(P[p].dropna().index) for p in PIPES]) & set(ph.index))
    ids = [i for i in ids if i not in AD.excluded(qc, PIPES, level)]
    y = ph.loc[ids, "AGE_AT_SCAN"].values
    Xd = {p: P[p].loc[ids].values for p in PIPES}
    acc = {t: {} for t in PIPES}
    for sd in SEEDS:
        r = E.run_mixed(Xd, y, sd, PIPES)
        for t in PIPES:
            for k, v in r[t].items():
                acc[t].setdefault(k, []).append(v)
    rows = []
    for t in PIPES:
        a = {k: np.mean(v, 0) for k, v in acc[t].items()}
        row = dict(target=t, N=len(ids), within=a["within"].mean(), xfer_raw_mean=a["xfer_raw"].mean(), xfer_raw_min=a["xfer_raw"].mean(1).min(),
                   xfer_std_mean=a["xfer_std"].mean(), mix_all=a["mix_all"].mean(), mix_loo=a["mix_loo"].mean(),
                   width_all=a["w_all"].mean(), width_loo=a["w_loo"].mean(), MAE_all=a["ae_all"].mean(), MAE_loo=a["ae_loo"].mean())
        for k in ("mix_all", "mix_loo"):
            row[f"{k}_lo"], row[f"{k}_hi"] = E.boot_mean(a[k])
        rows.append(row); say(f"[mixed] target={t} N={len(ids)} within={row['within']:.3f} xfer_raw={row['xfer_raw_mean']:.3f} xfer_std={row['xfer_std_mean']:.3f} mix_all={row['mix_all']:.3f} mix_loo={row['mix_loo']:.3f}")
    pd.DataFrame(rows).to_csv(f"{OUT}/mixed_pipeline.csv", index=False)


# ---------------------------------------------------------------- leave-one-site-out
def stage_loso(level=1, top=6):
    ids = sorted(set.intersection(*[set(P[p].dropna().index) for p in PIPES]) & set(ph.index))
    ids = [i for i in ids if i not in AD.excluded(qc, PIPES, level)]
    sizes = ph.loc[ids, "SITE_ID"].value_counts()
    sites6 = list(sizes.index[:top])
    say("[loso] six largest sites (all-pipeline QC1 set):", sizes.head(top).to_dict())
    rows = []
    for s, t in PAIRS:
        ids_, Xs, Xt, y, sites, dx = get(s, t, level)
        for H in sites6:
            te = np.where(sites == H)[0]
            if len(te) < 20:
                continue
            rest = np.where(sites != H)[0]
            acc = {k: [] for k in ("src", "raw", "std", "src_site", "sitestd")}
            for sd in SEEDS:
                rng = np.random.default_rng(sd); perm = rng.permutation(rest); k = int(0.75 * len(perm)); tr, ca = perm[:k], perm[k:]
                sc_s, sc_t = E.StandardScaler().fit(Xs[tr]), E.StandardScaler().fit(Xt[tr])
                rg = E.RidgeCV(alphas=E.ALPHAS).fit(sc_s.transform(Xs[tr]), y[tr])
                q = E.cq(np.abs(y[ca] - rg.predict(sc_s.transform(Xs[ca]))))
                for nm, p in (("src", rg.predict(sc_s.transform(Xs[te]))), ("raw", rg.predict(sc_s.transform(Xt[te]))), ("std", rg.predict(sc_t.transform(Xt[te])))):
                    acc[nm].append(np.abs(y[te] - p) <= q)
                ref = np.r_[tr, te]  # held-out site's own UNLABELLED features give its site statistics
                Zs, Zt = E.site_z(Xs, sites, ref), E.site_z(Xt, sites, ref)
                rg2 = E.RidgeCV(alphas=E.ALPHAS).fit(Zs[tr], y[tr])
                q2 = E.cq(np.abs(y[ca] - rg2.predict(Zs[ca])))
                acc["src_site"].append(np.abs(y[te] - rg2.predict(Zs[te])) <= q2)
                acc["sitestd"].append(np.abs(y[te] - rg2.predict(Zt[te])) <= q2)
            r = dict(src=s, tgt=t, heldout=H, n=len(te))
            for k in acc: r[f"cov_{k}"] = np.mean(acc[k])
            rows.append(r)
        say(f"[loso] {s}->{t} done")
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/loso_by_site.csv", index=False)
    g = []
    for (s, t), d in df.groupby(["src", "tgt"]):
        w = d.n.values
        g.append(dict(src=s, tgt=t, n_total=int(w.sum()), **{f"cov_{k}": float(np.average(d[f"cov_{k}"], weights=w)) for k in ("src", "raw", "std", "src_site", "sitestd")},
                      min_site_cov_raw=float(d.cov_raw.min()), min_site_cov_std=float(d.cov_std.min())))
    pd.DataFrame(g).to_csv(f"{OUT}/loso_pooled.csv", index=False)


if __name__ == "__main__":
    stages = sys.argv[3].split(",") if len(sys.argv) > 3 else ["matchedN", "main", "strata", "atlas", "qc", "lgbm", "monitor", "mixed", "loso"]
    t0 = time.time()
    if "matchedN" in stages: stage_matched_n()
    if "main" in stages: stage_main(1, True, PAIRS, "main"); stage_strata("main") if "strata" in stages else None
    if "atlas" in stages: stage_main(1, False, ATLAS_PAIRS, "atlas")
    if "qc" in stages:
        stage_main(0, False, PAIRS, "qc0"); stage_main(2, False, PAIRS, "qc2")
    if "lgbm" in stages: stage_lgbm(1)
    if "monitor" in stages:
        for b in (25, 50, 100): stage_monitor(1, b)
    if "mixed" in stages: stage_mixed(1)
    if "loso" in stages: stage_loso(1)
    say(f"done in {time.time()-t0:.0f}s")
