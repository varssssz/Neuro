"""
Phase 2: recalibration repair. Writes results_v2/recalibration/.
Usage: python -I code_v2/run_recalibration_v2.py DATA_DIR OUT_DIR [n_procs]
Designs (each: 20 transitions x 5 seeds x 5 folds x 100 draws per k):
  D0_v1_reproduction : v1 split (60/20/20), pool = train+cal, v1 rank      -> must equal original recal_main.csv
  D1_v1split_trca_tb : v1 split, pool = train+cal, textbook rank           -> effect of the quantile rule alone
  D2_v1split_ca_v1   : v1 split, pool = source-cal only, v1 rank           -> effect of removing training subjects (pool overlaps source cal)
  D3_v1split_ca_tb   : v1 split, pool = source-cal only, textbook rank
  D4_v2split_v1rank  : v2 four-way split (40/16/24/20), disjoint pool, v1 rank
  D5_v2split_tb      : v2 four-way split, disjoint pool, textbook rank     -> PRIMARY corrected analysis
"""
import sys, os
import numpy as np
import pandas as pd
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import abide_data as AD
import swap_engine_v2 as V2

DATA, OUT = sys.argv[1], sys.argv[2]
NP = int(sys.argv[3]) if len(sys.argv) > 3 else 2
DESIGNS = {
    "D0_v1_reproduction": dict(split="v1", pool="tr+ca", rule="v1"),
    "D1_v1split_trca_tb": dict(split="v1", pool="tr+ca", rule="textbook"),
    "D2_v1split_ca_v1": dict(split="v1", pool="ca", rule="v1"),
    "D3_v1split_ca_tb": dict(split="v1", pool="ca", rule="textbook"),
    "D4_v2split_v1rank": dict(split="v2", pool="recal", rule="v1"),
    "D5_v2split_tb": dict(split="v2", pool="recal", rule="textbook"),
}
CONTROL = {("FS60", "FS53"), ("FS53", "FS60")}
P, ph, common, qc = AD.build(DATA)
PAIRS = [(s, t) for s in AD.PIPES5 for t in AD.PIPES5 if s != t]


def one(pair):
    s, t = pair
    ids = AD.matched_ids(P, ph, qc, s, t, 1)
    Xs, Xt, y = P[s].loc[ids].values, P[t].loc[ids].values, ph.loc[ids, "AGE_AT_SCAN"].values
    rows, brows = [], []
    for dn, d in DESIGNS.items():
        acc = {k: dict(cov=[], width=[], cov_h=[], width_h=[]) for k in V2.KS}
        b = dict(cov_src=[], cov_raw=[], ae_src=[], ae_raw=[], width=[], n_train=[], n_cal=[], n_pool=[], n_test=[], overlap=[])
        for sd in range(5):
            out, base = V2.run_recal(Xs, Xt, y, sd, d)
            for k in V2.KS:
                for kk in acc[k]: acc[k][kk] += out[k][kk]
            for kk in ("cov_src", "cov_raw", "ae_src", "ae_raw", "width"): b[kk].append(base[kk].mean())
            for kk in ("n_train", "n_cal", "n_pool", "n_test", "overlap"): b[kk] += base[kk]
        brows.append(dict(design=dn, src=s, tgt=t, N=len(ids), control=(s, t) in CONTROL, cov_src=np.mean(b["cov_src"]),
                          cov_raw=np.mean(b["cov_raw"]), MAE_src=np.mean(b["ae_src"]), MAE_raw=np.mean(b["ae_raw"]),
                          width_src=np.mean(b["width"]), n_train_mean=np.mean(b["n_train"]), n_cal_mean=np.mean(b["n_cal"]),
                          n_pool_min=int(np.min(b["n_pool"])), n_pool_mean=np.mean(b["n_pool"]), n_test_mean=np.mean(b["n_test"]),
                          pool_test_overlap_max=int(np.max(b["overlap"]))))
        for k in V2.KS:
            for variant, ck, wk in (("reform_raw", "cov", "width"), ("std_then_reform", "cov_h", "width_h")):
                c, w = np.array(acc[k][ck]), np.array(acc[k][wk])
                if len(c) == 0:
                    rows.append(dict(design=dn, src=s, tgt=t, control=(s, t) in CONTROL, variant=variant, k=k, n_draws=0)); continue
                fin = np.isfinite(w)
                rows.append(dict(design=dn, src=s, tgt=t, control=(s, t) in CONTROL, variant=variant, k=k, n_draws=len(c),
                                 mean_cov=c.mean(), median_cov=np.median(c), cov_p025=np.percentile(c, 2.5), cov_p975=np.percentile(c, 97.5),
                                 frac_below_088=float((c < 0.88).mean()), frac_inf_width=float((~fin).mean()),
                                 mean_width_finite=float(w[fin].mean()) if fin.any() else np.inf,
                                 median_width=float(np.median(w)), expected_cov_theory=V2.expected_recal_coverage(k, rule=d["rule"])))
    print(s, t, "done", flush=True)
    return rows, brows


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    with Pool(NP) as pool:
        res = pool.map(one, PAIRS)
    R = pd.DataFrame([r for rr, _ in res for r in rr]); B = pd.DataFrame([b for _, bb in res for b in bb])
    R.to_csv(f"{OUT}/recal_by_design_transition_k.csv", index=False)
    B.to_csv(f"{OUT}/base_model_by_design.csv", index=False)
    print("written")
