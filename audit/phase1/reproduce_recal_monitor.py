"""Phase 1: exact re-execution of v1 recalibration (one transition, all seeds/k) and monitor (b=50) code paths.
Usage: python -I reproduce_recal_monitor.py ORIGINAL_SOURCE_DIR ORIGINAL_RESULTS_DIR DATA_DIR"""
import sys, numpy as np, pandas as pd
SRC, RES, DATA = sys.argv[1:4]; sys.path.insert(0, SRC)
import abide_data as AD, swap_engine as E
P, ph, common, qc = AD.build(DATA)
for s, t in [("FS60", "ANTS"), ("CIVET", "FS51")]:
    ids = AD.matched_ids(P, ph, qc, s, t, 1)
    Xs, Xt, y, sites = P[s].loc[ids].values, P[t].loc[ids].values, ph.loc[ids, "AGE_AT_SCAN"].values, ph.loc[ids, "SITE_ID"].values
    rc = {k: [] for k in E.KS}
    for sd in range(5):
        _, _, _, _, rec, _ = E.run_transition(Xs, Xt, y, sites, sd, rrep=100, do_wcp=True, do_reform=True, do_site=True)
        for k in E.KS: rc[k] += rec[k]["cov"]
    st = pd.read_csv(f"{RES}/recal_main.csv"); st = st[(st.src == s) & (st.tgt == t) & (st.variant == "reform_raw")].set_index("k")
    print(s, t, "recal max|diff| =", max(abs(np.mean(rc[k]) - st.loc[k, "mean_cov"]) for k in E.KS))
    m = E.run_monitor(Xs, Xt, y, seed=0, bsz=50)
    mo = pd.read_csv(f"{RES}/monitors_b50.csv"); mo = mo[(mo.src == s) & (mo.tgt == t)].iloc[0]
    for cond in ("raw", "std"):
        thr = np.quantile(m[cond]["mmd_cal"], 0.95)
        det = float((np.array(m[cond]["mmd_alt"]) > thr).mean())
        print(f"  monitor mmd {cond}: detect repro={det:.3f} stored={mo[f'mmd_{cond}_detect']:.3f}; batches per pool = {len(m[cond]['mmd_alt'])}")
