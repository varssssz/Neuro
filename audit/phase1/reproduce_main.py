"""
Phase 1: re-execute the v1 code path for transitions_main.csv (src/raw/std/sitestd/wcp columns; recalibration
draws skipped, rrep=0) and compare to the stored v1 file. Uses run_all.summarize() unmodified.
Usage: python -I reproduce_main.py ORIGINAL_SOURCE_DIR ORIGINAL_RESULTS_DIR DATA_DIR OUT_DIR
"""
import sys, types
import numpy as np
import pandas as pd

SRC, RES, DATA, OUT = sys.argv[1:5]
sys.path.insert(0, SRC)
import abide_data as AD
import swap_engine as E

# import summarize() from run_all.py without executing its module-level data loading / CLI
code = open(f"{SRC}/run_all.py").read()
start = code.index("def summarize("); end = code.index("def stage_main(")
ns = {"np": np, "E": E}
exec(code[start:end], ns)
summarize = ns["summarize"]

P, ph, common, qc = AD.build(DATA)
PIPES = AD.PIPES5
stored = pd.read_csv(f"{RES}/transitions_main.csv")
rows = []
for s, t in [(a, b) for a in PIPES for b in PIPES if a != b]:
    ids = AD.matched_ids(P, ph, qc, s, t, 1)
    Xs, Xt = P[s].loc[ids].values, P[t].loc[ids].values
    y, sites = ph.loc[ids, "AGE_AT_SCAN"].values, ph.loc[ids, "SITE_ID"].values
    accs = {k: [] for k in ["src", "raw", "std", "src_site", "sitestd", "wcp", "wcp_reg"]}
    aes = {k: [] for k in ["src", "raw", "std", "sitestd"]}
    qvs, wcps, wcprs, qss = [], [], [], []
    fcs = {"src": [], "raw": [], "std": []}
    for sd in range(5):
        cov, ae, qv, wcpw, rec, ex = E.run_transition(Xs, Xt, y, sites, sd, rrep=0, do_wcp=True, do_reform=False, do_site=True)
        for k in accs: accs[k].append(cov[k])
        for k in aes: aes[k].append(ae[k])
        qvs.append(qv); wcps.append(wcpw); wcprs.append(ex["wcp_reg_w"]); qss.append(ex["qsite"])
        for k in fcs: fcs[k] += ex["foldcov"][k]
    cov = {k: np.array(v) for k, v in accs.items()}
    ae = {k: np.array(v).mean(0) for k, v in aes.items()}
    sdv = {k: float(np.std(cov[k].mean(1))) for k in ("src", "raw", "std")}
    out, _ = summarize(s, t, 1, cov, ae, np.mean(qvs, 0), np.concatenate(wcps), True, seed_sd=sdv,
                       wcpw_reg=np.concatenate(wcprs), qsite=np.mean(qss, 0), foldcov=fcs)
    rows.append(out)
    print(s, t, "done", flush=True)
new = pd.DataFrame(rows)
new.to_csv(f"{OUT}/transitions_main_reproduced.csv", index=False)
m = stored.merge(new, on=["src", "tgt"], suffixes=("_stored", "_repro"))
cols = [c for c in new.columns if c not in ("src", "tgt", "qc_level") and f"{c}_stored" in m]
diff = pd.DataFrame([dict(column=c, max_abs_diff=float(np.nanmax(np.abs(m[f"{c}_stored"].astype(float) - m[f"{c}_repro"].astype(float)))))
                     for c in cols])
diff.to_csv(f"{OUT}/transitions_main_repro_diff.csv", index=False)
print(diff.to_string())
