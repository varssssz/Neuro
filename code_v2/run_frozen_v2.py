"""
Phases 3-5: frozen-model conditions with per-seed, per-subject outputs, under the v1 rank (reproduction) and the
textbook rank (v2 primary). Same 60/20/20 design and test subjects as v1.
Usage: python -I code_v2/run_frozen_v2.py DATA_DIR OUT_DIR [n_procs]
Writes OUT/per_subject_seed.csv.gz and OUT/widths_by_transition.csv
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
P, ph, common, qc = AD.build(DATA)
PAIRS = [(s, t) for s in AD.PIPES5 for t in AD.PIPES5 if s != t]
ATLAS = [("FS60DK", "FS60"), ("FS60", "FS60DK"), ("FS60DK", "FS53"), ("FS53", "FS60DK"), ("FS60DK", "FS51"), ("FS51", "FS60DK")]


def one(job):
    s, t = job
    arm = "atlas" if (s, t) in ATLAS else "main"
    ids = AD.matched_ids(P, ph, qc, s, t, 1)
    Xs, Xt = P[s].loc[ids].values, P[t].loc[ids].values
    y, sites, dx = ph.loc[ids, "AGE_AT_SCAN"].values, ph.loc[ids, "SITE_ID"].values, ph.loc[ids, "DX_GROUP"].values
    frames, wrows = [], []
    for rule in ("v1", "textbook"):
        for sd in range(5):
            r = V2.run_frozen(Xs, Xt, y, sites, sd, rule=rule, do_site=(arm == "main"))
            d = dict(src=s, tgt=t, arm=arm, rule=rule, seed=sd, sid=ids, site=sites, age=y, dx=dx, width=r["width"], width_site=r["width_site"])
            for k, v in r["cov"].items(): d[f"cov_{k}"] = v.astype(np.int8)
            for k, v in r["err"].items(): d[f"err_{k}"] = np.round(v, 5)
            frames.append(pd.DataFrame(d))
            wrows.append(dict(src=s, tgt=t, arm=arm, rule=rule, seed=sd, width_src=r["width"].mean(), width_site=r["width_site"].mean()))
    print(s, t, "done", flush=True)
    return pd.concat(frames), wrows


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    with Pool(NP) as pool:
        res = pool.map(one, PAIRS + ATLAS)
    pd.concat([f for f, _ in res]).to_csv(f"{OUT}/per_subject_seed.csv.gz", index=False)
    pd.DataFrame([w for _, ww in res for w in ww]).to_csv(f"{OUT}/widths_by_transition.csv", index=False)
    print("written")
