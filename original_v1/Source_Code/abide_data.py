"""
Data layer: ABIDE I compare-surf-tools derivative tables, QC flags, phenotype.

Loaders for FS5.1, FS5.3, FS6.0 (DKT), ANTs and CIVET 2.1 follow the user's pilot script
(pilot_pipeline_swap.py). Added here: FS6.0 Desikan-Killiany (DK) table ("FS60DK") so that
atlas and software version can be separated, and QC flags from the repo's data/QC folder.

QC definitions (decided before looking at any coverage result; follow the repo notebook
notebooks/Outlier_QC_analysis.ipynb):
  QC1 (primary)  : exclude a subject if ANY of
                     - manual lenient fail   : QC_maarten != 0 in master_QC_table.csv
                     - auto low-dim outlier  : outlier_ind == 1 for the subject in either pipeline of the pair
                     - auto high-dim outlier : OD_high_dim == 1
                   Subjects absent from a QC table are NOT excluded by that table (no rating = no flag).
  QC2 (stringent): QC1 AND explicit manual stringent pass (QC_gleb in {'1','-+1'}); unrated subjects dropped.
  QC0            : no exclusions.
"""
import re
import numpy as np
import pandas as pd

PIPES5 = ["FS51", "FS53", "FS60", "ANTS", "CIVET"]
QC_NAME = {"FS51": "fs51", "FS53": "fs53", "FS60": "fs60", "FS60DK": "fs60", "ANTS": "ants", "CIVET": "civet"}


def sid(x):
    m = re.findall(r"\d{5,}", str(x))
    return int(m[-1]) if m else np.nan


def _long_names(df, pattern, hemi_fn):
    out = {}
    for c in df.columns:
        m = re.match(pattern, str(c))
        if m:
            h, r = hemi_fn(m)
            out[c] = f"{h}_{r.lower()}"
    return df[list(out)].rename(columns=out)


def load_fs5(path):
    d = pd.read_csv(path)
    X = _long_names(d, r"^([LR])_(\w+)_thickavg$", lambda m: (m[1].lower(), m[2]))
    X.index = d["SubjID"].map(sid)
    return X


def load_fs6_dkt(dirpath):
    parts = []
    for h in ["lh", "rh"]:
        d = pd.read_csv(f"{dirpath}/DKTatlas_thickness_{h}.csv", sep="\t")
        idx = d.iloc[:, 0].map(sid)
        X = _long_names(d, rf"^{h}_(\w+?)_thickness$", lambda m, h=h: (h[0], m[1]))
        X = X.drop(columns=[c for c in X.columns if "meanthickness" in c], errors="ignore")
        X.index = idx
        parts.append(X)
    return pd.concat(parts, axis=1)


def load_fs6_dk(dirpath):
    """FreeSurfer 6.0 Desikan-Killiany (aparc) thickness: same software/version as FS60 DKT, different atlas."""
    parts = []
    for h in ["lh", "rh"]:
        d = pd.read_csv(f"{dirpath}/{h}.aparc.thickness.table.csv", sep="\t")
        idx = d.iloc[:, 0].map(sid)
        X = _long_names(d, rf"^{h}_(\w+?)_thickness$", lambda m, h=h: (h[0], m[1]))
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


def load_qc(data):
    qc = f"{data}/QC"
    man = pd.read_csv(f"{qc}/manual/master_QC_table.csv", dtype={"SubjID": str, "QC_gleb": str})
    man["sid"] = man["SubjID"].map(sid)
    man = man.drop_duplicates("sid").set_index("sid")
    manual_fail = set(man.index[man["QC_maarten"].astype(float) != 0])
    stringent_pass = set(man.index[man["QC_gleb"].isin(["1", "-+1"])])
    rated = set(man.index)
    sw = pd.read_csv(f"{qc}/auto/software/software_outlier_data.csv")
    auto_low = {p: set(sw.loc[(sw.pipeline == p) & (sw.outlier_ind == 1), "SubjID"].map(sid)) for p in sw.pipeline.unique()}
    hd = pd.read_csv(f"{qc}/auto/high-dim/high_dim_outlier_data.csv")
    hd["sid"] = hd["SubjID"].map(sid)
    high_flag = set(hd.loc[hd["OD_high_dim"] == 1, "sid"])
    return dict(manual_fail=manual_fail, stringent_pass=stringent_pass, rated=rated,
                auto_low=auto_low, high_flag=high_flag)


def excluded(qc, pipes, level):
    """Set of subject IDs excluded for a comparison involving `pipes` at QC level 0/1/2."""
    if level == 0:
        return set()
    ex = set(qc["manual_fail"]) | set(qc["high_flag"])
    for p in pipes:
        ex |= qc["auto_low"].get(QC_NAME[p], set())
    return ex


def stringent_ok(qc, level):
    return qc["stringent_pass"] if level == 2 else None


def build(data, with_dk=True):
    P = {
        "FS51": load_fs5(f"{data}/FS5.1/cortical_fs5.1_measuresenigma_thickavg.csv"),
        "FS53": load_fs5(f"{data}/FS5.3/ABIDE_fs5.3_thickness.csv"),
        "FS60": load_fs6_dkt(f"{data}/FS6.0"),
        "ANTS": load_ants(f"{data}/ANTs/ABIDE_ants_thickness_data.csv"),
        "CIVET": load_civet(f"{data}/CIVET/ABIDE_civet2.1_thickness.csv", f"{data}/CIVET/DKT/DKT_parcel_map_FS_CIVET.csv"),
    }
    if with_dk:
        P["FS60DK"] = load_fs6_dk(f"{data}/FS6.0")
    for k, v in P.items():
        v = v[~v.index.isna()]
        v = v[~v.index.duplicated()]
        P[k] = v.astype(float)
    common = sorted(set.intersection(*[set(P[k].columns) for k in PIPES5]))
    if with_dk:
        common = sorted(set(common) & set(P["FS60DK"].columns))
    ph = pd.read_csv(f"{data}/ABIDE_Phenotype.csv", low_memory=False)
    ph = ph.set_index("SUB_ID")[["AGE_AT_SCAN", "DX_GROUP", "SITE_ID", "SEX"]].dropna(subset=["AGE_AT_SCAN", "DX_GROUP", "SITE_ID"])
    P = {k: v[common] for k, v in P.items()}
    return P, ph, common, load_qc(data)


def matched_ids(P, ph, qc, s, t, level):
    ids = set(P[s].dropna().index) & set(P[t].dropna().index) & set(ph.index)
    ids -= excluded(qc, [s, t], level)
    ok = stringent_ok(qc, level)
    if ok is not None:
        ids &= ok
    return sorted(ids)


if __name__ == "__main__":
    import sys
    D = sys.argv[1]
    P, ph, common, qc = build(D)
    print("common features:", len(common), common[:6])
    print({k: v.shape for k, v in P.items()})
    print("QC sizes: manual_fail", len(qc["manual_fail"]), "high", len(qc["high_flag"]),
          {k: len(v) for k, v in qc["auto_low"].items()}, "stringent_pass", len(qc["stringent_pass"]))
    for lvl in (0, 1, 2):
        print("level", lvl, {f"{s}-{t}": len(matched_ids(P, ph, qc, s, t, lvl)) for s, t in [("FS60", "CIVET"), ("FS51", "ANTS"), ("FS53", "FS60")]})
