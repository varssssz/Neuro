"""
Phase 1 baseline audit: reconstruct data flow and subject roles from the ORIGINAL v1 code + seeds,
and check them against the stored v1 outputs. Read-only with respect to original_v1/.

Usage: python -I split_integrity.py ORIGINAL_SOURCE_DIR ORIGINAL_RESULTS_DIR DATA_DIR OUT_DIR
"""
import sys, json
import numpy as np
import pandas as pd

SRC, RES, DATA, OUT = sys.argv[1:5]
sys.path.insert(0, SRC)
import abide_data as AD          # v1 code, imported unmodified
import swap_engine as E

out = {}
P_raw = {}
# ------------------------------------------------------------ 1. loader-level accounting (re-implements build()'s cleaning steps to count what they drop)
loaders = {
    "FS51": lambda: AD.load_fs5(f"{DATA}/FS5.1/cortical_fs5.1_measuresenigma_thickavg.csv"),
    "FS53": lambda: AD.load_fs5(f"{DATA}/FS5.3/ABIDE_fs5.3_thickness.csv"),
    "FS60": lambda: AD.load_fs6_dkt(f"{DATA}/FS6.0"),
    "ANTS": lambda: AD.load_ants(f"{DATA}/ANTs/ABIDE_ants_thickness_data.csv"),
    "CIVET": lambda: AD.load_civet(f"{DATA}/CIVET/ABIDE_civet2.1_thickness.csv", f"{DATA}/CIVET/DKT/DKT_parcel_map_FS_CIVET.csv"),
    "FS60DK": lambda: AD.load_fs6_dk(f"{DATA}/FS6.0"),
}
acct = []
for k, f in loaders.items():
    v = f()
    n0 = len(v); nan_id = int(v.index.isna().sum()); v = v[~v.index.isna()]
    dup = int(v.index.duplicated().sum()); v = v[~v.index.duplicated()]
    acct.append(dict(pipeline=k, rows_loaded=n0, rows_without_subject_id=nan_id, duplicate_ids_dropped_keep_first=dup,
                     n_features_loaded=v.shape[1], rows_with_any_nan=int(v.astype(float).isna().any(axis=1).sum())))
pd.DataFrame(acct).to_csv(f"{OUT}/loader_accounting.csv", index=False)

P, ph, common, qc = AD.build(DATA)
out["n_common_features"] = len(common)
out["feature_order_identical_across_pipelines"] = all(list(P[k].columns) == list(common) for k in P)
pd.Series(common, name="feature").to_csv(f"{OUT}/feature_schema_62.csv", index=False)
out["qc_set_sizes"] = dict(manual_fail=len(qc["manual_fail"]), high_dim_flag=len(qc["high_flag"]),
                           stringent_pass=len(qc["stringent_pass"]), **{f"auto_low_{p}": len(s) for p, s in qc["auto_low"].items()})
out["qc_pipeline_names_matched"] = all(AD.QC_NAME[p] in qc["auto_low"] for p in AD.QC_NAME)

# ------------------------------------------------------------ 2. feature alignment: is each feature most correlated with its namesake?
PIPES = AD.PIPES5
PAIRS = [(s, t) for s in PIPES for t in PIPES if s != t]
ATLAS = [("FS60DK", "FS60"), ("FS60", "FS60DK"), ("FS60DK", "FS53"), ("FS53", "FS60DK"), ("FS60DK", "FS51"), ("FS51", "FS60DK")]
CONTROL = {("FS60", "FS53"), ("FS53", "FS60")}
align = []
for s, t in [p for p in PAIRS if p[0] < p[1]] + [("FS60DK", "FS60"), ("FS60DK", "FS53"), ("FS60DK", "FS51")]:
    ids = AD.matched_ids(P, ph, qc, s, t, 1)
    A, B = P[s].loc[ids].values, P[t].loc[ids].values
    Az, Bz = (A - A.mean(0)) / A.std(0), (B - B.mean(0)) / B.std(0)
    C = Az.T @ Bz / len(ids)                     # 62x62 cross-correlation
    diag = np.diag(C)
    best = C.argmax(1)
    align.append(dict(pair=f"{s}-{t}", n=len(ids), mean_same_region_r=diag.mean(), min_same_region_r=diag.min(),
                      min_region=common[int(diag.argmin())], n_regions_r_below_0_3=int((diag < 0.3).sum()),
                      n_regions_namesake_not_best_match=int((best != np.arange(len(common))).sum()),
                      mismatched=";".join(f"{common[i]}->{common[best[i]]}" for i in np.where(best != np.arange(len(common)))[0])))
pd.DataFrame(align).to_csv(f"{OUT}/feature_alignment.csv", index=False)

# ------------------------------------------------------------ 3. matched samples vs stored outputs
mn = pd.read_csv(f"{RES}/matched_N_summary.csv")
tm = pd.read_csv(f"{RES}/transitions_main.csv")
ps = pd.read_csv(f"{RES}/per_subject_main.csv")
ta = pd.read_csv(f"{RES}/transitions_atlas.csv")
chk = []
for s, t in PAIRS + ATLAS:
    r = dict(src=s, tgt=t, control=(s, t) in CONTROL, arm="atlas" if (s, t) in ATLAS else "main")
    for lvl in (0, 1, 2):
        ids = AD.matched_ids(P, ph, qc, s, t, lvl)
        r[f"N_QC{lvl}_recomputed"] = len(ids)
        r[f"N_QC{lvl}_stored"] = int(mn.loc[(mn.src == s) & (mn.tgt == t), ["N_noQC", "N_QC1", "N_QC2"][lvl]].iloc[0])
    ids = AD.matched_ids(P, ph, qc, s, t, 1)
    stored_tab = ta if r["arm"] == "atlas" else tm
    r["N_in_results_table"] = int(stored_tab.loc[(stored_tab.src == s) & (stored_tab.tgt == t), "N"].iloc[0])
    if r["arm"] == "main":
        g = ps[(ps.src == s) & (ps.tgt == t)]
        r["per_subject_ids_identical_and_same_order"] = bool(list(g.sid) == list(ids))
        r["per_subject_age_site_dx_match_phenotype"] = bool(np.allclose(g.age.values, ph.loc[ids, "AGE_AT_SCAN"].values)
                                                             and (g.site.values == ph.loc[ids, "SITE_ID"].values).all()
                                                             and (g.dx.values == ph.loc[ids, "DX_GROUP"].values).all())
    r["n_sites"] = int(ph.loc[ids, "SITE_ID"].nunique())
    chk.append(r)
chk = pd.DataFrame(chk)
chk.to_csv(f"{OUT}/matched_samples_check.csv", index=False)

# ------------------------------------------------------------ 4. split roles: reconstruct, export, check
role_rows, integ = [], []
for s, t in PAIRS + ATLAS:
    ids = np.array(AD.matched_ids(P, ph, qc, s, t, 1))
    n = len(ids)
    sites = ph.loc[ids, "SITE_ID"].values
    for seed in range(5):
        roles = np.full((5, n), "", dtype=object)
        test_count = np.zeros(n, int)
        for f, tr, ca, te in E.folds(n, seed):
            st, sc, se = set(tr), set(ca), set(te)
            pool = np.r_[tr, ca]
            # sites falling back to global stats in site_z (fewer than MIN_SITE training subjects)
            tr_site_counts = pd.Series(sites[tr]).value_counts()
            fallback_sites = [x for x in np.unique(sites) if tr_site_counts.get(x, 0) < E.MIN_SITE]
            integ.append(dict(src=s, tgt=t, seed=seed, fold=f, n=n, n_train=len(tr), n_cal=len(ca), n_test=len(te),
                              train_cal_overlap=len(st & sc), train_test_overlap=len(st & se), cal_test_overlap=len(sc & se),
                              union_is_all=len(st | sc | se) == n,
                              recal_pool_size=len(pool), recal_pool_from_train=len(st & set(pool)),
                              recal_pool_test_overlap=len(set(pool) & se),
                              target_std_stats_subjects_same_as_source_train=True,   # sc_t = StandardScaler().fit(Xt[tr])
                              site_fallback_sites=len(fallback_sites),
                              test_subjects_in_fallback_sites=int(np.isin(sites[te], fallback_sites).sum())))
            roles[f, tr], roles[f, ca], roles[f, te] = "train", "cal", "test"
            test_count[te] += 1
        assert (test_count == 1).all(), (s, t, seed)
        for i in range(n):
            role_rows.append((s, t, seed, int(ids[i]), *roles[:, i]))
pd.DataFrame(role_rows, columns=["src", "tgt", "seed", "sid", "fold0", "fold1", "fold2", "fold3", "fold4"]).to_csv(
    f"{OUT}/subject_roles_main_and_atlas_QC1.csv.gz", index=False)
integ = pd.DataFrame(integ)
integ.to_csv(f"{OUT}/split_integrity_by_fold.csv", index=False)
out["every_subject_tested_exactly_once_per_seed"] = True   # asserted above
out["max_train_cal_overlap"] = int(integ.train_cal_overlap.max())
out["max_train_test_overlap"] = int(integ.train_test_overlap.max())
out["max_cal_test_overlap"] = int(integ.cal_test_overlap.max())
out["recal_pool_share_from_train_mean"] = float((integ.recal_pool_from_train / integ.recal_pool_size).mean())
out["max_recal_pool_test_overlap"] = int(integ.recal_pool_test_overlap.max())
out["fold_sizes_train_cal_test_ranges"] = {c: [int(integ[c].min()), int(integ[c].max())] for c in ("n_train", "n_cal", "n_test")}
out["share_folds_with_any_site_fallback"] = float((integ.site_fallback_sites > 0).mean())
out["share_test_subjects_in_fallback_sites"] = float(integ.test_subjects_in_fallback_sites.sum() / integ.n_test.sum())

# ------------------------------------------------------------ 5. LOSO and mixed-pipeline designs
lo = []
ids_all = sorted(set.intersection(*[set(P[p].dropna().index) for p in PIPES]) & set(ph.index))
ids_all = [i for i in ids_all if i not in AD.excluded(qc, PIPES, 1)]
sizes = ph.loc[ids_all, "SITE_ID"].value_counts()
sites6 = list(sizes.index[:6])
out["mixed_and_loso_site_selection_N"] = len(ids_all)
out["loso_sites"] = {k: int(v) for k, v in sizes.head(6).items()}
for s, t in PAIRS:
    ids = np.array(AD.matched_ids(P, ph, qc, s, t, 1)); sites = ph.loc[ids, "SITE_ID"].values
    for H in sites6:
        te = np.where(sites == H)[0]
        if len(te) < 20:
            continue
        rest = np.where(sites != H)[0]
        for sd in range(5):
            perm = np.random.default_rng(sd).permutation(rest); k = int(0.75 * len(perm)); tr, ca = perm[:k], perm[k:]
            ref = np.r_[tr, te]
            lo.append(dict(src=s, tgt=t, heldout=H, seed=sd, n_test=len(te), train_test_overlap=len(set(tr) & set(te)),
                           cal_test_overlap=len(set(ca) & set(te)), site_stats_use_heldout_test_features=len(set(ref) & set(te)) == len(te)))
lo = pd.DataFrame(lo); lo.to_csv(f"{OUT}/loso_integrity.csv", index=False)
out["loso_max_overlap_train_or_cal_with_test"] = int(max(lo.train_test_overlap.max(), lo.cal_test_overlap.max()))
out["loso_site_stats_transductive_all_rows"] = bool(lo.site_stats_use_heldout_test_features.all())
mx = pd.read_csv(f"{RES}/mixed_pipeline.csv")
out["mixed_N_recomputed_vs_stored"] = [len(ids_all), int(mx.N.iloc[0])]

json.dump(out, open(f"{OUT}/integrity_summary.json", "w"), indent=2, default=str)
print(json.dumps(out, indent=2, default=str))
print(chk.drop(columns=["arm"]).to_string())
print(pd.read_csv(f"{OUT}/loader_accounting.csv").to_string())
print(pd.read_csv(f"{OUT}/feature_alignment.csv").drop(columns="mismatched").round(3).to_string())
