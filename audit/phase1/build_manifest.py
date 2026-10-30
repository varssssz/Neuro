"""Phase 1: results provenance manifest (file level + headline-number level). Usage: python -I build_manifest.py ORIGINAL_V1_DIR OUT_DIR"""
import sys, hashlib, csv
V1, OUT = sys.argv[1:3]
R = f"{V1}/Experimental_Results"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]
files = [
 ("transitions_main.csv","run_all.stage_main(level=1, extras=True, PAIRS, 'main')","swap_engine.run_transition + run_all.summarize","0-4","QC1","20 main transitions","REPRODUCED exactly (all 43 numeric cols; Phase 1 reproduce_main.py)"),
 ("per_subject_main.csv","stage_main(...'main')","run_transition (seed-mean indicators)","0-4","QC1","20","IDs/order/age/site/dx verified vs phenotype; values implied by transitions_main reproduction"),
 ("recal_main.csv","stage_main(...'main') do_reform","run_transition do_reform block (pool=np.r_[tr,ca])","0-4 (+100 draws/fold)","QC1","20 x k{10,20,40,80,160} x 2 variants","REPRODUCED exactly for FS60->ANTS and CIVET->FS51 (reform_raw); protocol issue: recal pool 75% training subjects"),
 ("strata_main.csv","run_all.stage_strata('main')","per_subject_main + swap_engine.wilson","-","QC1","20 x strata n>=20","NOT re-executed; derived deterministically from per_subject_main; manuscript Table 2 recomputed from per_subject_main and matches"),
 ("transitions_atlas.csv / per_subject_atlas.csv","stage_main(1, extras=False, ATLAS_PAIRS, 'atlas')","run_transition (src/raw/std only)","0-4","QC1","6 atlas pairs","N verified; values not re-executed"),
 ("transitions_qc0.csv / per_subject_qc0.csv","stage_main(0, False, PAIRS, 'qc0')","run_transition","0-4","QC0","20","N verified; values not re-executed"),
 ("transitions_qc2.csv / per_subject_qc2.csv","stage_main(2, False, PAIRS, 'qc2')","run_transition","0-4","QC2","20","N verified; values not re-executed"),
 ("lgbm.csv","run_all.stage_lgbm(1)","swap_engine.run_lgbm","0-2","QC1","20","NOT re-executable here (lightgbm unavailable: PyPI blocked)"),
 ("monitors_b25/b50/b100.csv","run_all.stage_monitor(1, b)","swap_engine.run_monitor","0","QC1","20","b50 REPRODUCED exactly for 2 transitions"),
 ("mixed_pipeline.csv","run_all.stage_mixed(1)","swap_engine.run_mixed","0-4","QC1 (all-5-pipeline intersection)","5 targets, N=669","N verified; values not re-executed"),
 ("loso_by_site.csv / loso_pooled.csv","run_all.stage_loso(1)","inline in stage_loso + swap_engine.site_z","0-4","QC1","20 x 6 sites","split disjointness verified; values not re-executed"),
 ("matched_N_by_site.csv / matched_N_summary.csv","run_all.stage_matched_n()","abide_data.matched_ids","-","QC0/1/2","26 pairs","REPRODUCED exactly (all 26 pairs x 3 QC levels)"),
 ("Figures_and_Tables/*","make_figs.py","reads transitions_main, recal_main, strata_main","-","-","-","not re-rendered in Phase 1"),
]
with open(f"{OUT}/RESULTS_PROVENANCE_MANIFEST.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["result_file","stage_call","code_path","seeds","qc_level","units","phase1_status","sha256_16"])
    for r in files:
        first=r[0].split(" / ")[0]; p=f"{V1}/{first}" if first.startswith("Figures") is False and "/" in first else f"{R}/{first}"
        try: h=sha(p)
        except Exception: h="(multiple)"
        w.writerow([*r,h])
heads = [
 ("Source coverage, mean 18 non-control","0.910","transitions_main.csv","cov_src","mean over control==False","18"),
 ("Raw transfer coverage, mean (range 0.001-0.876)","0.489","transitions_main.csv","cov_raw","mean/min/max, control==False","18"),
 ("Per-pipeline std coverage, mean; 12/18 < 0.88; worst 0.757","0.853","transitions_main.csv","cov_std","mean; count<0.88","18"),
 ("Per-site std coverage, mean (0.875-0.923); width 14.9 vs 12.1","0.897","transitions_main.csv","cov_sitestd, width_src_site, width_src","mean","18"),
 ("Near-identity control raw coverage","0.899 / 0.908","transitions_main.csv","cov_raw","control==True","2"),
 ("McNemar std: 16/18 p<0.05, 11 after Bonferroni","16 / 11","transitions_main.csv","mcnemar_p_std","seed 0 only; threshold 0.05/18 applied outside code","18"),
 ("WCP infinite share 90% / reg 79%","0.90 / 0.79","transitions_main.csv","wcp_frac_inf, wcp_reg_frac_inf","mean","18 (verified: 0.896 / 0.794)"),
 ("Recalibration coverage k=10..160","0.813/0.900/0.948/0.922/0.909","recal_main.csv","mean_cov (variant=reform_raw)","mean over control==False","18"),
 ("Age-bin coverage (Table 2)","see table","per_subject_main.csv","cov_src/raw/std/sitestd by age","subject-weighted pooled, control==False","18"),
 ("LightGBM raw/std mean","0.716 / 0.854","lgbm.csv","lgbm_cov_raw/std","mean over ALL","20 (inconsistent)"),
 ("Ridge comparator in LGBM sentence","0.531 / 0.859","transitions_main.csv","cov_raw/std","mean over ALL","20 (inconsistent)"),
 ("LOSO src/raw/std/sitestd/refit","0.867/0.504/0.818/0.910/0.919","loso_pooled.csv","cov_*","mean over ALL","20 (inconsistent; 18 gives 0.865/0.462/0.810/0.908/0.917)"),
 ("LOSO worst-site std","0.672","loso_pooled.csv","min_site_cov_std","mean over ALL","20 (18 gives 0.666)"),
 ("MMD std detection, non-control","0.139","monitors_b50.csv","mmd_std_detect","mean, control==False","18"),
 ("Mixed-pipeline mix_loo range","0.216-0.906","mixed_pipeline.csv","mix_loo","per target","5 targets"),
 ("QC0 / QC2 std range","0.771-0.923 / 0.719-0.957","transitions_qc0/qc2.csv","cov_std","min/max, control==False","18"),
]
with open(f"{OUT}/HEADLINE_NUMBER_PROVENANCE.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["manuscript_claim","value","file","column","aggregation","denominator"]); w.writerows(heads)
print("ok")
