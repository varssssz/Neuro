"""Render manuscript_v2 from the template, numbers.json and result tables. Fails on any unresolved placeholder.
Usage: python -I code_v2/render_manuscript.py RESULTS_V2 MANUSCRIPT_DIR"""
import sys, re, json
import numpy as np, pandas as pd

RV2, MD = sys.argv[1:3]
N = json.load(open(f"{MD}/numbers.json"))
NAME = {"FS51": "FS5.1", "FS53": "FS5.3", "FS60": "FS6.0", "ANTS": "ANTs", "CIVET": "CIVET", "FS60DK": "FS6.0-DK"}
f3 = lambda x: f"{x:.3f}"
fp = lambda p: "<0.001" if p < 0.001 else f"{p:.3f}"

T = pd.read_csv(f"{RV2}/inference/per_transition_inference.csv"); t = T[(T.rule == "textbook") & (T.arm == "main")]
rows = ["| Source → target | N | Sites | Source | Raw | Std [95% CI] | Std − source [95% CI] | Holm p | Site-std (refit) |", "|---|---|---|---|---|---|---|---|---|"]
for _, r in t.iterrows():
    star = "†" if r.control else ""
    rows.append(f"| {NAME[r.src]} → {NAME[r.tgt]}{star} | {r.N} | {r.n_sites} | {f3(r.cov_src)} | {f3(r.cov_raw)} | {f3(r.cov_std)} [{f3(r.cov_std_clu_lo)}, {f3(r.cov_std_clu_hi)}] | "
                f"{r.d_std_vs_src:+.3f} [{r.d_std_vs_src_clu_lo:+.3f}, {r.d_std_vs_src_clu_hi:+.3f}] | {'—' if r.control else fp(r.p_std_vs_src_holm)} | {f3(r.cov_sitestd)} |")
TABLE1 = ("**Table 1.** Coverage of the 90% interval for all 20 directed transitions (QC1; textbook conformal rank; seed-averaged). "
          "Std = per-pipeline standardisation (frozen model). CIs: site-cluster bootstrap. Holm p: cluster-robust test of std − source within the family of 18 non-control transitions. "
          "Site-std = per-site standardisation with model refitting (class C). † FreeSurfer 5.3↔6.0 low-impact pair, reported separately.\n\n" + "\n".join(rows))

C = pd.read_csv(f"{RV2}/inference/category_sensitivity.csv")
order = ["atlas only (same FS6.0 run, DK vs DKT)", "FS version, same DK atlas (6.0 vs 5.3)", "within FreeSurfer: 5.3<->6.0 (designated control)", "FS version, same DK atlas (6.0 vs 5.1)",
         "within FreeSurfer: involving 5.1", "FreeSurfer<->ANTs", "FreeSurfer<->CIVET", "ANTs<->CIVET"]
lab = {order[0]: "Atlas only: FS6.0-DK ↔ FS6.0 (negative control)", order[1]: "FreeSurfer version, same atlas: FS6.0-DK ↔ FS5.3", order[2]: "FreeSurfer 5.3 ↔ 6.0 (low-impact pair)",
       order[3]: "FreeSurfer version, same atlas: FS6.0-DK ↔ FS5.1", order[4]: "Within FreeSurfer, involving 5.1", order[5]: "FreeSurfer ↔ ANTs", order[6]: "FreeSurfer ↔ CIVET", order[7]: "ANTs ↔ CIVET"}
rows = ["| Category | Transitions | Mean same-region r | Raw − source [95% CI] | Std − source [95% CI] |", "|---|---|---|---|---|"]
for c in order:
    g = C[C.category == c].set_index("contrast")
    rr, ss = g.loc["raw_vs_src"], g.loc["std_vs_src"]
    rows.append(f"| {lab[c]} | {int(rr.n_transitions)} | {rr.mean_feature_r:.2f} | {rr.mean_d:+.3f} [{rr.lo:+.3f}, {rr.hi:+.3f}] | {ss.mean_d:+.3f} [{ss.lo:+.3f}, {ss.hi:+.3f}] |")
TABLE2 = "**Table 2.** Coverage change by software-change category (mean over transitions in the category; joint site-cluster 95% CI). Atlas-arm rows are exploratory.\n\n" + "\n".join(rows)

S = pd.read_csv(f"{RV2}/recalibration/summary_by_design_k_18noncontrol.csv")
rows = ["| k | Expected coverage (theory) | Coverage, raw [range over 18] | Width, raw (y) | Coverage, std→recal | Width, std→recal (y) | Initial analysis: coverage / width |", "|---|---|---|---|---|---|---|"]
for k in (10, 20, 40, 80, 160):
    a = S[(S.design == "D5_v2split_tb") & (S.variant == "reform_raw") & (S.k == k)].iloc[0]; h = S[(S.design == "D5_v2split_tb") & (S.variant == "std_then_reform") & (S.k == k)].iloc[0]
    v = S[(S.design == "D0_v1_reproduction") & (S.variant == "reform_raw") & (S.k == k)].iloc[0]
    if not np.isfinite(a.mean_width):
        rows.append(f"| {k} | 1 (infinite interval) | — | ∞ | — | ∞ | {f3(v.mean_cov)} / {v.mean_width:.1f} |")
    else:
        rows.append(f"| {k} | {f3(a.expected_cov_theory)} | {f3(a.mean_cov)} [{f3(a.min_cov)}, {f3(a.max_cov)}] | {a.mean_width:.1f} | {f3(h.mean_cov)} | {h.mean_width:.1f} | {f3(v.mean_cov)} / {v.mean_width:.1f} |")
TABLE3 = (f"**Table 3.** Labelled target recalibration in the corrected four-way design (source width {N['rc_width_src']['value']} years). "
          "Initial analysis: recalibration pool included training subjects and the quantile was one rank higher; shown for transparency, superseded.\n\n" + "\n".join(rows))

A = pd.read_csv(f"{RV2}/age/age_bin_coverage.csv"); D = pd.read_csv(f"{RV2}/age/bin_description.csv").set_index("bin")
rows = ["| Age bin (years) | Unique subjects (sites) | Per-transition n (median) | Source | Raw | Std | Std − source [95% CI] | Mean signed error at source (y) |", "|---|---|---|---|---|---|---|---|"]
for b in ["<12", "12-15", "15-18", "18-25", ">=25"]:
    g = A[A.bin == b].set_index("quantity")
    q = lambda k: g.loc[k, "mean_over_transitions"]
    rows.append(f"| {b.replace('>=', '≥')} | {int(D.loc[b, 'n_unique_subjects'])} ({int(D.loc[b, 'n_sites'])}) | {int(D.loc[b, 'n_per_transition_median'])} | "
                f"{f3(q('cov_src'))} [{f3(g.loc['cov_src', 'lo'])}, {f3(g.loc['cov_src', 'hi'])}] | {f3(q('cov_raw'))} | {f3(q('cov_std'))} | "
                f"{q('d_std'):+.3f} [{g.loc['d_std', 'lo']:+.3f}, {g.loc['d_std', 'hi']:+.3f}] | {q('mean_signed_error_src_years'):+.1f} |")
TABLE4 = "**Table 4.** Coverage by age bin (mean over 18 non-control transitions; joint site-cluster 95% CI). The ≥25 bin spans 25–35 years.\n\n" + "\n".join(rows)

tmpl = open(f"{MD}/manuscript_v2.md.tmpl").read()
tabs = dict(TABLE1=TABLE1, TABLE2=TABLE2, TABLE3=TABLE3, TABLE4=TABLE4)
used = set()
def sub(m):
    k = m.group(1)
    if k in tabs: return tabs[k]
    if k not in N: raise KeyError(f"unresolved placeholder {{{k}}}")
    used.add(k); return N[k]["value"]
out = re.sub(r"\{([A-Za-z0-9_]+)\}", sub, tmpl)
open(f"{MD}/manuscript_v2.md", "w").write(out)
json.dump(sorted(used), open(f"{MD}/numbers_used.json", "w"), indent=0)
print("rendered;", len(used), "numbers used")

# ---- markdown hygiene: pandoc needs a blank line before a list that follows a paragraph
lines = out.split("\n"); fixed = []
for i, l in enumerate(lines):
    if re.match(r"^(- |\d+\. )", l) and fixed and fixed[-1].strip() and not re.match(r"^(\s*- |\d+\. |\s+)", fixed[-1]):
        fixed.append("")
    fixed.append(l)
out = "\n".join(fixed)
open(f"{MD}/manuscript_v2.md", "w").write(out)

# ---- docx build + table formatting
import subprocess, docx
from docx.shared import Pt, Cm
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
ref = f"{MD}/../original_v1/Manuscript/Same_scan_different_pipeline_manuscript.docx"
dx = f"{MD}/Same_scan_different_pipeline_v2.docx"
subprocess.run(["pandoc", f"{MD}/manuscript_v2.md", "-o", dx, f"--resource-path={MD}"], check=True)
d = docx.Document(dx)
sec = d.sections[0]
if sec.page_width is None:
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.left_margin = sec.right_margin = Cm(2.5); sec.top_margin = sec.bottom_margin = Cm(2.5)
usable = sec.page_width - sec.left_margin - sec.right_margin
for tb in d.tables:
    tbl = tb._tbl; pr = tbl.tblPr
    for tag in ("w:tblW", "w:tblLayout", "w:tblBorders"):
        for e in pr.findall(qn(tag)): pr.remove(e)
    w = OxmlElement("w:tblW"); w.set(qn("w:w"), str(int(usable / 635))); w.set(qn("w:type"), "dxa"); pr.append(w)
    lay = OxmlElement("w:tblLayout"); lay.set(qn("w:type"), "fixed"); pr.append(lay)
    bd = OxmlElement("w:tblBorders")
    for side in ("top", "bottom", "insideH"):
        e = OxmlElement(f"w:{side}"); e.set(qn("w:val"), "single"); e.set(qn("w:sz"), "4"); e.set(qn("w:color"), "999999"); bd.append(e)
    pr.append(bd)
    ncol = len(tb.columns)
    texts = [[c.text for c in r.cells] for r in tb.rows]
    lens = np.array([max(len(t[j]) for t in texts[1:]) for j in range(ncol)], float); lens = np.clip(lens, 6, 30)
    widths = usable * lens / lens.sum()
    grid = tbl.tblGrid
    for j, gc in enumerate(grid.findall(qn("w:gridCol"))): gc.set(qn("w:w"), str(int(widths[j] / 635)))
    for r in tb.rows:
        for j, c in enumerate(r.cells):
            c.width = int(widths[j])
            for p in c.paragraphs:
                p.paragraph_format.space_after = Pt(1)
                for run in p.runs: run.font.size = Pt(7.5)
d.save(dx)
print("docx built with", len(d.tables), "formatted tables")
