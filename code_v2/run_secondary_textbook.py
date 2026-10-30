"""
Re-execute ORIGINAL v1 secondary stages (QC0/QC2 sensitivity, mixed-pipeline, LOSO) with ONE change:
swap_engine.cq is replaced by the textbook conformal quantile (code_v2/swap_engine_v2.conformal_quantile).
Usage: python -I code_v2/run_secondary_textbook.py ORIGINAL_SOURCE_DIR DATA_DIR OUT_DIR
"""
import sys, os, runpy
SRC, DATA, OUT = sys.argv[1:4]
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SRC); sys.path.insert(1, HERE)
import swap_engine as E          # v1 engine
import swap_engine_v2 as V2
E.cq = lambda res, alpha=E.ALPHA: V2.conformal_quantile(res, 1 - alpha)
sys.argv = ["run_all.py", DATA, OUT, "qc,mixed,loso"]
runpy.run_path(f"{SRC}/run_all.py", run_name="__main__")
