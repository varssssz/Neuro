"""Optional: re-run the v1 LightGBM check with the textbook conformal rank (needs lightgbm==4.7.0).
Usage: python code_v2/run_lgbm_textbook.py ORIGINAL_SOURCE_DIR DATA_DIR OUT_DIR"""
import sys, os, runpy
SRC, DATA, OUT = sys.argv[1:4]
sys.path[:0] = [SRC, os.path.dirname(os.path.abspath(__file__))]
import swap_engine as E, swap_engine_v2 as V2
E.cq = lambda res, alpha=E.ALPHA: V2.conformal_quantile(res, 1 - alpha)
sys.argv = ["run_all.py", DATA, OUT, "lgbm"]
runpy.run_path(f"{SRC}/run_all.py", run_name="__main__")
