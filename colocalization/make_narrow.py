# -*- coding: utf-8 -*-
"""从 step5_coloc.R 生成窄窗版 step5c_coloc_narrow.R (避免 PowerShell 编码破坏)"""
src = open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis\step5_coloc.R", encoding="utf-8").read()
out = src.replace("coloc_inputs", "coloc_inputs_narrow").replace("coloc_results", "coloc_results_narrow")
open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis\step5c_coloc_narrow.R", "w", encoding="utf-8", newline="\n").write(out)
print("written")
