# -*- coding: utf-8 -*-
"""QC: prot-a-2076 (ANP) vs prot-a-2078 (NT-proBNP) 是否为不同数据"""
import json
WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
OUT = open(WORK + r"\step9e_output.txt", "w", encoding="utf-8")
a = json.load(open(WORK + r"\nppa__prot-a-2076_assoc.json"))
b = json.load(open(WORK + r"\nppa__prot-a-2078_assoc.json"))
shared = sorted(set(a) & set(b))
diff = sum(1 for r in shared if abs((a[r]["beta"] or 0) - (b[r]["beta"] or 0)) > 1e-12)
OUT.write(f"shared={len(shared)}  beta-differ={diff} ({diff/len(shared)*100:.1f}%)\n")
for r in shared[:5]:
    OUT.write(f"{r}: 2076 beta={a[r]['beta']} se={a[r]['se']} | 2078 beta={b[r]['beta']} se={b[r]['se']}\n")
OUT.close()
