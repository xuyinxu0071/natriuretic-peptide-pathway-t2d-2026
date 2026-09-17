# -*- coding: utf-8 -*-
"""step11b: 排查 rs1421811 在 NPR3 各 assoc 文件中的状态"""
import json, os

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
OUT = []
def p(*a): OUT.append(" ".join(str(x) for x in a))

v3 = {x["rsid"]: x for x in json.load(open(os.path.join(WORK, "npr3_locus_variants.json")))}
p("variant list rs1421811:", json.dumps(v3.get("rs1421811")))
# 也查 rs3828589 是否在列表
p("variant list rs3828589:", json.dumps(v3.get("rs3828589")))

for ds in ["eqtl-a-ENSG00000113389", "ieu-b-38", "ieu-a-7", "ebi-a-GCST005195",
           "ebi-a-GCST011364", "finn-b-I9_MI"]:
    f = os.path.join(WORK, f"npr3__{ds}_assoc.json")
    if not os.path.exists(f):
        p(ds, ": file missing"); continue
    a = json.load(open(f))
    p(f"{ds}: rs1421811 -> {json.dumps(a.get('rs1421811'))}; rs3828589 -> {json.dumps(a.get('rs3828589'))}")

# LD 矩阵中两者的位置关系
lv = json.load(open(os.path.join(WORK, "ld_variants_npr3.json")))
p("ld_variants rs1421811:", json.dumps(lv.get("rs1421811")))
p("ld_variants rs3828589:", json.dumps(lv.get("rs3828589")))

with open(os.path.join(WORK, "step11b_output.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT))
