# -*- coding: utf-8 -*-
"""Step 11 前置：探测数据结构（LD 矩阵 / 变体列表 / assoc 文件）"""
import json, os

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
OUT = []

def p(*a):
    OUT.append(" ".join(str(x) for x in a))

# 1. locus variants 结构
v = json.load(open(os.path.join(WORK, "nppa_locus_variants.json")))
p("nppa_locus_variants:", len(v), "entries; first:", json.dumps(v[0], ensure_ascii=False))
keysnps = [x for x in v if x.get("rsid") in ("rs5068", "rs5066", "rs198411")]
p("keysnp entries:", json.dumps(keysnps, ensure_ascii=False))

v3 = json.load(open(os.path.join(WORK, "npr3_locus_variants.json")))
k3 = [x for x in v3 if x.get("rsid") == "rs1421811"]
p("npr3 rs1421811:", json.dumps(k3, ensure_ascii=False))

# 2. LD variants 结构
lv = json.load(open(os.path.join(WORK, "ld_variants_nppa.json")))
if isinstance(lv, list):
    p("ld_variants_nppa: list", len(lv), "; first:", json.dumps(lv[0], ensure_ascii=False))
else:
    p("ld_variants_nppa: dict keys", list(lv.keys())[:5])
    p(json.dumps(lv, ensure_ascii=False)[:400])

# 3. LD csv 结构
with open(os.path.join(WORK, "ld_nppa.csv")) as f:
    head = [next(f) for _ in range(2)]
p("ld_nppa.csv first 2 lines (truncated 300 chars):")
for h in head:
    p("  ", h[:300].rstrip())
p("ld_nppa.csv total lines:", sum(1 for _ in open(os.path.join(WORK, "ld_nppa.csv"))))

# 4. assoc 文件结构
a = json.load(open(os.path.join(WORK, "nppa__ebi-a-GCST006061_assoc.json")))
p("AF assoc entries:", len(a), "; rs5068:", json.dumps(a.get("rs5068")), "; rs5066:", json.dumps(a.get("rs5066")))

a2 = json.load(open(os.path.join(WORK, "nppa__ebi-a-GCST90012082_assoc.json")))
p("SCALLOP assoc rs5068:", json.dumps(a2.get("rs5068")))

# 5. coloc_inputs 目录
d = os.path.join(WORK, "coloc_inputs")
if os.path.isdir(d):
    p("coloc_inputs:", sorted(os.listdir(d))[:10])
    fp = os.path.join(d, sorted(os.listdir(d))[0])
    p("first input file head:", open(fp).read()[:200].replace("\n", " | "))

with open(os.path.join(WORK, "step11a_output.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT))
