# -*- coding: utf-8 -*-
"""检查关键SNP是否在变体列表, 缺失则用OpenGWAS已知等位补录(多等位基因被双等位过滤跳过)"""
import json

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
OUT = open(WORK + r"\check_keysnps.txt", "w", encoding="utf-8")

# OpenGWAS 报告的等位(GRCh37): rs5068 G/A af~0.054-0.059; rs1421811 G/C af~0.35-0.374
EXTRA = {
    "nppa": dict(rsid="rs5068", chrom="1", pos=11905974, ref="A", alt="G", af=0.0589),
    "npr3": dict(rsid="rs1421811", chrom="5", pos=32714270, ref="C", alt="G", af=0.374),
}

for name in ("nppa", "npr3"):
    v = json.load(open(WORK + f"\\{name}_locus_variants.json"))
    tgt = EXTRA[name]["rsid"]
    hit = [x for x in v if x["rsid"] == tgt]
    if hit:
        OUT.write(f"{tgt}: present {hit[0]}\n")
    else:
        v.append(EXTRA[name])
        v.sort(key=lambda x: x["pos"])
        json.dump(v, open(WORK + f"\\{name}_locus_variants.json", "w"))
        OUT.write(f"{tgt}: ABSENT -> added manually {EXTRA[name]}\n")
    OUT.write(f"{name}: total {len(v)} variants\n")
OUT.close()
print("done")
