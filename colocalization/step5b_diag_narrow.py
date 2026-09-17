# -*- coding: utf-8 -*-
"""Step 5b: 诊断 eQTL lead 变体 + 窄窗(±100kb)敏感性重跑"""
import json, os, subprocess, sys

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
INDIR = os.path.join(WORK, "coloc_inputs")
NARROW = os.path.join(WORK, "coloc_inputs_narrow")
os.makedirs(NARROW, exist_ok=True)

CENTERS = {"nppa": 11905974, "npr3": 32714270}

out = open(os.path.join(WORK, "step5b_log.txt"), "w", encoding="utf-8")
def P(s):
    print(s, flush=True); out.write(s + "\n"); out.flush()

# ---- 1. eQTL lead 诊断 ----
for pair in ("nppa_eqtl_af", "nppa_eqtl_stroke", "npr3_sbp_cad_vdh"):
    de = open(os.path.join(INDIR, pair + "_exposure.csv")).read().strip().split("\n")
    do = open(os.path.join(INDIR, pair + "_outcome.csv")).read().strip().split("\n")
    def top(lines, k=5):
        hdr = lines[0].split(",")
        idx = {c: i for i, c in enumerate(hdr)}
        rows = [l.split(",") for l in lines[1:] if l.strip()]
        rows = [r for r in rows if r[idx["se"]] not in ("", "nan")]
        import math
        rows = [(r[0], float(r[2]), float(r[3]), abs(float(r[2])/max(float(r[3]),1e-12))) for r in rows]
        rows.sort(key=lambda x: -x[3])
        return rows[:k]
    P(f"\n[{pair}] exposure top-5 by |z|:")
    for rsid, b, se, z in top(de):
        P(f"   {rsid}: beta={b:.4g} se={se:.4g} |z|={z:.1f}")
    P(f"[{pair}] outcome top-5 by |z|:")
    for rsid, b, se, z in top(do):
        P(f"   {rsid}: beta={b:.4g} se={se:.4g} |z|={z:.1f}")

# rs5068 / rs1421811 的 eQTL z 值
for pair, target in (("nppa_eqtl_af", "rs5068"), ("npr3_sbp_cad_vdh", "rs1421811")):
    lines = open(os.path.join(INDIR, pair + "_exposure.csv")).read().strip().split("\n")
    hdr = lines[0].split(","); idx = {c: i for i, c in enumerate(hdr)}
    for l in lines[1:]:
        r = l.split(",")
        if r[0] == target:
            P(f"\n{target} in {pair} exposure: beta={r[idx['beta']]} se={r[idx['se']]} pos={r[idx['pos']]}")

# ---- 2. 窄窗输入生成 ----
for locus, center in CENTERS.items():
    lo, hi = center - 100000, center + 100000
    variants = json.load(open(os.path.join(WORK, f"{locus}_locus_variants.json")))
    keep = {v["rsid"] for v in variants if lo <= v["pos"] <= hi}
    P(f"\n{locus} narrow window {lo}-{hi}: {len(keep)} variants")
    for f in os.listdir(INDIR):
        if f.endswith(".csv") and (f.startswith(locus + "_")):
            lines = open(os.path.join(INDIR, f)).read().strip().split("\n")
            hdr = lines[0]; rows = [l for l in lines[1:] if l.split(",")[0] in keep]
            open(os.path.join(NARROW, f), "w").write(hdr + "\n" + "\n".join(rows) + "\n")
    P(f"{locus}: narrow inputs written")

out.close()
print("DONE")
