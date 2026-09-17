# -*- coding: utf-8 -*-
"""
Step 9b: 蛋白层暴露 × 已有结局 的等位协调 + rs5068/rs5066 直接查表
================================================================================
新 coloc 对 (暴露 = NT-proBNP/ANP pQTL, 结局 = 已有 6 个 CVD 数据集)
另输出: rs5068 / rs5066 在新数据集中的 beta/se/p/eaf 直接查表 (含 BBJ 东亚)
"""
import json, os, sys

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
sys.path.insert(0, WORK)
from step4_harmonize import harmonize  # 复用 step4 的协调规则

INDIR = os.path.join(WORK, "coloc_inputs")
LOG = open(WORK + r"\step9b_output.txt", "w", encoding="utf-8")
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()

variants = json.load(open(WORK + r"\nppa_locus_variants.json"))
OUTCOMES = {
    "cad_vdh": "ebi-a-GCST005195", "cad_cardio": "ieu-a-7",
    "af": "ebi-a-GCST006061", "hf": "ebi-a-GCST009541",
    "stroke": "ebi-a-GCST005838", "ces": "ebi-a-GCST006910",
}
NEW_EXP = {
    "pqtl": "ebi-a-GCST90012082",   # NT-proBNP SCALLOP n=21,758
    "pqtlr": "prot-a-2078",         # NT-proBNP SomaLogic n=3,301
    "anp":  "prot-a-2076",          # ANP SomaLogic n=3,301
}

for tag, exp_ds in NEW_EXP.items():
    exp_f = WORK + f"\\nppa__{exp_ds}_assoc.json"
    if not os.path.exists(exp_f):
        P(f"[skip] {exp_ds} 无关联文件"); continue
    exp_h = harmonize(variants, json.load(open(exp_f)))
    for outkey, out_ds in OUTCOMES.items():
        out_f = WORK + f"\\nppa__{out_ds}_assoc.json"
        out_h = harmonize(variants, json.load(open(out_f)))
        shared = sorted(set(exp_h) & set(out_h), key=lambda r: exp_h[r][0])
        pair = f"nppa_{tag}_{outkey}"
        def rows(h, keys):
            return [",".join(str(x) for x in (r, h[r][0], h[r][1], h[r][2], h[r][3], h[r][4] or "", h[r][5], h[r][6])) for r in keys]
        with open(os.path.join(INDIR, f"{pair}_exposure.csv"), "w") as f:
            f.write("rsid,pos,beta,se,eaf,n,ea,nea\n" + "\n".join(rows(exp_h, shared)) + "\n")
        with open(os.path.join(INDIR, f"{pair}_outcome.csv"), "w") as f:
            f.write("rsid,pos,beta,se,eaf,n,ea,nea\n" + "\n".join(rows(out_h, shared)) + "\n")
        P(f"{pair}: exp {len(exp_h)} / out {len(out_h)} -> shared {len(shared)}")

# ---- rs5068 / rs5066 直接查表 ----
P("\n=== rs5068 / rs5066 直接查表 (协调到 alt 方向) ===")
KEY = ["rs5068", "rs5066", "rs198411", "rs1421811"]
LOOKUP_DS = {
    "ebi-a-GCST90012082": "NT-proBNP pQTL (SCALLOP, EUR)",
    "prot-a-2078": "NT-proBNP pQTL (SomaLogic, EUR)",
    "prot-a-2076": "ANP pQTL (SomaLogic, EUR)",
    "bbj-a-159": "CAD (BBJ, EA)",
    "bbj-a-109": "CHF (BBJ, EA)",
    "bbj-a-129": "Ischemic stroke (BBJ, EA)",
    "bbj-a-71": "AF (BBJ, EA)",
    "eqtl-a-ENSG00000175206": "NPPA eQTL (eQTLGen, EUR)",
    "ebi-a-GCST005195": "CAD (van der Harst, EUR)",
    "ebi-a-GCST009541": "HF (HERMES, EUR)",
    "ebi-a-GCST005838": "Stroke (MEGASTROKE, EUR)",
}
for ds, label in LOOKUP_DS.items():
    f = WORK + f"\\nppa__{ds}_assoc.json"
    if not os.path.exists(f):
        continue
    h = harmonize(variants, json.load(open(f)))
    for k in KEY:
        if k in h:
            pos, beta, se, eaf, n, ea, nea = h[k]
            z = beta / se if se else float("nan")
            P(f"  {label:35s} {k:9s} beta={beta:+.4f} se={se:.4f} z={z:+.2f} eaf={eaf if eaf is None else round(eaf,4)} ea={ea}")
        else:
            P(f"  {label:35s} {k:9s} -- 不在协调后集合")
P("DONE")
LOG.close()
