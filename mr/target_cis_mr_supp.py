# -*- coding: utf-8 -*-
"""补充分析：cis-MR 工具变量位置验证 + Cochran Q 异质性 + 逐 SNP Wald（leave-one-out 等价）+ NPPA 信号 LD 一致性"""
import json, math, csv

DIR = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"
inst = json.load(open(DIR + r"\target_cis_instruments.json", encoding="utf-8"))
assoc = json.load(open(DIR + r"\target_outcome_assoc.json", encoding="utf-8"))
assoc = {tuple(k.split("|")): v for k, v in assoc.items()}

print("===== 1. cis 工具变量位置验证（基因座位核对）=====")
GENE_POS = {"MME": ("3", 155024124), "CORIN": ("4", 47593999), "NPPA": ("1", 11845701),
            "IL6R": ("1", 154405141), "CRP": ("1", 159711438)}
for name, iv in inst.items():
    if "cis" not in iv: continue
    g = iv["gene"]
    for c in iv["cis"]:
        print(f"{name:15s} {g:6s} {c['rsid']:12s} chr{c['chr']}:{c['pos']:,} (基因中心 chr{GENE_POS[g][0]}:{GENE_POS[g][1]:,}, 距离 {abs(c['pos']-GENE_POS[g][1])//1000} kb) p={c['p']:.2g} F={c['F']:.1f}")

print("\n===== 2. Cochran Q + 逐 SNP Wald（关键结果）=====")
def wald(bx, sx, by, sy):
    if not bx or not by: return None, None
    b = by / bx
    se = abs(b) * math.sqrt((sy/by)**2 + (sx/bx)**2)
    return b, se

def z_p(z): return math.erfc(abs(z)/math.sqrt(2))

KEY = [("IL6R_SCALLOP", "ieu-a-7", "CHD CARDIoGRAM"), ("IL6R_SCALLOP", "ebi-a-GCST005838", "Any stroke"),
       ("IL6R_INTERVAL", "ieu-a-7", "CHD CARDIoGRAM"),
       ("MME_eQTL", "ieu-a-7", "CHD CARDIoGRAM"), ("MME_eQTL", "ebi-a-GCST005838", "Any stroke"),
       ("CRP_Ligthart", "ieu-a-7", "CHD CARDIoGRAM"), ("CRP_Barton", "ieu-a-7", "CHD CARDIoGRAM"),
       ("CORIN_eQTL", "finn-b-I9_CHD", "CHD FinnGen")]
for name, oid, label in KEY:
    iv = inst[name]
    walds = []
    for c in iv["cis"]:
        a = assoc.get((c["rsid"], oid))
        if not a: continue
        ea, nea, oea, onea = c.get("ea"), c.get("nea"), a.get("ea"), a.get("nea")
        if {ea, nea} == {oea, onea}:
            sign = 1 if ea == oea else -1
        else:
            continue
        b, se = wald(c["beta"], c["se"], sign * a["beta"], a["se"])
        if b is not None:
            walds.append((c["rsid"], b, se))
    if len(walds) >= 2:
        # Cochran Q
        bs = [w[1] for w in walds]; ses = [w[2] for w in walds]
        ws = [1/s**2 for s in ses]
        B = sum(w*b for w,b in zip(ws,bs))/sum(ws)
        Q = sum(w*(b-B)**2 for w,b in zip(ws,bs))
        Qp = 1 - 0.5*(1+math.erf(Q/2/math.sqrt(2)))  # chi2 df=k-1 近似（k-1=1 或 2 时精确到误差函数）
        print(f"{name} -> {label}: Q={Q:.2f} p={Qp:.3f}")
        for rsid, b, se in walds:
            print(f"   {rsid}: OR {math.exp(b):.3f} ({math.exp(b-1.96*se):.3f}-{math.exp(b+1.96*se):.3f}) p={z_p(b/se):.4g}")

print("\n===== 3. NPPA 表达-卒中信号 LD 一致性（coloc-lite）=====")
# 主仪器 rs198364；查它在 MEGASTROKE any stroke 的关联 + 它与主稿蛋白仪器 rs5068 的关系
# rs5068 与 rs198364 均在 NPPA 区域（chr1:11.85Mb vs rs5068 chr1:11.86Mb）
lead = assoc.get(("rs198364", "ebi-a-GCST005838"))
print("rs198364 (NPPA eQTL lead) -> Any stroke MEGASTROKE:", json.dumps(lead, ensure_ascii=False))
# NPPA eQTL 仪器详情
for c in inst["NPPA_eQTL"]["cis"]:
    print("NPPA eQTL IV:", json.dumps(c, ensure_ascii=False))

# 主稿蛋白层信号（供报告引用一致性）
print("\nrs198364 与 rs5068 是否同一信号待 LD 矩阵确认（Ensembl 若可用）")
try:
    import urllib.request
    req = urllib.request.Request(
        "https://rest.ensembl.org/ld/human/rs198364/rs5068?population_name=1000GENOMES:phase_3:EUR",
        headers={"Content-Type": "application/json"})
    r = json.load(urllib.request.urlopen(req, timeout=30))
    if isinstance(r, list) and r:
        print("LD rs198364-rs5068: r2 =", r[0].get("r2"), "Dprime =", r[0].get("d_prime"))
    else:
        print("LD 查询无结果:", r)
except Exception as e:
    print("Ensembl LD 查询失败（稍后重试）:", e)
