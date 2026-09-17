# -*- coding: utf-8 -*-
"""
rs5068 全数据集再审计：CAD/MI/HF 单 SNP Wald + 随机效应 meta + 异质性 + 重叠矩阵
================================================================================
目的：裁决 A 篇（rs5068→CAD/MI meta OR 0.752 强保护）与 B 篇（NT-proBNP→CAD/MI/HF
全 null）的矛盾——用 rs5068 单 SNP 在全部可得结局 GWAS 上逐个计算 Wald ratio，
做 DerSimonian-Laird 随机效应 meta，报告 Cochran Q / I2，并附样本重叠矩阵。

暴露：rs5068 G 等位基因 → NT-proBNP（SCALLOP, N=21,758）：beta=0.1452, se=0.0338
（与原两稿一致；F=18.5）。

结局原始关联来源（真实 API 落盘数据，非模拟）：
  - replication_assoc_v2.json（vdH / GeneBANK / Hartiala / FinnGen / AF / 卒中等）
  - mr_v2_results_full.csv（ieu-a-7 / ieu-a-798 / HERMES / FinnGen HF 的 rs5068 Wald，
    由 mr_ntprobnp_cvd_v2.py 自 mr_*_assoc.json 原始关联计算）
排除：ebi-a-GCST90038610（Sakaue MI）——beta/se 量纲异常（原始比例尺度），与原稿
  a priori 排除规则一致，仅作 p 值方向参考。
"""
import json, math

DIR = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"

BX, SX = 0.1452, 0.0338  # rs5068 -> NT-proBNP (SCALLOP)

def wald(by, sy):
    r = by / BX
    se = sy / abs(BX)
    z = r / se
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return r, se, p

# ---- 从 replication_assoc_v2.json 取 rs5068 结局原始 beta/se 并算 Wald ----
rep = json.load(open(DIR + r"\replication_assoc_v2.json", encoding="utf-8"))

def from_rep(key):
    a = rep[key]
    assert a["rsid"] == "rs5068" and a["ea"] == "G", key
    r, se, p = wald(a["beta"], a["se"])
    return r, se, p

# (数据集ID, 结局, 表型家族, logOR, SE, p, N, 纳入主meta?)
rows = []

def add(oid, name, family, logOR, se, p, n, include):
    rows.append(dict(id=oid, outcome=name, family=family, logOR=logOR, se=se,
                     OR=math.exp(logOR), lo=math.exp(logOR - 1.96 * se),
                     hi=math.exp(logOR + 1.96 * se), p=p, n=n, include=include))

# B 篇三大大数据集（来自 mr_v2_results_full.csv 已算好的 rs5068 Wald）
add("ieu-a-7",          "CAD (CARDIoGRAMplusC4D)",      "CAD/MI", -0.0317424, 0.1482183, 0.8304, 184305, True)
add("ieu-a-798",        "MI (UKB+C4D)",                 "CAD/MI",  0.0865427, 0.1631116, 0.5957, None,  True)
add("ebi-a-GCST009541", "HF (HERMES)",                  "HF",      0.1838843, 0.1198347, 0.1249, 977323, True)
add("finn-b-I9_HEARTFAIL","HF (FinnGen R9)",            "HF",      0.0564738, 0.1735537, 0.7449, None,  True)

# A 篇四个数据集 + FinnGen CAD/MI（从原始关联重算）
r, se, p = from_rep("rs5068|ebi-a-GCST005195");    add("ebi-a-GCST005195",   "CAD (van der Harst 2018)",  "CAD/MI", r, se, p, 547261, True)
r, se, p = from_rep("rs5068|ebi-a-GCST90013868");  add("ebi-a-GCST90013868", "CAD (GeneBANK, SPA)",       "CAD/MI", r, se, p, 352063, True)
r, se, p = from_rep("rs5068|ebi-a-GCST011364");    add("ebi-a-GCST011364",   "MI (Hartiala 2021)",        "CAD/MI", r, se, p, 17505,  True)
r, se, p = from_rep("rs5068|finn-b-I9_CHD");       add("finn-b-I9_CHD",      "CHD (FinnGen R9)",          "CAD/MI", r, se, p, None,   True)
r, se, p = from_rep("rs5068|finn-b-I9_MI");        add("finn-b-I9_MI",       "MI (FinnGen R9)",           "CAD/MI", r, se, p, None,   True)

# 参照：AF 与卒中（非主 meta，展示轴异质性）
r, se, p = from_rep("rs5068|ebi-a-GCST006061");    add("ebi-a-GCST006061",   "AF (Nielsen 2018)",         "AF",     r, se, p, 537409, False)
r, se, p = from_rep("rs5068|ebi-a-GCST006414");    add("ebi-a-GCST006414",   "AF (Roselli 2018)",         "AF",     r, se, p, 1030836, False)
r, se, p = from_rep("rs5068|ebi-a-GCST005838");    add("ebi-a-GCST005838",   "Stroke (MEGASTROKE)",       "Stroke", r, se, p, 446696, False)
r, se, p = from_rep("rs5068|finn-b-C_STROKE");     add("finn-b-C_STROKE",    "Stroke broad (FinnGen R9)", "Stroke", r, se, p, None,   False)
r, se, p = from_rep("rs5068|ebi-a-GCST006910");    add("ebi-a-GCST006910",   "Cardioembolic stroke (MEGASTROKE)", "Stroke", r, se, p, 413304, False)

# ---- DerSimonian-Laird 随机效应 + 固定效应 meta ----
def meta(sub):
    b = [d["logOR"] for d in sub]; v = [d["se"] ** 2 for d in sub]
    w = [1 / x for x in v]
    k = len(b)
    bfe = sum(wi * bi for wi, bi in zip(w, b)) / sum(w)
    sefe = math.sqrt(1 / sum(w))
    Q = sum(wi * (bi - bfe) ** 2 for wi, bi in zip(w, b))
    df = k - 1
    # chi2 p ( Wilson-Hilferty 近似 )
    z = ((Q / df) ** (1 / 3) - (1 - 2 / (9 * df))) / math.sqrt(2 / (9 * df)) if df > 0 else 0
    Qp = 1 - 0.5 * (1 + math.erf(z / math.sqrt(2)))
    C = sum(w) - sum(wi ** 2 for wi in w) / sum(w)
    tau2 = max(0.0, (Q - df) / C)
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    wr = [1 / (vi + tau2) for vi in v]
    bre = sum(wi * bi for wi, bi in zip(wr, b)) / sum(wr)
    sere = math.sqrt(1 / sum(wr))
    def pval(bb, ss):
        z2 = bb / ss
        return 2 * (1 - 0.5 * (1 + math.erf(abs(z2) / math.sqrt(2))))
    return dict(k=k, fe=(bfe, sefe, pval(bfe, sefe)), re=(bre, sere, pval(bre, sere)),
                Q=Q, Qp=Qp, tau2=tau2, I2=I2)

def show(title, sub):
    m = meta(sub)
    print(f"\n### {title} (k={m['k']})")
    for d in sub:
        flag = "  " if d["include"] else "  "
        print(f"  {d['id']:20s} {d['outcome']:38s} OR={d['OR']:6.3f} ({d['lo']:.3f}-{d['hi']:.3f}) p={d['p']:.4g}")
    bfe, sefe, pfe = m["fe"]; bre, sere, pre_ = m["re"]
    print(f"  固定效应:  OR={math.exp(bfe):.3f} ({math.exp(bfe-1.96*sefe):.3f}-{math.exp(bfe+1.96*sefe):.3f}) p={pfe:.4g}")
    print(f"  随机效应:  OR={math.exp(bre):.3f} ({math.exp(bre-1.96*sere):.3f}-{math.exp(bre+1.96*sere):.3f}) p={pre_:.4g}")
    print(f"  异质性:    Q={m['Q']:.1f} (df={m['k']-1}) p={m['Qp']:.3g}  I2={m['I2']:.0f}%  tau2={m['tau2']:.4f}")
    return m

print("=" * 100)
print("rs5068 (NP 升高等位基因 G) 单 SNP Wald：全数据集再审计")
print("=" * 100)

cadm_i = [d for d in rows if d["family"] == "CAD/MI"]
hf     = [d for d in rows if d["family"] == "HF"]
af     = [d for d in rows if d["family"] == "AF"]
stk    = [d for d in rows if d["family"] == "Stroke"]

m1 = show("CAD/MI 家族（7 数据集，主裁决）", cadm_i)
m2 = show("HF 家族（2 数据集）", hf)
m3 = show("AF 家族（2 数据集，参照）", af)
m4 = show("Stroke 家族（3 数据集，参照）", stk)

# 亚组：仅 A 篇四数据集 vs 仅 B 篇三大大数据集
a4 = [d for d in cadm_i if d["id"] in ("ebi-a-GCST005195", "ebi-a-GCST90013868", "ebi-a-GCST011364", "finn-b-I9_CHD")]
b3 = [d for d in rows if d["id"] in ("ieu-a-7", "ieu-a-798", "ebi-a-GCST009541")]
show("A 篇原四数据集（对照复核：应复现 OR~0.75）", a4)
show("B 篇三大大数据集（CAD/MI/HF 混合家族，仅示意）", b3)

# ---- 输出 CSV ----
import csv
out = DIR + r"\rs5068_reaudit\rs5068_wald_all_datasets.csv"
import os
os.makedirs(DIR + r"\rs5068_reaudit", exist_ok=True)
with open(out, "w", newline="", encoding="utf-8-sig") as f:
    wcsv = csv.DictWriter(f, fieldnames=["id", "outcome", "family", "logOR", "se", "OR", "lo", "hi", "p", "n", "include"])
    wcsv.writeheader()
    for d in rows:
        wcsv.writerow(d)
print(f"\n已保存: {out}")
