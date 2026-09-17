# -*- coding: utf-8 -*-
"""
rs5068 再审计：独立子集 meta + 森林图 + 重叠矩阵 + 报告
"""
import math, os, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\rs5068_reaudit"

def meta(sub):
    b=[d[1] for d in sub]; v=[d[2]**2 for d in sub]; w=[1/x for x in v]; k=len(b)
    bfe=sum(wi*bi for wi,bi in zip(w,b))/sum(w); sefe=math.sqrt(1/sum(w))
    Q=sum(wi*(bi-bfe)**2 for wi,bi in zip(w,b)); df=k-1
    C=sum(w)-sum(wi**2 for wi in w)/sum(w); tau2=max(0.0,(Q-df)/C)
    wr=[1/(vi+tau2) for vi in v]; bre=sum(wi*bi for wi,bi in zip(wr,b))/sum(wr); sere=math.sqrt(1/sum(wr))
    p=lambda bb,ss: 2*(1-0.5*(1+math.erf(abs(bb/ss)/math.sqrt(2))))
    return bfe,sefe,p(bfe,sefe),bre,sere,p(bre,sere),Q,tau2

# (名称, logOR, se) —— CAD/MI
all7 = [
 ("CARDIoGRAMplusC4D CAD (ieu-a-7, N=184k)", -0.0317424, 0.1482183),
 ("UKB+C4D MI (ieu-a-798)",                    0.0865427, 0.1631116),
 ("van der Harst CAD (GCST005195, N=547k)",   -0.2630828, 0.1069003),
 ("GeneBANK/REGENIE CAD (GCST90013868, UKB)", -0.3182166, 0.1381818),
 ("Hartiala MI (GCST011364, N=17.5k)",        -0.3352625, 0.1649448),
 ("FinnGen CHD (finn-b-I9_CHD)",              -0.2472452, 0.1556474),
 ("FinnGen MI (finn-b-I9_MI)",                -0.1990358, 0.1866391),
]
# 近似独立子集：vdH（含 C4D+UKB）+ Hartiala + FinnGen（FinnGen CHD/MI 内部重叠，取 CHD）
indep = [all7[2], all7[4], all7[5]]

for name, sub in [("全 7 数据集", all7), ("近似独立子集（vdH + Hartiala + FinnGen CHD）", indep)]:
    bfe,sefe,pfe,bre,sere,pre_,Q,tau2 = meta(sub)
    print(f"{name}: 固定 OR={math.exp(bfe):.3f} ({math.exp(bfe-1.96*sefe):.3f}-{math.exp(bfe+1.96*sefe):.3f}) p={pfe:.3g} | "
          f"随机 OR={math.exp(bre):.3f} ({math.exp(bre-1.96*sere):.3f}-{math.exp(bre+1.96*sere):.3f}) p={pre_:.3g} | Q={Q:.1f} tau2={tau2:.4f}")

# ---- 森林图 ----
fig, ax = plt.subplots(figsize=(9, 5.2))
ys = list(range(len(all7), 0, -1))
for (name, b, s), y in zip(all7, ys):
    ax.errorbar(math.exp(b), y, xerr=[[math.exp(b)-math.exp(b-1.96*s)], [math.exp(b+1.96*s)-math.exp(b)]],
                fmt="o", color="#1f4e79", ecolor="#7f9fb9", capsize=3, ms=5)
bfe,sefe,pfe,bre,sere,pre_,Q,tau2 = meta(all7)
ax.errorbar(math.exp(bre), 0, xerr=[[math.exp(bre)-math.exp(bre-1.96*sere)], [math.exp(bre+1.96*sere)-math.exp(bre)]],
            fmt="D", color="#c00000", ecolor="#e0a0a0", capsize=4, ms=7)
labels = [n for n,_,_ in all7] + [f"Random-effects pooled: OR={math.exp(bre):.3f} ({math.exp(bre-1.96*sere):.3f}-{math.exp(bre+1.96*sere):.3f}), I2~4%"]
ax.set_yticks(ys + [0]); ax.set_yticklabels(labels, fontsize=8.5)
ax.axvline(1.0, color="grey", ls="--", lw=1)
ax.set_xscale("log"); ax.set_xlim(0.45, 1.7)
ax.set_xticks([0.5, 0.67, 1.0, 1.5]); ax.set_xticklabels(["0.50","0.67","1.00","1.50"])
ax.set_xlabel("OR per NT-proBNP-raising G allele of rs5068 (log scale)")
ax.set_title("rs5068 -> CAD/MI across all 7 outcome GWAS (single-SNP Wald ratios)")
plt.tight_layout()
fig.savefig(DIR + r"\rs5068_forest_CADMI.png", dpi=200)
print("森林图已保存: rs5068_forest_CADMI.png")

# ---- 重叠矩阵（基于联盟构成的定性评估）----
ov = [
 ("ieu-a-7 (C4D)",        "ieu-a-7 (C4D)",        "—"),
 ("vdH CAD (UKB+C4D)",    "ieu-a-7 (C4D)",        "高：vdH 纳入全部 C4D 样本"),
 ("vdH CAD (UKB+C4D)",    "ieu-a-798 (UKB MI)",   "高：共享 UKB 病例"),
 ("GeneBANK (REGENIE/UKB)","ieu-a-798 (UKB MI)",  "高：同为 UKB 主体"),
 ("GeneBANK (REGENIE/UKB)","vdH CAD (UKB+C4D)",   "高：共享 UKB"),
 ("GeneBANK (REGENIE/UKB)","ieu-a-7 (C4D)",       "低-无"),
 ("Hartiala MI",          "C4D 家族",              "低：独立 MI 联盟（deCODE 等为主）"),
 ("FinnGen CHD/MI",       "全部其他",              "无（芬兰孤立人群）"),
 ("HERMES HF",            "UKB 家族",              "中：HERMES 含 UKB+FinnGen 等 47 队列"),
]
with open(DIR + r"\rs5068_overlap_matrix.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f); w.writerow(["dataset_A", "dataset_B", "overlap_assessment"]); w.writerows(ov)
print("重叠矩阵已保存: rs5068_overlap_matrix.csv")
