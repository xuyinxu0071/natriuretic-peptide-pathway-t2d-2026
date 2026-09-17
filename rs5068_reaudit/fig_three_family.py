# -*- coding: utf-8 -*-
"""合并版核心图：rs5068 三家族分裂森林图（CAD/MI 保护 | HF null | AF/CE 卒中 有害）"""
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\rs5068_reaudit"

panels = [
    ("A. CAD / MI (atherosclerotic)", [
        ("CARDIoGRAMplusC4D CAD (ieu-a-7)", 0.969, 0.725, 1.295),
        ("UKB+C4D MI (ieu-a-798)", 1.090, 0.792, 1.501),
        ("van der Harst CAD (GCST005195)", 0.769, 0.623, 0.948),
        ("GeneBANK/REGENIE CAD (GCST90013868)", 0.727, 0.555, 0.954),
        ("Hartiala MI (GCST011364)", 0.715, 0.518, 0.988),
        ("FinnGen CHD (R9)", 0.781, 0.576, 1.060),
        ("FinnGen MI (R9)", 0.820, 0.568, 1.181),
    ], ("Random effects (k=7): OR 0.819 (0.733-0.914), p=3.8e-4, I2=4%",
        0.819, 0.733, 0.914)),
    ("B. Heart failure", [
        ("HERMES HF (GCST009541)", 1.202, 0.950, 1.520),
        ("FinnGen HF (R9)", 1.058, 0.753, 1.487),
    ], ("Fixed effect (k=2): OR 1.153 (0.951-1.399), p=0.15, I2=0%",
        1.153, 0.951, 1.399)),
    ("C. Atrial fibrillation & cardioembolic stroke", [
        ("Nielsen AF (GCST006061)", 1.487, 1.188, 1.860),
        ("Roselli AF (GCST006414)", 1.302, 1.066, 1.590),
        ("Cardioembolic stroke (MEGASTROKE)", 2.019, 1.124, 3.627),
        ("Stroke, broad (FinnGen R9)", 0.630, 0.468, 0.848),
        ("Any stroke (MEGASTROKE)", 0.900, 0.693, 1.169),
    ], ("AF fixed effect (k=2): OR 1.381 (1.189-1.603), p=2.2e-5, I2=0%",
        1.381, 1.189, 1.603)),
]

fig, ax = plt.subplots(figsize=(9.5, 8.2))
y = 0
yticks, ylabels = [], []
for title, rows, (pool_lab, po, plo, phi) in panels:
    ax.text(0.42, y, title, fontsize=10, fontweight="bold", va="center")
    y -= 1
    for name, orr, lo, hi in rows:
        ax.errorbar(orr, y, xerr=[[orr - lo], [hi - orr]], fmt="o", color="#1f4e79",
                    ecolor="#8aa9c4", capsize=3, ms=5)
        yticks.append(y); ylabels.append(name)
        y -= 1
    ax.errorbar(po, y, xerr=[[po - plo], [phi - po]], fmt="D", color="#c00000",
                ecolor="#e0a0a0", capsize=4, ms=7)
    yticks.append(y); ylabels.append(pool_lab)
    y -= 1.5

ax.set_yticks(yticks); ax.set_yticklabels(ylabels, fontsize=8.5)
ax.axvline(1.0, color="grey", ls="--", lw=1)
ax.set_xscale("log"); ax.set_xlim(0.4, 4.0)
ax.set_xticks([0.5, 0.67, 1.0, 1.5, 2.0, 3.0])
ax.set_xticklabels(["0.50", "0.67", "1.00", "1.50", "2.00", "3.00"])
ax.set_xlabel("OR per NT-proBNP-raising G allele of rs5068 (log scale)")
ax.set_title("Phenotype-family-specific effects of the natriuretic-peptide secretion variant rs5068\n(single-SNP Wald ratios; exposure: SCALLOP NT-proBNP, N=21,758)")
plt.tight_layout()
fig.savefig(DIR + r"\Fig_three_family_forest.png", dpi=220)
fig.savefig(DIR + r"\Fig_three_family_forest.pdf")
print("saved Fig_three_family_forest.png/.pdf")
