# -*- coding: utf-8 -*-
"""make_companion_figures.py — companion 稿全部图表（真实数据）"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import pandas as pd, numpy as np, os

OUT = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\Companion_Submission"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 8.5, "axes.linewidth": 0.8})

# ============ Fig 1 研究设计 ============
fig, ax = plt.subplots(figsize=(8.5, 5.2))
ax.set_xlim(0, 100); ax.set_ylim(0, 62); ax.axis("off")
def box(x, y, w, h, title, body, fc, fs=8.0, bold=False):
    ax.add_patch(mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4",
                 fc=fc, ec="#555555", lw=0.8))
    ax.text(x + w/2, y + h - 1.7, title, ha="center", va="top", fontsize=fs,
            fontweight="bold" if bold else "normal")
    ax.text(x + w/2, y + h - 4.2, body, ha="center", va="top", fontsize=fs - 0.7, linespacing=1.35)
box(1, 44, 30, 17, "STEP 1 — cis-MR target screen",
    "8 a priori targets of the natriuretic\npathway (NPPA/NPPB/NPR1-3/MME/CORIN/IL6R)\n\nProtein & expression GWAS instruments\n→ 10 cardiometabolic outcomes\nDual-control validation:\nIL6R→CHD 0.950 (positive)\nCRP null ×9 (negative)", "#d6eaf8", bold=True)
box(35, 44, 30, 17, "STEP 2 — Independent replication",
    "19 outcome GWAS across ancestries\n(MEGASTROKE · Sakaue FinnGen+UKB · BBJ\n· Mbatchou · van der Harst · FinnGen)\n\nFixed-effect meta-analysis:\nNPPA→stroke 0.920, p=0.0014\nNPR3→MI 0.984, p=0.00083\nrs5068→CAD/MI 0.752, p=0.00018", "#d5f5e3", bold=True)
box(69, 44, 30, 17, "STEP 3 — Mechanistic dissection",
    "Two-step mediation MR:\nNPR3 rs1421811 → SBP → MI/CAD\n\nBP-mediated fraction: 5–14%\n→ majority BP-independent\n(cGMP vascular signalling)\nLocus heterogeneity:\nrs5068 (atrial) vs NPR3 (clearance)", "#fdebd0", bold=True)
ax.annotate("", xy=(34.6, 52), xytext=(31.5, 52), arrowprops=dict(arrowstyle="-|>", lw=1.4, color="#333"))
ax.annotate("", xy=(68.6, 52), xytext=(65.5, 52), arrowprops=dict(arrowstyle="-|>", lw=1.4, color="#333"))
box(1, 30, 98, 11, "PRINCIPLE — dual-control anchoring before target prioritization",
    "The analytic pipeline was validated by an established drug target (IL6R blockade → CHD OR 0.950, 0.932–0.968; replicates Georgakis et al. JAMA 2021)\nand a non-druggable biomarker (CRP → 9 outcomes, all null) before candidate targets were ranked.", "#f2f2f2", fs=8.2)
box(1, 2, 47.5, 25, "FINDINGS — target prioritization",
    "1. NPPA/natriuretic pathway augmentation\n   Expression (rs198364) & protein (rs5068) layers,\n   LD-independent, directionally concordant\n   → strongest three-layer evidence\n\n2. NPR3 (clearance receptor) inhibition\n   rs1421811: MI 0.962 · CHD 0.973 · AF 0.983\n   Bonferroni-passing meta (p=0.00083)\n\n3. IL6R blockade (repurposing benchmark)\n   OR 0.950, CANTOS/tocilizumab-consistent", "#eaf2f8", fs=8.0)
box(51.5, 2, 47.5, 25, "TRANSLATION — within the CKM framework",
    "• Natriuretic pathway: existing pharmacology\n  (ARNI class) with genetic concordance\n\n• NPR3: genetically-supported, BP-lowering,\n  but 86–95% of protection is BP-independent\n  → direct vascular signalling, not just\n  blood pressure\n\n• Aligns biomarker-guided stratification\n  (companion manuscript) with pathway\n  intervention targets", "#fef9e7", fs=8.0)
fig.suptitle("From risk barometer to therapeutic target: a three-step Mendelian randomization workflow", fontsize=10.5, fontweight="bold", y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT + r"\Fig1_design.png", dpi=300); fig.savefig(OUT + r"\Fig1_design.pdf")
plt.close(fig)
print("Fig1 done")

# ============ Fig 2 cis-MR 全靶点森林图 ============
t = pd.read_csv(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\target_cis_mr_results.csv", on_bad_lines="skip")
t = t[t["OR"].notna() & t["n_iv"].notna()].copy()
rep = pd.read_csv(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\replication_results_v2.csv")
PICK = {
    ("NPPA_eQTL", "Any stroke (MEGASTROKE)"), ("NPPA_eQTL", "Ischemic stroke (MEGASTROKE)"),
    ("NPPA_eQTL", "MI (FinnGen)"), ("NPPA_eQTL", "Heart failure (HERMES)"),
    ("MME_eQTL", "CHD (FinnGen)"), ("MME_eQTL", "MI (FinnGen)"),
    ("CORIN_eQTL", "Any stroke (MEGASTROKE)"), ("CORIN_eQTL", "MI (FinnGen)"),
    ("IL6R_INTERVAL", "CHD (CARDIoGRAMplusC4D)"), ("IL6R_SCALLOP", "CHD (CARDIoGRAMplusC4D)"),
    ("CRP_Barton", "Any stroke (MEGASTROKE)"), ("CRP_Ligthart", "CHD (CARDIoGRAMplusC4D)"),
}
LBL = {"NPPA_eQTL": "NPPA expression (cis-eQTL)", "NPR3_funcvar": "NPR3 (functional rs1421811)",
       "MME_eQTL": "MME/neprilysin expression", "CORIN_eQTL": "CORIN expression",
       "IL6R_INTERVAL": "IL6R protein (pos. ctrl, INTERVAL)", "IL6R_SCALLOP": "IL6R protein (pos. ctrl, SCALLOP)",
       "CRP_Barton": "CRP protein (neg. ctrl, Barton)", "CRP_Ligthart": "CRP protein (neg. ctrl, Ligthart)"}
rows_f2 = []
for _, r in t.iterrows():
    if (r["target"], r["outcome"]) in PICK:
        rows_f2.append({"lbl": LBL[r["target"]], "oc": r["outcome"], "OR": r["OR"], "lo": r["lo"], "hi": r["hi"], "p": r["p"],
                        "grp": "candidate" if r["target"] in ("NPPA_eQTL", "MME_eQTL", "CORIN_eQTL") else ("pos" if "IL6R" in r["target"] else "neg")})
# NPR3 功能变异行（来自复制结果文件）
for oc, show in [("MI (FinnGen) [DISCOVERY]", "MI (FinnGen)"), ("CHD (FinnGen) [PRIOR]", "CHD (FinnGen)"), ("CAD (vanderHarst UKB+Leipzig) [REP]", "CAD (vanderHarst)")]:
    rr = rep[(rep["gene"] == "NPR3-functional") & (rep["outcome"] == oc)]
    if len(rr):
        r0 = rr.iloc[0]
        rows_f2.append({"lbl": LBL["NPR3_funcvar"], "oc": show, "OR": r0["OR"], "lo": r0["lo"], "hi": r0["hi"], "p": r0["p"], "grp": "candidate"})
sel = pd.DataFrame(rows_f2)
sel["_o"] = sel["grp"].map({"candidate": 0, "pos": 1, "neg": 2})
sel = sel.sort_values(["_o", "lbl", "oc"]).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(7.8, 0.38 * len(sel) + 1.7))
ys = np.arange(len(sel))[::-1]
prev_lbl = None
for i, (_, r) in enumerate(sel.iterrows()):
    y = ys[i]
    col = "#1f6f43" if r["p"] < 0.05 else "#7f8c8d"
    ax.plot([np.log10(r["lo"]), np.log10(r["hi"])], [y, y], color=col, lw=1.2)
    ax.plot(np.log10(r["OR"]), y, "s", color=col, ms=4.5)
    ax.text(0.06, y, f'{r["lbl"]} → {r["oc"]}', va="center", fontsize=7.4, transform=ax.get_yaxis_transform())
    ptxt = f'{r["p"]:.3g}' if r["p"] >= 0.001 else "<0.001"
    ax.text(1.01, y, f'{r["OR"]:.2f} ({r["lo"]:.2f}–{r["hi"]:.2f})   p={ptxt}', va="center", fontsize=7.0, transform=ax.get_yaxis_transform())
ax.axvline(0, color="black", lw=0.8, ls="--")
ax.set_xlim(np.log10(0.55), np.log10(1.7))
ax.set_xticks([np.log10(x) for x in [0.6, 0.8, 1.0, 1.2, 1.6]])
ax.set_xticklabels(["0.6", "0.8", "1.0", "1.2", "1.6"])
ax.set_yticks([]); ax.set_ylim(-1, len(sel))
ax.set_xlabel("Odds ratio (log scale) per genetically predicted increment of target level", fontsize=8.5)
ax.spines[["left", "top", "right"]].set_visible(False)
fig.suptitle("Fig 2 | cis-Mendelian randomization across the natriuretic pathway and controls", fontsize=9.5, fontweight="bold", x=0.06, ha="left", y=0.99)
fig.text(0.06, 0.905, "Green = p<0.05; grey = null. Instruments: cis-eQTL (eQTLGen whole blood), cis-pQTL (INTERVAL/SCALLOP/Barton/Ligthart), functional variants (rs1421811). NPR1/NPR2: no genome-wide-significant cis instruments available (not shown).", fontsize=6.6, color="#444")
fig.tight_layout(rect=[0.06, 0.02, 0.70, 0.90])
fig.savefig(OUT + r"\Fig2_cisMR_forest.png", dpi=300); fig.savefig(OUT + r"\Fig2_cisMR_forest.pdf")
plt.close(fig)
print("Fig2 done")

# ============ Fig 3 复制 + meta 森林图 ============
r = pd.read_csv(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\replication_results_v2.csv")
r = r[~r["scale_flag"].fillna(False)].copy()
panels = [
    ("A | NPPA expression (rs198364) → any/ischemic stroke", "NPPA-expression", None, 0.920, 0.874, 0.968, "Fixed-effect meta: OR 0.920 (0.874–0.968), p=0.0014"),
    ("B | NPR3 functional variant (rs1421811) → MI/CAD", "NPR3-functional", None, 0.984, 0.975, 0.993, "Fixed-effect meta (k=4): OR 0.984 (0.975–0.993), p=0.00083"),
    ("C | NPPA protein (rs5068) → CAD/MI", "NPPA-protein", None, 0.752, 0.648, 0.873, "Fixed-effect meta: OR 0.752 (0.648–0.873), p=0.00018"),
]
fig, axes = plt.subplots(1, 3, figsize=(10.5, 4.4))
for ax, (title, gene, _, mor, mlo, mhi, mtxt) in zip(axes, panels):
    sub = r[r["gene"] == gene]
    if gene == "NPPA-expression":
        sub = sub[sub["outcome"].str.contains("stroke|Stroke", case=False, na=False)]
    elif gene == "NPR3-functional":
        sub = sub[sub["outcome"].str.contains("MI |myocardial|CHD|coronary|CAD", case=False, na=False)]
    else:
        sub = sub[sub["outcome"].str.contains("CAD|MI|coronary|CHD", case=False, na=False)]
    sub = sub.head(6)
    ys = np.arange(len(sub) + 1)[::-1]
    for i, (_, row) in enumerate(sub.iterrows()):
        y = ys[i]
        col = "#1f6f43" if row["p"] < 0.05 else "#7f8c8d"
        ax.plot([np.log(row["lo"]), np.log(row["hi"])], [y, y], color=col, lw=1.2)
        ax.plot(np.log(row["OR"]), y, "s", color=col, ms=4.5)
        ax.text(0.03, y, row["outcome"].replace(" [DISCOVERY]", " *").replace(" [REP]", "").replace(" [XANC]", " †").replace(" [SCLAE]", ""), va="center", fontsize=6.4, transform=ax.get_yaxis_transform())
    ym = ys[-1]
    ax.plot([np.log(mlo), np.log(mhi)], [ym, ym], color="#8e44ad", lw=2.0)
    ax.plot(np.log(mor), ym, "D", color="#8e44ad", ms=5.5)
    ax.text(0.03, ym - 0.75, "Meta-analysis", va="center", fontsize=6.8, fontweight="bold", color="#8e44ad", transform=ax.get_yaxis_transform())
    ax.axvline(0, color="black", lw=0.7, ls="--")
    ax.set_xlim(np.log(0.55), np.log(1.55)); ax.set_ylim(-1.6, len(sub) + 0.8)
    ax.set_yticks([])
    ax.set_xlabel("OR (log scale)", fontsize=7.5)
    ax.spines[["left", "top", "right"]].set_visible(False)
    ax.set_title(title, fontsize=7.8, fontweight="bold", loc="left")
    ax.text(0.5, -0.24, mtxt, transform=ax.transAxes, ha="center", fontsize=6.6, color="#8e44ad", fontweight="bold")
fig.suptitle("Fig 3 | Independent replication and fixed-effect meta-analysis of the two prioritized signals", fontsize=9.5, fontweight="bold", x=0.01, ha="left")
fig.text(0.01, 0.015, "* discovery dataset;  † East-Asian cross-ancestry. Scale-flagged datasets (Dönertaş/Neale UKB series) excluded a priori.", fontsize=6.5, color="#444")
fig.tight_layout(rect=[0, 0.045, 1, 0.93])
fig.savefig(OUT + r"\Fig3_replication_forest.png", dpi=300); fig.savefig(OUT + r"\Fig3_replication_forest.pdf")
plt.close(fig)
print("Fig3 done")

# ============ Fig 4 中介路径图 ============
fig, ax = plt.subplots(figsize=(7.6, 4.0))
ax.set_xlim(0, 100); ax.set_ylim(0, 56); ax.axis("off")
def node(x, y, w, h, txt, fc):
    ax.add_patch(mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.5", fc=fc, ec="#333", lw=1.0))
    ax.text(x + w/2, y + h/2, txt, ha="center", va="center", fontsize=8.0, linespacing=1.4)
node(2, 20, 24, 16, "NPR3 rs1421811\nG allele\n(clearance receptor\nfunction)", "#d6eaf8")
node(38, 32, 24, 16, "Systolic blood pressure\n−0.57 mmHg per allele\n(Evangelou, N=757,601)", "#fdebd0")
node(74, 20, 24, 16, "MI / CAD\nMI: OR 0.962 (0.937–0.988)\nCHD: OR 0.973 (0.947–0.999)", "#d5f5e3")
ax.annotate("", xy=(37.5, 40), xytext=(26.5, 32), arrowprops=dict(arrowstyle="-|>", lw=1.3, color="#b9770e"))
ax.text(32.0, 40.5, "Step 1\np=4.6×10^-75", fontsize=7.2, ha="center", color="#b9770e")
ax.annotate("", xy=(74, 36), xytext=(62.5, 40), arrowprops=dict(arrowstyle="-|>", lw=1.3, color="#b9770e"))
ax.text(68.2, 44.5, "Step 2 (IVW, 27 strong IVs)\nper 10 mmHg: OR 1.03–1.07", fontsize=7.2, ha="center", color="#b9770e")
ax.annotate("", xy=(74, 28), xytext=(26.5, 28), arrowprops=dict(arrowstyle="-|>", lw=1.6, color="#1f6f43"))
ax.text(50, 30.0, "Total effect", fontsize=8.0, ha="center", fontweight="bold", color="#1f6f43")
ax.annotate("", xy=(50, 31.6), xytext=(50, 33.4), arrowprops=dict(arrowstyle="-", lw=0.7, color="#999"))
# 中介比例标注
ax.add_patch(mpatches.FancyBboxPatch((30, 2), 40, 12, boxstyle="round,pad=0.5", fc="#fef9e7", ec="#b7950b", lw=1.0))
ax.text(50, 8, "Mediated fraction (two-step MR): 4.7% (MI) · 13.3% (CAD) · 14.4% (CHD)\n→ 86–95% of NPR3 protection is blood-pressure-independent\nConsistent with direct vascular cGMP signalling; estimator sensitivity in Supplement",
        ha="center", va="center", fontsize=7.4, linespacing=1.5)
fig.suptitle("Fig 4 | Two-step Mendelian randomization: dissecting the NPR3 → myocardial infarction pathway", fontsize=9.5, fontweight="bold", x=0.01, ha="left", y=0.98)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(OUT + r"\Fig4_mediation.png", dpi=300); fig.savefig(OUT + r"\Fig4_mediation.pdf")
plt.close(fig)
print("Fig4 done")

# ============ Graphical abstract ============
fig, ax = plt.subplots(figsize=(8.5, 4.6))
ax.set_xlim(0, 100); ax.set_ylim(0, 54); ax.axis("off")
box2 = lambda x, y, w, h, t, b, c: (ax.add_patch(mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4", fc=c, ec="#555", lw=0.8)),
    ax.text(x + w/2, y + h - 1.5, t, ha="center", va="top", fontsize=8.5, fontweight="bold"),
    ax.text(x + w/2, y + h - 4.5, b, ha="center", va="top", fontsize=7.6, linespacing=1.4))
box2(1, 30, 31, 22, "VALIDATED PIPELINE", "cis-MR across natriuretic\npathway with dual controls\nIL6R → CHD 0.950 (replicated)\nCRP → 9 outcomes null", "#d6eaf8")
box2(34.5, 30, 31, 22, "REPLICATED TARGETS", "NPPA pathway → stroke\nmeta 0.920, p=0.0014\nNPR3 → MI meta 0.984\np=0.00083 (Bonferroni-pass)", "#d5f5e3")
box2(68, 30, 31, 22, "MECHANISM", "NPR3 protection only 5–14%\nBP-mediated → direct\nvascular cGMP signalling\nsupports natriuretic strategy", "#fdebd0")
ax.annotate("", xy=(34.1, 41), xytext=(32.5, 41), arrowprops=dict(arrowstyle="-|>", lw=1.5, color="#333"))
ax.annotate("", xy=(67.6, 41), xytext=(66, 41), arrowprops=dict(arrowstyle="-|>", lw=1.5, color="#333"))
box2(1, 2, 98, 24, "FROM BAROMETER TO TARGET — genetic validation of natriuretic-pathway intervention in cardiometabolic disease",
    "Step 1: 8 a priori targets screened by cis-MR (protein + expression instruments)\nStep 2: independent replication across ancestries (MEGASTROKE · FinnGen+UKB · BBJ · GeneBANK · van der Harst)\nStep 3: two-step mediation MR quantifies the blood-pressure-independent component of protection\nRanking: NPPA pathway (3-layer evidence) > NPR3 inhibition > IL6R blockade (benchmark)", "#f2f2f2")
fig.savefig(OUT + r"\Graphical_abstract.png", dpi=300)
plt.close(fig)
print("Graphical abstract done")
print("ALL FIGURES DONE ->", OUT)
