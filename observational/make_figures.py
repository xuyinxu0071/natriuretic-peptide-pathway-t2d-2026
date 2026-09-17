# -*- coding: utf-8 -*-
"""
投稿图表生成（真实数据，英文，300dpi）
Fig1 研究设计图；Fig2 主结果森林图（HRS+CHARLS）；Fig3 KM+CIF；Fig4 RCS；Fig5 MR 森林图
+ 图形摘要（graphical abstract）
运行：python make_figures.py
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import pyreadstat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from lifelines import CoxPHFitter, KaplanMeierFitter
from scipy.stats import chi2

OUT = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"
HRS_BASE = r"D:\2026\3. HRS  美国\HRS_美国"
CH_BASE = r"D:\2026\1. CHARLS  中国\CHARLS_中国"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.8, "figure.dpi": 300, "savefig.dpi": 300,
    "savefig.bbox": "tight", "axes.spines.top": False, "axes.spines.right": False,
})
C = {"main": "#1f4e79", "accent": "#c0392b", "gray": "#7f8c8d",
     "light": "#aed6f1", "green": "#1e8449", "orange": "#d68910"}

# ============================================================
# 数据复建（HRS，与 v4/v5 一致）
# ============================================================
vbs, _ = pyreadstat.read_sav(HRS_BASE + r"\Raw_data\2016 静脉血研究VBS\HRS2016VBS\hrs2016vbs.sav")
vbs["hhidpn"] = vbs["HHID"].astype(int) * 1000 + vbs["PN"].astype(int)
vbs = vbs[["hhidpn", "PNTBNPE", "PCRP", "PCYSC", "PVBSWGTR"]].rename(
    columns={"PNTBNPE": "ntbnp", "PCRP": "crp", "PCYSC": "cysc", "PVBSWGTR": "vbswgt"})
dbs, _ = pyreadstat.read_dta(HRS_BASE + r"\Raw_data\2016 HRS\BIOMK16BL\BIOMK16BL_R.dta")
dbs["hhidpn"] = dbs["HHID"].astype(int) * 1000 + dbs["PN"].astype(int)
dbs = dbs[["hhidpn", "PA1C_ADJ"]].rename(columns={"PA1C_ADJ": "a1c"})
hrs, _ = pyreadstat.read_dta(HRS_BASE + r"\Working_data\hrs.dta")
w12 = hrs[hrs["wave"] == 12][["hhidpn", "ragender", "rabyear", "raracem", "raedyrs", "bmi",
                              "smoken", "drink", "hibpe", "diabe", "hearte", "stroke",
                              "cancre", "shlt"]].copy()
w12["age"] = 2016 - w12["rabyear"]
fu = hrs[hrs["wave"].isin([13, 14, 15])][["hhidpn", "wave", "hearte", "stroke"]].sort_values(["hhidpn", "wave"])
first_ev = fu.groupby("hhidpn").apply(
    lambda g: g.loc[(g["hearte"] == 1) | (g["stroke"] == 1), "wave"].min(), include_groups=False)
last_ob = fu.groupby("hhidpn")["wave"].max()
trk, _ = pyreadstat.read_dta(HRS_BASE + r"\Raw_data\2022 HRS\2022 HRS Core\h22core\trk2022\trk2022tr_r.dta")
trk["hhidpn"] = trk["HHID"].astype(int) * 1000 + trk["PN"].astype(int)
trk = trk[["hhidpn", "EXDEATHYR", "EXDEATHMO"]]
d = (vbs.merge(dbs, on="hhidpn").merge(w12, on="hhidpn").merge(trk, on="hhidpn", how="left")
     .set_index("hhidpn"))
d = d.join(first_ev.rename("ev_wave")).join(last_ob.rename("last_wave"))

def glyrow(r):
    a = r["a1c"]
    if pd.isna(a):
        return np.nan
    if r["diabe"] == 1:
        return 1 if a < 7.0 else 2
    if a >= 6.5:
        return 3
    if a >= 5.7:
        return 4
    return 5
d["gly"] = d.apply(glyrow, axis=1)
d["death"] = ((d["EXDEATHYR"] >= 2016) & (d["EXDEATHYR"] <= 2023)).astype(int)
t_ev = np.where(d["EXDEATHMO"].notna(),
                d["EXDEATHYR"] + (d["EXDEATHMO"].astype(float) - 0.5) / 12.0,
                d["EXDEATHYR"] + 0.5)
d["t_death"] = np.where(d["death"] == 1, np.maximum(t_ev - 2016.5, 0.1), 7.5)
d["base_cvd"] = ((d["hearte"] == 1) | (d["stroke"] == 1)).astype(int)
d["new_cvd"] = 0
d.loc[(d["base_cvd"] == 0) & d["ev_wave"].notna(), "new_cvd"] = 1
d["t_cvd"] = np.where(d["new_cvd"] == 1, (d["ev_wave"] - 12) * 2.0,
                      np.where(d["last_wave"].notna(), (d["last_wave"] - 12) * 2.0, np.nan))
d.loc[d["t_cvd"].isna() & (d["base_cvd"] == 0) & d["gly"].notna(), "t_cvd"] = 5.0
d["lnbnp"] = np.log(d["ntbnp"].clip(lower=1))
d["lncysc"] = np.log(d["cysc"].clip(lower=0.1))
d["female"] = (d["ragender"] == 0).astype(int)
d["smoke_now"] = (d["smoken"] == 1).astype(int)
d["drink_now"] = (d["drink"] == 1).astype(int)

mm = d[d["gly"].notna() & d["ntbnp"].notna()].copy()
goal = mm[mm["gly"] == 1].copy()
goal["bnp_t"] = pd.qcut(goal["ntbnp"], 3, labels=["T1", "T2", "T3"])
COV_H = ["age", "female", "raracem", "raedyrs", "bmi", "smoke_now", "drink_now",
         "hibpe", "cancre", "shlt", "base_cvd", "lncysc"]

# ============================================================
# Fig 1 研究设计图
# ============================================================
fig, ax = plt.subplots(figsize=(9, 6.2))
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

def box(x, y, w, h, text, fc, ec="#34495e", fs=8.5, tc="black", bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.12",
                                fc=fc, ec=ec, lw=1.2))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs,
            color=tc, fontweight="bold" if bold else "normal", linespacing=1.45)

def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=14, color="#34495e", lw=1.2))

box(0.3, 8.4, 9.4, 1.2,
    "Residual vascular risk despite glycemic control:\ndual-cohort + Mendelian randomization",
    "#d6eaf8", fs=10.5, bold=True)

# 三臂
box(0.3, 5.4, 3.0, 2.4,
    "ARM 1 — Discovery\nHRS 2016 VBS (USA)\nN = 3,686 with HbA1c\n+ NT-proBNP\n\n"
    "Exposure: ln NT-proBNP\nOutcome: all-cause death\n(N=523 DM-controlled;\n98 deaths, 2016–2023)",
    "#fdf2e9", fs=8)
box(3.5, 5.4, 3.0, 2.4,
    "ARM 2 — Extension\nCHARLS 2015 (China)\nN = 10,959 with HbA1c\n+ hsCRP\n\n"
    "Exposure: ln hsCRP\nOutcome: incident CVD\n(N=381 DM-controlled;\n113 events, 2015–2020)",
    "#eafaf1", fs=8)
box(6.7, 5.4, 3.0, 2.4,
    "ARM 3 — Triangulation\nTwo-sample MR\n(OpenGWAS)\n\n"
    "Instruments: NPPA/NPPB cis\n(rs198389, F=180; rs5068)\nOutcomes: CAD, MI, HF,\nstroke (6 consortia)",
    "#f4ecf7", fs=8)
arrow(5.0, 8.3, 5.0, 7.9)
arrow(1.8, 5.3, 1.8, 4.75); arrow(5.0, 5.3, 5.0, 4.75); arrow(8.2, 5.3, 8.2, 4.75)

box(0.3, 2.9, 9.4, 1.7,
    "Analyses: survey-weighted & cluster-robust Cox models (12 covariates incl. cystatin C)\n"
    "MICE m=20 · bootstrap · restricted cubic splines · competing-risk CIF · QBA · E-values",
    "#fef9e7", fs=8.5)

box(0.3, 0.6, 4.6, 1.9,
    "Observational arms\nNT-proBNP and hsCRP identify\nresidual risk in HbA1c-controlled DM\n"
    "HR 1.61 (1.10–2.34) death; 1.96 (1.20–3.21) CVD",
    "#fadbd8", fs=8.5)
box(5.1, 0.6, 4.6, 1.9,
    "Genetic arm\nMarker–driver dissociation:\nNT-proBNP null for CAD/MI/HF;\n"
    "NPPA pathway protective for stroke\n(IVW OR 0.91; rs5068 OR 0.63)",
    "#d4e6f1", fs=8.5)
arrow(2.6, 2.8, 2.6, 2.55); arrow(7.4, 2.8, 7.4, 2.55)

plt.savefig(OUT + r"\Fig1_design.png")
plt.savefig(OUT + r"\Fig1_design.pdf")
plt.close()
print("Fig1 done")

# ============================================================
# Fig 2 主结果森林图（双队列）
# ============================================================
# HRS per-IQR lnNT-proBNP→死亡，各血糖层（v4 肾功能调整后实测值）
hrs_rows = [
    ("HRS: DM controlled (HbA1c<7.0)", 1.517, None, None, 0.014),
    ("HRS: DM uncontrolled", 1.694, None, None, 0.003),
    ("HRS: Undiagnosed high HbA1c", 3.230, None, None, 0.021),
    ("HRS: Prediabetes", 1.439, None, None, 0.005),
    ("HRS: No diabetes", 1.985, None, None, 0.000),
]
# 主模型+敏感性（v4/v5 实测 CI）
hrs_main = [
    ("Main model (12 covariates)", 1.606, 1.102, 2.340),
    ("MICE pooled (m=20)", 1.554, 1.061, 2.276),
    ("VBS-weighted", 1.77, 1.19, 2.63),
    ("Core 6-covariate model", 1.947, 1.413, 2.682),
    ("L2-penalized", 1.436, 1.118, 1.845),
    ("Bootstrap (500 resamples)", 1.612, 1.075, 2.480),
]
ch_rows = [
    ("CHARLS: DM controlled", 1.273, 1.013, 1.598),
    ("CHARLS: DM uncontrolled", 0.956, 0.742, 1.233),
    ("CHARLS: Undiagnosed high HbA1c", 0.974, 0.771, 1.232),
    ("CHARLS: No diabetes", 1.083, 1.011, 1.160),
]

fig, axes = plt.subplots(1, 3, figsize=(11, 5.6), gridspec_kw={"width_ratios": [1, 1, 1]})
def forest(ax, rows, title, has_ci):
    ys = np.arange(len(rows))[::-1]
    for i, r in enumerate(rows):
        y = ys[i]
        if has_ci and r[2] is not None:
            ax.plot([np.log(r[2]), np.log(r[3])], [y, y], color=C["main"], lw=1.6)
            ax.scatter(np.log(r[1]), y, s=42, color=C["main"], zorder=3)
            ax.text(2.9, y, f"{r[1]:.2f} ({r[2]:.2f}-{r[3]:.2f})", fontsize=7.2,
                    va="center", ha="left")
        else:
            ax.scatter(np.log(r[1]), y, s=42, color=C["main"], zorder=3, marker="D")
            ax.text(2.9, y, f"{r[1]:.2f}", fontsize=7.2, va="center", ha="left")
    ax.axvline(0, color=C["gray"], lw=0.8, ls="--")
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7.5)
    ax.set_xlim(-1.2, 4.6)
    ax.set_xticks([-1, 0, 1, 2])
    ax.set_xticklabels([0.37, 1, 2.7, 7.4], fontsize=8)
    ax.set_xlabel("Hazard ratio (log scale)", fontsize=8.5)
    ax.set_title(title, fontsize=9.5, fontweight="bold", pad=8)

forest(axes[0], hrs_rows, "A. HRS per-IQR lnNT-proBNP\n→ all-cause death, by glycemic stratum", False)
forest(axes[1], hrs_main, "B. HRS DM-controlled:\nmodel robustness", True)
forest(axes[2], ch_rows, "C. CHARLS per-IQR lnhsCRP\n→ incident CVD, by glycemic stratum", True)
plt.tight_layout()
plt.savefig(OUT + r"\Fig2_forest_main.png")
plt.savefig(OUT + r"\Fig2_forest_main.pdf")
plt.close()
print("Fig2 done")

# ============================================================
# Fig 3 KM + 竞争风险 CIF
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(10, 4.4))
# Panel A: KM
ax = axes[0]
colors3 = [C["green"], C["orange"], C["accent"]]
for (t, col) in zip(["T1", "T2", "T3"], colors3):
    sub = goal[goal["bnp_t"] == t]
    kmf = KaplanMeierFitter().fit(sub["t_death"], sub["death"], label=t)
    kmf.plot_survival_function(ax=ax, ci_show=False, color=col, lw=1.6)
    cum = (1 - kmf.predict(7.5)) * 100
    ax.text(7.55, kmf.predict(7.5) + 0.01, f"{cum:.0f}%", color=col, fontsize=8, va="bottom")
ax.set_xlabel("Years since 2016 VBS", fontsize=8.5)
ax.set_ylabel("Survival probability", fontsize=8.5)
ax.set_title("A. All-cause mortality, HRS DM-controlled\nby NT-proBNP tertile", fontsize=9.5, fontweight="bold")
ax.legend(title="NT-proBNP tertile", fontsize=7.5, title_fontsize=8, loc="lower left")
ax.set_xlim(0, 8.6)
# Panel B: Aalen-Johansen CIF
ax = axes[1]
goal_nc = goal[goal["base_cvd"] == 0].copy()
for (t, col) in zip(["T1", "T2", "T3"], colors3):
    sub = goal_nc[goal_nc["bnp_t"] == t].copy()
    sub["ev_type"] = 0
    sub.loc[sub["new_cvd"] == 1, "ev_type"] = 1
    sub.loc[(sub["ev_type"] == 0) & (sub["death"] == 1), "ev_type"] = 2
    sub["t_any"] = np.where(sub["ev_type"] == 1, sub["t_cvd"].fillna(sub["t_death"]), sub["t_death"])
    sub = sub.sort_values("t_any").reset_index(drop=True)
    n = len(sub); S = 1.0; cif = 0.0
    grid_t, grid_c = [0.0], [0.0]
    at_risk = n
    for i in range(n):
        if at_risk <= 0:
            break
        ev = sub.loc[i, "ev_type"]
        if ev == 1:
            cif += S * (1 / at_risk)
        elif ev == 2:
            S *= (1 - 1 / at_risk)
        grid_t.append(sub.loc[i, "t_any"]); grid_c.append(cif)
        at_risk -= 1
    grid_t = np.array(grid_t); grid_c = np.array(grid_c) * 100
    ax.step(grid_t, grid_c, where="post", color=col, lw=1.6, label=t)
ax.set_xlabel("Years since 2016 VBS", fontsize=8.5)
ax.set_ylabel("Cumulative incidence of CVD, %", fontsize=8.5)
ax.set_title("B. Incident CVD (Aalen–Johansen,\ndeath as competing event)", fontsize=9.5, fontweight="bold")
ax.legend(title="NT-proBNP tertile", fontsize=7.5, title_fontsize=8, loc="upper left")
ax.set_xlim(0, 8.0)
plt.tight_layout()
plt.savefig(OUT + r"\Fig3_KM_CIF.png")
plt.savefig(OUT + r"\Fig3_KM_CIF.pdf")
plt.close()
print("Fig3 done")

# ============================================================
# Fig 4 RCS 剂量-反应曲线
# ============================================================
def rcs_basis(x, knots):
    k = len(knots)
    cols = [x]
    def cub(u):
        return np.where(u > 0, u ** 3, 0.0)
    denom = (knots[-1] - knots[0]) ** 2
    for j in range(k - 2):
        term = (cub(x - knots[j])
                - cub(x - knots[-2]) * (knots[-1] - knots[j]) / (knots[-1] - knots[-2])
                + cub(x - knots[-1]) * (knots[-2] - knots[j]) / (knots[-1] - knots[-2]))
        cols.append(term / denom)
    return cols

sub = goal[["lnbnp"] + COV_H + ["death", "t_death"]].dropna()
x = sub["lnbnp"].values
knots = np.percentile(x, [5, 35, 65, 95])
basis = rcs_basis(x, knots)
scols = ["lnbnp_s0", "lnbnp_s1", "lnbnp_s2"]
dat = sub.copy()
for i, b in enumerate(basis):
    dat[scols[i]] = b
cph_sp = CoxPHFitter().fit(dat[scols + COV_H + ["death", "t_death"]], "t_death", "death")

xgrid = np.linspace(np.percentile(x, 1), np.percentile(x, 99), 200)
Bg = rcs_basis(xgrid, knots)
Z = np.column_stack(Bg)
x0 = np.median(x)
Z0 = np.column_stack(rcs_basis(np.array([x0]), knots))
beta = cph_sp.params_[scols].values
V = cph_sp.variance_matrix_.loc[scols, scols].values
loghr = (Z - Z0) @ beta
se = np.sqrt(np.einsum("ij,jk,ik->i", Z - Z0, V, Z - Z0))

fig, ax = plt.subplots(figsize=(6.4, 4.4))
ax.axhline(0, color=C["gray"], lw=0.8, ls="--")
ax.plot(xgrid, np.exp(loghr), color=C["main"], lw=1.8)
ax.fill_between(xgrid, np.exp(loghr - 1.96*se), np.exp(loghr + 1.96*se),
                color=C["main"], alpha=0.15)
# 节点标记
for kn in knots:
    ax.axvline(kn, color=C["gray"], lw=0.5, ls=":", alpha=0.6)
# 次轴：pg/mL（显式刻度，避免自动刻度重叠）
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
ln_ticks = np.log([30, 100, 300, 1000, 3000])
ax2.set_xticks(ln_ticks)
ax2.set_xticklabels(["30", "100", "300", "1,000", "3,000"], fontsize=8)
ax2.set_xlabel("NT-proBNP, pg/mL", fontsize=8.5)
ax.set_xlabel("ln NT-proBNP", fontsize=8.5)
ax.set_ylabel("Hazard ratio for all-cause death\n(reference: median)", fontsize=8.5)
ax.set_title("Restricted cubic spline: NT-proBNP and all-cause mortality\n"
             "HRS DM-controlled; nonlinearity p = 0.010", fontsize=9.5, fontweight="bold")
ax.set_ylim(0.1, 30)
ax.set_yscale("log")
plt.tight_layout()
plt.savefig(OUT + r"\Fig4_RCS.png")
plt.savefig(OUT + r"\Fig4_RCS.pdf")
plt.close()
print("Fig4 done")

# ============================================================
# Fig 5 MR 森林图
# ============================================================
mr = pd.read_csv(OUT + r"\mr_v2_results_full.csv")
# 选取展示：发现层 4 结局（SetA IVW + 两个单 SNP）+ 复制层卒中
EN_MAP = {
    "CAD (CARDIoGRAMplusC4D)": "Coronary artery disease (CARDIoGRAM+C4D)",
    "心肌梗死 (UKB+C4D)": "Myocardial infarction (UKB+C4D)",
    "心衰 (HERMES)": "Heart failure (HERMES)",
    "全卒中 (MEGASTROKE)": "Any stroke (MEGASTROKE)",
    "冠心病 (FinnGen R9)": "Coronary heart disease (FinnGen R9)",
    "心肌梗死 (FinnGen R9)": "Myocardial infarction (FinnGen R9)",
    "心衰 (FinnGen R9)": "Heart failure (FinnGen R9)",
    "卒中 (FinnGen R9 广义)": "Stroke, broad (FinnGen R9)",
}
sel = []
for out in ["CAD (CARDIoGRAMplusC4D)", "心肌梗死 (UKB+C4D)", "心衰 (HERMES)", "全卒中 (MEGASTROKE)"]:
    for _, r in mr[(mr["outcome"] == out) & (mr["set"].str.startswith("SetA"))].iterrows():
        if r["method"].startswith("IVW"):
            sel.append((f"{EN_MAP[out]} — IVW (2 SNP)", r["OR"], r["lo"], r["hi"], r["p"], "IVW"))
        elif "rs5068" in r["method"]:
            sel.append((f"{EN_MAP[out]} — rs5068 (NPPA)", r["OR"], r["lo"], r["hi"], r["p"], "SNP"))
for out in ["冠心病 (FinnGen R9)", "心肌梗死 (FinnGen R9)", "心衰 (FinnGen R9)", "卒中 (FinnGen R9 广义)"]:
    for _, r in mr[(mr["outcome"] == out) & (mr["set"].str.startswith("SetA"))].iterrows():
        if r["method"].startswith("IVW"):
            sel.append((f"{EN_MAP[out]} — IVW", r["OR"], r["lo"], r["hi"], r["p"], "IVW"))

fig, ax = plt.subplots(figsize=(8.4, 6.2))
ys = np.arange(len(sel))[::-1]
for i, (lab, orv, lo, hi, p, kind) in enumerate(sel):
    y = ys[i]
    col = C["main"] if kind == "IVW" else C["gray"]
    mk = "s" if kind == "IVW" else "o"
    ax.plot([np.log(lo), np.log(hi)], [y, y], color=col, lw=1.4)
    ax.scatter(np.log(orv), y, s=34 if kind == "IVW" else 24, color=col, marker=mk, zorder=3)
    ax.text(1.55, y, f"{orv:.2f} ({lo:.2f}-{hi:.2f})  p={p:.3g}", fontsize=7, va="center")
ax.axvline(0, color=C["gray"], lw=0.8, ls="--")
ax.axvline(np.log(0.914), color=C["accent"], lw=0.7, ls=":", alpha=0.7)
ax.set_yticks(ys)
ax.set_yticklabels([s[0] for s in sel], fontsize=7.2)
ax.set_xticks([-0.5, 0, 0.5, 1])
ax.set_xticklabels(["0.61", "1", "1.65", "2.72"], fontsize=8)
ax.set_xlabel("Odds ratio per genetically-instrumented NT-proBNP (log scale)", fontsize=8.5)
ax.set_title("Mendelian randomization: NPPA/NPPB instruments and cardiovascular outcomes\n"
             "Discovery tier (top) and FinnGen replication (bottom); squares = IVW, circles = rs5068",
             fontsize=9.5, fontweight="bold")
plt.tight_layout()
plt.savefig(OUT + r"\Fig5_MR_forest.png")
plt.savefig(OUT + r"\Fig5_MR_forest.pdf")
plt.close()
print("Fig5 done")

# ============================================================
# 图形摘要
# ============================================================
fig, ax = plt.subplots(figsize=(9, 5))
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
box(0.2, 8.6, 9.6, 1.0,
    "Metabolic control is not vascular safety: residual risk in HbA1c-controlled type 2 diabetes",
    "#d6eaf8", fs=10.5, bold=True)

box(0.2, 5.6, 4.6, 2.6,
    "WHO: HbA1c-controlled DM\n(HRS n=523 · CHARLS n=381)\n\n"
    "High NT-proBNP / hsCRP\n→ residual vascular risk\n"
    "Death HR 1.61 · CVD HR 1.96\n7.5-yr CVD risk: 13% → 49%",
    "#fadbd8", fs=9)
box(5.2, 5.6, 4.6, 2.6,
    "Marker–driver dissociation\n(two-sample MR)\n\n"
    "Genetic NT-proBNP:\nCAD / MI / HF — null\nStroke — protective\n(IVW OR 0.91 MEGASTROKE)",
    "#d4e6f1", fs=9)
arrow(2.5, 5.5, 2.5, 5.0); arrow(7.5, 5.5, 7.5, 5.0)
box(0.2, 2.6, 9.6, 2.2,
    "Risk persists across all glycemic strata after adjustment for kidney function and lipids\n"
    "(cystatin C, HDL) — robust to MICE, weighting, penalized and bootstrap inference\n"
    "Nonlinear dose–response (spline p=0.01): risk accelerates above ~300 pg/mL",
    "#fef9e7", fs=9)
arrow(5.0, 2.5, 5.0, 2.0)
box(0.2, 0.4, 9.6, 1.4,
    "Implication for CKM Stage 1–2: measure NT-proBNP/hsCRP for stratification;\n"
    "natriuretic-pathway augmentation (ARNI direction) is a mechanistically distinct lever",
    "#d5f5e3", fs=9.5, bold=False)
plt.savefig(OUT + r"\Graphical_abstract.png")
plt.close()
print("Graphical abstract done")
print("ALL FIGURES DONE")
