# -*- coding: utf-8 -*-
"""
协变量审计敏感性分析（2026-09-11）
文献惯例对照后补跑：
  (A) HRS 主模型 2b + ln(cystatin C)   —— 肾功能调整（NT-proBNP 经肾清除，审稿人必问）
  (B) HRS 主模型 2b + 血脂(HDL)        —— 代谢性协变量
  (C) CHARLS 主模型 3a + ln(cystatin C) + cancre —— 肾功能 + 癌症
  (D) CHARLS 主模型 3a + ln(TG/HDL)    —— 血脂轴
运行：python covariate_audit_sensitivity.py
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import pyreadstat
from lifelines import CoxPHFitter

OUT = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"
HRS_BASE = r"D:\2026\3. HRS  美国\HRS_美国"
CH_BASE = r"D:\2026\1. CHARLS  中国\CHARLS_中国"

log_lines = []
def log(msg):
    print(msg)
    log_lines.append(str(msg))

def run_cox(df, covs, event, time, cluster=None, label=""):
    cols = covs + [event, time] + ([cluster] if cluster else [])
    sub = df[cols].dropna().copy()
    cph = CoxPHFitter()
    cph.fit(sub, duration_col=time, event_col=event, robust=True,
            cluster_col=cluster if cluster else None)
    s = cph.summary
    log(f"\n--- {label} (n={len(sub)}, events={int(sub[event].sum())}) ---")
    for v in covs[:4]:
        if v in s.index:
            log(f"  {v:12s} HR={s.loc[v,'exp(coef)']:.3f} "
                f"({s.loc[v,'exp(coef) lower 95%']:.3f}-{s.loc[v,'exp(coef) upper 95%']:.3f}) p={s.loc[v,'p']:.4f}")
    return cph, s

# ============================================================
# HRS 臂：复建主分析数据（同 v3 逻辑）
# ============================================================
vbs, _ = pyreadstat.read_sav(HRS_BASE + r"\Raw_data\2016 静脉血研究VBS\HRS2016VBS\hrs2016vbs.sav")
vbs["hhidpn"] = vbs["HHID"].astype(int) * 1000 + vbs["PN"].astype(int)
vbs = vbs[["hhidpn", "PNTBNPE", "PCRP", "PCYSC", "PHDLD", "PVBSWGTR"]].rename(
    columns={"PNTBNPE": "ntbnp", "PCRP": "crp", "PCYSC": "cysc", "PHDLD": "hdl", "PVBSWGTR": "vbswgt"})

dbs, _ = pyreadstat.read_dta(HRS_BASE + r"\Raw_data\2016 HRS\BIOMK16BL\BIOMK16BL_R.dta")
dbs["hhidpn"] = dbs["HHID"].astype(int) * 1000 + dbs["PN"].astype(int)
dbs = dbs[["hhidpn", "PA1C_ADJ"]].rename(columns={"PA1C_ADJ": "a1c"})

hrs, _ = pyreadstat.read_dta(HRS_BASE + r"\Working_data\hrs.dta")
w12 = hrs[hrs["wave"] == 12][["hhidpn", "ragender", "rabyear", "raracem", "raedyrs", "bmi",
                              "smoken", "drink", "hibpe", "diabe", "hearte", "stroke",
                              "cancre", "shlt"]].copy()
w12["age"] = 2016 - w12["rabyear"]

trk, _ = pyreadstat.read_dta(HRS_BASE + r"\Raw_data\2022 HRS\2022 HRS Core\h22core\trk2022\trk2022tr_r.dta")
trk["hhidpn"] = trk["HHID"].astype(int) * 1000 + trk["PN"].astype(int)
trk = trk[["hhidpn", "EXDEATHYR", "EXDEATHMO"]]

d = (vbs.merge(dbs, on="hhidpn").merge(w12, on="hhidpn").merge(trk, on="hhidpn", how="left")
     .set_index("hhidpn"))

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

d["lnbnp"] = np.log(d["ntbnp"].clip(lower=1))
d["lncysc"] = np.log(d["cysc"].clip(lower=0.1))
d["lnhdl"] = np.log(d["hdl"].clip(lower=1))
d["female"] = (d["ragender"] == 0).astype(int)
d["smoke_now"] = (d["smoken"] == 1).astype(int)
d["drink_now"] = (d["drink"] == 1).astype(int)
d["base_cvd"] = ((d["hearte"] == 1) | (d["stroke"] == 1)).astype(int)

COV_H = ["age", "female", "raracem", "raedyrs", "bmi", "smoke_now", "drink_now",
         "hibpe", "cancre", "shlt", "base_cvd"]

goal = d[d["gly"] == 1].copy()
goal["bnp_iqr"] = goal["lnbnp"] / (goal["lnbnp"].quantile(0.75) - goal["lnbnp"].quantile(0.25))
log(f"HRS DM达标组 N={len(goal)}，死亡={int(goal['death'].sum())}")
log(f"  cysc 非缺失 {goal['lncysc'].notna().sum()}/{len(goal)}"
    f"（{100*goal['lncysc'].notna().mean():.1f}%）；HDL 非缺失 {goal['lnhdl'].notna().sum()}")

log("\n=============== HRS 敏感性（协变量审计补跑） ===============")
# 基准（复现 v3 主模型）
run_cox(goal, ["bnp_iqr"] + COV_H, "death", "t_death",
        label="HRS 基准：v3 主模型 2b（对照）")
# A. + 肾功能
run_cox(goal, ["bnp_iqr", "lncysc"] + COV_H, "death", "t_death",
        label="A. 主模型 + ln(cystatin C)（肾功能调整）")
# B. + 血脂
run_cox(goal, ["bnp_iqr", "lnhdl"] + COV_H, "death", "t_death",
        label="B. 主模型 + ln(HDL)（血脂调整）")
# A+B 联合
run_cox(goal, ["bnp_iqr", "lncysc", "lnhdl"] + COV_H, "death", "t_death",
        label="A+B. 主模型 + cysc + HDL")

# ============================================================
# CHARLS 臂：复建（同 v3 逻辑）+ 补 cysc/cancre/血脂
# ============================================================
ch, _ = pyreadstat.read_dta(CH_BASE + r"\Working_data\charls.dta")
w3 = ch[ch["wave"] == 3][["ID", "ragender", "age", "bmi", "raeducl", "hibpe", "diabe",
                          "hearte", "stroke", "bl_crp", "bl_hbalc", "bl_cysc",
                          "bl_tg", "bl_hdl", "cancre",
                          "communityID", "bloodweight", "smoken", "drinkl"]].copy()
fu4 = ch[ch["wave"].isin([4, 5])][["ID", "wave", "hearte", "stroke"]].sort_values(["ID", "wave"])
first_ev4 = fu4.groupby("ID").apply(
    lambda g: g.loc[(g["hearte"] == 1) | (g["stroke"] == 1), "wave"].min(), include_groups=False)
last_ob4 = fu4.groupby("ID")["wave"].max()
c = w3.set_index("ID").join(first_ev4.rename("ev_wave")).join(last_ob4.rename("last_wave"))

def glyrow_c(r):
    a = r["bl_hbalc"]
    if pd.isna(a):
        return np.nan
    if r["diabe"] == 1:
        return 1 if a < 7.0 else 2
    if a >= 6.5:
        return 3
    return 5
c["gly"] = c.apply(glyrow_c, axis=1)
c["base_cvd"] = ((c["hearte"] == 1) | (c["stroke"] == 1)).astype(int)
c["new_cvd"] = 0
c.loc[(c["base_cvd"] == 0) & c["ev_wave"].notna(), "new_cvd"] = 1
c["t_cvd"] = np.where(c["new_cvd"] == 1, (c["ev_wave"] - 3) * 2.5,
                      np.where(c["last_wave"].notna(), (c["last_wave"] - 3) * 2.5, np.nan))
c.loc[c["t_cvd"].isna() & (c["base_cvd"] == 0) & c["gly"].notna(), "t_cvd"] = 5.0

cm = c[c["gly"].notna() & c["bl_crp"].notna() & (c["base_cvd"] == 0)].copy()
cm["lncrp"] = np.log(cm["bl_crp"].clip(lower=0.1))
cm["crp_iqr"] = cm["lncrp"] / (cm["lncrp"].quantile(0.75) - cm["lncrp"].quantile(0.25))
cm["lncysc_c"] = np.log(cm["bl_cysc"].clip(lower=0.1))
cm["lnhdl_c"] = np.log(cm["bl_hdl"].clip(lower=1))
cm["lntg_c"] = np.log(cm["bl_tg"].clip(lower=1))
cm["female"] = (cm["ragender"] == 0).astype(int)
cm["smoke_now"] = (cm["smoken"] == 1).astype(int)
cm["drink_now"] = (cm["drinkl"] == 1).astype(int)

COV_C = ["age", "female", "bmi", "raeducl", "hibpe", "smoke_now", "drink_now"]
goalc = cm[cm["gly"] == 1].copy()
log(f"\nCHARLS DM达标组 N={len(goalc)}，新发CVD={int(goalc['new_cvd'].sum())}")
log(f"  cysc 非缺失 {goalc['lncysc_c'].notna().sum()}；cancre 非缺失 {goalc['cancre'].notna().sum()}；"
    f"HDL 非缺失 {goalc['lnhdl_c'].notna().sum()}")

log("\n=============== CHARLS 敏感性（协变量审计补跑） ===============")
run_cox(goalc, ["crp_iqr"] + COV_C, "new_cvd", "t_cvd",
        cluster="communityID", label="CHARLS 基准：v3 主模型 3a（对照）")
# C. + 肾功能 + 癌症
run_cox(goalc, ["crp_iqr", "lncysc_c", "cancre"] + COV_C, "new_cvd", "t_cvd",
        cluster="communityID", label="C. 主模型 + ln(cystatin C) + cancre")
# D. + 血脂
run_cox(goalc, ["crp_iqr", "lntg_c", "lnhdl_c"] + COV_C, "new_cvd", "t_cvd",
        cluster="communityID", label="D. 主模型 + ln(TG) + ln(HDL)")
# C+D 联合
run_cox(goalc, ["crp_iqr", "lncysc_c", "cancre", "lntg_c", "lnhdl_c"] + COV_C, "new_cvd", "t_cvd",
        cluster="communityID", label="C+D. 全协变量联合模型")

with open(OUT + r"\covariate_audit_sensitivity_results.txt", "w", encoding="utf-8") as f:
    f.write("协变量审计敏感性分析（文献惯例对照补跑）2026-09-11\n")
    f.write("=" * 60 + "\n\n")
    f.write("\n".join(log_lines))
log("\n完成。结果已存 covariate_audit_sensitivity_results.txt")
