# -*- coding: utf-8 -*-
"""
投稿前专家驱动终版补充分析 v6（基于 v4 数据结构）
四个目标：
  (1) 逆转因果敏感性：剔除随访 12/24 个月内死亡
  (2) 剔除基线 CVD 的死亡敏感性（排除"已知疾病标志物"解释）
  (3) BNP 与 CRP 相互调整（心脏应激轴 vs 炎症轴独立性）
  (4) 达标组 lnNT-proBNP IQR 真实数值核定（手稿 2.2 数字来源）
运行：python pre_submission_v6.py
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import pyreadstat
from lifelines import CoxPHFitter

OUT = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"
HRS_BASE = r"D:\2026\3. HRS  美国\HRS_美国"

results = []
def log(msg):
    print(msg)

# ---------- 数据构建（同 v4 PART 1）----------
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
d["t_death"] = np.where(d["death"] == 1,
                        np.maximum(t_ev - 2016.5, 0.1), 7.5)

d["base_cvd"] = ((d["hearte"] == 1) | (d["stroke"] == 1)).astype(int)
d["lnbnp"] = np.log(d["ntbnp"].clip(lower=1))
d["lncrp"] = np.log(d["crp"].clip(lower=0.1))
d["lncysc"] = np.log(d["cysc"].clip(lower=0.1))
d["female"] = (d["ragender"] == 0).astype(int)
d["smoke_now"] = (d["smoken"] == 1).astype(int)
d["drink_now"] = (d["drink"] == 1).astype(int)

COV_H = ["age", "female", "raracem", "raedyrs", "bmi", "smoke_now", "drink_now",
         "hibpe", "cancre", "shlt", "base_cvd", "lncysc"]

mm = d[d["gly"].notna() & d["ntbnp"].notna()].copy()
goal = mm[mm["gly"] == 1].copy()
goal["bnp_iqr"] = goal["lnbnp"] / (goal["lnbnp"].quantile(0.75) - goal["lnbnp"].quantile(0.25))

def run_cox(df, covs, event, time, label=""):
    cols = covs + [event, time]
    sub = df[cols].dropna().copy()
    cph = CoxPHFitter()
    cph.fit(sub, duration_col=time, event_col=event, robust=True)
    s = cph.summary
    log(f"\n--- {label} (n={len(sub)}, events={int(sub[event].sum())}) ---")
    log(s.loc[covs[0], ["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].round(3).to_string())
    return s

# ---------- (4) IQR 核定 ----------
q25, q75 = goal["lnbnp"].quantile(0.25), goal["lnbnp"].quantile(0.75)
iqr_ln = q75 - q25
log(f"\n[IQR核定] 达标组 lnNT-proBNP: Q1={q25:.2f}, Q3={q75:.2f}, IQR={iqr_ln:.2f} ln-units")
log(f"           对应 NT-proBNP(pg/mL): Q1={np.exp(q25):.0f}, Q3={np.exp(q75):.0f}, IQR倍数={np.exp(iqr_ln):.2f}x")
results.append(("IQR_核定", {"Q1_ln": round(q25, 2), "Q3_ln": round(q75, 2),
                             "IQR_ln": round(iqr_ln, 2),
                             "Q1_pg": round(float(np.exp(q25)), 0),
                             "Q3_pg": round(float(np.exp(q75)), 0)}))

# ---------- (1) 逆转因果敏感性 ----------
log("\n================= (1) 逆转因果敏感性：剔除早期死亡 =================")
for cut_y in [1.0, 2.0]:
    g = goal.copy()
    early = (g["death"] == 1) & (g["t_death"] < cut_y)
    n_early = int(early.sum())
    g2 = g[~early].copy()  # 剔除早期死亡者（事件），保留其余
    # 也可做 landmark（t 起点平移）；此处用剔除法（conservative，等同 exclude early deaths）
    s = run_cox(g2, ["bnp_iqr"] + COV_H, "death", "t_death",
                label=f"剔除 <{cut_y:.0f} 年死亡（早期死亡 {n_early} 例）后主模型")
    results.append((f"逆转因果_剔除{cut_y:.0f}年", 
                    s.loc["bnp_iqr", ["exp(coef)", "exp(coef) lower 95%",
                                      "exp(coef) upper 95%", "p"]].round(3).to_dict()))

# ---------- (2) 剔除基线 CVD 的死亡敏感性 ----------
log("\n================= (2) 剔除基线 CVD：BNP→死亡（基线无 CVD 达标组） =================")
g_ncvd = goal[goal["base_cvd"] == 0].copy()
s = run_cox(g_ncvd, ["bnp_iqr"] + [c for c in COV_H if c != "base_cvd"], "death", "t_death",
              label="基线无CVD 达标组 per-IQR BNP→全因死亡（+cysc）")
results.append(("剔除基线CVD_死亡", 
                s.loc["bnp_iqr", ["exp(coef)", "exp(coef) lower 95%",
                                  "exp(coef) upper 95%", "p"]].round(3).to_dict()))

# ---------- (3) BNP 与 CRP 相互调整 ----------
log("\n================= (3) 心脏应激轴 vs 炎症轴：相互调整 =================")
goal["crp_iqr"] = goal["lncrp"] / (goal["lncrp"].quantile(0.75) - goal["lncrp"].quantile(0.25))
# CRP 单独（达标组，死亡）——此前未正式报告
s = run_cox(goal, ["crp_iqr"] + COV_H, "death", "t_death",
              label="达标组 per-IQR lnCRP→全因死亡（单暴露，+cysc）")
results.append(("CRP单独_死亡", 
                s.loc["crp_iqr", ["exp(coef)", "exp(coef) lower 95%",
                                  "exp(coef) upper 95%", "p"]].round(3).to_dict()))
# 双标志物模型
s = run_cox(goal, ["bnp_iqr", "crp_iqr"] + COV_H, "death", "t_death",
              label="双标志物模型：BNP 与 CRP 相互调整（+cysc）")
results.append(("双标志物_BNP", 
                s.loc["bnp_iqr", ["exp(coef)", "exp(coef) lower 95%",
                                  "exp(coef) upper 95%", "p"]].round(3).to_dict()))
results.append(("双标志物_CRP", 
                s.loc["crp_iqr", ["exp(coef)", "exp(coef) lower 95%",
                                  "exp(coef) upper 95%", "p"]].round(3).to_dict()))

# ---------- 保存 ----------
with open(OUT + r"\pre_submission_v6_results.txt", "w", encoding="utf-8") as f:
    f.write("投稿前专家驱动终版补充分析 v6\n"
            "(1) 逆转因果剔除 (2) 剔除基线CVD (3) BNP/CRP 相互调整 (4) IQR 核定\n"
            "========================\n\n")
    for k, v in results:
        f.write(f"[{k}]\n{v}\n\n")
log("\n完成。结果已存 pre_submission_v6_results.txt")
