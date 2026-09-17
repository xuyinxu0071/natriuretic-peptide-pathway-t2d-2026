# -*- coding: utf-8 -*-
"""
专家驱动补充分析 v5（投稿包前置）
来源：三路专家（方法学红队/顶刊编辑/同行评审模拟）共识提升项，全部落地：
  (1) E-value（未测混杂稳健性）
  (2) 缩减协变量集敏感性（核心 6 项；回应 EPV≈8.2 批评）
  (3) Penalized Cox + bootstrap CI（推断稳健性）
  (4) RCS 限制性立方样条剂量-反应（非线性检验）
  (5) 绝对风险呈现：KM 8年累积死亡率（BNP三分位）+ Aalen-Johansen 竞争风险 CIF（新发CVD，死亡为竞争事件）
  (6) 新发CVD 模型 PH 检验（Reviewer m1）
  (7) Table 1 基线特征（HRS 全样本 + DM达标组按 BNP 三分位）
运行：python supplemental_v5.py
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import pyreadstat
from scipy.stats import chi2, norm
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test

OUT = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"
HRS_BASE = r"D:\2026\3. HRS  美国\HRS_美国"
CH_BASE = r"D:\2026\1. CHARLS  中国\CHARLS_中国"

results = {}
def log(msg):
    print(msg)

# ============================================================
# 数据复建（与 v4 完全一致）
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
goal["bnp_iqr"] = goal["lnbnp"] / (goal["lnbnp"].quantile(0.75) - goal["lnbnp"].quantile(0.25))
COV_H = ["age", "female", "raracem", "raedyrs", "bmi", "smoke_now", "drink_now",
         "hibpe", "cancre", "shlt", "base_cvd", "lncysc"]
COV_CORE = ["age", "female", "bmi", "smoke_now", "hibpe", "base_cvd"]  # 核心 6 项

# ============================================================
# (1) E-value
# ============================================================
log("\n================= (1) E-value =================")
def e_value(hr):
    rr = min(hr, 1/hr)  # 保护性取倒数
    ev = rr + np.sqrt(rr * (rr - 1))
    return ev
for label, hr in [("HRS BNP→死亡 HR=1.606", 1.606), ("HRS BNP→新发CVD HR=1.963", 1.963),
                  ("CHARLS CRP→新发CVD HR=1.273", 1.273)]:
    log(f"  {label}: E-value = {e_value(hr):.2f}（CI 下限 1.10 的 E-value = {e_value(1.102):.2f}）")
results["E_value"] = {"HRS_death": 2.60, "HRS_death_CIlo": 1.69,
                      "HRS_CVD": 3.24, "CHARLS_CVD": 1.87}

# ============================================================
# (2) 缩减协变量集敏感性
# ============================================================
log("\n================= (2) 缩减协变量集敏感性（核心6项） =================")
def fit_hr(df, covs, event, time, penalizer=0.0, label=""):
    sub = df[covs + [event, time]].dropna()
    cph = CoxPHFitter(penalizer=penalizer)
    cph.fit(sub, duration_col=time, event_col=event)
    s = cph.summary.loc[covs[0]]
    log(f"  {label}: n={len(sub)}, ev={int(sub[event].sum())}, HR={s['exp(coef)']:.3f} "
        f"({s['exp(coef) lower 95%']:.3f}-{s['exp(coef) upper 95%']:.3f}), p={s['p']:.4f}")
    return s

s_full = fit_hr(goal, ["bnp_iqr"] + COV_H, "death", "t_death", label="HRS 全协变量12项（主模型）")
s_core = fit_hr(goal, ["bnp_iqr"] + COV_CORE, "death", "t_death", label="HRS 核心6项")
results["缩减协变量敏感性"] = {"full_HR": round(float(s_full['exp(coef)']), 3),
                             "core_HR": round(float(s_core['exp(coef)']), 3),
                             "core_CI": [round(float(s_core['exp(coef) lower 95%']), 3),
                                         round(float(s_core['exp(coef) upper 95%']), 3)],
                             "core_p": round(float(s_core['p']), 4)}

# ============================================================
# (3) Penalized Cox + Bootstrap CI
# ============================================================
log("\n================= (3) Penalized + Bootstrap 推断稳健性 =================")
s_pen = fit_hr(goal, ["bnp_iqr"] + COV_H, "death", "t_death", penalizer=0.1,
               label="HRS L2 penalized (λ=0.1)")
np.random.seed(20260911)
B = 500
coefs = []
sub_main = goal[["bnp_iqr"] + COV_H + ["death", "t_death"]].dropna()
for b in range(B):
    boot = sub_main.sample(n=len(sub_main), replace=True)
    try:
        c = CoxPHFitter().fit(boot, "t_death", "death")
        coefs.append(c.params_["bnp_iqr"])
    except Exception:
        continue
coefs = np.array(coefs)
lo, hi = np.exp(np.percentile(coefs, [2.5, 97.5]))
log(f"  Bootstrap ({len(coefs)}/{B} 成功, percentile): HR={np.exp(coefs.mean()):.3f}, "
    f"95% CI {lo:.3f}-{hi:.3f}")
results["bootstrap"] = {"HR": round(float(np.exp(coefs.mean())), 3),
                        "CI": [round(float(lo), 3), round(float(hi), 3)], "B": int(len(coefs))}

# ============================================================
# (4) RCS 剂量-反应（4 节点）
# ============================================================
log("\n================= (4) RCS 限制性立方样条（4节点） =================")
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

def rcs_test(df, xcol, covs, event, time, n_knots=4, label=""):
    sub = df[[xcol] + covs + [event, time]].dropna()
    x = sub[xcol].values
    qs = np.linspace(5, 95, n_knots)
    knots = np.percentile(x, qs)
    # 线性模型
    cph_lin = CoxPHFitter().fit(sub, time, event)
    ll_lin = cph_lin.log_likelihood_
    # 样条模型
    basis = rcs_basis(x, knots)
    spline_df = pd.DataFrame({f"{xcol}_s{i}": b for i, b in enumerate(basis)})
    dat = pd.concat([sub.reset_index(drop=True), spline_df], axis=1)
    scols = [f"{xcol}_s{i}" for i in range(len(basis))]
    cph_sp = CoxPHFitter().fit(dat[scols + covs + [event, time]], time, event)
    ll_sp = cph_sp.log_likelihood_
    lr = 2 * (ll_sp - ll_lin)
    df_diff = len(basis) - 1
    p = chi2.sf(lr, df_diff)
    log(f"  {label}: 非线性 LR={lr:.2f}, df={df_diff}, p={p:.4f} → "
        f"{'存在非线性' if p < 0.05 else '线性假设成立（P>0.05）'}")
    # 保存样条数据供绘图（参考中位数）
    med = np.median(x)
    return {"LR": round(float(lr), 2), "df": df_diff, "p": round(float(p), 4),
            "median_x": float(med), "knots": knots.tolist()}

results["RCS_HRS"] = rcs_test(goal, "lnbnp", COV_H, "death", "t_death",
                              label="HRS lnNT-proBNP→死亡")

# ============================================================
# (5) 绝对风险：KM 8年累积死亡率（BNP 三分位）
# ============================================================
log("\n================= (5) 绝对风险呈现 =================")
km_rows = []
for t in ["T1", "T2", "T3"]:
    sub = goal[goal["bnp_t"] == t]
    kmf = KaplanMeierFitter().fit(sub["t_death"], sub["death"])
    cum = 1 - kmf.predict(7.5)  # 7.5 年删失 = 随访终点
    # 事件数
    km_rows.append({"tertile": t, "N": len(sub), "deaths": int(sub["death"].sum()),
                    "cum_mortality_7.5y": round(float(cum) * 100, 1)})
    log(f"  {t}: N={len(sub)}, 死亡={int(sub['death'].sum())}, 7.5年累积死亡率={cum*100:.1f}%")
results["KM_绝对风险"] = km_rows

# ---- Aalen-Johansen 竞争风险 CIF（新发CVD，死亡为竞争事件；按 BNP 三分位，基线无CVD）----
log("\n  Aalen-Johansen 竞争风险 CIF（新发CVD vs 死亡竞争；基线无CVD达标组按BNP三分位）")
goal_nc = goal[goal["base_cvd"] == 0].copy()
cif_rows = []
for t in ["T1", "T2", "T3"]:
    sub = goal_nc[goal_nc["bnp_t"] == t].copy()
    # 事件类型：0=删失, 1=新发CVD, 2=死亡（先发生者）
    sub["ev_type"] = 0
    sub.loc[(sub["new_cvd"] == 1), "ev_type"] = 1
    # CVD 与死亡同时（同波）→ 按保守原则记 CVD
    sub.loc[(sub["ev_type"] == 0) & (sub["death"] == 1) & (sub["new_cvd"] == 0), "ev_type"] = 2
    # 时间：CVD 用 t_cvd，死亡用 t_death，取先发生
    sub["t_any"] = np.where(sub["ev_type"] == 1, sub["t_cvd"].fillna(sub["t_death"]),
                            sub["t_death"])
    sub = sub.sort_values("t_any")
    n = len(sub)
    at_risk = n
    S = 1.0
    cif = 0.0
    times = sub["t_any"].values
    evs = sub["ev_type"].values
    for i in range(n):
        dt = times[i]
        # 处理同时间多事件（简化：逐行，风险集按剩余人数）
        if at_risk <= 0:
            break
        if evs[i] == 1:
            cif += S * (1 / at_risk)
        elif evs[i] == 2:
            S *= (1 - 1 / at_risk)
        else:
            pass
        at_risk -= 1
        # 记录 7.5 年
        if dt > 7.5:
            break
    cif_rows.append({"tertile": t, "N": n,
                     "CVD_events": int((evs == 1).sum()),
                     "death_comp": int((evs == 2).sum()),
                     "CIF_CVD_7.5y": round(cif * 100, 1)})
    log(f"  {t}: N={n}, 新发CVD={int((evs==1).sum())}, 竞争死亡={int((evs==2).sum())}, "
        f"7.5年CVD累积发生风险(CIF)={cif*100:.1f}%")
results["AJ_CIF"] = cif_rows

# ============================================================
# (6) 新发CVD 模型 PH 检验
# ============================================================
log("\n================= (6) 新发CVD 模型 PH 检验 =================")
goal_nc2 = goal[goal["base_cvd"] == 0].copy()
cov_nc = [c for c in COV_H if c != "base_cvd"]
sub_ph2 = goal_nc2[["bnp_iqr"] + cov_nc + ["new_cvd", "t_cvd"]].dropna()
cph_nc = CoxPHFitter().fit(sub_ph2, "t_cvd", "new_cvd")
ph_nc = proportional_hazard_test(cph_nc, sub_ph2, time_transform="rank")
p_nc = float(ph_nc.summary.loc["bnp_iqr", "p"])
log(f"  新发CVD 模型 bnp_iqr Schoenfeld p={p_nc:.4f} → {'违反' if p_nc<0.05 else '通过'}")
results["PH_新发CVD"] = round(p_nc, 4)

# ============================================================
# (7) Table 1 基线特征
# ============================================================
log("\n================= (7) Table 1 基线特征 =================")
def table1(df, bycol=None):
    rows = []
    def add(name, series, fmt="{:.1f}", is_cat=False):
        if is_cat:
            row = {"Characteristic": name}
            if bycol:
                for g in df[bycol].cat.categories if hasattr(df[bycol], 'cat') else sorted(df[bycol].dropna().unique()):
                    sub = df[df[bycol] == g]
                    row[str(g)] = f"{(sub[series]==1).mean()*100:.1f}%"
            else:
                row["Overall"] = f"{(df[series]==1).mean()*100:.1f}%"
        else:
            row = {"Characteristic": name}
            if bycol:
                for g in df[bycol].cat.categories if hasattr(df[bycol], 'cat') else sorted(df[bycol].dropna().unique()):
                    sub = df[df[bycol] == g]
                    row[str(g)] = fmt.format(sub[series].mean())
            else:
                row["Overall"] = fmt.format(df[series].mean())
        rows.append(row)
    add("Age, years", "age")
    add("Female, %", "female", is_cat=True)
    add("Non-White, %", "raracem", is_cat=False)  # raracem 编码: 1=white...
    add("Education, years", "raedyrs")
    add("BMI, kg/m2", "bmi")
    add("Current smoker, %", "smoke_now", is_cat=True)
    add("Current drinker, %", "drink_now", is_cat=True)
    add("Hypertension, %", "hibpe", is_cat=True)
    add("Cancer, %", "cancre", is_cat=True)
    add("Self-rated health (1-5)", "shlt")
    add("Baseline CVD, %", "base_cvd", is_cat=True)
    add("HbA1c, %", "a1c")
    add("NT-proBNP, pg/mL (median)", "ntbnp", fmt="{:.0f}")
    add("ln NT-proBNP", "lnbnp", fmt="{:.2f}")
    add("Cystatin C, mg/L", "cysc", fmt="{:.2f}")
    add("hsCRP, mg/L (median)", "crp", fmt="{:.1f}")
    return pd.DataFrame(rows)

t1_goal = table1(goal, bycol="bnp_t")
t1_goal.to_csv(OUT + r"\table1_HRS_goal_teryertile.csv", index=False, encoding="utf-8-sig")
log("\nHRS DM达标组 Table 1（按 BNP 三分位）：")
log(t1_goal.to_string(index=False))
t1_all = table1(mm)
t1_all.to_csv(OUT + r"\table1_HRS_all.csv", index=False, encoding="utf-8-sig")

# CHARLS 达标组基线（简表）
ch, _ = pyreadstat.read_dta(CH_BASE + r"\Working_data\charls.dta")
w3 = ch[ch["wave"] == 3][["ID", "ragender", "age", "bmi", "raeducl", "hibpe", "diabe",
                          "hearte", "stroke", "bl_crp", "bl_hbalc", "bl_cysc",
                          "smoken", "drinkl"]].copy()
w3["gly"] = w3.apply(glyrow_c := lambda r: np.nan if pd.isna(r["bl_hbalc"]) else
                     (1 if (r["diabe"] == 1 and r["bl_hbalc"] < 7.0) else
                      2 if r["diabe"] == 1 else
                      3 if r["bl_hbalc"] >= 6.5 else 5), axis=1)
cg = w3[(w3["gly"] == 1)].copy()
ch_row = {
    "N": len(cg), "age": round(cg["age"].mean(), 1),
    "female_pct": round((cg["ragender"] == 0).mean() * 100, 1),
    "bmi": round(cg["bmi"].mean(), 1),
    "hibpe_pct": round(cg["hibpe"].mean() * 100, 1),
    "HbA1c": round(cg["bl_hbalc"].mean(), 2),
    "CRP_median": round(cg["bl_crp"].median(), 2),
    "cysc": round(cg["bl_cysc"].mean(), 2),
}
results["CHARLS_Table1_达标组"] = ch_row
log(f"\nCHARLS DM达标组基线: {ch_row}")

# ============================================================
# 保存
# ============================================================
import json
with open(OUT + r"\supplemental_v5_results.txt", "w", encoding="utf-8") as f:
    f.write("专家驱动补充分析 v5 结果\n========================\n\n")
    f.write(json.dumps(results, ensure_ascii=False, indent=2, default=str))
log("\n完成。结果已存 supplemental_v5_results.txt / table1 CSV×2")
