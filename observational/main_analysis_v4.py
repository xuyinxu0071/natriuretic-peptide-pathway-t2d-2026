# -*- coding: utf-8 -*-
"""
CD 特刊主分析 v4（协变量处置终版）
v3 → v4 变更（协变量科学处置，其余方法学沿用 v3）：
  (1) HRS 臂：主模型协变量集 COV_H 全程加入 ln(cystatin C)（肾功能调整）
      —— NT-proBNP 经肾清除，属 DAG 混杂因子；审稿人必问项，主动调整为主结果。
      未调整版降级为敏感性（2b-sens）。
      波及模型：2a 三分位 / 2b 主效应 / 2c 新发CVD / 2d 交互 / 2e 分层 /
                2f SABV / 2g 加权 / 2h PH / 2i MICE / 2k CVD死亡。
  (2) CHARLS 臂：主模型协变量集 COV_C 维持 v3 的 7 项不变（DAG 最小充分调整集；
      cysc/cancre 与血脂均按预先定位为敏感性）。
      本版新增正式输出两项敏感性供 Supplement：
      3s1: + ln(cysc) + cancre（肾功能+癌症，完全稳健）
      3s2: + ln(TG) + ln(HDL)（血脂轴，点估计稳定、显著性边缘）
  (3) 血脂在两臂均不入主模型（非混杂：与暴露同源/无因果箭头，实测信息量近零）。
运行：python main_analysis_v4.py
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import pyreadstat
from scipy.stats import chi2, norm
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from statsmodels.imputation.mice import MICEData

OUT = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"
HRS_BASE = r"D:\2026\3. HRS  美国\HRS_美国"
CH_BASE = r"D:\2026\1. CHARLS  中国\CHARLS_中国"

results = []
def log(msg):
    print(msg)

# ============================================================
# PART 1 —— HRS 臂：建分析数据集（同 v3，月粒度死亡时间）
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

fu = hrs[hrs["wave"].isin([13, 14, 15])][["hhidpn", "wave", "hearte", "stroke"]]
fu = fu.sort_values(["hhidpn", "wave"])
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
GLY = {1: "DM达标", 2: "DM不达标", 3: "未诊断高A1C", 4: "前期", 5: "无DM"}

d["death"] = ((d["EXDEATHYR"] >= 2016) & (d["EXDEATHYR"] <= 2023)).astype(int)
t_ev = np.where(d["EXDEATHMO"].notna(),
                d["EXDEATHYR"] + (d["EXDEATHMO"].astype(float) - 0.5) / 12.0,
                d["EXDEATHYR"] + 0.5)
d["t_death"] = np.where(d["death"] == 1,
                        np.maximum(t_ev - 2016.5, 0.1), 7.5)

d["base_cvd"] = ((d["hearte"] == 1) | (d["stroke"] == 1)).astype(int)
d["new_cvd"] = 0
d.loc[(d["base_cvd"] == 0) & d["ev_wave"].notna(), "new_cvd"] = 1
d["t_cvd"] = np.where(d["new_cvd"] == 1,
                      (d["ev_wave"] - 12) * 2.0,
                      np.where(d["last_wave"].notna(), (d["last_wave"] - 12) * 2.0, np.nan))
d.loc[d["t_cvd"].isna() & (d["base_cvd"] == 0) & d["gly"].notna(), "t_cvd"] = 5.0

d["lnbnp"] = np.log(d["ntbnp"].clip(lower=1))
d["lncrp"] = np.log(d["crp"].clip(lower=0.1))
d["lncysc"] = np.log(d["cysc"].clip(lower=0.1))
d["female"] = (d["ragender"] == 0).astype(int)
d["smoke_now"] = (d["smoken"] == 1).astype(int)
d["drink_now"] = (d["drink"] == 1).astype(int)

log("\n================= HRS 臂：样本清点（v4 = v3 数据结构不变） =================")
mm = d[d["gly"].notna() & d["ntbnp"].notna()].copy()
log(f"分析总样本 N={len(mm)}（A1C+NT-proBNP 齐备）；cysc 非缺失 {mm['lncysc'].notna().sum()}/{len(mm)}"
    f"（{100*mm['lncysc'].notna().mean():.1f}%）")
tab = mm.groupby("gly").agg(N=("death", "size"), 死亡=("death", "sum"),
                            新发CVD=("new_cvd", "sum")).rename(index=GLY)
log("\n" + tab.to_string())

goal = mm[mm["gly"] == 1].copy()
goal["bnp_t"] = pd.qcut(goal["ntbnp"], 3, labels=["T1", "T2", "T3"])
log("\nDM达标组 NT-proBNP 三分位：")
log(goal.groupby("bnp_t", observed=True).agg(N=("death", "size"), 死亡=("death", "sum"),
                                             新发CVD=("new_cvd", "sum")).to_string())
results.append(("HRS_达标组_BNP三分位", goal.groupby("bnp_t", observed=True)
               .agg(N=("death", "size"), 死亡=("death", "sum"), 新发CVD=("new_cvd", "sum"))
               .to_csv()))

# ---- v4 核心：协变量集升级 ----
# 主模型协变量集（v4）：在 v3 的 11 项上加入 lncysc（肾功能调整）
COV_H = ["age", "female", "raracem", "raedyrs", "bmi", "smoke_now", "drink_now",
         "hibpe", "cancre", "shlt", "base_cvd", "lncysc"]
# 敏感性协变量集（= v3 主模型，未调整肾功能）
COV_H_NO肾 = [c for c in COV_H if c != "lncysc"]
log(f"\n【v4 主模型协变量集】{len(COV_H)} 项（v3 + ln(cystatin C)）；"
    f"未调整版降级为敏感性")

# ============================================================
# PART 2 —— HRS 臂：主分析（v4：全程肾功能调整）
# ============================================================
def run_cox(df, covs, event, time, weights=None, label="", cluster=None):
    cols = covs + [event, time] + ([weights] if weights else []) + ([cluster] if cluster else [])
    sub = df[cols].dropna().copy()
    cph = CoxPHFitter()
    cph.fit(sub, duration_col=time, event_col=event,
            weights_col=weights, robust=True,
            cluster_col=cluster if cluster else None)
    s = cph.summary
    log(f"\n--- {label} (n={len(sub)}, events={int(sub[event].sum())}) ---")
    log(s[["coef", "exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].round(3).to_string())
    return cph, s

log("\n================= HRS 臂：主分析（v4，肾功能调整） =================")
goal2 = goal.copy()
goal2["bnp_T2"] = (goal2["bnp_t"] == "T2").astype(int)
goal2["bnp_T3"] = (goal2["bnp_t"] == "T3").astype(int)
run_cox(goal2, ["bnp_T3", "bnp_T2"] + COV_H, "death", "t_death",
        label="2a. DM达标组 NT-proBNP 三分位→全因死亡（+cysc）")

goal2["bnp_iqr"] = goal2["lnbnp"] / (goal2["lnbnp"].quantile(0.75) - goal2["lnbnp"].quantile(0.25))
cph2b, s2b = run_cox(goal2, ["bnp_iqr"] + COV_H, "death", "t_death",
                     label="2b. 主效应：per-IQR lnNT-proBNP→全因死亡（+cysc，月粒度）【v4 主结果】")
results.append(("HRS_2b_主效应_肾功能调整", s2b.loc["bnp_iqr", ["exp(coef)", "exp(coef) lower 95%",
                "exp(coef) upper 95%", "p"]].round(3).to_dict()))

# 2b-sens. 未调整肾功能（v3 主模型）→ 降级为敏感性
_, s2bs = run_cox(goal2, ["bnp_iqr"] + COV_H_NO肾, "death", "t_death",
                  label="2b-sens. 同模型不含 cysc（v3 主模型，降级为敏感性）")
results.append(("HRS_2b_sens_未调整肾", s2bs.loc["bnp_iqr", ["exp(coef)", "exp(coef) lower 95%",
                "exp(coef) upper 95%", "p"]].round(3).to_dict()))

goal3 = goal2[goal2["base_cvd"] == 0].copy()
_, s2c = run_cox(goal3, ["bnp_iqr"] + [c for c in COV_H if c != "base_cvd"], "new_cvd", "t_cvd",
                 label="2c. DM达标组(基线无CVD) per-IQR lnNT-proBNP→新发CVD(至2022, +cysc)")
results.append(("HRS_2c_新发CVD_肾功能调整", s2c.loc["bnp_iqr", ["exp(coef)", "exp(coef) lower 95%",
                "exp(coef) upper 95%", "p"]].round(3).to_dict()))

# 2d. 全人群交互（嵌套 LR；交互模型同样 +cysc）
full = mm.copy()
full["bnp_iqr"] = full["lnbnp"] / (full["lnbnp"].quantile(0.75) - full["lnbnp"].quantile(0.25))
full = pd.get_dummies(full, columns=["gly"], prefix="gly", drop_first=True, dtype=int)
gcols = [c for c in full.columns if c.startswith("gly_")]
inter_cols = [f"{g}_x" for g in gcols]
for g in gcols:
    full[f"{g}_x"] = full[g] * full["bnp_iqr"]
base_covs = ["bnp_iqr"] + gcols + COV_H
m_full, _ = run_cox(full, base_covs + inter_cols, "death", "t_death",
                    label="2d. 全人群 gly×BNP 交互模型（+cysc）")
m_base, _ = run_cox(full, base_covs, "death", "t_death", label="2d-nested")
lr = 2 * (m_full.log_likelihood_ - m_base.log_likelihood_)
log(f"\n嵌套LR检验（交互项联合）: LR={lr:.2f}, df={len(inter_cols)}, p={chi2.sf(lr, len(inter_cols)):.4f}")
results.append(("HRS_交互LR检验", f"LR={lr:.2f}, df={len(inter_cols)}, p={chi2.sf(lr, len(inter_cols)):.5f}"))

# 2e. 分层 HR（+cysc）
log("\n--- 2e. 各血糖层 per-IQR lnNT-proBNP→死亡 HR（+cysc，月粒度） ---")
for g in sorted(mm["gly"].unique()):
    sub = mm[mm["gly"] == g].copy()
    sub["bnp_iqr"] = sub["lnbnp"] / (mm["lnbnp"].quantile(0.75) - mm["lnbnp"].quantile(0.25))
    try:
        cph, s = run_cox(sub, ["bnp_iqr"] + COV_H, "death", "t_death", label=f"  {GLY[g]}")
        results.append((f"HR_{GLY[g]}_肾功能调整", s.loc["bnp_iqr", ["exp(coef)", "p"]].round(3).to_dict()))
    except Exception as e:
        log(f"  {GLY[g]}: 拟合失败 {e}")

# 2f. SABV（+cysc）
log("\n--- 2f. SABV 分性别 per-IQR BNP→死亡（+cysc，月粒度） ---")
for gen, name in [(0, "女性"), (1, "男性")]:
    sub = goal2[goal2["ragender"] == gen].copy()
    try:
        cph, s = run_cox(sub, ["bnp_iqr"] + [c for c in COV_H if c != "female"],
                         "death", "t_death", label=f"  {name}")
        results.append((f"SABV_{name}_肾功能调整", s.loc["bnp_iqr", ["exp(coef)", "p"]].round(3).to_dict()))
    except Exception as e:
        log(f"  {name}: 拟合失败 {e}")

# 2g. 加权敏感性（+cysc）
sub = goal2[goal2["vbswgt"].notna() & (goal2["vbswgt"] > 0)].copy()
run_cox(sub, ["bnp_iqr"] + COV_H, "death", "t_death", weights="vbswgt",
        label="2g. 达标组 加权敏感性(VBS权重, +cysc)")

# 2h. PH 检验（v4 主模型 2b）
log("\n================= 2h. PH 检验（v4 主模型 2b） =================")
sub_ph = goal2[["bnp_iqr"] + COV_H + ["death", "t_death"]].dropna()
ph = proportional_hazard_test(cph2b, sub_ph, time_transform="rank")
ph_s = ph.summary
log("\nSchoenfeld 全局/逐变量 p 值（rank 变换）：")
log(ph_s[["test_statistic", "p"]].round(4).to_string())
results.append(("HRS_PH检验_bnp_iqr", float(ph_s.loc["bnp_iqr", "p"])))
if ph_s.loc["bnp_iqr", "p"] < 0.05:
    log("⚠️ bnp_iqr 违反 PH 假设 → 补充运行时间分段敏感性")
    goal2["late"] = (goal2["t_death"] > 2).astype(int)
    goal2["bnp_late"] = goal2["bnp_iqr"] * goal2["late"]
    try:
        run_cox(goal2, ["bnp_iqr", "bnp_late", "late"] + COV_H, "death", "t_death",
                label="2h-sens. 时间交互（0-2 年 vs >2 年）")
    except Exception as e:
        log(f"时间交互拟合失败: {e}")
else:
    log("bnp_iqr PH 检验通过，主模型无需时间交互修正")

# 2i. MICE（m=20；v4 插补变量集 + lncysc）
log("\n================= 2i. MICE 多重插补 m=20（v4 主模型，含 cysc） =================")
mice_cols = ["bnp_iqr", "death", "age", "female", "raracem", "raedyrs", "bmi",
             "smoke_now", "drink_now", "hibpe", "cancre", "shlt", "base_cvd", "lncysc"]
md = goal2[mice_cols].copy()
md["lt"] = np.log(pd.to_numeric(goal2["t_death"], errors="coerce").clip(lower=0.1))
n_cc = goal2[mice_cols].dropna().shape[0]
log(f"达标组 N={len(goal2)}，完全病例 n={n_cc}（插补前缺失 {len(goal2)-n_cc} 人）")

M = 20
imp = MICEData(md, perturbation_method="boot")
imp.update_all(10)
ests, vars_ = [], []
bin_cols = ["death", "female", "smoke_now", "drink_now", "hibpe", "cancre", "base_cvd"]
for i in range(M):
    imp.update_all(1)
    comp = imp.data.copy()
    for b in bin_cols:
        comp[b] = (comp[b] > 0.5).astype(int)
    comp["t"] = np.exp(comp["lt"])
    try:
        cph_m = CoxPHFitter()
        cph_m.fit(comp[mice_cols + ["t"]], "t", "death")
        ests.append(cph_m.params_["bnp_iqr"])
        vars_.append(cph_m.standard_errors_["bnp_iqr"] ** 2)
    except Exception:
        continue

if len(ests) >= 10:
    ests = np.array(ests); vars_ = np.array(vars_)
    Qbar = ests.mean(); Ubar = vars_.mean(); B = ests.var(ddof=1)
    Tvar = Ubar + (1 + 1 / len(ests)) * B
    se = np.sqrt(Tvar)
    z = Qbar / se
    pval = 2 * norm.sf(abs(z))
    hr = np.exp(Qbar)
    lo, hi = np.exp(Qbar - 1.96 * se), np.exp(Qbar + 1.96 * se)
    gamma = (1 + 1 / len(ests)) * B / Tvar
    log(f"\nMICE 合并（m={len(ests)}）: per-IQR lnNT-proBNP→死亡（肾功能调整）")
    log(f"  HR = {hr:.2f} (95% CI {lo:.2f}-{hi:.2f}), p = {pval:.4f}")
    log(f"  组内方差 Ubar={Ubar:.4f}, 组间方差 B={B:.4f}, 缺失信息比例 γ={gamma:.2f}")
    results.append(("HRS_MICE合并_肾功能调整", {"HR": round(hr, 3), "CI_lo": round(lo, 3),
                                     "CI_hi": round(hi, 3), "p": round(pval, 4),
                                     "m": len(ests), "gamma": round(gamma, 2)}))
    log(f"  与完全病例 2b（HR={s2b.loc['bnp_iqr','exp(coef)']:.2f}）一致性："
        f"{'方向一致' if (Qbar > 0) == (s2b.loc['bnp_iqr','coef'] > 0) else '方向不一致'}")
else:
    log("MICE 收敛失败次数过多，报告完全病例结果为准")

# 2j. QBA（结局误分类；与协变量无关，沿用 v3 逻辑）
log("\n================= 2j. QBA：新发 CVD 自报结局非差异化误分类校正 =================")
log("注：HRS 主终点（全因死亡，Tracker/NDI 链接）为登记式确定，无需校正")

def qba_table(hr_obs, se_obs, label):
    rows = []
    for se_ in [0.65, 0.75, 0.85, 0.90]:
        for sp_ in [0.85, 0.90, 0.95, 0.98]:
            f = se_ + sp_ - 1
            if f <= 0.3:
                continue
            lhr = np.log(hr_obs) / f
            lse = se_obs / f
            lo, hi = np.exp(lhr - 1.96 * lse), np.exp(lhr + 1.96 * lse)
            rows.append({"Se": se_, "Sp": sp_, "衰减因子": round(f, 2),
                         "校正HR": round(np.exp(lhr), 2),
                         "95%CI": f"{lo:.2f}-{hi:.2f}"})
    t = pd.DataFrame(rows)
    log(f"\n--- {label}（观测 HR={hr_obs:.2f}） ---")
    log(t.to_string(index=False))
    return t

hr_c = float(s2c.loc["bnp_iqr", "exp(coef)"])
se_c = float(s2c.loc["bnp_iqr", "se(coef)"]) if "se(coef)" in s2c.columns else float(
    (np.log(hr_c) - np.log(float(s2c.loc["bnp_iqr", "exp(coef) lower 95%"]))) / 1.96)
t_qba = qba_table(hr_c, se_c, "HRS 2c：达标组 per-IQR BNP→自报新发 CVD（+cysc）")
results.append(("HRS_QBA表", t_qba.to_csv(index=False)))

# 2k. CVD 死亡敏感性（+cysc）
log("\n================= 2k. CVD 死亡敏感性（Exit 死因 121-129） =================")
def load_exit(path, var):
    df, _ = pyreadstat.read_dta(path)
    df["hhidpn"] = df["HHID"].astype(int) * 1000 + df["PN"].astype(int)
    return df[["hhidpn", var]].rename(columns={var: "cod"})

exit_parts = []
for path, var, tag in [
    (HRS_BASE + r"\Raw_data\2018 HRS\2018 HRS Exit\x18exit\X18A_R.dta", "XQA133M1M", "x18"),
    (HRS_BASE + r"\Raw_data\2020 HRS\2020 HRS Exit\x20exit\X20A_R.dta", "XRA133M1M", "x20"),
    (HRS_BASE + r"\Raw_data\2022 HRS\2022 HRS Core\h22core\x22exit\X22A_R.dta", "XSA133M1M", "x22")]:
    try:
        e = load_exit(path, var)
        n_valid = e["cod"].notna().sum()
        log(f"  {tag}: N={len(e)}, 死因非缺失={n_valid}")
        if n_valid > 0:
            exit_parts.append(e)
    except Exception as ex:
        log(f"  {tag}: 读取失败 {str(ex)[:60]}")

if exit_parts:
    cause = pd.concat(exit_parts).drop_duplicates("hhidpn").set_index("hhidpn")["cod"]
    d2 = mm.copy()
    d2["cod"] = cause.reindex(d2.index)
    deaths = d2[d2["death"] == 1]
    cov_cause = deaths["cod"].notna().sum()
    cvd = deaths["cod"].between(121, 129).sum()
    log(f"\n死亡者中死因可及 {cov_cause}/{len(deaths)}；其中 CVD 死因 {cvd} 例")
    d2["t_cvdd"] = np.minimum(d2["t_death"], 4.25)
    d2["cvd_death"] = ((d2["death"] == 1) & d2["cod"].between(121, 129)
                       & (d2["EXDEATHYR"] <= 2020)).astype(int)
    g3 = d2[d2["gly"] == 1].copy()
    g3["bnp_iqr"] = g3["lnbnp"] / (mm["lnbnp"].quantile(0.75) - mm["lnbnp"].quantile(0.25))
    _, s2k = run_cox(g3, ["bnp_iqr"] + COV_H, "cvd_death", "t_cvdd",
                     label="2k. 达标组 per-IQR BNP→CVD死亡（cause-specific, +cysc）")
    results.append(("HRS_2k_CVD死亡_肾功能调整", s2k.loc["bnp_iqr", ["exp(coef)", "exp(coef) lower 95%",
                    "exp(coef) upper 95%", "p"]].round(3).to_dict()))
    log(f"达标组 CVD 死亡事件数：{int(g3['cvd_death'].sum())} / N={len(g3)}")
    results.append(("HRS_2k_事件数", {"N": len(g3), "cvd_death": int(g3['cvd_death'].sum())}))
else:
    log("Exit 死因数据不可用，跳过")

# ============================================================
# PART 3 —— CHARLS 臂（v4：主模型不变 + 两项正式敏感性）
# ============================================================
log("\n================= CHARLS 臂：验证分析（v4） =================")
ch, _ = pyreadstat.read_dta(CH_BASE + r"\Working_data\charls.dta")
w3 = ch[ch["wave"] == 3][["ID", "ragender", "age", "bmi", "raeducl", "hibpe", "diabe",
                          "hearte", "stroke", "bl_crp", "bl_hbalc", "bl_cysc",
                          "bl_tg", "bl_hdl", "cancre",
                          "communityID", "bloodweight",
                          "smoken", "drinkl"]].copy()
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

log(f"CHARLS 分析样本（基线无CVD, CRP+HbA1c齐备）N={len(cm)}，新发CVD={int(cm['new_cvd'].sum())}")
tabc = cm.groupby("gly").agg(N=("new_cvd", "size"), 新发CVD=("new_cvd", "sum")).rename(index=GLY)
log("\n" + tabc.to_string())
results.append(("CHARLS_分层事件", tabc.to_csv()))

# v4：主模型协变量集维持 7 项不变（DAG 最小充分调整集）
COV_C = ["age", "female", "bmi", "raeducl", "hibpe", "smoke_now", "drink_now"]
log(f"\n【CHARLS v4 主模型协变量集】维持 {len(COV_C)} 项不变；"
    f"肾功能/癌症与血脂为 Supplement 敏感性")

goalc = cm[cm["gly"] == 1].copy()
goalc["crp_t"] = pd.qcut(goalc["lncrp"], 3, labels=["T1", "T2", "T3"])
log("\nCHARLS DM达标组 CRP 三分位：")
log(goalc.groupby("crp_t", observed=True).agg(N=("new_cvd", "size"), 新发CVD=("new_cvd", "sum")).to_string())

# 3a. 主模型 + 聚类稳健 SE
_, s3a = run_cox(goalc, ["crp_iqr"] + COV_C, "new_cvd", "t_cvd",
                 cluster="communityID",
                 label="3a. CHARLS DM达标组 per-IQR lnCRP→新发CVD（主模型不变, 聚类SE）【v4 主结果】")
results.append(("CHARLS_3a_主效应", s3a.loc["crp_iqr", ["exp(coef)", "exp(coef) lower 95%",
                "exp(coef) upper 95%", "p"]].round(3).to_dict()))

# 3a-sens. 加血检权重 + 聚类
sub_w = goalc[goalc["bloodweight"].notna() & (goalc["bloodweight"] > 0)].copy()
run_cox(sub_w, ["crp_iqr"] + COV_C, "new_cvd", "t_cvd",
         weights="bloodweight", cluster="communityID",
         label="3a-sens. 加权(bloodweight)+聚类")

# ---- v4 新增：两项正式 Supplement 敏感性 ----
# 3s1. + 肾功能 + 癌症
_, s3s1 = run_cox(goalc, ["crp_iqr", "lncysc_c", "cancre"] + COV_C, "new_cvd", "t_cvd",
                  cluster="communityID",
                  label="3s1. Supplement敏感性：+ ln(cysc) + cancre（肾功能+癌症）")
results.append(("CHARLS_3s1_肾癌症敏感性", s3s1.loc["crp_iqr", ["exp(coef)", "exp(coef) lower 95%",
                "exp(coef) upper 95%", "p"]].round(3).to_dict()))

# 3s2. + 血脂（TG + HDL）
_, s3s2 = run_cox(goalc, ["crp_iqr", "lntg_c", "lnhdl_c"] + COV_C, "new_cvd", "t_cvd",
                  cluster="communityID",
                  label="3s2. Supplement敏感性：+ ln(TG) + ln(HDL)（血脂轴）")
results.append(("CHARLS_3s2_血脂敏感性", s3s2.loc["crp_iqr", ["exp(coef)", "exp(coef) lower 95%",
                "exp(coef) upper 95%", "p"]].round(3).to_dict()))

# 3b. 全样本交互
fullc = cm.copy()
fullc = pd.get_dummies(fullc, columns=["gly"], prefix="gly", drop_first=True, dtype=int)
gcols_c = [x for x in fullc.columns if x.startswith("gly_")]
inter_c = [f"{g}_x" for g in gcols_c]
for g in gcols_c:
    fullc[f"{g}_x"] = fullc[g] * fullc["crp_iqr"]
mf, _ = run_cox(fullc, ["crp_iqr"] + gcols_c + COV_C + inter_c, "new_cvd", "t_cvd",
                cluster="communityID", label="3b. CHARLS 全样本 gly×CRP 交互（聚类SE）")
mb, _ = run_cox(fullc, ["crp_iqr"] + gcols_c + COV_C, "new_cvd", "t_cvd",
                cluster="communityID", label="3b-nested")
lrc = 2 * (mf.log_likelihood_ - mb.log_likelihood_)
log(f"\nCHARLS 嵌套LR: LR={lrc:.2f}, df={len(inter_c)}, p={chi2.sf(lrc, len(inter_c)):.4f}")
results.append(("CHARLS_交互LR", f"LR={lrc:.2f}, p={chi2.sf(lrc, len(inter_c)):.5f}"))

# 3c. 各层
log("\n--- 3c. CHARLS 各血糖层 per-IQR CRP→新发CVD（聚类SE） ---")
for g in sorted(cm["gly"].unique()):
    sub = cm[cm["gly"] == g].copy()
    try:
        cph, s = run_cox(sub, ["crp_iqr"] + COV_C, "new_cvd", "t_cvd",
                         cluster="communityID", label=f"  {GLY.get(g, g)}")
        results.append((f"CHARLS_HR_{GLY.get(g, g)}",
                        s.loc["crp_iqr", ["exp(coef)", "p"]].round(3).to_dict()))
    except Exception as e:
        log(f"  {GLY.get(g, g)}: 拟合失败 {e}")

# 3d. CHARLS QBA（主模型口径）
hr_cc = float(s3a.loc["crp_iqr", "exp(coef)"])
se_cc = (np.log(hr_cc) - np.log(float(s3a.loc["crp_iqr", "exp(coef) lower 95%"]))) / 1.96
t_qba_c = qba_table(hr_cc, se_cc, "CHARLS 3a：达标组 per-IQR CRP→自报新发 CVD")
results.append(("CHARLS_QBA表", t_qba_c.to_csv(index=False)))

# ============================================================
# PART 4 —— 保存
# ============================================================
with open(OUT + r"\main_analysis_v4_results.txt", "w", encoding="utf-8") as f:
    f.write("CD特刊主分析 v4（协变量科学处置终版）\n"
            "HRS 主模型全程 +ln(cystatin C)；CHARLS 主模型不变 + 两项 Supplement 敏感性\n"
            "========================\n\n")
    for k, v in results:
        f.write(f"[{k}]\n{v}\n\n")
log("\n完成。结果摘要已存 main_analysis_v4_results.txt")
