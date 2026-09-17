# -*- coding: utf-8 -*-
"""
CD 特刊 MR 三角臂：遗传预测 NT-proBNP → 心血管结局（OpenGWAS API）
暴露：
  主：ebi-a-GCST90012082（NT-proBNP，N=21,758，Folkersen 2020）——rs198389（NPPA/NPPB 簇）
  敏感性：prot-a-2078（NT-proBNP，N=3,301，Sun 2018）——rs198389 + rs2175388（双 SNP IVW）
结局（10）：
  发现集：ieu-a-7 CAD｜ieu-a-798 MI｜ebi-a-GCST009541 HF(HERMES)｜
          ebi-a-GCST005838 卒中(MEGASTROKE)｜ebi-a-GCST006908 缺血性卒中
  复制集：finn-b-I9_CHD｜finn-b-I9_MI｜finn-b-I9_HEARTFAIL｜finn-b-I9_STR_SAH｜finn-b-C_STROKE
方法：Wald ratio（单 SNP）+ IVW（双 SNP）+ Cochran Q 异质性
数据源：本地缓存 mr_instruments.json / mr_outcome_assoc.json（OpenGWAS API 抓取，可复现）
"""
import json
import numpy as np
from scipy.stats import norm, chi2

DIR = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"

# ---------- 暴露工具变量 ----------
# 主暴露 GCST90012082（N=21,758）
exp_main = {"rs198389": {"beta": 0.2122, "se": 0.0158, "ea": "G", "nea": "A", "n": 21758}}
# 敏感性暴露 prot-a-2078（N=3,301，Sun 2018 plasma proteome）
exp_sens = {
    "rs198389":   {"beta": 0.2934, "se": 0.0253, "ea": "G", "nea": "A", "n": 3301},
    "rs2175388":  {"beta": 0.1830, "se": 0.0330, "ea": "T", "nea": "C", "n": 3301},
}

# ---------- 结局 ----------
OUTCOMES = [
    ("ieu-a-7",            "CAD (CARDIoGRAMplusC4D)",      "发现"),
    ("ieu-a-798",          "心肌梗死 (UKB+C4D)",           "发现"),
    ("ebi-a-GCST009541",   "心衰 (HERMES)",                "发现"),
    ("ebi-a-GCST005838",   "全卒中 (MEGASTROKE)",          "发现"),
    ("ebi-a-GCST006908",   "缺血性卒中 (MEGASTROKE)",      "发现"),
    ("finn-b-I9_CHD",      "冠心病 (FinnGen R9)",          "复制"),
    ("finn-b-I9_MI",       "心肌梗死 (FinnGen R9)",        "复制"),
    ("finn-b-I9_HEARTFAIL","心衰 (FinnGen R9)",            "复制"),
    ("finn-b-I9_STR_SAH",  "卒中 (FinnGen R9)",            "复制"),
    ("finn-b-C_STROKE",    "卒中 (FinnGen R9 广义)",       "复制"),
]

# ---------- 载入结局关联 ----------
assoc = json.load(open(DIR + r"\mr_outcome_assoc.json", encoding="utf-8"))
assoc_map = {(a["rsid"], a["id"]): a for a in assoc}

# ---------- 工具变量强度 ----------
print("==================== 工具变量强度 ====================")
for k, v in exp_main.items():
    F = (v["beta"] / v["se"]) ** 2
    print(f"主暴露 {k}: beta={v['beta']}, se={v['se']} → F = {F:.1f}")
for k, v in exp_sens.items():
    F = (v["beta"] / v["se"]) ** 2
    print(f"敏感性暴露 {k}: beta={v['beta']}, se={v['se']} → F = {F:.1f}")

def wald(beta_exp, se_exp, beta_out, se_out):
    """Wald ratio + delta-method SE"""
    r = beta_out / beta_exp
    se = se_out / abs(beta_exp)
    z = r / se
    return r, se, 2 * norm.sf(abs(z))

def ivw(snps, beta_exp, se_out_map):
    """固定效应 IVW（多元加权回归过原点）"""
    X = np.array([beta_exp[s] for s in snps])
    Y = np.array([beta_out[s] for s in snps])
    W = np.array([1.0 / se_out_map[s] ** 2 for s in snps])
    # 加权最小二乘过原点
    beta_ivw = np.sum(W * X * Y) / np.sum(W * X * X)
    se_ivw = np.sqrt(1.0 / np.sum(W * X * X))
    z = beta_ivw / se_ivw
    # Cochran Q
    Q = np.sum(W * (Y - beta_ivw * X) ** 2)
    dof = len(snps) - 1
    Qp = chi2.sf(Q, dof) if dof > 0 else np.nan
    return beta_ivw, se_ivw, 2 * norm.sf(abs(z)), Q, Qp

print("\n==================== 主分析：Wald ratio（rs198389, GCST90012082 主暴露） ====================")
rows = []
beta_exp = {s: v["beta"] for s, v in exp_main.items()}
se_exp = {s: v["se"] for s, v in exp_main.items()}
for oid, oname, tier in OUTCOMES:
    a = assoc_map.get(("rs198389", oid))
    if a is None or a.get("beta") is None:
        print(f"  {oname}: rs198389 结局关联不可得，跳过")
        continue
    # 效应等位基因协调
    b_out, se_out = a["beta"], a["se"]
    sign = 1.0
    if a["ea"] != exp_main["rs198389"]["ea"]:
        if a["ea"] == exp_main["rs198389"]["nea"] and a["nea"] == exp_main["rs198389"]["ea"]:
            sign = -1.0
        else:
            print(f"  {oname}: 等位基因不匹配（{a['ea']}/{a['nea']}），跳过")
            continue
    b_out *= sign
    r, se, p = wald(beta_exp["rs198389"], se_exp["rs198389"], b_out, se_out)
    lo, hi = np.exp(r - 1.96 * se), np.exp(r + 1.96 * se)
    flag = "*" if p < 0.05 else " "
    print(f"  {flag}{oname:28s} OR={np.exp(r):6.3f} ({lo:.3f}-{hi:.3f})  p={p:.4g}   [结局β={b_out:+.4f}, SE={se_out}]")
    rows.append({"outcome": oname, "id": oid, "tier": tier, "or": np.exp(r),
                 "lo": lo, "hi": hi, "p": p, "b_out": b_out, "se_out": se_out, "se_exp": se_out})

# ---------- 敏感性：双 SNP IVW（prot-a-2078） ----------
print("\n==================== 敏感性：双 SNP IVW（rs198389+rs2175388, Sun 2018） ====================")
beta_exp2 = {s: v["beta"] for s, v in exp_sens.items()}
se_exp2 = {s: v["se"] for s, v in exp_sens.items()}
for oid, oname, tier in OUTCOMES:
    snps, beta_out_d, se_out_d = [], {}, {}
    for s in exp_sens:
        a = assoc_map.get((s, oid))
        if a is None or a.get("beta") is None:
            continue
        sign = 1.0
        if a["ea"] != exp_sens[s]["ea"]:
            if a["ea"] == exp_sens[s]["nea"] and a["nea"] == exp_sens[s]["ea"]:
                sign = -1.0
            else:
                continue
        snps.append(s)
        beta_out_d[s] = sign * a["beta"]
        se_out_d[s] = a["se"]
    beta_out = beta_out_d  # 供 ivw 使用
    if len(snps) < 1:
        print(f"  {oname}: 无可用 SNP")
        continue
    if len(snps) == 1:
        r, se, p = wald(beta_exp2[snps[0]], se_exp2[snps[0]], beta_out_d[snps[0]], se_out_d[snps[0]])
        print(f"  {oname:28s} Wald(1SNP) OR={np.exp(r):6.3f} ({np.exp(r-1.96*se):.3f}-{np.exp(r+1.96*se):.3f}) p={p:.4g}")
    else:
        b, se, p, Q, Qp = ivw(snps, beta_exp2, se_out_d)
        print(f"  {oname:28s} IVW OR={np.exp(b):6.3f} ({np.exp(b-1.96*se):.3f}-{np.exp(b+1.96*se):.3f}) p={p:.4g}  Q={Q:.2f}(df={len(snps)-1}, p={Qp:.3f})")

# ---------- 保存 ----------
import csv
with open(DIR + r"\mr_results_main.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["outcome", "id", "tier", "or", "lo", "hi", "p"])
    w.writeheader()
    for r in rows:
        w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()
                    if k in ("outcome", "id", "tier", "or", "lo", "hi", "p")})
print("\n主分析表已存 mr_results_main.csv")
