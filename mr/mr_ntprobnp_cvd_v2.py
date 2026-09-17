# -*- coding: utf-8 -*-
"""
NT-proBNP/BNP → CVD 结局 MR 分析 v2：扩展工具变量
================================================
工具变量扩展策略（对 v1 单 IV 的升级）：
  Set A（主, SCALLOP NT-proBNP, N=21,758, ebi-a-GCST90012082）：
      rs198389 (NPPB 5', p=4.4e-41, F=180.4) + rs5068 (NPPA 3'UTR 功能位点,
      p=1.75e-5, F=18.5) —— 两两 LD r2=0.083 (1000G EUR, Ensembl)，独立。
      rs5068 为文献已知功能变异（miRNA 结合位点，升高利钠肽），按先验知识入选。
  Set B（复制暴露, Sun 2018 NT-proBNP, N=3,301, prot-a-2078）：
      rs198389 + rs2175388 (p=2.8e-8, OpenGWAS clump r2<0.001 判独立)
  Set C（BNP 暴露, Folkersen 2018, N=3,394, prot-b-80）：
      rs198379 (cis NPPB, p=9.6e-25) + rs6557662 (trans chr8, p=1.5e-08)
      —— 含一个远端 trans 工具，可分离顺式区域多效性。

已知 LD（Ensembl, 1000G EUR）: rs198389-rs198379 r2=0.984; rs198389-rs549596 r2=0.898;
  rs198379-rs632793 r2=0.872（同一 LD 块，各集合只取一个代表）；rs198389-rs5068 r2=0.083。

方法：Wald ratio（单 SNP）/ 固定效应 IVW（2 SNP）+ Cochran Q；回文 SNP 用 EAF 校验。
数据：OpenGWAS API 真实关联（已落盘 JSON）。
"""
import json, math
import numpy as np
import pandas as pd

DIR = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"

# ---------- 载入已落盘的关联数据 ----------
cand = json.load(open(DIR + r"\mr_cand_snp_assoc.json", encoding="utf-8"))
newiv = json.load(open(DIR + r"\mr_new_iv_outcome_assoc.json", encoding="utf-8"))
oldiv = json.load(open(DIR + r"\mr_outcome_assoc.json", encoding="utf-8"))
assoc_all = cand + newiv + oldiv
amap = {}
for a in assoc_all:
    if a.get("beta") is not None and a.get("se") is not None:
        amap[(a["rsid"], a["id"])] = a

# ---------- 暴露侧工具变量 ----------
EXPOSURES = {
    "SetA_SCALLOP_NTproBNP_N21758": {
        "ebi-a-GCST90012082": {
            "rs198389": {"beta": 0.2122, "se": 0.0158, "ea": "G", "nea": "A", "n": 21758, "eaf": None},
            "rs5068":   {"beta": 0.1452, "se": 0.0338, "ea": "G", "nea": "A", "n": 21758, "eaf": None},
        }
    },
    "SetB_Sun2018_NTproBNP_N3301": {
        "prot-a-2078": {
            "rs198389":  {"beta": 0.2934, "se": 0.0253, "ea": "G", "nea": "A", "n": 3301},
            "rs2175388": {"beta": 0.1830, "se": 0.0330, "ea": "T", "nea": "C", "n": 3301},
        }
    },
    "SetC_Folkersen2018_BNP_N3394": {
        "prot-b-80": {
            "rs198379":  {"beta": 0.2560, "se": 0.0247, "ea": "C", "nea": "T", "n": 3394},
            "rs6557662": {"beta": -0.2320, "se": 0.0409, "ea": "G", "nea": "A", "n": 3394},
        }
    },
}

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

def complement(a):
    return {"A": "T", "T": "A", "G": "C", "C": "G"}.get(a, a)

def harmonize(exp_snp, out_assoc):
    """结局关联对齐到暴露效应等位基因。返回 (beta, se, ok, note)"""
    if out_assoc is None:
        return None, None, False, "结局无该 SNP"
    ea_o, nea_o = out_assoc["ea"], out_assoc["nea"]
    ea_e, nea_e = exp_snp["ea"], exp_snp["nea"]
    beta_o, se_o = out_assoc["beta"], out_assoc["se"]
    # 直接匹配
    if {ea_o, nea_o} == {ea_e, nea_e}:
        if ea_o == ea_e:
            return beta_o, se_o, True, ""
        else:
            return -beta_o, se_o, True, "翻转"
    # 互补匹配
    if {complement(ea_o), complement(nea_o)} == {ea_e, nea_e}:
        if complement(ea_o) == ea_e:
            return beta_o, se_o, True, "互补同向"
        else:
            return -beta_o, se_o, True, "互补翻转"
    return None, None, False, f"等位基因不匹配 {ea_o}/{nea_o}"

def palindrome_check(exp_snp, out_assoc, tol=0.02):
    """回文 SNP：比较暴露与结局 EAF 是否同侧（都<0.5 或都>0.5）。"""
    ea_e, ea_o = exp_snp["ea"], out_assoc.get("ea")
    eaf_o = out_assoc.get("eaf")
    eaf_e = exp_snp.get("eaf")
    # 暴露 eaf 优先用 API 返回（cand json 有）；否则跳过检查（返回 True）
    if eaf_o is None or eaf_e is None:
        return True, "eaf 缺失,未校验"
    # 对齐效应等位基因后再比较
    if ea_o != ea_e:
        eaf_o = 1 - eaf_o
    same_side = (eaf_e < 0.5) == (eaf_o < 0.5)
    diff = abs(eaf_e - eaf_o)
    if same_side and diff < 0.5 - tol:
        return True, f"回文 eaf 一致 ({eaf_e:.3f} vs {eaf_o:.3f})"
    return False, f"回文 eaf 冲突 ({eaf_e:.3f} vs {eaf_o:.3f}) → 剔除"

# 从 cand json 提取暴露 eaf（若可得）
for (rsid, gwasid), a in amap.items():
    for setname, d in EXPOSURES.items():
        for gid, snps in d.items():
            if gid == gwasid and rsid in snps and snps[rsid].get("eaf") in (None,):
                snps[rsid]["eaf"] = a.get("eaf")

def wald(bx, sx, by, sy):
    r = by / bx
    se = sy / abs(bx)
    z = r / se
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return r, se, p

def ivw2(pairs):
    """固定效应 IVW（2 SNP, 一阶矩权重, 过原点）+ Cochran Q"""
    bxs = np.array([p[0] for p in pairs]); sxs = np.array([p[1] for p in pairs])
    bys = np.array([p[2] for p in pairs]); sys_ = np.array([p[3] for p in pairs])
    # 比率及其 SE
    ratios = bys / bxs
    se_ratios = sys_ / np.abs(bxs)
    w = 1 / se_ratios ** 2
    b = np.sum(w * ratios) / np.sum(w)
    se = math.sqrt(1 / np.sum(w))
    z = b / se
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    # Cochran Q
    Q = np.sum(w * (ratios - b) ** 2)
    from scipy.stats import chi2
    Qp = 1 - chi2.cdf(Q, df=len(pairs) - 1)
    return b, se, p, Q, Qp

# ---------- 主分析 ----------
rows = []
print("=" * 110)
print("扩展工具变量 MR：NT-proBNP/BNP → 10 心血管结局（3 工具集 × Wald/IVW）")
print("=" * 110)
for oid, oname, tier in OUTCOMES:
    print(f"\n--- {oname} [{oid}, {tier}] ---")
    for setname, d in EXPOSURES.items():
        gid = list(d.keys())[0]; snps = d[gid]
        pairs, singles = [], []
        notes = []
        for rsid, esnp in snps.items():
            oa = amap.get((rsid, oid))
            if oa is None:
                notes.append(f"{rsid}:结局缺失"); continue
            by, sy, ok, hnote = harmonize(esnp, oa)
            if not ok:
                notes.append(f"{rsid}:{hnote}"); continue
            # 回文校验
            palin = set(esnp["ea"]) == {complement(esnp["ea"])} or \
                    (esnp["ea"] == complement(esnp["nea"]) and esnp["nea"] == complement(esnp["ea"]))
            if esnp["ea"] in "AT" and esnp["nea"] in "AT" and esnp["ea"] != esnp["nea"] and set(esnp["ea"] + esnp["nea"]) <= {"A", "T"}:
                okp, pnote = palindrome_check(esnp, oa)
                if not okp:
                    notes.append(f"{rsid}:{pnote}"); continue
            b, s, p = wald(esnp["beta"], esnp["se"], by, sy)
            singles.append((rsid, b, s, p))
            pairs.append((esnp["beta"], esnp["se"], by, sy))
        if not pairs:
            print(f"  [{setname}] 无可用工具 ({'; '.join(notes)})"); continue
        for rsid, b, s, p in singles:
            rows.append(dict(outcome=oname, outcome_id=oid, tier=tier, set=setname, method=f"Wald({rsid})",
                             snp=rsid, logOR=b, SE=s, OR=math.exp(b), lo=math.exp(b - 1.96 * s), hi=math.exp(b + 1.96 * s), p=p, Q=np.nan, Qp=np.nan))
        if len(pairs) >= 2:
            b, s, p, Q, Qp = ivw2(pairs)
            rows.append(dict(outcome=oname, outcome_id=oid, tier=tier, set=setname, method="IVW(2SNP)",
                             snp="+".join(snps.keys()), logOR=b, SE=s, OR=math.exp(b), lo=math.exp(b - 1.96 * s), hi=math.exp(b + 1.96 * s), p=p, Q=Q, Qp=Qp))
            print(f"  [{setname:32s}] IVW  OR={math.exp(b):6.3f} ({math.exp(b-1.96*s):.3f}-{math.exp(b+1.96*s):.3f}) p={p:.4g} | Q={Q:.2f} Qp={Qp:.3f}")
        for rsid, b, s, p in singles:
            print(f"    · Wald {rsid:10s} OR={math.exp(b):6.3f} ({math.exp(b-1.96*s):.3f}-{math.exp(b+1.96*s):.3f}) p={p:.4g}")
        if notes:
            print(f"    剔除/警告: {'; '.join(notes)}")

df = pd.DataFrame(rows)
df.to_csv(DIR + r"\mr_v2_results_full.csv", index=False, encoding="utf-8-sig")
print("\n" + "=" * 110)
print("关键结局汇总（全卒中）")
print("=" * 110)
sub = df[df["outcome"].str.contains("卒中|STROKE", case=False)]
print(sub[["outcome", "set", "method", "OR", "lo", "hi", "p", "Qp"]].to_string(index=False))
print("\n结果已保存: mr_v2_results_full.csv")
