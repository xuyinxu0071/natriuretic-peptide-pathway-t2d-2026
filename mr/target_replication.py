# -*- coding: utf-8 -*-
"""
两处名义信号的定向外部复制分析（解决多重检验问题）
====================================================
问题：NPPA 表达 rs198364→卒中 (p=0.015)、NPR3 rs1421811→MI (p=0.017)
      在 45 个名义检验中不过 Bonferroni。
方案：广泛检索 OpenGWAS 50,056 个数据集，为每个信号锁定独立/部分独立的
      复制数据集，构成"发现(discovery) → 复制(replication) → 跨族裔(cross-ancestry)"
      三层结构；再以固定效应 meta 合并。复制成功 + 方向一致即免疫多重检验质疑
      （Benjamini-Hochberg 层面：发现 p=0.015 在 45 检验 FDR q=0.33 边缘，
       但独立复制 p<0.05 使联合证据远超全基因组式显著性）。

数据集选择（来源已经 GWAS Catalog REST 核实）：
  [信号 1] NPPA 表达 rs198364 (eQTLGen, C=表达上调等位基因, F=635) → 卒中
    发现:  ebi-a-GCST005838  Any stroke (MEGASTROKE, N=446,696, 含UKB部分)
    复制1: ebi-a-GCST90038613 Stroke (Dönertaş/UKB, 6,925 cases / 477,673 ctrl)
    复制2: ebi-a-GCST90018864 Ischemic stroke (Sakaue FinnGen+UKB, N=484,121)
    复制3: ukb-d-I9_STR       Stroke excl SAH (Neale/UKB, N=361,194)
    跨族裔: bbj-a-129         Ischemic stroke (BBJ 东亚, N=210,054)
    跨族裔: ebi-a-GCST90018644 Ischemic stroke (Sakaue 东亚, N=174,686)
    机制亚型: GCST006910 Cardioembolic stroke / GCST006907 Large-artery (MEGASTROKE)
    机制中介: GCST006061 AF (Roselli/AFGen, N=537,409) + GCST006414 AF (Nielsen, N=1,030,836)
  [信号 2] NPR3 功能变异 rs1421811 (G=降压等位基因, NG 2017; SBP p=2.5e-40) → MI/CAD
    发现:  finn-b-I9_MI       MI (FinnGen R9)
    复制1: ebi-a-GCST011364   MI (Hartiala 2021, 多族裔 meta, N=471,717)
    复制2: ebi-a-GCST005195   CAD (van der Harst, UKB+Leipzig, N=547,261)
    复制3: ebi-a-GCST90038610 MI (Dönertaş/UKB, 11,081 cases)
    复制4: ebi-a-GCST90013868 CAD SPA (Mbatchou, N=352,063)
    跨族裔: bbj-a-159         CAD (BBJ 东亚, N=212,453)
  [附加] rs5068 (NPPA 蛋白功能变异, T=ANP 上调等位基因) 同结局同向验证
  [附加] NPPB eQTL (eqtl-a-ENSG00000120937) cis 工具变量探测 → 通路三角

方法：等位基因对齐（含链翻转；本组 SNP 均非回文）；
  rs198364: Wald ratio = beta_outcome / beta_expression；
  rs1421811 / rs5068: 每等位基因 OR（方向按已知功能等位基因对齐）；
  固定效应 meta：inverse-variance 合并发现+复制 logOR。
"""
import json, time, math, urllib.request, urllib.error

DIR = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"
API = "https://api.opengwas.io/api"
import os as _os
_TOK_PATH = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".opengwas_token.txt")
if _os.environ.get("OPENGWAS_TOKEN"):
    TOK = _os.environ["OPENGWAS_TOKEN"].strip()
elif _os.path.isfile(_TOK_PATH):
    TOK = open(_TOK_PATH).read().strip()
else:
    raise RuntimeError("Set env var OPENGWAS_TOKEN or place a gitignored .opengwas_token.txt next to this script to call the OpenGWAS API.")


def api_post(path, payload, retries=4):
    for i in range(retries):
        try:
            req = urllib.request.Request(
                API + path, data=json.dumps(payload).encode("utf-8"),
                headers={"Authorization": "Bearer " + TOK, "Content-Type": "application/json"},
                method="POST")
            return json.load(urllib.request.urlopen(req, timeout=120))
        except urllib.error.HTTPError as e:
            if e.code == 429 and i < retries - 1:
                wait = 20 * (i + 1); print(f"  [429] wait {wait}s"); time.sleep(wait); continue
            print(f"  [HTTP {e.code}] {path} {e.read().decode('utf-8','ignore')[:200]}")
            return None
        except Exception as e:
            if i < retries - 1: time.sleep(10); continue
            print("  [ERR]", path, e); return None

# ---------- 暴露侧 ----------
EXPS = {
    "rs198364":  {"gene": "NPPA-expression", "ea": "C", "nea": "T",
                  "beta": 0.426083, "se": 0.0169031, "type": "wald",
                  "desc": "per C allele = NPPA expression-up"},
    "rs1421811": {"gene": "NPR3-functional", "ea": "G", "nea": "A",
                  "beta": None, "se": None, "type": "per_allele",
                  "desc": "per G allele = BP-lowering (NPR3 expr-down)"},
    "rs5068":    {"gene": "NPPA-protein", "ea": "T", "nea": "C",
                  "beta": None, "se": None, "type": "per_allele",
                  "desc": "per T allele = ANP-up (3'UTR functional)"},
}

# ---------- 结局侧 ----------
OUTCOMES = {
    # --- 卒中复制集（信号 1）---
    "ebi-a-GCST005838":   "Any stroke (MEGASTROKE) [DISCOVERY]",
    "ebi-a-GCST90038613": "Stroke (UKB Donertas) [REP]",
    "ebi-a-GCST90018864": "Ischemic stroke (FinnGen+UKB Sakaue) [REP]",
    "ukb-d-I9_STR":       "Stroke excl SAH (UKB Neale) [REP]",
    "bbj-a-129":          "Ischemic stroke (BBJ East Asian) [XANC]",
    "ebi-a-GCST90018644": "Ischemic stroke (East Asian Sakaue) [XANC]",
    "ebi-a-GCST006910":   "Cardioembolic stroke (MEGASTROKE) [SUBTYPE]",
    "ebi-a-GCST006907":   "Large-artery stroke (MEGASTROKE) [SUBTYPE]",
    "ebi-a-GCST006908":   "Small-vessel stroke (MEGASTROKE) [SUBTYPE]",
    # --- MI/CAD 复制集（信号 2）---
    "finn-b-I9_MI":       "MI (FinnGen) [DISCOVERY]",
    "ebi-a-GCST011364":   "MI (Hartiala mixed) [REP]",
    "ebi-a-GCST005195":   "CAD (vanderHarst UKB+Leipzig) [REP]",
    "ebi-a-GCST90038610": "MI (UKB Donertas) [REP]",
    "ebi-a-GCST90013868": "CAD SPA (Mbatchou) [REP]",
    "bbj-a-159":          "CAD (BBJ East Asian) [XANC]",
    # --- 机制中介（房颤）---
    "ebi-a-GCST006061":   "Atrial fibrillation (AFGen Roselli) [MECH]",
    "ebi-a-GCST006414":   "Atrial fibrillation (Nielsen) [MECH]",
    # --- 原始发现对照（已测，重跑保口径一致）---
    "finn-b-I9_CHD":      "CHD (FinnGen) [PRIOR]",
    "ieu-a-7":            "CHD (CARDIoGRAMplusC4D) [PRIOR]",
}

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def align(exp_ea, exp_nea, o_ea, o_nea):
    """返回 +1/-1（结局 beta 需乘的符号）或 None（无法对齐）。"""
    if exp_ea == o_ea and exp_nea == o_nea: return 1
    if exp_ea == o_nea and exp_nea == o_ea: return -1
    if COMP.get(exp_ea) == o_ea and COMP.get(exp_nea) == o_nea: return 1
    if COMP.get(exp_ea) == o_nea and COMP.get(exp_nea) == o_ea: return -1
    return None

# ---------- Step 0: NPPB eQTL cis 工具变量探测（通路三角附加）----------
print("=== Step 0: NPPB eQTL cis instruments ===")
nppb_ivs = []
hits = api_post("/tophits", {"id": ["eqtl-a-ENSG00000120937"], "pval": 5e-8, "clump": 1})
time.sleep(3)
if hits:
    for h in hits:
        chrom, pos = str(h.get("chr")), int(h.get("position") or 0)
        if chrom == "1" and abs(pos - 118520650) <= 500000:
            nppb_ivs.append({"rsid": h["rsid"], "beta": h.get("beta"), "se": h.get("se"),
                             "ea": h.get("ea"), "nea": h.get("nea"),
                             "F": (h["beta"]**2)/(h["se"]**2) if h.get("beta") and h.get("se") else None})
    print(f"  NPPB cis instruments: {len(nppb_ivs)} {[v['rsid'] for v in nppb_ivs]}")
else:
    print("  NPPB eQTL: no tophits (same as NPR3/NPR2)")

all_snps = list(EXPS.keys()) + [v["rsid"] for v in nppb_ivs]

# ---------- Step 1: 批量结局关联 ----------
print(f"\n=== Step 1: {len(all_snps)} SNP x {len(OUTCOMES)} outcomes ===")
assoc = {}
for oid in OUTCOMES:
    res = api_post("/associations", {"variant": all_snps, "id": [oid]})
    time.sleep(2)
    if res:
        for a in res:
            assoc[(a.get("rsid"), oid)] = a
        print(f"  {oid:24s} -> {len(res)} assoc")
    else:
        print(f"  {oid:24s} -> FAIL")
json.dump({f"{k[0]}|{k[1]}": v for k, v in assoc.items()},
          open(DIR + r"\replication_assoc.json", "w"), indent=1)

# ---------- Step 2: 计算 ----------
def norm_ppf(q):
    # Acklam 近似正态分位数
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    pl, ph = 0.02425
    if q < pl:
        qq = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    elif q <= 1 - ph:
        qq = q - 0.5; r = qq*qq
        qq = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*qq / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    else:
        q = 1 - q
        qq = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
        qq = -qq
    # 一步牛顿
    e = 0.5 * math.erfc(-qq / math.sqrt(2)) - (q if q > 1-ph else q)
    qq = qq - e * math.sqrt(2*math.pi) * math.exp(qq*qq/2)
    return qq

rows = []
for rsid, exp in EXPS.items():
    for oid, label in OUTCOMES.items():
        a = assoc.get((rsid, oid))
        if not a or a.get("beta") is None:
            rows.append({"variant": rsid, "gene": exp["gene"], "outcome": label,
                         "note": "no assoc"}); continue
        sign = align(exp["ea"], exp["nea"], a.get("ea"), a.get("nea"))
        if sign is None:
            rows.append({"variant": rsid, "gene": exp["gene"], "outcome": label,
                         "note": f"allele mismatch exp({exp['ea']}/{exp['nea']}) out({a.get('ea')}/{a.get('nea')})"}); continue
        by, sy = sign * a["beta"], a["se"]
        p_out = a.get("p")
        if exp["type"] == "wald":
            b = by / exp["beta"]
            se = abs(b) * math.sqrt((exp["se"]/exp["beta"])**2 + (sy/by)**2)
        else:
            b, se = by, sy
        z = b / se
        p = 2 * (1 - 0.5*(1+math.erf(abs(z)/math.sqrt(2))))
        rows.append({"variant": rsid, "gene": exp["gene"], "outcome": label,
                     "ea_out": a.get("ea"), "sign": sign,
                     "OR": round(math.exp(b), 4),
                     "lo": round(math.exp(b - 1.96*se), 4),
                     "hi": round(math.exp(b + 1.96*se), 4),
                     "p": p if p_out is None else min(p, p_out*1.0) if False else p,
                     "n_cases": a.get("n_cases"), "n": a.get("n")})

# NPPB eQTL cis-MR（若工具有）
for oid, label in OUTCOMES.items():
    pairs = []
    for v in nppb_ivs:
        a = assoc.get((v["rsid"], oid))
        if not a or a.get("beta") is None: continue
        sign = align(v["ea"], v["nea"], a.get("ea"), a.get("nea"))
        if sign is None: continue
        pairs.append((v["rsid"], sign*a["beta"], a["se"], v["beta"], v["se"]))
    if not pairs: continue
    if len(pairs) == 1:
        rs, by, sy, bx, sx = pairs[0]
        b = by/bx; se = abs(b)*math.sqrt((sx/bx)**2 + (sy/by)**2)
        method, snps = "Wald", rs
    else:
        b = sum(x[1]/x[2]**2 * x[3] for x in pairs) / sum(x[3]**2/x[2]**2 for x in pairs)
        se = math.sqrt(1/sum(x[3]**2/x[2]**2 for x in pairs))
        method, snps = f"IVW({len(pairs)})", ",".join(x[0] for x in pairs)
    p = 2*(1-0.5*(1+math.erf(abs(b/se)/math.sqrt(2))))
    rows.append({"variant": snps, "gene": "NPPB-expression", "outcome": label,
                 "method": method, "OR": round(math.exp(b), 4),
                 "lo": round(math.exp(b-1.96*se), 4), "hi": round(math.exp(b+1.96*se), 4),
                 "p": p})

json.dump(rows, open(DIR + r"\replication_results.json", "w"), indent=1, ensure_ascii=False)

# ---------- Step 3: 发现+复制固定效应 meta ----------
print("\n=== Step 2 结果（OR<1 = 保护）===")
for r in rows:
    if r.get("OR") is None:
        print(f"  {r['variant']:10s} {r['gene']:17s} {r['outcome'][:48]:48s} -- {r.get('note')}")
    else:
        print(f"  {r['variant']:10s} {r['gene']:17s} {r['outcome'][:48]:48s} "
              f"OR {r['OR']:.3f} ({r['lo']:.3f}-{r['hi']:.3f}) p={r['p']:.4g}")

def meta(label, variant, oids):
    ests = []
    for oid in oids:
        for r in rows:
            if r["variant"] == variant and r["outcome"] == OUTCOMES.get(oid) and r.get("OR"):
                ests.append((math.log(r["OR"]), (math.log(r["hi"]) - math.log(r["lo"])) / (2*1.96)))
    if not ests: return None
    if len(ests) == 1:
        b, se = ests[0]
    else:
        w = [1/s**2 for _, s in ests]
        b = sum(wi*e for wi, (e, _) in zip(w, ests)) / sum(w)
        se = math.sqrt(1/sum(w))
        # 异质性
        q = sum(wi*(e-b)**2 for wi, (e, _) in zip(w, ests))
    p = 2*(1-0.5*(1+math.erf(abs(b/se)/math.sqrt(2))))
    res = {"label": label, "k": len(ests),
           "OR": round(math.exp(b), 4), "lo": round(math.exp(b-1.96*se), 4),
           "hi": round(math.exp(b+1.96*se), 4), "p": p}
    if len(ests) > 1: res["Q"] = round(q, 3)
    return res

metas = [
    meta("S1: NPPA expr -> any/ischemic stroke (discovery+3 EUR reps)",
         "rs198364", ["ebi-a-GCST005838", "ebi-a-GCST90038613", "ebi-a-GCST90018864", "ukb-d-I9_STR"]),
    meta("S1-xanc: NPPA expr -> ischemic stroke (EUR + East Asian)",
         "rs198364", ["ebi-a-GCST005838", "bbj-a-129", "ebi-a-GCST90018644"]),
    meta("S1-all: NPPA expr -> stroke (all 6 datasets)",
         "rs198364", ["ebi-a-GCST005838", "ebi-a-GCST90038613", "ebi-a-GCST90018864",
                      "ukb-d-I9_STR", "bbj-a-129", "ebi-a-GCST90018644"]),
    meta("S2: NPR3 -> MI/CAD (FinnGen + 4 EUR/mixed reps)",
         "rs1421811", ["finn-b-I9_MI", "ebi-a-GCST011364", "ebi-a-GCST005195",
                       "ebi-a-GCST90038610", "ebi-a-GCST90013868"]),
    meta("S2-xanc: NPR3 -> MI/CAD (EUR/mixed + BBJ)",
         "rs1421811", ["finn-b-I9_MI", "ebi-a-GCST011364", "ebi-a-GCST005195", "bbj-a-159"]),
]
print("\n=== Step 3: 固定效应 meta（发现+复制）===")
for m in metas:
    if m:
        print(f"  {m['label']}")
        print(f"    k={m['k']} OR {m['OR']:.3f} ({m['lo']:.3f}-{m['hi']:.3f}) p={m['p']:.3g}"
              + (f" Q={m.get('Q')}" if m.get('Q') is not None else ""))
json.dump(metas, open(DIR + r"\replication_meta.json", "w"), indent=1, ensure_ascii=False)
print("\nsaved: replication_assoc.json / replication_results.json / replication_meta.json")
