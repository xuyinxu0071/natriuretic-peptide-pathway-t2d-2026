# -*- coding: utf-8 -*-
"""
两信号定向复制分析 v2（修正版）
================================
v1 修正项：
  1. rs1421811 等位基因 = G/C（v1 误写 G/A），G = 降压等位基因（NG 2017）。
  2. rs5068 等位基因 = G/A，G = 利钠肽上调等位基因（SCALLOP NT-proBNP beta=+0.1452，
     与主稿 MR Set A 完全同口径）；v1 的 T 是链向错误，方向已全部反转。
  3. 剔除尺度异常数据集（beta/se 缩放 ~1/45，z/p 一致但 OR 不可用）：
     ebi-a-GCST90038613/38610 (Dönertaş UKB)、ukb-d-I9_STR (Neale)——仅作定性参考。
  4. 补充 finn-b-C_STROKE / finn-b-I9_STR_SAH（主稿 rs5068 信号的原始发现数据集）。
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
            print(f"  [HTTP {e.code}] {path}")
            return None
        except Exception as e:
            if i < retries - 1: time.sleep(10); continue
            print("  [ERR]", path, e); return None

# ---------- 暴露侧（修正） ----------
EXPS = {
    "rs198364":  {"gene": "NPPA-expression", "ea": "C", "nea": "T",
                  "beta": 0.426083, "se": 0.0169031, "type": "wald",
                  "desc": "per C = NPPA expression-up (eQTLGen, F=635)"},
    "rs1421811": {"gene": "NPR3-functional", "ea": "G", "nea": "C",
                  "beta": None, "se": None, "type": "per_allele",
                  "desc": "per G = BP-lowering / NPR3-expression-down (NG 2017)"},
    "rs5068":    {"gene": "NPPA-protein", "ea": "G", "nea": "A",
                  "beta": 0.1452, "se": 0.0338, "type": "wald",
                  "desc": "per G = NT-proBNP-up (SCALLOP, main-MR Set A)"},
}

OUTCOMES = {
    # 卒中
    "ebi-a-GCST005838":   "Any stroke (MEGASTROKE) [DISCOVERY]",
    "ebi-a-GCST90018864": "Ischemic stroke (FinnGen+UKB Sakaue) [REP]",
    "bbj-a-129":          "Ischemic stroke (BBJ East Asian) [XANC]",
    "finn-b-C_STROKE":    "Stroke (FinnGen broad) [FG]",
    "finn-b-I9_STR_SAH":  "Stroke incl SAH (FinnGen) [FG]",
    "ebi-a-GCST006910":   "Cardioembolic stroke (MEGASTROKE) [SUBTYPE]",
    "ebi-a-GCST006907":   "Large-artery stroke (MEGASTROKE) [SUBTYPE]",
    "ebi-a-GCST006908":   "Small-vessel stroke (MEGASTROKE) [SUBTYPE]",
    # MI/CAD
    "finn-b-I9_MI":       "MI (FinnGen) [DISCOVERY]",
    "finn-b-I9_CHD":      "CHD (FinnGen) [PRIOR]",
    "ebi-a-GCST011364":   "MI (Hartiala mixed-ancestry) [REP]",
    "ebi-a-GCST005195":   "CAD (vanderHarst UKB+Leipzig) [REP]",
    "ebi-a-GCST90013868": "CAD SPA (Mbatchou GeneBANK) [REP]",
    "bbj-a-159":          "CAD (BBJ East Asian) [XANC]",
    # AF（机制）
    "ebi-a-GCST006061":   "Atrial fibrillation (AFGen Roselli) [MECH]",
    "ebi-a-GCST006414":   "Atrial fibrillation (Nielsen) [MECH]",
    # 尺度异常（仅定性）
    "ebi-a-GCST90038613": "Stroke (UKB Donertas) [SCALE-FLAG]",
    "ebi-a-GCST90038610": "MI (UKB Donertas) [SCALE-FLAG]",
    "ukb-d-I9_STR":       "Stroke excl SAH (Neale UKB) [SCALE-FLAG]",
}
SCALE_FLAG = {"ebi-a-GCST90038613", "ebi-a-GCST90038610", "ukb-d-I9_STR"}

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
def align(e1, e2, o1, o2):
    if e1 == o1 and e2 == o2: return 1
    if e1 == o2 and e2 == o1: return -1
    if COMP.get(e1) == o1 and COMP.get(e2) == o2: return 1
    if COMP.get(e1) == o2 and COMP.get(e2) == o1: return -1
    return None

# ---------- 关联数据：缓存 + 补两个 FinnGen ----------
cache = json.load(open(DIR + r"\replication_assoc.json"))
assoc = {}
for k, v in cache.items():
    rsid, oid = k.split("|", 1)
    assoc[(rsid, oid)] = v

# 缺失补全：任何 (variant, outcome) 缺失即重新请求该结局
for oid in list(OUTCOMES.keys()):
    missing = [rs for rs in EXPS if (rs, oid) not in assoc]
    if missing:
        res = api_post("/associations", {"variant": missing, "id": [oid]})
        time.sleep(2)
        if res:
            for a in res:
                assoc[(a.get("rsid"), oid)] = a
            print(f"refetch {oid}: {len(res)} assoc (was missing {len(missing)})")
        else:
            print(f"refetch FAIL {oid}")
    elif oid in ["finn-b-C_STROKE", "finn-b-I9_STR_SAH"]:
        pass

json.dump({f"{k[0]}|{k[1]}": v for k, v in assoc.items()},
          open(DIR + r"\replication_assoc_v2.json", "w"), indent=1)
print(f"assoc cache complete: {len(assoc)} records")

# ---------- 计算 ----------
def pval(z):
    return 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))

rows = []
for rsid, exp in EXPS.items():
    for oid, label in OUTCOMES.items():
        a = assoc.get((rsid, oid))
        if not a or a.get("beta") is None:
            rows.append({"variant": rsid, "gene": exp["gene"], "oid": oid,
                         "outcome": label, "note": "no assoc"}); continue
        sign = align(exp["ea"], exp["nea"], a.get("ea"), a.get("nea"))
        if sign is None:
            rows.append({"variant": rsid, "gene": exp["gene"], "oid": oid,
                         "outcome": label, "note": "allele mismatch"}); continue
        by, sy = sign * a["beta"], a["se"]
        if exp["type"] == "wald":
            b = by / exp["beta"]
            se = abs(b) * math.sqrt((exp["se"] / exp["beta"]) ** 2 + (sy / by) ** 2)
        else:
            b, se = by, sy
        p = pval(b / se)
        r = {"variant": rsid, "gene": exp["gene"], "oid": oid, "outcome": label,
             "OR": round(math.exp(b), 4), "lo": round(math.exp(b - 1.96 * se), 4),
             "hi": round(math.exp(b + 1.96 * se), 4), "p": p,
             "beta_out": by, "se_out": sy, "scale_flag": oid in SCALE_FLAG}
        rows.append(r)

print("\n===== v2 全部结果（OR<1 = 暴露上调→风险降低）=====")
for r in rows:
    if r.get("OR") is None:
        print(f"  {r['variant']:10s} {r['gene']:17s} {r['outcome'][:48]:48s} -- {r.get('note')}")
    else:
        flag = " [SCALE-FLAG]" if r["scale_flag"] else ""
        print(f"  {r['variant']:10s} {r['gene']:17s} {r['outcome'][:48]:48s} "
              f"OR {r['OR']:.3f} ({r['lo']:.3f}-{r['hi']:.3f}) p={r['p']:.4g}{flag}")

# ---------- 固定效应 meta（剔除 SCALE-FLAG）----------
def meta(label, variant, oids, direction="protective-if-OR<1"):
    ests = []
    for oid in oids:
        if oid in SCALE_FLAG: continue
        for r in rows:
            if r["variant"] == variant and r["oid"] == oid and r.get("OR"):
                ests.append((math.log(r["OR"]),
                             (math.log(r["hi"]) - math.log(r["lo"])) / (2 * 1.96)))
    if not ests: return None
    w = [1 / s ** 2 for _, s in ests]
    b = sum(wi * e for wi, (e, _) in zip(w, ests)) / sum(w)
    se = math.sqrt(1 / sum(w))
    q = sum(wi * (e - b) ** 2 for wi, (e, _) in zip(w, ests))
    df = len(ests) - 1
    # 方向一致性
    same = all((e < 0) == (ests[0][0] < 0) for e, _ in ests)
    return {"label": label, "k": len(ests), "direction_consistent": same,
            "OR": round(math.exp(b), 4), "lo": round(math.exp(b - 1.96 * se), 4),
            "hi": round(math.exp(b + 1.96 * se), 4), "p": pval(b / se),
            "Q": round(q, 3), "df": df}

metas = [
    # 信号 1：NPPA 表达 → 卒中（发现 + 独立复制 + 跨族裔）
    meta("S1a: NPPA-expr -> any stroke [MEGASTROKE only]", "rs198364", ["ebi-a-GCST005838"]),
    meta("S1b: NPPA-expr -> ischemic stroke [Sakaue FinnGen+UKB only]", "rs198364", ["ebi-a-GCST90018864"]),
    meta("S1c: NPPA-expr -> stroke META (MEGASTROKE + Sakaue + BBJ)",
         "rs198364", ["ebi-a-GCST005838", "ebi-a-GCST90018864", "bbj-a-129"]),
    # 信号 2：NPR3 功能变异 → MI/CAD
    meta("S2a: NPR3 -> MI/CAD META (FinnGen + Hartiala + vanderHarst + Mbatchou)",
         "rs1421811", ["finn-b-I9_MI", "ebi-a-GCST011364", "ebi-a-GCST005195", "ebi-a-GCST90013868"]),
    meta("S2b: NPR3 -> MI/CAD + 跨族裔 (加 BBJ)",
         "rs1421811", ["finn-b-I9_MI", "ebi-a-GCST011364", "ebi-a-GCST005195", "ebi-a-GCST90013868", "bbj-a-159"]),
    meta("S2c: NPR3 -> AF META (Nielsen + Roselli)", "rs1421811",
         ["ebi-a-GCST006414", "ebi-a-GCST006061"]),
    # 信号 3（新发现）：NPPA 蛋白功能变异 rs5068 → CAD/MI（与 NPR3 方向互证）
    meta("S3a: NPPA-protein -> CAD/MI META (Hartiala + vanderHarst + Mbatchou + FinnGen CHD)",
         "rs5068", ["ebi-a-GCST011364", "ebi-a-GCST005195", "ebi-a-GCST90013868", "finn-b-I9_CHD"]),
    meta("S3b: NPPA-protein -> stroke (FinnGen broad + SAH, 主稿信号复制)",
         "rs5068", ["finn-b-C_STROKE", "finn-b-I9_STR_SAH"]),
    meta("S3c: NPPA-protein -> AF META (Nielsen + Roselli)", "rs5068",
         ["ebi-a-GCST006414", "ebi-a-GCST006061"]),
]
print("\n===== 固定效应 meta（已剔除 SCALE-FLAG 数据集）=====")
for m in metas:
    if m:
        print(f"  {m['label']}")
        print(f"    k={m['k']} 方向一致={m['direction_consistent']} "
              f"OR {m['OR']:.3f} ({m['lo']:.3f}-{m['hi']:.3f}) p={m['p']:.3g} Q={m['Q']} (df={m['df']})")

json.dump(rows, open(DIR + r"\replication_results_v2.json", "w"), indent=1, ensure_ascii=False)
json.dump(metas, open(DIR + r"\replication_meta_v2.json", "w"), indent=1, ensure_ascii=False)

# ---------- CSV 输出 ----------
import csv
with open(DIR + r"\replication_results_v2.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["variant", "gene", "outcome", "OR", "lo", "hi", "p", "scale_flag", "note"])
    for r in rows:
        w.writerow([r["variant"], r["gene"], r["outcome"], r.get("OR", ""), r.get("lo", ""),
                    r.get("hi", ""), r.get("p", ""), r.get("scale_flag", ""), r.get("note", "")])
print("\nsaved: replication_assoc_v2.json / replication_results_v2.{json,csv} / replication_meta_v2.json")
