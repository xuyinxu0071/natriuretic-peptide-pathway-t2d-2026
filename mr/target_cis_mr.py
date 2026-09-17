# -*- coding: utf-8 -*-
"""
治疗靶点 cis-MR 筛选：利钠肽通路 + 对照锚定
==================================================
设计（先验驱动靶点面板，全部真实 OpenGWAS 数据）：
  通路候选（表达 pQTL/eQTL，eQTLGen via eqtl-a-*，N=31,684）：
      NPR3  eqtl-a-ENSG00000113389  chr5:32.69Mb   头号候选（清除受体，抑制→升利钠肽）
      NPR2  eqtl-a-ENSG00000159899  chr9:35.79Mb   GC-B 受体
      MME   eqtl-a-ENSG00000196549  chr3:155.02Mb  neprilysin（阳性对照：ARNI 已验证）
      CORIN eqtl-a-ENSG00000145244  chr4:47.59Mb   pro-ANP 转化酶
      NPPA  eqtl-a-ENSG00000175206  chr1:11.85Mb   利钠肽配体（表达层）
  阳性对照（蛋白 pQTL）：
      IL6R  ebi-a-GCST90012025 (SCALLOP, N=21,758) + prot-a-1540 (INTERVAL, N=3,301)
  阴性对照（蛋白 pQTL）：
      CRP   ieu-b-35 (Ligthart, N=204,402) + ebi-a-GCST90025959 (Barton, N=436,939)
  结局（9 个，与主稿 MR 臂一致）：
      ieu-a-7                 CHD (CARDIoGRAMplusC4D, N=184,305)
      finn-b-I9_CHD           CHD (FinnGen R9)
      finn-b-I9_MI            MI (FinnGen R9)
      ebi-a-GCST009541        Heart failure (HERMES, N=977,323)
      finn-b-I9_HEARTFAIL     HF strict (FinnGen R9)
      ebi-a-GCST005838        Any stroke (MEGASTROKE, N=446,696)
      ebi-a-GCST006908        Ischemic stroke (MEGASTROKE, N=440,328)
      finn-b-C_STROKE         Stroke (FinnGen)
      finn-b-I9_STR_SAH       Stroke incl SAH (FinnGen)
方法：
  cis 窗口 = 基因 TSS ±500 kb；工具变量 = tophits (p<5e-8, clump r2<0.001) 后取 cis 内；
  F = beta^2/se^2；Wald ratio（单 SNP）/ 固定效应 IVW（多 SNP）；回文 SNP 用 EAF 校验。
  位置构建：eqtl-a-*/prot-a-*/ieu-b-35/ebi-a-* 均为 GRCh37（gwasinfo build 字段）。
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
                API + path,
                data=json.dumps(payload).encode("utf-8"),
                headers={"Authorization": "Bearer " + TOK, "Content-Type": "application/json"},
                method="POST")
            return json.load(urllib.request.urlopen(req, timeout=120))
        except urllib.error.HTTPError as e:
            if e.code == 429 and i < retries - 1:
                wait = 20 * (i + 1)
                print(f"  [429] throttled, wait {wait}s"); time.sleep(wait); continue
            body = e.read().decode("utf-8", "ignore")[:300]
            print(f"  [HTTP {e.code}] {path} {body}")
            return None
        except Exception as e:
            if i < retries - 1:
                time.sleep(10); continue
            print("  [ERR]", path, e); return None

# ---------- 暴露侧定义 ----------
# NPR3/NPR2 eQTLGen 数据集存在但无 p<5e-8 工具变量（实测 0 hits）；
# NPR3 改用文献锚定功能变异 rs1421811（NG 2017：NPR3 表达降低等位基因 → 降血压；
# phewas 实证 SBP p=2.5e-40 / DBP / MAP / PP / 降压药使用全部命中）。
EXPOSURES = {
    "MME_eQTL":    {"id": "eqtl-a-ENSG00000196549", "gene": "MME",   "chrom": "3", "center": 155024124, "role": "positive-ctrl"},
    "CORIN_eQTL":  {"id": "eqtl-a-ENSG00000145244", "gene": "CORIN", "chrom": "4", "center": 47593999, "role": "candidate"},
    "NPPA_eQTL":   {"id": "eqtl-a-ENSG00000175206", "gene": "NPPA",  "chrom": "1", "center": 11845701, "role": "candidate"},
    "IL6R_SCALLOP":  {"id": "ebi-a-GCST90012025", "gene": "IL6R", "chrom": "1", "center": 154405141, "role": "positive-ctrl"},
    "IL6R_INTERVAL": {"id": "prot-a-1540",         "gene": "IL6R", "chrom": "1", "center": 154405141, "role": "positive-ctrl"},
    "CRP_Ligthart":  {"id": "ieu-b-35",            "gene": "CRP",  "chrom": "1", "center": 159711438, "role": "negative-ctrl"},
    "CRP_Barton":    {"id": "ebi-a-GCST90025959",  "gene": "CRP",  "chrom": "1", "center": 159711438, "role": "negative-ctrl"},
}
OUTCOMES = {
    "ieu-a-7": "CHD (CARDIoGRAMplusC4D)",
    "finn-b-I9_CHD": "CHD (FinnGen)",
    "finn-b-I9_MI": "MI (FinnGen)",
    "ebi-a-GCST009541": "Heart failure (HERMES)",
    "finn-b-I9_HEARTFAIL": "HF strict (FinnGen)",
    "ebi-a-GCST005838": "Any stroke (MEGASTROKE)",
    "ebi-a-GCST006908": "Ischemic stroke (MEGASTROKE)",
    "finn-b-C_STROKE": "Stroke (FinnGen)",
    "finn-b-I9_STR_SAH": "Stroke incl SAH (FinnGen)",
}
CIS_KB = 500

# ---------- Step 1: 提取 cis 工具变量 ----------
instruments = {}
for name, ex in EXPOSURES.items():
    print(f"[tophits] {name} ({ex['id']}) ...")
    hits = api_post("/tophits", {"id": [ex["id"]], "pval": 5e-8, "clump": 1})
    time.sleep(3)
    if not hits:
        instruments[name] = {"error": "tophits failed"}; continue
    cis, trans = [], []
    for h in hits:
        try:
            chrom = str(h.get("chr")); pos = int(h.get("position"))
        except (TypeError, ValueError):
            continue
        rec = {"rsid": h["rsid"], "chr": chrom, "pos": pos,
               "beta": h.get("beta"), "se": h.get("se"), "p": h.get("p"),
               "ea": h.get("ea"), "nea": h.get("nea"), "eaf": h.get("eaf"), "n": h.get("n")}
        if chrom == ex["chrom"] and abs(pos - ex["center"]) <= CIS_KB * 1000:
            cis.append(rec)
        else:
            trans.append(rec)
    # F 统计量
    for c in cis:
        c["F"] = (c["beta"] ** 2) / (c["se"] ** 2) if c.get("beta") is not None and c.get("se") else None
    instruments[name] = {"dataset": ex["id"], "gene": ex["gene"], "role": ex["role"],
                         "cis": cis, "trans_n": len(trans)}
    print(f"   tophits={len(hits)} cis={len(cis)} trans={len(trans)}; cis SNPs: {[c['rsid'] for c in cis][:6]}")

json.dump(instruments, open(DIR + r"\target_cis_instruments.json", "w"), indent=1)

# ---------- Step 2: 收集全部 cis SNP，批量结局关联 ----------
all_snps = sorted({c["rsid"] for v in instruments.values() for c in v.get("cis", [])})
print(f"\n[associations] 共 {len(all_snps)} 个 cis SNP × {len(OUTCOMES)} 个结局 ...")

outcome_assoc = {}  # (rsid, outcome_id) -> rec
for oid in OUTCOMES:
    res = api_post("/associations", {"variant": all_snps, "id": [oid]})
    time.sleep(3)
    if not res:
        print(f"   {oid}: FAIL"); continue
    cnt = 0
    for a in res:
        if a.get("beta") is None or a.get("se") is None:
            continue
        outcome_assoc[(a["rsid"], oid)] = {
            "beta": a["beta"], "se": a["se"], "p": a.get("p"),
            "ea": a.get("ea"), "nea": a.get("nea"), "eaf": a.get("eaf"),
            "n": a.get("n"), "trait": OUTCOMES[oid]}
        cnt += 1
    print(f"   {oid}: {cnt} assoc")

json.dump({f"{k[0]}|{k[1]}": v for k, v in outcome_assoc.items()},
          open(DIR + r"\target_outcome_assoc.json", "w"), indent=1)

# ---------- Step 3: Wald / IVW ----------
def wald(bx, sx, by, sy):
    if not bx or not by:
        return None, None
    b = by / bx
    se = abs(b) * math.sqrt((sy / by) ** 2 + (sx / bx) ** 2)  # delta 法
    return b, se

def ivw(pairs):
    # pairs: [(bx,sx,by,sy)]
    bs, ws = [], []
    for bx, sx, by, sy in pairs:
        b, se = by / bx, (sy / abs(bx))
        bs.append(b); ws.append(1 / se ** 2)
    B = sum(w * b for w, b in zip(ws, bs)) / sum(ws)
    SE = math.sqrt(1 / sum(ws))
    return B, SE

def z_p(z):
    return math.erfc(abs(z) / math.sqrt(2))

rows = []
for name, iv in instruments.items():
    if "cis" not in iv or not iv["cis"]:
        rows.append({"target": name, "role": iv.get("role"), "n_iv": 0, "note": "no cis IV"}); continue
    for oid, olabel in OUTCOMES.items():
        pairs, singles = [], []
        palin_warn = []
        for c in iv["cis"]:
            a = outcome_assoc.get((c["rsid"], oid))
            if not a:
                continue
            # 等位基因对齐（简化回文逻辑）
            ea, nea = c.get("ea"), c.get("nea")
            oea, onea = a.get("ea"), a.get("nea")
            if not ea or not nea or not oea or not onea:
                continue
            sign = 1
            if {ea, nea} == {oea, onea}:
                sign = 1 if ea == oea else -1
            else:
                comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
                if comp.get(ea) == oea and comp.get(nea) == onea:
                    fe, fo = c.get("eaf"), a.get("eaf")
                    if fe is not None and fo is not None and abs(fe - fo) < 0.02:
                        sign = 1  # 回文且 EAF 一致，视为同链
                    else:
                        palin_warn.append(c["rsid"]); continue  # 回文且无法定链，跳过
                else:
                    continue  # 无法对齐
            bx, sx = c["beta"], c["se"]
            by, sy = sign * a["beta"], a["se"]
            pairs.append((bx, sx, by, sy))
            singles.append((c["rsid"], wald(bx, sx, by, sy), c.get("F")))
        if not pairs:
            rows.append({"target": name, "role": iv["role"], "outcome": olabel, "n_iv": 0, "note": "no aligned assoc"}); continue
        if len(pairs) == 1:
            b, se = singles[0][1]
            if b is None or se is None or se == 0:
                rows.append({"target": name, "role": iv["role"], "outcome": olabel, "n_iv": 1, "note": "wald undefined (zero beta)"}); continue
            method = "Wald"
            snps = singles[0][0]; F = singles[0][2]
        else:
            b, se = ivw(pairs)
            method = f"IVW({len(pairs)})"
            snps = ",".join(s[0] for s in singles)
            Fs = [s[2] for s in singles if s[2]]
            F = min(Fs) if Fs else None
        z = b / se if se else 0
        p = z_p(z)
        orr = math.exp(b)
        lo, hi = math.exp(b - 1.96 * se), math.exp(b + 1.96 * se)
        rows.append({"target": name, "role": iv["role"], "outcome": olabel, "n_iv": len(pairs),
                     "method": method, "snps": snps, "min_F": round(F, 1) if F else None,
                     "OR": round(orr, 3), "lo": round(lo, 3), "hi": round(hi, 3),
                     "p": round(p, 5), "palindrome_warn": palin_warn})

import csv
keys = ["target", "role", "outcome", "n_iv", "method", "snps", "min_F", "OR", "lo", "hi", "p", "palindrome_warn", "note"]
with open(DIR + r"\target_cis_mr_results.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=keys)
    w.writeheader()
    for r in rows:
        w.writerow({k: r.get(k, "") for k in keys})

print("\n===== cis-MR 主结果（OR<1 = 暴露升高→结局风险降低）=====")
for r in rows:
    if r.get("n_iv"):
        print(f"{r['target']:16s} {r['role']:14s} {r['outcome']:28s} {r['method']:10s} "
              f"OR {r['OR']:.3f} ({r['lo']:.3f}-{r['hi']:.3f}) p={r['p']:.4g} minF={r['min_F']}")
print("\nsaved: target_cis_mr_results.csv")

# ---------- Step 4: NPR3 功能变异 rs1421811（文献锚定仪器，NG 2017）----------
print("\n===== NPR3 功能变异 rs1421811（表达降低等位基因 = 降压等位基因, NG 2017）=====")
# 先取 SBP (ieu-b-38, Evangelou 2018 UKB) 定义降压等位基因
npr3_rows = []
snp = "rs1421811"
res_bp = api_post("/associations", {"variant": [snp], "id": ["ieu-b-38"]})
time.sleep(3)
bp_ea, bp_beta = None, None
if res_bp:
    a0 = res_bp[0]
    bp_ea, bp_beta = a0.get("ea"), a0.get("beta")
    print(f"SBP (ieu-b-38): ea={bp_ea} beta={bp_beta:.4f} (降压等位基因 = {bp_ea})")
for oid, olabel in OUTCOMES.items():
    res = api_post("/associations", {"variant": [snp], "id": [oid]})
    time.sleep(3)
    if not res:
        continue
    for a in res:
        if a.get("beta") is None or a.get("se") is None:
            continue
        sign = 1 if a.get("ea") == bp_ea else -1
        b, se = sign * a["beta"], a["se"]
        z = b / se
        p = z_p(z)
        npr3_rows.append({
            "target": "NPR3_funcVar", "role": "candidate", "outcome": olabel,
            "variant": snp, "bp_lowering_allele": bp_ea,
            "OR": round(math.exp(b), 3), "lo": round(math.exp(b - 1.96 * se), 3),
            "hi": round(math.exp(b + 1.96 * se), 3), "p": round(p, 5),
            "interpretation": "OR per BP-lowering (= NPR3-expression-lowering) allele"})
        print(f"   {olabel:28s} OR {npr3_rows[-1]['OR']:.3f} ({npr3_rows[-1]['lo']:.3f}-{npr3_rows[-1]['hi']:.3f}) p={p:.4g}")

# 追加写入 CSV
with open(DIR + r"\target_cis_mr_results.csv", "a", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["target", "role", "outcome", "n_iv", "method", "snps", "min_F",
                                      "OR", "lo", "hi", "p", "palindrome_warn", "note",
                                      "variant", "bp_lowering_allele", "interpretation"])
    for r in npr3_rows:
        w.writerow(r)
json.dump(npr3_rows, open(DIR + r"\npr3_funcvar_results.json", "w"), indent=1)
print("\nsaved: npr3_funcvar_results.json (appended to target_cis_mr_results.csv)")
