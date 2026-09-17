# -*- coding: utf-8 -*-
"""
Step 8: 数据集重筛可行性实测
================================================================================
Q: 是否存在比"全血 eQTL"更合适的层级（蛋白 pQTL / 心脏组织 eQTL / 东亚结局）
   的公开汇总统计量，可用于对 NPPA/NPR3 位点做先验正当的 coloc 重筛？
路线:
  A. GWAS Catalog: rs5068 / rs1421811 的所有关联 -> 归属研究 + 是否有全量 sumstats
  B. OpenGWAS: 候选 pQTL GCST 是否已有 ebi-a 镜像（可直接复用流水线）
  C. 已知 pQTL 资源盘点 (UKB-PPP prot-c / deCODE / Fenland / INTERVAL prot-b)
"""
import json, urllib.request, urllib.error, time

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
DIR  = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"
import os as _os
_TOK_PATH = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".opengwas_token.txt")
if _os.environ.get("OPENGWAS_TOKEN"):
    TOK = _os.environ["OPENGWAS_TOKEN"].strip()
elif _os.path.isfile(_TOK_PATH):
    TOK = open(_TOK_PATH).read().strip()
else:
    raise RuntimeError("Set env var OPENGWAS_TOKEN or place a gitignored .opengwas_token.txt next to this script to call the OpenGWAS API.")

OUT  = open(WORK + r"\step8_output.txt", "w", encoding="utf-8")

def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

def get(url, headers=None, retries=4):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers or {"User-Agent": "Mozilla/5.0"})
            return json.load(urllib.request.urlopen(req, timeout=120))
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and i < retries - 1:
                time.sleep(10 * (i + 1)); continue
            log(f"  [HTTP {e.code}] {url[:110]}"); return None
        except Exception as e:
            if i < retries - 1:
                time.sleep(8); continue
            log(f"  [ERR] {e} {url[:110]}"); return None

# ---------- A. GWAS Catalog: 两个 SNP 的关联 -> 研究 ----------
GWAS = "https://www.ebi.ac.uk/gwas/rest/api"
for snp in ["rs5068", "rs1421811"]:
    log(f"\n===== A. GWAS Catalog associations for {snp} =====")
    url = f"{GWAS}/singleNucleotidePolymorphisms/{snp}/associations?size=500"
    d = get(url)
    if not d: continue
    assoc = d.get("_embedded", {}).get("associations", [])
    log(f"  associations: {len(assoc)}")
    studies = {}
    for a in assoc:
        sl = a.get("_links", {}).get("study", {}).get("href")
        if sl: studies[sl] = None
    log(f"  unique studies: {len(studies)}")
    for sl in sorted(studies):
        st = get(sl)
        if not st: continue
        acc = st.get("accessionId", "?")
        pub = (st.get("publicationInfo") or {}).get("title", "")[:70]
        trait = st.get("diseaseTrait", {}).get("trait", "?") if st.get("diseaseTrait") else "?"
        full = st.get("fullPvalueSet", False)
        snp_study = st.get("initialSampleSize", "")[:40]
        studies[sl] = acc
        log(f"  {acc}  full={full}  {trait[:45]:45s}  {pub}")
        time.sleep(0.3)

# ---------- B/C. OpenGWAS: 候选 pQTL / 结局数据集是否可用 ----------
OG = "https://api.opengwas.io/api"
def ogwas_post(path, payload, retries=4):
    for i in range(retries):
        try:
            req = urllib.request.Request(OG + path, data=json.dumps(payload).encode(),
                headers={"Authorization": "Bearer " + TOK, "Content-Type": "application/json"}, method="POST")
            return json.load(urllib.request.urlopen(req, timeout=180))
        except urllib.error.HTTPError as e:
            if e.code == 401:
                log("  [OpenGWAS 401] token 过期，需刷新"); return None
            if e.code in (429, 500, 502, 503, 504) and i < retries - 1:
                time.sleep(10 * (i + 1)); continue
            log(f"  [OpenGWAS HTTP {e.code}]"); return None
        except Exception as e:
            if i < retries - 1:
                time.sleep(8); continue
            log(f"  [OpenGWAS ERR] {e}"); return None

# 候选: deCODE pQTL / UKB-PPP / Fenland / INTERVAL / BBJ / GIGASTROKE
candidates = [
    # deCODE pQTL Ferkingstad 2021 (SomaScan, 4,907 proteins) — NPPA/NPPB/NPR3 探针
    "ebi-a-GCST90012062",  # 占位: 需以 A 部分结果为准, 先探测命名空间可用性
    # INTERVAL SomaLogic Sun 2018 (prot-b 命名空间前缀样本)
    "prot-b-1", "prot-b-100",
    # UKB-PPP Olink (prot-c 命名空间前缀样本)
    "prot-c-1", "prot-c-100",
    # BBJ 东亚结局 (Kanai 2018 / Ishigaki 2020)
    "bbj-a-159", "bbj-a-109", "bbj-a-151",
    # GIGASTROKE 2022 (更大卒中)
    "ebi-a-GCST90018695",
]
log("\n===== B. OpenGWAS 命名空间探测 =====")
r = ogwas_post("/gwasinfo", {"id": candidates})
if r:
    found = {x["id"]: x for x in r}
    for cid in candidates:
        if cid in found:
            m = found[cid]
            log(f"  [OK] {cid}: {str(m.get('trait'))[:60]} n={m.get('sample_size')} year={m.get('year')}")
        else:
            log(f"  [--] {cid}: 不存在")

OUT.close()
log("\nDONE -> step8_output.txt")
