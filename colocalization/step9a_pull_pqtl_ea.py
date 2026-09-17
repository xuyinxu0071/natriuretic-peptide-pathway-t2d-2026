# -*- coding: utf-8 -*-
"""
Step 9a: 拉取蛋白层 pQTL + 东亚结局 (NPPA 位点)
================================================================================
新增暴露 (蛋白层 — rs5068 真正作用的层级):
  ebi-a-GCST90012082  NT-proBNP levels, SCALLOP Olink, n=21,758  [主]
  prot-a-2078         NT-proBNP, SomaLogic, n=3,301              [重复]
  prot-a-2076         ANF (ANP), SomaLogic, n=3,301              [NPPA 基因产物]
新增结局 (东亚, 与 CHARLS 人群匹配):
  bbj-a-159           BBJ CAD, n=212,453
  bbj-a-109           BBJ CHF
  + 索引扫描发现的 BBJ stroke/AF/MI
"""
import json, os, time, urllib.request, urllib.error

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
DIR  = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"
API  = "https://api.opengwas.io/api"
import os as _os
_TOK_PATH = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".opengwas_token.txt")
if _os.environ.get("OPENGWAS_TOKEN"):
    TOK = _os.environ["OPENGWAS_TOKEN"].strip()
elif _os.path.isfile(_TOK_PATH):
    TOK = open(_TOK_PATH).read().strip()
else:
    raise RuntimeError("Set env var OPENGWAS_TOKEN or place a gitignored .opengwas_token.txt next to this script to call the OpenGWAS API.")

LOG  = open(WORK + r"\step9a_output.txt", "w", encoding="utf-8")

def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()

def api_post(path, payload, retries=5):
    for i in range(retries):
        try:
            req = urllib.request.Request(API + path, data=json.dumps(payload).encode(),
                headers={"Authorization": "Bearer " + TOK, "Content-Type": "application/json"}, method="POST")
            return json.load(urllib.request.urlopen(req, timeout=300))
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504, 401, 403) and i < retries - 1:
                w = 20 * (i + 1) if e.code == 429 else 10
                log(f"  [HTTP {e.code}] retry {w}s"); time.sleep(w); continue
            log(f"  [HTTP {e.code}] {e.read().decode('utf-8','ignore')[:150]}"); return None
        except Exception as e:
            if i < retries - 1:
                time.sleep(15); continue
            log(f"  [ERR] {e}"); return None

# ---- 1. 扫描 BBJ 心血管结局 ----
idx = json.load(open(WORK + r"\gwasinfo_full.json", encoding="utf-8"))
log("=== BBJ 心血管结局扫描 ===")
EA_OUT = {}
for e in idx.values():
    i = str(e.get("id", ""))
    if not i.startswith("bbj-a"):
        continue
    t = str(e.get("trait", "")).lower()
    if any(k in t for k in ["coronary", "myocardial", "heart failure", "stroke",
                            "atrial fibrillation", "ischaemic", "ischemic", "angina"]):
        EA_OUT[i] = (e.get("trait"), e.get("sample_size"), e.get("ncase"))
        log(f"  {i}: {e.get('trait')} n={e.get('sample_size')} ncase={e.get('ncase')}")

# ---- 2. 拉取 ----
NEW = ["ebi-a-GCST90012082", "prot-a-2078", "prot-a-2076", "bbj-a-159", "bbj-a-109"]
variants = json.load(open(WORK + r"\nppa_locus_variants.json"))
rsids = [v["rsid"] for v in variants]
log(f"\n=== NPPA locus: {len(rsids)} variants x {len(NEW)} datasets ===")

def pull_dataset(dsid, batch=64):
    out, done = {}, 0
    while done < len(rsids):
        chunk = rsids[done:done + batch]
        res = api_post("/associations", {"variant": chunk, "id": [dsid]})
        if res is None:
            log(f"    batch@{done} FAILED — abort"); return None
        for a in res:
            if a.get("id") == dsid and a.get("beta") is not None and a.get("se") is not None:
                out[a["rsid"]] = {k: a.get(k) for k in ("beta", "se", "p", "ea", "nea", "eaf", "n")}
        done += batch
        if (done // batch) % 10 == 0:
            log(f"    {done}/{len(rsids)} (have {len(out)})")
        time.sleep(0.6)
    return out

for dsid in NEW:
    dest = WORK + f"\\nppa__{dsid}_assoc.json"
    if os.path.exists(dest):
        log(f"  {dsid}: exists, skip"); continue
    log(f"  [{dsid}] pulling ...")
    res = pull_dataset(dsid)
    if res is None:
        log(f"  {dsid}: FAILED"); continue
    json.dump(res, open(dest, "w"))
    log(f"  {dsid}: {len(res)} variants saved")
    time.sleep(3)
log("\nALL DONE")
LOG.close()
