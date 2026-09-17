# -*- coding: utf-8 -*-
"""
Step 3: OpenGWAS 批量拉取两位点 × 全部暴露/结局数据集的区域关联
================================================================================
产出: {locus}__{dataset}_assoc.json  (rsid -> {beta,se,p,ea,nea,eaf,n})
"""
import json, os, time, urllib.request, urllib.error

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
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


PULLS = {
    "nppa": [
        "eqtl-a-ENSG00000175206",   # NPPA eQTLGen (暴露)
        "ieu-a-7",                  # CAD CARDIoGRAMplusC4D
        "ebi-a-GCST005195",         # CAD van der Harst 2018
        "ebi-a-GCST006061",         # AF Nielsen 2018
        "ebi-a-GCST009541",         # HF HERMES
        "ebi-a-GCST005838",         # Any stroke MEGASTROKE
        "ebi-a-GCST006910",         # Cardioembolic stroke MEGASTROKE
    ],
    "npr3": [
        "eqtl-a-ENSG00000113389",   # NPR3 eQTLGen (暴露)
        "ieu-b-38",                 # SBP Evangelou 2018 (暴露)
        "ieu-a-7",                  # CAD
        "ebi-a-GCST005195",         # CAD
        "ebi-a-GCST011364",         # MI Hartiala 2021
        "finn-b-I9_MI",             # MI FinnGen R9
    ],
}

def api_post(path, payload, retries=5):
    for i in range(retries):
        try:
            req = urllib.request.Request(
                API + path,
                data=json.dumps(payload).encode("utf-8"),
                headers={"Authorization": "Bearer " + TOK, "Content-Type": "application/json"},
                method="POST")
            return json.load(urllib.request.urlopen(req, timeout=300))
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504, 401, 403) and i < retries - 1:
                wait = 20 * (i + 1) if e.code == 429 else 10
                print(f"  [HTTP {e.code}] retry in {wait}s", flush=True)
                time.sleep(wait); continue
            body = e.read().decode("utf-8", "ignore")[:200]
            print(f"  [HTTP {e.code}] {body}", flush=True)
            return None
        except Exception as e:
            if i < retries - 1:
                time.sleep(15); continue
            print("  [ERR]", e, flush=True); return None

def pull_dataset(locus, dsid, rsids, batch=64):
    out = {}
    done = 0
    while done < len(rsids):
        chunk = rsids[done:done + batch]
        res = api_post("/associations", {"variant": chunk, "id": [dsid]})
        if res is None:
            print(f"    batch at {done} FAILED — aborting dataset", flush=True)
            return None
        for a in res:
            if a.get("id") == dsid and a.get("beta") is not None and a.get("se") is not None:
                out[a["rsid"]] = {k: a.get(k) for k in ("beta", "se", "p", "ea", "nea", "eaf", "n")}
        done += batch
        if (done // batch) % 10 == 0:
            print(f"    {done}/{len(rsids)} (have {len(out)})", flush=True)
        time.sleep(0.6)
    return out

if __name__ == "__main__":
    # 先实测批量上限 (N(id)*N(variant) <= 64)
    test = api_post("/associations", {"variant": [f"rs{i}" for i in range(100000, 100064)],
                                      "id": ["ieu-a-7"]})
    print(f"batch-1000 probe: {'OK ' + str(len(test)) if test is not None else 'FAILED'}", flush=True)
    if test is None:
        exit(1)

    for locus, datasets in PULLS.items():
        vpath = os.path.join(WORK, f"{locus}_locus_variants.json")
        variants = json.load(open(vpath))
        rsids = [v["rsid"] for v in variants]
        print(f"\n=== {locus} locus: {len(rsids)} variants x {len(datasets)} datasets ===", flush=True)
        for dsid in datasets:
            dest = os.path.join(WORK, f"{locus}__{dsid}_assoc.json")
            if os.path.exists(dest):
                print(f"  {dsid}: exists, skip", flush=True)
                continue
            print(f"  [{dsid}] pulling ...", flush=True)
            res = pull_dataset(locus, dsid, rsids)
            if res is None:
                print(f"  {dsid}: FAILED", flush=True)
                continue
            json.dump(res, open(dest, "w"))
            print(f"  {dsid}: {len(res)} variants saved", flush=True)
            time.sleep(3)
    print("\nALL DONE", flush=True)
