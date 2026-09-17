# -*- coding: utf-8 -*-
"""Step 9d: 补拉 BBJ ischemic stroke / AF (NPPA 位点)"""
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

LOG  = open(WORK + r"\step9d_output.txt", "w", encoding="utf-8")

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

variants = json.load(open(WORK + r"\nppa_locus_variants.json"))
rsids = [v["rsid"] for v in variants]
log(f"NPPA locus: {len(rsids)} variants x 2 datasets")

for dsid in ["bbj-a-129", "bbj-a-71"]:
    dest = WORK + f"\\nppa__{dsid}_assoc.json"
    if os.path.exists(dest):
        log(f"  {dsid}: exists, skip"); continue
    log(f"  [{dsid}] pulling ...")
    out, done = {}, 0
    ok = True
    while done < len(rsids):
        chunk = rsids[done:done + 64]
        res = api_post("/associations", {"variant": chunk, "id": [dsid]})
        if res is None:
            log(f"    batch@{done} FAILED"); ok = False; break
        for a in res:
            if a.get("id") == dsid and a.get("beta") is not None and a.get("se") is not None:
                out[a["rsid"]] = {k: a.get(k) for k in ("beta", "se", "p", "ea", "nea", "eaf", "n")}
        done += 64
        time.sleep(0.6)
    if ok:
        json.dump(out, open(dest, "w"))
        log(f"  {dsid}: {len(out)} variants saved")
log("ALL DONE")
LOG.close()
