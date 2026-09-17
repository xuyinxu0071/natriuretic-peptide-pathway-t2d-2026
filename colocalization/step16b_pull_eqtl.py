# -*- coding: utf-8 -*-
"""
Step 16b: 拉取 NP 级联暴露 eQTL 关联 + FURIN 蛋白 + 既有位点 T2D
================================================================================
- corin/furin/mme/dpp4/npr2 locus variants x 各自 eqtl-a-ENSG... (eQTLGen)
- furin variants x prot-a-1150 (INTERVAL Furin protein)
- nppa/npr3 locus variants x ebi-a-GCST007515 (T2D Mahajan 2018)
"""
import json, os, time, urllib.error, urllib.request

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
TOKEN = open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\.opengwas_token.txt").read().strip()
API = "https://api.opengwas.io/api"

log = open(os.path.join(WORK, "step16b_log.txt"), "w", encoding="utf-8")
def P(s):
    print(s, flush=True); log.write(s + "\n"); log.flush()

def api_post(path, payload, retries=6):
    for i in range(retries):
        try:
            req = urllib.request.Request(
                API + path, data=json.dumps(payload).encode("utf-8"),
                headers={"Authorization": "Bearer " + TOKEN,
                         "Content-Type": "application/json"})
            return json.load(urllib.request.urlopen(req, timeout=300))
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and i < retries - 1:
                wait = 20 * (i + 1) if e.code == 429 else 10
                P("  [HTTP %d] retry in %ds" % (e.code, wait))
                time.sleep(wait); continue
            body = e.read().decode("utf-8", "ignore")[:200]
            P("  [HTTP %d] %s" % (e.code, body))
            return None
        except Exception as e:
            if i < retries - 1:
                time.sleep(15); continue
            P("  [ERR] %s" % e); return None

def pull_dataset(locus, dsid, rsids, batch=64):
    out = {}
    done = 0
    while done < len(rsids):
        chunk = rsids[done:done + batch]
        res = api_post("/associations", {"variant": chunk, "id": [dsid]})
        if res is None:
            P("    batch at %d FAILED - aborting dataset" % done)
            return None
        for a in res:
            if a.get("id") == dsid and a.get("beta") is not None and a.get("se") is not None:
                out[a["rsid"]] = {k: a.get(k) for k in ("beta", "se", "p", "ea", "nea", "eaf", "n")}
        done += batch
        if (done // batch) % 20 == 0:
            P("    %d/%d (have %d)" % (done, len(rsids), len(out)))
        time.sleep(0.6)
    return out

PULLS = [
    ("corin", "eqtl-a-ENSG00000145244"),
    ("furin", "eqtl-a-ENSG00000140564"),
    ("furin", "prot-a-1150"),           # INTERVAL Furin protein (SomaScan)
    ("mme",   "eqtl-a-ENSG00000196549"),
    ("dpp4",  "eqtl-a-ENSG00000197635"),
    ("npr2",  "eqtl-a-ENSG00000159899"),
    ("nppa",  "ebi-a-GCST007515"),      # T2D Mahajan 2018
    ("npr3",  "ebi-a-GCST007515"),
]

if __name__ == "__main__":
    for locus, dsid in PULLS:
        dest = os.path.join(WORK, "%s__%s_assoc.json" % (locus, dsid))
        if os.path.exists(dest):
            P("%s x %s: exists, skip" % (locus, dsid))
            continue
        rsids = [v["rsid"] for v in json.load(
            open(os.path.join(WORK, "%s_locus_variants.json" % locus)))]
        P("=== %s (%d variants) x %s ===" % (locus, len(rsids), dsid))
        res = pull_dataset(locus, dsid, rsids)
        if res is None:
            P("  %s: FAILED" % dsid)
            continue
        json.dump(res, open(dest, "w"))
        P("  %s x %s: %d variants saved" % (locus, dsid, len(res)))
        time.sleep(3)
    P("ALL DONE")
    log.close()
