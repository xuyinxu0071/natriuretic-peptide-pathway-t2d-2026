# -*- coding: utf-8 -*-
"""
Step 16c: M2 中介拉取
================================================================================
A) nppa locus variants x ieu-b-38 (SBP) / ieu-b-39 (DBP) / ieu-b-40 (BMI) /
   ebi-a-GCST007344 (eGFR) —— NP 暴露 -> 中介 的关联
B) 各中介 tophits (POST /tophits, clumped p<5e-8)
C) 中介 tophits x AF (ebi-a-GCST006061) —— 中介 -> 结局 的关联
"""
import json, os, time, urllib.error, urllib.request

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
TOKEN = open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\.opengwas_token.txt").read().strip()
API = "https://api.opengwas.io/api"

log = open(os.path.join(WORK, "step16c_log.txt"), "w", encoding="utf-8")
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

def pull_assoc(dest, dsid, rsids, batch=64):
    if os.path.exists(dest):
        P("  exists, skip: %s" % os.path.basename(dest))
        return True
    out = {}
    done = 0
    while done < len(rsids):
        chunk = rsids[done:done + batch]
        res = api_post("/associations", {"variant": chunk, "id": [dsid]})
        if res is None:
            P("    batch at %d FAILED - aborting" % done)
            return False
        for a in res:
            if a.get("id") == dsid and a.get("beta") is not None and a.get("se") is not None:
                out[a["rsid"]] = {k: a.get(k) for k in ("beta", "se", "p", "ea", "nea", "eaf", "n")}
        done += batch
        if (done // batch) % 20 == 0:
            P("    %d/%d (have %d)" % (done, len(rsids), len(out)))
        time.sleep(0.6)
    json.dump(out, open(dest, "w"))
    P("  saved %d -> %s" % (len(out), os.path.basename(dest)))
    return True

MEDIATORS = {
    "SBP":  "ieu-b-38",
    "DBP":  "ieu-b-39",
    "BMI":  "ieu-b-40",
    "eGFR": "ebi-a-GCST007344",
    "T2D":  "ebi-a-GCST007515",
}
AF = "ebi-a-GCST006061"

if __name__ == "__main__":
    # --- A) nppa locus x mediator GWAS ---
    nppa_rsids = [v["rsid"] for v in json.load(
        open(os.path.join(WORK, "nppa_locus_variants.json")))]
    P("nppa locus: %d variants" % len(nppa_rsids))
    for mname, mid in MEDIATORS.items():
        if mname == "T2D":
            continue  # T2D 已在 step16b 拉取
        dest = os.path.join(WORK, "nppa__%s_assoc.json" % mid)
        P("=== A: nppa x %s (%s) ===" % (mname, mid))
        pull_assoc(dest, mid, nppa_rsids)
        time.sleep(2)

    # --- B) mediator tophits ---
    for mname, mid in MEDIATORS.items():
        dest = os.path.join(WORK, "step16c_tophits_%s.json" % mname)
        if os.path.exists(dest):
            P("=== B: %s tophits exists, skip ===" % mname)
            continue
        P("=== B: tophits %s (%s) ===" % (mname, mid))
        res = api_post("/tophits", {"id": [mid], "clump": 1})
        if res is None:
            P("  FAILED")
            continue
        json.dump(res, open(dest, "w"))
        P("  %d tophits saved" % len(res))
        time.sleep(2)

    # --- C) AF x mediator tophits ---
    for mname, mid in MEDIATORS.items():
        tp = os.path.join(WORK, "step16c_tophits_%s.json" % mname)
        if not os.path.exists(tp):
            P("=== C: %s no tophits, skip ===" % mname)
            continue
        hits = json.load(open(tp))
        rsids = [h.get("rsid") for h in hits if h.get("rsid")]
        dest = os.path.join(WORK, "step16c_AF__%s.json" % mname)
        P("=== C: AF x %s tophits (%d snps) ===" % (mname, len(rsids)))
        pull_assoc(dest, AF, rsids)
        time.sleep(2)

    P("ALL DONE")
    log.close()
