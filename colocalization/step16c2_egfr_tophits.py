# -*- coding: utf-8 -*-
"""step16c2: 补拉 eGFR tophits (GCST90026654 Yu 2021 n=1.16M) + AF 关联"""
import json, os, time, urllib.error, urllib.request

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
TOKEN = open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\.opengwas_token.txt").read().strip()
API = "https://api.opengwas.io/api"
out = []

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
                time.sleep(10 * (i + 1)); continue
            out.append("[HTTP %d] %s" % (e.code, e.read().decode("utf-8", "ignore")[:150]))
            return None
        except Exception as e:
            if i < retries - 1:
                time.sleep(15); continue
            out.append("[ERR] %s" % e); return None

for dsid, label in [("ebi-a-GCST90026654", "eGFR")]:
    dest = os.path.join(WORK, "step16c_tophits_%s.json" % label)
    if not os.path.exists(dest):
        res = api_post("/tophits", {"id": [dsid], "clump": 1})
        if res is None:
            out.append("%s tophits FAILED" % dsid)
            continue
        json.dump(res, open(dest, "w"))
        out.append("%s: %d tophits" % (dsid, len(res)))
        time.sleep(2)
    hits = json.load(open(dest))
    rsids = [h.get("rsid") for h in hits if h.get("rsid")]
    afd = os.path.join(WORK, "step16c_AF__%s.json" % label)
    res2 = {}
    done = 0
    while done < len(rsids):
        r = api_post("/associations",
                     {"variant": rsids[done:done + 64],
                      "id": ["ebi-a-GCST006061"]})
        if r is None:
            break
        for a in r:
            if a.get("id") == "ebi-a-GCST006061" and a.get("beta") is not None:
                res2[a["rsid"]] = a
        done += 64
        time.sleep(0.6)
    json.dump(res2, open(afd, "w"))
    out.append("AF at %s tophits: %d/%d" % (label, len(res2), len(rsids)))

open(os.path.join(WORK, "step16c2_output.txt"), "w", encoding="utf-8").write("\n".join(out))
print("done")
