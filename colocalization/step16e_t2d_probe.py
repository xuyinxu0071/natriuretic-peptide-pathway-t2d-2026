# -*- coding: utf-8 -*-
"""step16e: 探测替代 T2D 数据集覆盖率 (finn-b-T2D / GCST90018926 / ieu-a-24)"""
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
                time.sleep(10 * (i + 1) if e.code == 429 else 10); continue
            out.append("[HTTP %d] %s" % (e.code, e.read().decode("utf-8", "ignore")[:150]))
            return None
        except Exception as e:
            if i < retries - 1:
                time.sleep(15); continue
            out.append("[ERR] %s" % e); return None

# 测试集: AF 数据集中命中的 nppa 前 128 个变体(代表常见覆盖)
af = json.load(open(os.path.join(WORK, "nppa__ebi-a-GCST006061_assoc.json")))
test = sorted(af.keys())[:128]

for dsid in ["finn-b-T2D", "ebi-a-GCST90018926", "ieu-a-24",
             "ebi-a-GCST007515"]:
    hits = 0
    done = 0
    while done < len(test):
        chunk = test[done:done + 64]
        r = api_post("/associations", {"variant": chunk, "id": [dsid]})
        if r is None:
            break
        hits += sum(1 for a in r if a.get("id") == dsid
                    and a.get("beta") is not None)
        done += 64
        time.sleep(0.6)
    out.append("%s: %d/%d hit" % (dsid, hits, len(test)))

open(os.path.join(WORK, "step16e_output.txt"), "w", encoding="utf-8").write("\n".join(out))
print("done")
