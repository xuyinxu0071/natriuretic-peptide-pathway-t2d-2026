# -*- coding: utf-8 -*-
"""step15a: M1 侦察 - 在 OpenGWAS 检索 NP 加工酶级联 pQTL 数据集 + T2D 结局数据集"""
import json, os, urllib.request

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
TOKEN = open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\.opengwas_token.txt").read().strip()
API = "https://api.opengwas.io/api"

def get(url):
    req = urllib.request.Request(url, headers={
        "Authorization": "Bearer " + TOKEN,
        "Accept": "application/json",
        "User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode("utf-8"))

out = []
# probe dataset listing endpoints
data = None
for ep in ["/datasets", "/datasets/?size=100000"]:
    try:
        data = get(API + ep)
        out.append("OK %s -> type=%s" % (ep, type(data).__name__))
        break
    except Exception as e:
        out.append("FAIL %s -> %s" % (ep, e))
if data is None:
    open(os.path.join(W, "_step15a_probe.txt"), "w", encoding="utf-8").write("\n".join(out))
    print("end")
    raise SystemExit

# normalize to list of dicts
if isinstance(data, dict):
    for k in ("results", "data", "datasets"):
        if k in data and isinstance(data[k], list):
            data = data[k]
            out.append("unwrapped key: %s, n=%d" % (k, len(data)))
            break
if isinstance(data, dict):
    keys = list(data.keys())
    out.append("dict keys sample: %s" % keys[:10])
    # values may be the records
    vals = list(data.values())
    if vals and isinstance(vals[0], dict):
        data = vals
        out.append("using dict values, n=%d" % len(data))

out.append("total datasets: %d" % len(data))
if data:
    out.append("sample record keys: %s" % sorted(data[0].keys()))

# candidate gene symbols to find in pQTL datasets
targets = ["CORIN", "FURIN", "MME", "DPP4", "NPR1", "NPR2", "NPR3",
           "NPPA", "NPPB", "NATRIURETIC", "NEPRILYSIN", "ATRIAL NATRIURETIC"]
# T2D outcome candidates
t2d_kw = ["type 2 diabetes", "type ii diabetes", "t2d", "fasting glucose",
          "glycated haemoglobin", "hba1c"]

def txt(rec):
    return json.dumps(rec, ensure_ascii=False).lower()

# focus on pQTL-like datasets (prot-, eqtl-, ebi-a-GCST90* SCALLOP, deCODE, UKB-PPP)
pqtl_hits = {}
t2d_hits = []
for rec in data:
    s = txt(rec)
    rid = rec.get("id", "?")
    trait = rec.get("trait", "")
    # pQTL search
    for t in targets:
        tl = t.lower()
        # require word-ish match on trait or note
        if tl in trait.lower() or (tl in s and ("protein" in s or "pqtl" in s or "plasma" in s or "olink" in s or "somascan" in s or "decode" in s)):
            pqtl_hits.setdefault(t, []).append((rid, trait[:90]))
    # T2D search
    for kw in t2d_kw:
        if kw in trait.lower():
            t2d_hits.append((rid, trait[:90], rec.get("sample_size", ""), rec.get("year", "")))
            break

out.append("")
out.append("=== pQTL candidate datasets per target ===")
for t in targets:
    hits = pqtl_hits.get(t, [])
    out.append("%s: %d candidates" % (t, len(hits)))
    for h in hits[:12]:
        out.append("   %s | %s" % h)

out.append("")
out.append("=== T2D outcome candidates (top 25 by sample_size str) ===")
seen = set()
n = 0
for rid, tr, ss, yr in sorted(t2d_hits, key=lambda x: str(x[2]), reverse=True):
    if rid in seen:
        continue
    seen.add(rid)
    out.append("   %s | n=%s | %s" % (rid, ss, tr))
    n += 1
    if n >= 25:
        break

open(os.path.join(W, "step15a_output.txt"), "w", encoding="utf-8").write("\n".join(out))
print("done")
