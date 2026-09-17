# -*- coding: utf-8 -*-
"""
Step 8b: 下载 OpenGWAS 全量 gwasinfo 索引并过滤利钠肽通路相关数据集
目标: 找出 trait 含 natriuretic / proBNP / NPPA / NPPB / NPR3 / CORIN / neprilysin /
      MME / atrial natriuretic / BNP 的全部数据集 (pQTL/eQTL/biomarker)
"""
import json, urllib.request, urllib.error, time, gzip, os

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

CACHE = WORK + r"\gwasinfo_full.json"
OUT  = open(WORK + r"\step8b_output.txt", "w", encoding="utf-8")

def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

# 1. 下载或读缓存
if os.path.exists(CACHE) and os.path.getsize(CACHE) > 10_000_000:
    log(f"[cache] {os.path.getsize(CACHE)/1e6:.1f} MB")
    data = json.load(open(CACHE, encoding="utf-8"))
else:
    log("[download] full gwasinfo ...")
    req = urllib.request.Request("https://api.opengwas.io/api/gwasinfo",
        headers={"Authorization": "Bearer " + TOK})
    buf = b""
    with urllib.request.urlopen(req, timeout=600) as r:
        while True:
            chunk = r.read(1 << 20)
            if not chunk: break
            buf += chunk
    log(f"[download] {len(buf)/1e6:.1f} MB")
    open(CACHE, "wb").write(buf)
    data = json.loads(buf.decode("utf-8"))
log(f"total datasets: {len(data)}")

# 2. 过滤
KEYS = ["natriuretic", "probnp", "nt-pro", "nppa", "nppb", "npr3", "corin",
        "neprilysin", "mme ", "atrial peptide", "bnp", "cnp", "guanylyl cyclase"]
hits = []
for d in data:
    t = str(d.get("trait", "")).lower()
    if any(k in t for k in KEYS):
        hits.append(d)
log(f"trait-matched: {len(hits)}")
for d in sorted(hits, key=lambda x: str(x.get("id"))):
    log(f"  {d.get('id')}  n={d.get('sample_size')}  {d.get('year')}  {str(d.get('population'))[:10]}  {str(d.get('trait'))[:75]}")

OUT.close()
