# -*- coding: utf-8 -*-
"""step17b: 探测 Heart Cell Atlas v2 cellxgene-explorer REST API"""
import json, os, urllib.request

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
out = []

def get(url, timeout=120):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0",
                                               "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")

BASE = "https://www.heartcellatlas.org/v2/global"
for ep in ["/api/v0.2/config", "/api/v0.2/schema",
           "/api/v0.2/annotations/obs?annotation-name=cell_type"]:
    try:
        s = get(BASE + ep)
        out.append("OK %s -> %d bytes" % (ep, len(s)))
        out.append(s[:1200])
    except Exception as e:
        out.append("FAIL %s -> %s" % (ep, e))
    out.append("---")

open(os.path.join(W, "_hca_api_probe.txt"), "w", encoding="utf-8").write("\n".join(out))
print("done")
