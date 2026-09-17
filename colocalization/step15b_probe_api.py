# -*- coding: utf-8 -*-
"""step15b: 探测 OpenGWAS API v3 现行数据集列表端点"""
import json, os, urllib.request

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
TOKEN = open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\.opengwas_token.txt").read().strip()
API = "https://api.opengwas.io/api"
out = []

def get(url, accept="application/json"):
    req = urllib.request.Request(url, headers={
        "Authorization": "Bearer " + TOKEN, "Accept": accept,
        "User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace")

def post(path, payload):
    req = urllib.request.Request(
        API + path, data=json.dumps(payload).encode("utf-8"),
        headers={"Authorization": "Bearer " + TOKEN,
                 "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode("utf-8", "replace")

# GET probes
for ep in ["/status", "/docs", "/openapi.json", "/swagger.json",
           "/gwasinfo", "/datasets/info", "/datasets/search?q=CORIN"]:
    try:
        s = get(API + ep)
        out.append("GET OK %s -> %d bytes" % (ep, len(s)))
        out.append(s[:400])
    except Exception as e:
        out.append("GET FAIL %s -> %s" % (ep, e))

# POST probes
for ep, pl in [("/datasets", {}), ("/datasets/search", {"query": "CORIN"}),
               ("/datasets", {"query": "CORIN"}),
               ("/gwasinfo", ["prot-a-2078"])]:
    try:
        s = post(ep, pl)
        out.append("POST OK %s %s -> %d bytes" % (ep, pl, len(s)))
        out.append(s[:500])
    except Exception as e:
        out.append("POST FAIL %s %s -> %s" % (ep, pl, e))

open(os.path.join(W, "_step15b_probe.txt"), "w", encoding="utf-8").write("\n".join(out))
print("done")
