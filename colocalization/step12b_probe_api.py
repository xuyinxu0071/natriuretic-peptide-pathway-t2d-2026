# -*- coding: utf-8 -*-
"""step12b: 探测 OpenGWAS tophits/clumps 正确端点"""
import json, os, time, urllib.request, urllib.error

DIR = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"
API = "https://api.opengwas.io/api"
TOK = open(os.path.join(DIR, ".opengwas_token.txt")).read().strip()

def try_req(method, path, payload=None):
    try:
        data = json.dumps(payload).encode() if payload is not None else None
        req = urllib.request.Request(API + path, data=data, method=method,
            headers={"Authorization": "Bearer " + TOK, "Content-Type": "application/json"})
        res = urllib.load(urllib.request.urlopen(req, timeout=120)) if hasattr(urllib, "load") else json.load(urllib.request.urlopen(req, timeout=120))
        return ("OK", res)
    except urllib.error.HTTPError as e:
        return (f"HTTP {e.code}", e.read().decode("utf-8", "ignore")[:120])
    except Exception as e:
        return ("ERR", str(e)[:120])

OUT = []
def p(s): OUT.append(str(s))

# 状态检查
p("status: " + str(try_req("GET", "/status"))[:200])

CANDIDATES = [
    ("GET",  "/tophits/ebi-a-GCST006061/", None),
    ("GET",  "/tophits/ebi-a-GCST006061?clump=1", None),
    ("GET",  "/clumps/ebi-a-GCST006061", None),
    ("POST", "/tophits", {"id": ["ebi-a-GCST006061"], "clump": 1}),
    ("POST", "/tophits", {"id": "ebi-a-GCST006061"}),
    ("GET",  "/association/ebi-a-GCST006061", None),
    ("POST", "/clump", {"id": "ebi-a-GCST006061", "p1": 5e-8}),
    ("GET",  "/docs", None),
    ("GET",  "/openapi.json", None),
]
for method, path, payload in CANDIDATES:
    code, body = try_req(method, path, payload)
    p(f"{method} {path} -> {code}: {str(body)[:200]}")
    time.sleep(0.8)

# openapi.json 若成功则解析端点列表
code, body = try_req("GET", "/openapi.json")
if code == "OK" and isinstance(body, dict):
    p("\n== openapi paths ==")
    for k in sorted(body.get("paths", {}).keys()):
        p("  " + k + "  " + ",".join(body["paths"][k].keys()))

with open(os.path.join(DIR, "coloc_analysis", "step12b_output.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT))
print("DONE")
