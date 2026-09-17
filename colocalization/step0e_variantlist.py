# -*- coding: utf-8 -*-
"""Step 0e: 变体列表获取路线实测——A) OpenGWAS tophits(p=1,clump=0); B) GTEx per-gene API"""
import json, time, urllib.request, urllib.error

DIR = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"
API = "https://api.opengwas.io/api"
import os as _os
_TOK_PATH = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".opengwas_token.txt")
if _os.environ.get("OPENGWAS_TOKEN"):
    TOK = _os.environ["OPENGWAS_TOKEN"].strip()
elif _os.path.isfile(_TOK_PATH):
    TOK = open(_TOK_PATH).read().strip()
else:
    raise RuntimeError("Set env var OPENGWAS_TOKEN or place a gitignored .opengwas_token.txt next to this script to call the OpenGWAS API.")


def api_post(path, payload, retries=3):
    for i in range(retries):
        try:
            req = urllib.request.Request(API + path,
                data=json.dumps(payload).encode("utf-8"),
                headers={"Authorization": "Bearer " + TOK, "Content-Type": "application/json"},
                method="POST")
            return json.load(urllib.request.urlopen(req, timeout=180))
        except urllib.error.HTTPError as e:
            if e.code == 429 and i < retries - 1:
                time.sleep(30); continue
            print(f"[HTTP {e.code}] {path} {e.read().decode('utf-8','ignore')[:300]}")
            return None
        except Exception as e:
            if i < retries - 1:
                time.sleep(10); continue
            print("[ERR]", e); return None

if __name__ == "__main__":
    print("=== A) tophits pval=1 clump=0 on NPPA eQTLGen ===")
    r = api_post("/tophits", {"id": ["eqtl-a-ENSG00000175206"], "pval": 1, "clump": 0})
    if r is not None:
        print(f"Returned variants: {len(r)}")
        if r:
            print("Sample rec keys:", sorted(r[0].keys()))
            with_eaf = sum(1 for x in r if x.get("eaf") is not None)
            print(f"with eaf: {with_eaf}")
            for x in sorted(r, key=lambda v: v.get("p", 1))[:5]:
                print(f"  {x['rsid']} chr{x.get('chr')}:{x.get('position')} beta={x.get('beta')} se={x.get('se')} p={x.get('p')} eaf={x.get('eaf')} n={x.get('n')}")

    print("\n=== A2) tophits pval=1 clump=0 on NPR3 eQTLGen ===")
    r2 = api_post("/tophits", {"id": ["eqtl-a-ENSG00000113389"], "pval": 1, "clump": 0})
    if r2 is not None:
        print(f"Returned variants: {len(r2)}")
        if r2:
            sig = [x for x in r2 if x.get("p", 1) < 5e-8]
            print(f"p<5e-8: {len(sig)}")
            for x in sorted(r2, key=lambda v: v.get("p", 1))[:5]:
                print(f"  {x['rsid']} chr{x.get('chr')}:{x.get('position')} p={x.get('p')}")

    print("\n=== B) GTEx v8 per-gene cis-eQTL (NPPA, Heart-Atrial Appendage) ===")
    try:
        req = urllib.request.Request(
            "https://gtexportal.org/api/v2/association/singleTissueEqtl",
            data=json.dumps({"gencodeId": ["ENSG00000175206.7"],
                             "tissueSiteDetailId": ["Heart_Atrial_Appendage"]}).encode(),
            headers={"Content-Type": "application/json"}, method="POST")
        g = json.load(urllib.request.urlopen(req, timeout=120))
        data = g.get("data", [])
        print(f"GTEx records: {len(data)}")
        if data:
            print("keys:", sorted(data[0].keys()))
    except Exception as e:
        print("GTEx failed:", e)
