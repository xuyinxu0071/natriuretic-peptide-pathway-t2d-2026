# -*- coding: utf-8 -*-
"""
Step 0: 前提实测——OpenGWAS 令牌活性 + 关键 SNP 在暴露/结局数据集的返回格式
================================================================================
"""
import json, sys, urllib.request, urllib.error

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
            req = urllib.request.Request(
                API + path,
                data=json.dumps(payload).encode("utf-8"),
                headers={"Authorization": "Bearer " + TOK, "Content-Type": "application/json"},
                method="POST")
            return json.load(urllib.request.urlopen(req, timeout=120))
        except urllib.error.HTTPError as e:
            if e.code == 429 and i < retries - 1:
                import time; time.sleep(20); continue
            print(f"[HTTP {e.code}] {path} {e.read().decode('utf-8','ignore')[:300]}")
            return None
        except Exception as e:
            if i < retries - 1:
                import time; time.sleep(10); continue
            print("[ERR]", e); return None

if __name__ == "__main__":
    # 测试 1: 关键 SNP 在两个 eQTL 数据集 + 主要结局数据集
    res = api_post("/associations", {
        "variant": ["rs5068", "rs1421811", "rs1734792"],
        "id": ["eqtl-a-ENSG00000175206",   # NPPA eQTLGen
               "eqtl-a-ENSG00000113389",   # NPR3 eQTLGen
               "eqtl-a-ENSG00000169174",   # NPPB eQTLGen(待验证存在性)
               "ieu-a-7",                  # CAD CARDIoGRAMplusC4D
               "ebi-a-GCST005195",         # CAD van der Harst
               "ebi-a-GCST006061"]})       # AF Nielsen
    if res is None:
        print("API FAILED"); sys.exit(1)
    for a in res:
        print(f"{a.get('rsid'):12s} {a.get('id'):24s} chr{a.get('chr')}:{a.get('position')} "
              f"beta={a.get('beta')} se={a.get('se')} p={a.get('p')} n={a.get('n')} "
              f"ea={a.get('ea')} nea={a.get('nea')} eaf={a.get('eaf')}")
    print(f"\nTotal records: {len(res)}")
