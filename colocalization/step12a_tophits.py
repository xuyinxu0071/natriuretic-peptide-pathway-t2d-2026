# -*- coding: utf-8 -*-
"""
Step 12a (P2): Reverse MR 第一步 — 拉取疾病 GWAS 的全基因组显著 IV (tophits)
================================================================================
暴露: AF (Nielsen), CAD (CARDIoGRAMplusC4D / van der Harst), MI (Hartiala), HF (HERMES)
结局: NT-proBNP pQTL (SCALLOP ebi-a-GCST90012082 主; INTERVAL prot-a-2078 敏感性)
输出: step12_tophits_{ds}.json
"""
import json, os, time, urllib.request, urllib.error

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
DIR = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"
API = "https://api.opengwas.io/api"
TOK = open(os.path.join(DIR, ".opengwas_token.txt")).read().strip()

EXPOSURES = {
    "ebi-a-GCST006061": "AF (Nielsen 2018, n~1M)",
    "ieu-a-7": "CAD (CARDIoGRAMplusC4D, n=184k)",
    "ebi-a-GCST005195": "CAD (van der Harst 2018, n=547k)",
    "ebi-a-GCST011364": "MI (Hartiala 2021)",
    "ebi-a-GCST009541": "HF (HERMES, n=977k)",
}

def api_post(path, payload, retries=5):
    for i in range(retries):
        try:
            req = urllib.request.Request(API + path,
                data=json.dumps(payload).encode("utf-8"),
                headers={"Authorization": "Bearer " + TOK, "Content-Type": "application/json"},
                method="POST")
            return json.load(urllib.request.urlopen(req, timeout=300))
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504, 401, 403) and i < retries - 1:
                wait = 20 * (i + 1) if e.code == 429 else 10
                time.sleep(wait); continue
            body = e.read().decode("utf-8", "ignore")[:300]
            return {"__error__": f"HTTP {e.code}: {body}"}
        except Exception as e:
            if i < retries - 1:
                time.sleep(15); continue
            return {"__error__": str(e)}

OUT = []
def p(s):
    OUT.append(str(s))

for dsid, label in EXPOSURES.items():
    dest = os.path.join(WORK, f"step12_tophits_{dsid}.json")
    if os.path.exists(dest):
        p(f"{dsid} [{label}]: 已存在, 跳过"); continue
    p(f"\n=== {dsid} [{label}] tophits ===")
    res = api_post("/tophits", {"id": [dsid], "clump": 1, "r2": 0.001, "kb": 10000, "p1": 5e-8})
    if isinstance(res, dict) and "__error__" in res:
        p(f"  GET /tophits 失败: {res['__error__']}")
        # 备选: POST /tophits
        continue
    if not res:
        p("  空结果")
        continue
    json.dump(res, open(dest, "w"))
    p(f"  保存 {len(res)} 个 IV; 首条: {json.dumps(res[0], ensure_ascii=False)[:300]}")
    time.sleep(1)

with open(os.path.join(WORK, "step12a_output.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT))
print("DONE")
