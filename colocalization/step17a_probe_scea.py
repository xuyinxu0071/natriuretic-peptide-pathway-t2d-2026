# -*- coding: utf-8 -*-
"""step17a: M3 侦察 - EBI Single Cell Expression Atlas 检索人类心脏单细胞数据集"""
import json, os, urllib.request

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
out = []

def get(url, timeout=90):
    req = urllib.request.Request(url, headers={
        "User-Agent": "Mozilla/5.0", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")

# 1) 实验列表
try:
    s = get("https://www.ebi.ac.uk/gxa/sc/json/experiments")
    j = json.loads(s)
    exps = j.get("experiments", [])
    out.append("total SCEA experiments: %d" % len(exps))
    # 找心脏相关
    heart = []
    for e in exps:
        blob = json.dumps(e).lower()
        if any(k in blob for k in ["heart", "cardiac", "cardiomyocyte", "myocard",
                                    "ventricle", "atrial"]):
            if "human" in blob or "homo" in blob or e.get("species", "").lower() in ("homo sapiens",):
                heart.append(e)
    out.append("heart-related (any species): %d" % len(heart))
    for e in heart[:20]:
        out.append("  %s | %s | species=%s" % (
            e.get("experimentAccession"), (e.get("experimentDescription") or "")[:70],
            e.get("species")))
except Exception as ex:
    out.append("experiments list FAIL: %s" % ex)

open(os.path.join(W, "step17a_output.txt"), "w", encoding="utf-8").write("\n".join(out))
print("done")
