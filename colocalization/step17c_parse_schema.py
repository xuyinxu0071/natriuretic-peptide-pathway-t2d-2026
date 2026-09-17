# -*- coding: utf-8 -*-
"""step17c: 解析 HCA v2 schema，定位细胞类型列，并取 obs 注释"""
import json, os, urllib.request

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
out = []
BASE = "https://www.heartcellatlas.org/v2/global"

def get(url, accept="application/json", timeout=300):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0",
                                               "Accept": accept})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()

# schema 已在探测时确认可用, 完整保存
raw = get(BASE + "/api/v0.2/schema").decode("utf-8")
sch = json.loads(raw)
open(os.path.join(W, "hca_v2_schema.json"), "w", encoding="utf-8").write(raw)

cols = sch["schema"]["annotations"]["obs"]["columns"]
out.append("obs columns: %s" % [c["name"] for c in cols])
for c in cols:
    if "categories" in c:
        cats = c["categories"]
        out.append("  %s (%d categories): %s" % (c["name"], len(cats), cats[:30]))

varcols = sch["schema"]["annotations"]["var"]["columns"]
out.append("var columns: %s" % [c["name"] for c in varcols])

# 尝试取细胞类型 obs 注释（多种列名/头组合）
for name in ["cell_type", "Cell_Type", "celltype", "major_celltype",
             "cell_state", "Subclass_cell_type"]:
    for acc in ["application/json", "*/*"]:
        try:
            data = get(BASE + "/api/v0.2/annotations/obs?annotation-name=" + name,
                       accept=acc).decode("utf-8")
            j = json.loads(data)
            out.append("")
            out.append("OK obs annotation %s (accept=%s): %s" % (name, acc, data[:500]))
            break
        except Exception as e:
            out.append("FAIL %s accept=%s -> %s" % (name, acc, e))

open(os.path.join(W, "step17c_output.txt"), "w", encoding="utf-8").write("\n".join(out))
print("done")
