# -*- coding: utf-8 -*-
"""
Step 17d: M3 - 从 HCA v2 Global_lognormalised.h5ad (8.7GB, 远程 Range 读取)
          提取 NP 通路基因的单细胞表达 + 细胞类型/区域注释
================================================================================
仅下载: /X 目标基因的 CSC 列、obs 编码列、var 基因名 —— 总传输量 < 100MB。
若 X 非 CSC, 打印诊断后退出(启用回退方案)。
产出: step17d_output.txt + hca_np_expression.csv + hca_np_region.csv
"""
import json, os
import numpy as np

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
URL = "https://cellgeni.cog.sanger.ac.uk/heartcellatlas/v2/Global_lognormalised.h5ad"
GENES = ["NPPA", "NPPB", "NPPC", "NPR1", "NPR2", "NPR3",
         "CORIN", "MME", "FURIN", "DPP4", "KLKB1", "NPPF"]

out = []
def P(s):
    print(s, flush=True); out.append(str(s))

import fsspec, h5py

P("opening remote h5ad via fsspec HTTP Range ...")
f = fsspec.open(URL, "rb", block_size=4 * 1024 * 1024).open()
h5 = h5py.File(f, "r")

P("root keys: %s" % list(h5.keys()))
enc = dict(h5["X"].attrs).get("encoding-type", "?")
P("X encoding-type: %s | keys: %s" % (enc, list(h5["X"].keys())))
P("obs keys: %s" % list(h5["obs"].keys()))
P("var keys: %s" % list(h5["var"].keys()))
P("obs attrs: %s" % dict(h5["obs"].attrs))
P("var attrs: %s" % dict(h5["var"].attrs))

if enc != "csc_matrix":
    P(">>> X is %s, NOT gene-major CSC - diagnostics only, aborting" % enc)
    open(os.path.join(W, "step17d_output.txt"), "w", encoding="utf-8").write("\n".join(out))
    raise SystemExit(0)

# ---------- var ----------
var = h5["var"]
vidx = var.attrs.get("_index")
vidx = vidx.decode() if isinstance(vidx, bytes) else str(vidx)
P("var index dataset: %s" % vidx)
gene_names = var[vidx][()].astype(str)
P("n_vars=%d, sample=%s" % (len(gene_names), list(gene_names[:5])))
targets = {}
for g in GENES:
    hits = np.where(gene_names == g)[0]
    if len(hits):
        targets[g] = int(hits[0])
        P("  %s -> var %d" % (g, hits[0]))
    else:
        P("  %s NOT FOUND" % g)

# ---------- obs ----------
obs = h5["obs"]
def read_cat(col):
    grp = obs[col]
    codes = grp["codes"][()]
    cats = grp["categories"][()].astype(str)
    return cats[codes]

labels = {}
for col in ["cell_type", "cell_state", "region"]:
    if col in obs:
        vals = read_cat(col)
        labels[col] = vals
        P("obs %s: %d cells" % (col, len(vals)))

# ---------- X CSC 目标列 ----------
X = h5["X"]
indptr = X["indptr"][()]
n_obs = indptr.shape[0] - 1
P("CSC indptr len=%d -> n_vars=%d ; first cells col indptr[0]=%d"
  % (indptr.shape[0], n_obs, indptr[0]))
# CSC: indptr 索引 = var(基因), 列内 = obs(细胞)
expr = {}   # gene -> dense vector (n_obs,)
for g, vi in targets.items():
    s, e = int(indptr[vi]), int(indptr[vi + 1])
    if e <= s:
        expr[g] = np.zeros(n_obs, dtype=np.float32)
        P("  %s: 0 nnz" % g)
        continue
    data = X["data"][s:e]
    rows = X["indices"][s:e]
    v = np.zeros(n_obs, dtype=np.float32)
    v[rows] = data
    expr[g] = v
    P("  %s: %d nnz, mean=%.4f" % (g, e - s, v.mean()))

# ---------- 汇总表 ----------
def summarize(colname):
    lab = labels.get(colname)
    if lab is None:
        return
    cats = sorted(set(lab.tolist()))
    P("")
    P("=== per-%s: mean lognorm | pct_expressing ===" % colname)
    hdr = ["%s" % colname, "n_cells"] + ["%s_mean" % g for g in expr] + \
          ["%s_pct" % g for g in expr]
    rows = [",".join(hdr)]
    for c in cats:
        m = lab == c
        n = int(m.sum())
        means = ["%.4f" % expr[g][m].mean() for g in expr]
        pcts = ["%.1f" % (100.0 * (expr[g][m] > 0).mean()) for g in expr]
        rows.append(",".join([c, str(n)] + means + pcts))
    fn = "hca_np_%s.csv" % colname
    open(os.path.join(W, fn), "w", encoding="utf-8").write("\n".join(rows))
    P("saved " + fn)
    # 控制台简表: 每基因 top3 细胞类型
    for g in expr:
        tops = sorted(cats, key=lambda c: -expr[g][lab == c].mean())[:3]
        P("  %-6s top: %s" % (g, "; ".join(
            "%s(%.3f)" % (t, expr[g][lab == t].mean()) for t in tops)))

summarize("cell_type")
summarize("cell_state")
summarize("region")

open(os.path.join(W, "step17d_output.txt"), "w", encoding="utf-8").write("\n".join(out))
print("DONE")
