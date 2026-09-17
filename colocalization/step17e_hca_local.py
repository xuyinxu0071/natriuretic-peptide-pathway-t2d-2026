# -*- coding: utf-8 -*-
"""
Step 17e: M3 本地提取 —— HCA v2 Global_lognormalised.h5ad (8.13GB 已下载)
          NP 通路基因 x 细胞类型/细胞状态/心脏区域 表达
================================================================================
CSR 单遍分块扫描: 每块读 indices/data, 命中目标基因列, 按注释累加 sum/n>0。
产出: hca_np_cell_type.csv / hca_np_cell_state.csv / hca_np_region.csv + 输出日志
"""
import json, os
import numpy as np
import h5py

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
H5 = os.path.join(W, "Global_lognormalised.h5ad")
GENES = ["NPPA", "NPPB", "NPPC", "NPR1", "NPR2", "NPR3",
         "CORIN", "MME", "FURIN", "DPP4", "KLKB1"]

out = []
def P(s):
    print(s, flush=True); out.append(str(s))

h5 = h5py.File(H5, "r")
P("root: %s" % list(h5.keys()))
X = h5["X"]
enc = dict(X.attrs).get("encoding-type", "?")
P("X encoding: %s shapes: data=%s indices=%s indptr=%s" %
  (enc, X["data"].shape, X["indices"].shape, X["indptr"].shape))
assert enc == "csr_matrix"

var = h5["var"]
vidx = var.attrs.get("_index")
vidx = vidx.decode() if isinstance(vidx, bytes) else str(vidx)
gene_names = var[vidx][()].astype(str)
P("n_vars=%d" % len(gene_names))
targets = {}
for g in GENES:
    hits = np.where(gene_names == g)[0]
    if len(hits):
        targets[int(hits[0])] = g
        P("  %s -> var %d" % (g, hits[0]))
    else:
        P("  %s NOT FOUND" % g)

obs = h5["obs"]
def read_cat(col):
    grp = obs[col]
    codes = grp["codes"][()]
    cats = grp["categories"][()].astype(str)
    return cats, cats[codes]

annots = {}
for col in ["cell_type", "cell_state", "region", "modality"]:
    if col in obs:
        cats, vals = read_cat(col)
        annots[col] = (cats, vals)
        P("obs %s: %d cats" % (col, len(cats)))

indptr = X["indptr"][()]
n_obs = indptr.shape[0] - 1
P("n_obs=%d, nnz=%d" % (n_obs, X["data"].shape[0]))

tgt_idx = np.array(sorted(targets.keys()), dtype=np.int64)
lookup = {int(t): targets[t] for t in tgt_idx}
tset = np.zeros(len(gene_names), dtype=bool)
tset[tgt_idx] = True

def summarize(colname):
    cats, vals = annots[colname]
    ncat = len(cats)
    # per-gene accumulators: mean over lognorm -> sum/n ; pct>0
    acc = {g: np.zeros(ncat) for g in targets.values()}
    cnt = {g: np.zeros(ncat) for g in targets.values()}
    celln = np.zeros(ncat, dtype=np.int64)
    cidx = np.searchsorted(np.arange(ncat), np.zeros(0))  # placeholder
    codes = np.array([np.where(cats == v)[0][0] for v in []])  # placeholder
    # 直接用 vals 的类别索引
    cat_to_i = {c: i for i, c in enumerate(cats)}
    vcode = np.array([cat_to_i[v] for v in vals.tolist()])
    CHUNK = 200000
    for r0 in range(0, n_obs, CHUNK):
        r1 = min(r0 + CHUNK, n_obs)
        s0, s1 = int(indptr[r0]), int(indptr[r1])
        idx = X["indices"][s0:s1]
        dat = X["data"][s0:s1]
        mask = tset[idx]
        if mask.any():
            rows_local = np.repeat(
                np.arange(r1 - r0),
                (indptr[r0 + 1:r1 + 1] - indptr[r0:r1]).astype(np.int64))
            hit_rows = rows_local[mask]
            hit_genes = idx[mask]
            hit_vals = dat[mask]
            celltype = vcode[r0 + hit_rows]
            for gi in tgt_idx:
                g = lookup[int(gi)]
                m = hit_genes == gi
                if not m.any():
                    continue
                np.add.at(acc[g], celltype[m], hit_vals[m])
                np.add.at(cnt[g], celltype[m], 1)
        celln += np.bincount(vcode[r0:r1], minlength=ncat)
    P("")
    P("=== per-%s: mean_lognorm | pct_expressing ===" % colname)
    rows = ["%s,n_cells,%s" % (colname, ",".join(
        "%s_mean,%s_pct" % (g, g) for g in targets.values()))]
    for i, c in enumerate(cats):
        if celln[i] == 0:
            continue
        vals_row = []
        for g in targets.values():
            mean = acc[g][i] / celln[i]
            pct = 100.0 * cnt[g][i] / celln[i]
            vals_row += ["%.4f" % mean, "%.2f" % pct]
        rows.append("%s,%d,%s" % (c, celln[i], ",".join(vals_row)))
    fn = os.path.join(W, "hca_np_%s.csv" % colname)
    open(fn, "w", encoding="utf-8").write("\n".join(rows))
    P("saved " + fn)
    for g in targets.values():
        means = acc[g] / np.maximum(celln, 1)
        order = np.argsort(-means)[:3]
        P("  %-6s top: %s" % (g, "; ".join(
            "%s(%.3f)" % (cats[i], means[i]) for i in order if celln[i] > 0)))

summarize("cell_type")
summarize("cell_state")
summarize("region")

open(os.path.join(W, "step17e_output.txt"), "w", encoding="utf-8").write("\n".join(out))
print("DONE")
