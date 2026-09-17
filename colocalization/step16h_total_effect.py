# -*- coding: utf-8 -*-
"""step16h: NP 暴露 -> AF 总效应 (IVW RE), 用于中介比例计算"""
import json, math, os
import numpy as np

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
out = []
def P(s):
    print(s, flush=True); out.append(str(s))

def harmonize_iv(iv, orec):
    ea_i, nea_i = iv["ea"].upper(), iv["nea"].upper()
    ea_o, nea_o = (orec.get("ea") or "").upper(), (orec.get("nea") or "").upper()
    if len(ea_o) != 1 or len(nea_o) != 1:
        return None
    if {ea_i, nea_i} in ({"A", "T"}, {"C", "G"}):
        return None
    if ea_o == ea_i and nea_o == nea_i:
        return orec["beta"], orec["se"]
    if ea_o == nea_i and nea_o == ea_i:
        return -orec["beta"], orec["se"]
    if ea_o == COMPLEMENT.get(ea_i) and nea_o == COMPLEMENT.get(nea_i):
        return orec["beta"], orec["se"]
    if ea_o == COMPLEMENT.get(nea_i) and nea_o == COMPLEMENT.get(ea_i):
        return -orec["beta"], orec["se"]
    return None

def p2(z):
    return math.erfc(abs(z) / math.sqrt(2))

NP_EXPOSURES = {
    "SCALLOP_NTproBNP": "ebi-a-GCST90012082",
    "INTERVAL_NTproBNP": "prot-a-2078",
    "NPPA_eQTL": "eqtl-a-ENSG00000175206",
}
af = json.load(open(os.path.join(WORK, "nppa__ebi-a-GCST006061_assoc.json")))

# 与 step16g 相同的工具选择
def load_ld(locus):
    path = os.path.join(WORK, "ld_%s.csv" % locus)
    with open(path) as f:
        hdr = f.readline().strip().split(",")
        idx = {r: i for i, r in enumerate(hdr)}
        R = {}
        for line in f:
            p = line.strip().split(",")
            R[p[0]] = [float(x) for x in p[1:]]
    return idx, R

idx, R = load_ld("nppa")
for name, dsid in NP_EXPOSURES.items():
    assoc = json.load(open(os.path.join(WORK, "nppa__%s_assoc.json" % dsid)))
    cands = [(rs, a) for rs, a in assoc.items()
             if a.get("p") is not None and a["p"] < 5e-8 and rs in idx]
    cands.sort(key=lambda x: x[1]["p"])
    kept = []
    for rs, a in cands:
        if all(R[rs][idx[k]] ** 2 <= 0.01 for k in kept):
            kept.append(rs)
    rows = []
    for rs in kept:
        orec = af.get(rs)
        if orec is None:
            continue
        h = harmonize_iv(assoc[rs], orec)
        if h is None or assoc[rs]["beta"] == 0:
            continue
        b_out, se_out = h
        b_exp = assoc[rs]["beta"]
        rows.append((b_out / b_exp, math.sqrt(
            (se_out ** 2 * b_exp ** 2 + assoc[rs]["se"] ** 2 * b_out ** 2)
            / b_exp ** 4)))
    if not rows:
        P("%s -> AF: no rows" % name); continue
    w = 1.0 / np.array([r[1] for r in rows]) ** 2
    b = np.array([r[0] for r in rows])
    bf = np.sum(w * b) / np.sum(w)
    Q = float(np.sum(w * (b - bf) ** 2))
    k = len(rows)
    if k > 1:
        tau2 = max(0.0, (Q - (k - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w)))
        w_re = 1.0 / (np.array([r[1] for r in rows]) ** 2 + tau2)
        br = np.sum(w_re * b) / np.sum(w_re)
        se_r = 1.0 / math.sqrt(np.sum(w_re))
    else:
        br, se_r = bf, rows[0][1]
    P("%s -> AF TOTAL: k=%d b=%.4g se=%.4g z=%.2f p=%.3g" %
      (name, k, br, se_r, br / se_r, p2(br / se_r)))

open(os.path.join(WORK, "step16h_output.txt"), "w", encoding="utf-8").write("\n".join(out))
print("DONE")
