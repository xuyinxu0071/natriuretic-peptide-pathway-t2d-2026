# -*- coding: utf-8 -*-
"""
Step 16f: M1 cis-MR 分析 —— NP 加工酶级联暴露 x 6 结局
================================================================================
每 SNP Wald ratio + IVW(FE/RE, DL tau2) + Q 异质性; harmonize 沿用 step12c 规则
(回文丢弃, 等位翻转/互补对齐)。结果 -> step16f_cascade_mr.csv + step16f_output.txt
"""
import csv, json, math, os
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

def chi2_sf(x, df):
    if x <= 0:
        return 1.0
    a, xx = df / 2.0, x / 2.0
    if xx < a + 1.0:
        term = 1.0 / a; s = term; n = 0
        while abs(term) > abs(s) * 1e-14 and n < 1000:
            n += 1; term *= xx / (a + n); s += term
        return 1.0 - s * math.exp(-xx + a * math.log(xx) - math.lgamma(a))
    tiny = 1e-300; b = xx + 1 - a; c = 1 / tiny; d = 1 / b; h = d
    for i in range(1, 1000):
        an = -i * (i - a); b += 2
        d = an * d + b; d = tiny if abs(d) < tiny else d
        c = b + an / c; c = tiny if abs(c) < tiny else c
        d = 1 / d; delta = d * c; h *= delta
        if abs(delta - 1) < 1e-14:
            break
    return math.exp(-xx + a * math.log(xx) - math.lgamma(a)) * h

def ivw(beta_xy, se_xy):
    k = len(beta_xy)
    w = 1.0 / np.asarray(se_xy) ** 2
    b = np.asarray(beta_xy)
    bf = np.sum(w * b) / np.sum(w)
    se_f = 1.0 / math.sqrt(np.sum(w))
    Q = float(np.sum(w * (b - bf) ** 2))
    if k > 1:
        tau2 = max(0.0, (Q - (k - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w)))
        w_re = 1.0 / (np.asarray(se_xy) ** 2 + tau2)
        br = np.sum(w_re * b) / np.sum(w_re)
        se_r = 1.0 / math.sqrt(np.sum(w_re))
    else:
        tau2, br, se_r = 0.0, bf, se_f
    return bf, se_f, br, se_r, Q, k - 1, tau2

def p2(z):
    return math.erfc(abs(z) / math.sqrt(2))

EXPOSURES = ["CORIN_eQTL", "FURIN_eQTL", "MME_eQTL", "DPP4_eQTL"]
OUTCOMES = ["AF", "CAD", "CAD2", "HF", "T2D", "T2D_FINN"]

inst = json.load(open(os.path.join(WORK, "step16d_instruments.json")))
rows_csv = []
results = []

for oname in OUTCOMES:
    oass = json.load(open(os.path.join(WORK, "step16d_outcome__%s.json" % oname)))
    for ename in EXPOSURES:
        e = inst[ename]
        ivs = e["ivs"]
        if not ivs:
            P("[%s -> %s] no instruments (underpowered exposure)" % (ename, oname))
            continue
        rows, dropped = [], 0
        for iv in ivs:
            orec = oass.get(iv["rsid"])
            if orec is None:
                dropped += 1; continue
            h = harmonize_iv(iv, orec)
            if h is None:
                dropped += 1; continue
            b_out, se_out = h
            if iv["beta"] == 0:
                dropped += 1; continue
            rows.append(dict(rsid=iv["rsid"], b_exp=iv["beta"], se_exp=iv["se"],
                             b_out=b_out, se_out=se_out,
                             ratio=b_out / iv["beta"],
                             se_ratio=math.sqrt(
                                 (se_out ** 2 * iv["beta"] ** 2 +
                                  iv["se"] ** 2 * b_out ** 2) /
                                 iv["beta"] ** 4)))
        k = len(rows)
        if k == 0:
            P("[%s -> %s] 0 harmonizable IVs" % (ename, oname)); continue
        # lead (最小暴露 p) Wald
        lead = min(rows, key=lambda r: abs(r["b_exp"]) / r["se_exp"])
        z_lead = lead["ratio"] / math.sqrt(
            (lead["se_out"] ** 2 * lead["b_exp"] ** 2 +
             lead["se_exp"] ** 2 * lead["b_out"] ** 2) / lead["b_exp"] ** 4)
        p_lead = p2(z_lead)
        if k >= 2:
            bf, se_f, br, se_r, Q, df_q, tau2 = ivw(
                [r["ratio"] for r in rows], [r["se_ratio"] for r in rows])
            z_ivw = br / se_r
            p_ivw = p2(z_ivw)
            pQ = chi2_sf(Q, df_q) if df_q > 0 else float("nan")
        else:
            br, se_r, z_ivw, p_ivw, Q, pQ = (
                lead["ratio"], math.nan, z_lead, p_lead, float("nan"), float("nan"))
        P("[%s -> %s] k=%d drop=%d | IVW_RE b=%.4g se=%.4g z=%.2f p=%.4g | "
          "lead %s b=%.4g p=%.4g | Q_p=%.3g" %
          (ename, oname, k, dropped, br, se_r, z_ivw, p_ivw,
           lead["rsid"], lead["ratio"], p_lead, pQ))
        results.append(dict(exposure=ename, outcome=oname, k=k, dropped=dropped,
                            b_ivw=br, se_ivw=se_r, z=z_ivw, p=p_ivw,
                            lead=lead["rsid"], b_lead=lead["ratio"], p_lead=p_lead,
                            Q_p=pQ))
        for r in rows:
            rows_csv.append(dict(exposure=ename, outcome=oname, rsid=r["rsid"],
                                 b_exp=r["b_exp"], se_exp=r["se_exp"],
                                 b_out=r["b_out"], se_out=r["se_out"],
                                 ratio=r["ratio"]))

with open(os.path.join(WORK, "step16f_cascade_mr.csv"), "w", newline="",
          encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["exposure", "outcome", "k", "dropped",
                                      "b_ivw", "se_ivw", "z", "p", "lead",
                                      "b_lead", "p_lead", "Q_p"])
    w.writeheader()
    for r in results:
        w.writerow({k2: (("%.6g" % v) if isinstance(v, float) else v)
                    for k2, v in r.items()})
with open(os.path.join(WORK, "step16f_cascade_snp.csv"), "w", newline="",
          encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["exposure", "outcome", "rsid", "b_exp",
                                      "se_exp", "b_out", "se_out", "ratio"])
    w.writeheader()
    for r in rows_csv:
        w.writerow({k2: (v if isinstance(v, str) else "%.6g" % v)
                    for k2, v in r.items()})

# 逐 SNP 明细(校验用)
P("")
P("=== per-SNP detail (exposure z, outcome z, ratio z) ===")
for oname in OUTCOMES:
    for ename in EXPOSURES:
        snps = [r for r in rows_csv
                if r["exposure"] == ename and r["outcome"] == oname]
        for r in snps:
            zx = r["b_exp"] / r["se_exp"]
            zy = r["b_out"] / r["se_out"]
            zr = r["ratio"] / math.sqrt(
                (r["se_out"] ** 2 * r["b_exp"] ** 2 +
                 r["se_exp"] ** 2 * r["b_out"] ** 2) / r["b_exp"] ** 4)
            P("  %-10s %-8s %-12s zx=%7.1f zy=%7.1f zratio=%6.2f ratio=%+.4g"
              % (ename, oname, r["rsid"], zx, zy, zr, r["ratio"]))

# 多重检验阈值
n_tests = len(results)
P("")
P("total MR tests: %d ; Bonferroni alpha = %.4g" % (n_tests, 0.05 / n_tests))

open(os.path.join(WORK, "step16f_output.txt"), "w", encoding="utf-8").write("\n".join(out))
print("DONE")
