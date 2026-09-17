# -*- coding: utf-8 -*-
"""
mediation_two_step_mr.py — 两步中介 MR（虚拟通路阻断）
=========================================================
Step 1: NPR3 rs1421811 (G = 降压等位基因) -> SBP (ieu-b-38, Evangelou 2018, N=757,601)  [缓存已有]
Step 2: SBP GWAS 工具变量 (ieu-b-38, p<5e-8, clumped) -> MI/CAD (IVW)
中介效应 = beta1 * beta2 ; 中介比例 = (beta1*beta2)/beta_total
beta_total: rs1421811 -> MI (finn-b-I9_MI lnOR=-0.0388) / CAD (finn-b-I9_CHD -0.0274)
数据源: OpenGWAS API 4.0 (真实数据)
"""
import urllib.request, urllib.error, json, time, math, csv

DIR = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"
import os as _os
_TOK_PATH = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".opengwas_token.txt")
if _os.environ.get("OPENGWAS_TOKEN"):
    TOK = _os.environ["OPENGWAS_TOKEN"].strip()
elif _os.path.isfile(_TOK_PATH):
    TOK = open(_TOK_PATH).read().strip()
else:
    raise RuntimeError("Set env var OPENGWAS_TOKEN or place a gitignored .opengwas_token.txt next to this script to call the OpenGWAS API.")


def post(path, payload, retries=5):
    for i in range(retries):
        try:
            req = urllib.request.Request("https://api.opengwas.io/api" + path,
                data=json.dumps(payload).encode(),
                headers={"Authorization": "Bearer " + TOK, "Content-Type": "application/json"},
                method="POST")
            return json.load(urllib.request.urlopen(req, timeout=120))
        except urllib.error.HTTPError as e:
            if e.code == 429:
                w = 30 * (i + 1); print(f"  429, wait {w}s"); time.sleep(w); continue
            print("  HTTP", e.code, e.read().decode("utf-8", "ignore")[:200]); return None
        except Exception as e:
            print("  err", e); time.sleep(8)
    return None

# ---------- Step 1 (缓存) ----------
phewas = json.load(open(DIR + r"\phewas_rs1421811.json"))
sbp = next(a for a in phewas if a["id"] == "ieu-b-38")
B1 = sbp["beta"]; SE1 = sbp["se"]  # mmHg per G allele
print(f"Step1  rs1421811(G)->SBP: beta={B1:.4f} mmHg (SE {SE1:.4f}), p={sbp['p']:.2e}")

# 总效应（缓存）
cache = json.load(open(DIR + r"\replication_assoc_v2.json"))
TOTAL = {
    "MI (FinnGen)":  (cache["rs1421811|finn-b-I9_MI"]["beta"],  cache["rs1421811|finn-b-I9_MI"]["se"]),
    "CHD (FinnGen)": (cache["rs1421811|finn-b-I9_CHD"]["beta"], cache["rs1421811|finn-b-I9_CHD"]["se"]),
}

# ---------- Step 2: SBP IVs ----------
print("\nFetching SBP tophits (ieu-b-38, p<5e-8, clumped)...")
ivs = post("/tophits", {"id": ["ieu-b-38"], "pval": 5e-8, "clump": 1})
print(f"SBP instruments: {len(ivs)} SNPs")
json.dump(ivs, open(DIR + r"\sbp_ivs.json", "w"))

OUTCOMES = {
    "MI (FinnGen)": "finn-b-I9_MI",
    "CAD (vanderHarst)": "ebi-a-GCST005195",
    "CHD (FinnGen)": "finn-b-I9_CHD",
}
snps = [iv["rsid"] for iv in ivs]

# 已有缓存？
try:
    assoc_cache = json.load(open(DIR + r"\mediation_assoc_cache.json"))
except FileNotFoundError:
    assoc_cache = {}

results = {}
for label, oid in OUTCOMES.items():
    pairs = []
    miss = [s for s in snps if f"{s}|{oid}" not in assoc_cache]
    if miss:
        for i in range(0, len(miss), 60):   # API 限制: N(id)*N(variant)<=64
            batch = miss[i:i+60]
            r = post("/associations", {"variant": batch, "id": [oid]})
            time.sleep(2)
            if r:
                for a in r:
                    assoc_cache[f"{a.get('rsid')}|{oid}"] = a
            print(f"  {label}: batch {i//60+1}/{(len(miss)+59)//60} fetched {len(r) if r else 0}")
        json.dump(assoc_cache, open(DIR + r"\mediation_assoc_cache.json", "w"))
    for iv in ivs:
        s = iv["rsid"]
        a = assoc_cache.get(f"{s}|{oid}")
        if not a or a.get("beta") is None:
            continue
        # 等位基因对齐（暴露 ieu-b-38 的 ea/nea vs 结局）
        ea, nea, eaf = iv.get("ea"), iv.get("nea"), iv.get("eaf")
        oea, onea, oeaf = a.get("ea"), a.get("nea"), a.get("eaf")
        if not oea or not onea:
            continue
        sign = 1
        if (ea, nea) == (oea, onea):
            sign = 1
        elif (ea, nea) == (onea, oea):
            sign = -1
        else:
            comp = {"A":"T","T":"A","C":"G","G":"C"}
            if comp.get(ea)==oea and comp.get(nea)==onea:
                # 回文 SNP：一律剔除（保守处理）
                continue
            else:
                continue
        bx = iv["beta"] * sign            # SBP (mmHg) per ea
        sx = iv.get("se")
        by = a["beta"] * sign              # lnOR per ea
        sy = a.get("se")
        if sx is None or sy is None or bx == 0:
            continue
        pairs.append((s, bx, sx, by, sy))
    # IVW (random-effects 退化为固定效应当无异质性)
    if len(pairs) >= 2:
        W = [1.0/(p[4]**2) for p in pairs]
        b_ivw = sum(w*p[3] for w,p in zip(W,pairs)) / sum(W)
        se_ivw = math.sqrt(1.0/sum(W))
        Q = sum(w*(p[3]-b_ivw)**2 for w,p in zip(W,pairs))
        dof = len(pairs)-1
        # 残差方差
        s2 = max(1.0, Q/dof) if dof>0 else 1.0
        se_re = se_ivw*math.sqrt(s2)
        z = b_ivw/se_re
        from math import erf
        pval = 2*(1-0.5*(1+erf(abs(z)/math.sqrt(2))))
        results[label] = {"n_iv": len(pairs), "b": b_ivw, "se": se_re, "p": pval,
                          "Q": Q, "dof": dof}
        print(f"Step2  SBP->{label}: lnOR={b_ivw:.6f}/mmHg (SE {se_re:.6f}) p={pval:.2e} n={len(pairs)} Q/dof={Q/dof:.2f}")
    json.dump(assoc_cache, open(DIR + r"\mediation_assoc_cache.json", "w"))

# ---------- 中介效应 ----------
print("\n===== 中介分解（乘积法）=====")
rows = []
for label in results:
    if label not in TOTAL:
        key2 = "CHD (FinnGen)" if "CHD" in label or "CAD" in label else label
    tot_label = "MI (FinnGen)" if "MI" in label else "CHD (FinnGen)"
    b_tot, se_tot = TOTAL[tot_label]
    r = results[label]
    b_med = B1 * r["b"]                      # lnOR 中介
    se_med = abs(b_med)*math.sqrt((SE1/B1)**2 + (r["se"]/r["b"])**2)  # delta
    prop = b_med / b_tot
    se_prop = abs(prop)*math.sqrt((se_med/b_med)**2 + (se_tot/b_tot)**2)
    rows.append({"outcome": label, "n_iv_step2": r["n_iv"],
        "beta1_rs1421811_SBP_mmHg": round(B1,4),
        "beta2_SBP_outcome_lnOR_per_mmHg": round(r["b"],6),
        "mediated_lnOR": round(b_med,5), "se_mediated": round(se_med,5),
        "total_lnOR": round(b_tot,5),
        "proportion_mediated": round(prop,3), "se_proportion": round(se_prop,3),
        "Q_step2": round(r["Q"],1), "Q_dof": r["dof"]})
    print(f"{label}: 中介 lnOR={b_med:.5f} ({b_med-1.96*se_med:.5f},{b_med+1.96*se_med:.5f}) | "
          f"总效应 {b_tot:.5f} | 中介比例={prop*100:.1f}% (SE {se_prop*100:.1f}%)")

with open(DIR + r"\mediation_results.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
json.dump(results, open(DIR + r"\mediation_step2_results.json", "w"), indent=1)
print("\nsaved: mediation_results.csv / mediation_step2_results.json / sbp_ivs.json")
