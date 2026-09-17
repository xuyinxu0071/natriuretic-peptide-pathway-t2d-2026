# -*- coding: utf-8 -*-
"""
Step 16g: M2 两步中介 MR —— NP(NT-proBNP/NPPA eQTL) -> 中介(SBP/DBP/BMI/eGFR/T2D) -> AF
================================================================================
Step A: NP 暴露 cis 工具(p<5e-8, clump r2<0.01) -> 中介 GWAS (nppa locus 已拉取)
        T2D 额外拉取工具 SNP 的 GCST90018926 关联
Step B: 中介 tophits(排除 MHC 与 NPPA/NPR3 位点 ±500kb) -> AF
间接效应 = b1*b2, delta 法 SE; 结果 -> step16g_mediation.csv + step16g_output.txt
"""
import json, math, os, time, urllib.error, urllib.request
import numpy as np

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
TOKEN = open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\.opengwas_token.txt").read().strip()
API = "https://api.opengwas.io/api"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

out = []
def P(s):
    print(s, flush=True); out.append(str(s))

def api_post(path, payload, retries=6):
    for i in range(retries):
        try:
            req = urllib.request.Request(
                API + path, data=json.dumps(payload).encode("utf-8"),
                headers={"Authorization": "Bearer " + TOKEN,
                         "Content-Type": "application/json"})
            return json.load(urllib.request.urlopen(req, timeout=300))
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and i < retries - 1:
                time.sleep(10 * (i + 1) if e.code == 429 else 10); continue
            P("[HTTP %d] %s" % (e.code, e.read().decode("utf-8", "ignore")[:150]))
            return None
        except Exception as e:
            if i < retries - 1:
                time.sleep(15); continue
            P("[ERR] %s" % e); return None

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

def ivw(betas, ses):
    k = len(betas)
    w = 1.0 / np.asarray(ses) ** 2
    b = np.asarray(betas)
    bf = np.sum(w * b) / np.sum(w)
    se_f = 1.0 / math.sqrt(np.sum(w))
    Q = float(np.sum(w * (b - bf) ** 2))
    if k > 1:
        tau2 = max(0.0, (Q - (k - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w)))
        w_re = 1.0 / (np.asarray(ses) ** 2 + tau2)
        br = np.sum(w_re * b) / np.sum(w_re)
        se_r = 1.0 / math.sqrt(np.sum(w_re))
    else:
        tau2, br, se_r = 0.0, bf, se_f
    return bf, se_f, br, se_r, Q, k - 1, tau2

def p2(z):
    return math.erfc(abs(z) / math.sqrt(2))

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

# ---------- NP 暴露工具 ----------
NP_EXPOSURES = {
    "SCALLOP_NTproBNP": ("nppa", "ebi-a-GCST90012082"),
    "INTERVAL_NTproBNP": ("nppa", "prot-a-2078"),
    "NPPA_eQTL": ("nppa", "eqtl-a-ENSG00000175206"),
}
GWS = 5e-8
CLUMP_R2 = 0.01

idx, R = load_ld("nppa")
np_ivs = {}
for name, (locus, dsid) in NP_EXPOSURES.items():
    assoc = json.load(open(os.path.join(WORK, "%s__%s_assoc.json" % (locus, dsid))))
    cands = [(rs, a) for rs, a in assoc.items()
             if a.get("p") is not None and a["p"] < GWS and rs in idx]
    cands.sort(key=lambda x: x[1]["p"])
    kept = []
    for rs, a in cands:
        if all(R[rs][idx[k]] ** 2 <= CLUMP_R2 for k in kept):
            kept.append(rs)
    np_ivs[name] = [dict(rsid=rs, **assoc[rs]) for rs in kept]
    P("NP exposure %s: %d instruments %s" %
      (name, len(kept), [x["rsid"] for x in np_ivs[name]]))

# ---------- Step A: NP -> mediator ----------
MEDIATORS = {
    "SBP":  ("nppa__ieu-b-38_assoc.json", "mmHg"),
    "DBP":  ("nppa__ieu-b-39_assoc.json", "mmHg"),
    "BMI":  ("nppa__ieu-b-40_assoc.json", "kg/m2"),
    "eGFR": ("nppa__ebi-a-GCST007344_assoc.json", "log eGFR"),
    "T2D":  (None, "log OR"),  # 特殊: 需在线拉取
}

def stepA(exposure, medname, medfile):
    ivs = np_ivs[exposure]
    if medfile:
        med = json.load(open(os.path.join(WORK, medfile)))
    else:
        med = {}
        # 拉取 T2D (GCST90018926) at instruments
        for i in range(0, len(ivs), 64):
            r = api_post("/associations",
                         {"variant": [x["rsid"] for x in ivs[i:i + 64]],
                          "id": ["ebi-a-GCST90018926"]})
            if r:
                for a in r:
                    if a.get("id") == "ebi-a-GCST90018926" and a.get("beta") is not None:
                        med[a["rsid"]] = a
            time.sleep(0.6)
    rows = []
    for iv in ivs:
        rec = med.get(iv["rsid"])
        if rec is None:
            continue
        h = harmonize_iv(iv, rec)
        if h is None:
            continue
        b_out, se_out = h
        if iv["beta"] == 0:
            continue
        rows.append(dict(rsid=iv["rsid"], b_exp=iv["beta"], se_exp=iv["se"],
                         b_out=b_out, se_out=se_out,
                         ratio=b_out / iv["beta"],
                         se_ratio=math.sqrt(
                             (se_out ** 2 * iv["beta"] ** 2 +
                              iv["se"] ** 2 * b_out ** 2) / iv["beta"] ** 4)))
    if len(rows) < 2:
        if len(rows) == 1:
            r0 = rows[0]
            return r0["ratio"], r0["se_ratio"], 1, r0
        return None, None, 0, None
    bf, se_f, br, se_r, Q, df_q, tau2 = ivw(
        [r["ratio"] for r in rows], [r["se_ratio"] for r in rows])
    return br, se_r, len(rows), dict(Q=Q, Q_p=chi2_sf(Q, df_q))

# ---------- Step B: mediator -> AF ----------
MHC = ("6", 25000000, 35000000)
EXCL_LOCI = [("1", 11905974), ("5", 32714270)]  # NPPA/NPR3 ±500kb

def stepB(medname):
    tpf = os.path.join(WORK, "step16c_tophits_%s.json" % medname)
    aff = os.path.join(WORK, "step16c_AF__%s.json" % medname)
    if not os.path.exists(tpf) or not os.path.exists(aff):
        return None, None, 0, "tophits/AF file missing"
    tp = json.load(open(tpf))
    af = json.load(open(aff))
    rows = []
    excluded = 0
    for iv in tp:
        c = str(iv.get("chr"))
        pos = iv.get("position") or iv.get("pos")
        if pos is None:
            excluded += 1; continue
        if c == MHC[0] and MHC[1] <= pos <= MHC[2]:
            excluded += 1; continue
        if any(c == lc[0] and abs(pos - lc[1]) <= 500000 for lc in EXCL_LOCI):
            excluded += 1; continue
        rec = af.get(iv["rsid"])
        if rec is None:
            continue
        h = harmonize_iv(iv, rec)
        if h is None:
            continue
        b_out, se_out = h
        if iv["beta"] == 0:
            continue
        rows.append(dict(rsid=iv["rsid"], b_exp=iv["beta"], se_exp=iv["se"],
                         b_out=b_out, se_out=se_out,
                         ratio=b_out / iv["beta"],
                         se_ratio=math.sqrt(
                             (se_out ** 2 * iv["beta"] ** 2 +
                              iv["se"] ** 2 * b_out ** 2) / iv["beta"] ** 4)))
    if len(rows) < 3:
        return None, None, len(rows), excluded
    bf, se_f, br, se_r, Q, df_q, tau2 = ivw(
        [r["ratio"] for r in rows], [r["se_ratio"] for r in rows])
    return br, se_r, len(rows), dict(excluded=excluded, Q=Q,
                                     Q_p=chi2_sf(Q, df_q))

results = []
P("")
P("=== Step B: mediator -> AF (IVW RE) ===")
B = {}
for medname in MEDIATORS:
    b, se, k, info = stepB(medname)
    B[medname] = (b, se)
    if b is None:
        P("  %s -> AF: k=%s insufficient (%s)" % (medname, k, info))
    else:
        extra = ""
        if isinstance(info, dict) and "Q_p" in info:
            extra = "(Q_p=%.3g, excl=%d)" % (info["Q_p"], info.get("excluded", 0))
        elif isinstance(info, str):
            extra = info
        else:
            extra = "(excl %d)" % info
        P("  %s -> AF: k=%d b=%.4g se=%.4g z=%.2f p=%.3g %s" %
          (medname, k, b, se, b / se, p2(b / se), extra))

P("")
P("=== Step A: NP -> mediator (IVW RE) ===")
for exposure in NP_EXPOSURES:
    for medname, (medfile, unit) in MEDIATORS.items():
        b1, se1, k, info = stepA(exposure, medname, medfile)
        if b1 is None:
            P("  %s -> %s: k=%s insufficient" % (exposure, medname, k))
            continue
        b2, se2 = B.get(medname, (None, None))
        P("  %s -> %s: k=%d b1=%.4g se=%.4g p=%.3g%s" %
          (exposure, medname, k, b1, se1, p2(b1 / se1),
           (" (Q_p=%.3g)" % info["Q_p"])
           if isinstance(info, dict) and info.get("Q_p") is not None else ""))
        if b2 is None:
            continue
        # 间接效应
        bind = b1 * b2
        se_ind = math.sqrt(b2 ** 2 * se1 ** 2 + b1 ** 2 * se2 ** 2)
        P("    INDIRECT %s via %s: b=%.5g se=%.5g z=%.2f p=%.3g" %
          (exposure, medname, bind, se_ind, bind / se_ind, p2(bind / se_ind)))
        results.append(dict(exposure=exposure, mediator=medname,
                            b1=b1, se1=se1, p1=p2(b1 / se1),
                            k_A=k,
                            b2=b2, se2=se2, p2_=p2(b2 / se2),
                            b_indirect=bind, se_indirect=se_ind,
                            p_indirect=p2(bind / se_ind)))

import csv
with open(os.path.join(WORK, "step16g_mediation.csv"), "w", newline="",
          encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["exposure", "mediator", "b1", "se1", "p1",
                                      "k_A", "b2", "se2", "p2_", "b_indirect",
                                      "se_indirect", "p_indirect"])
    w.writeheader()
    for r in results:
        w.writerow({k2: (("%.6g" % v) if isinstance(v, float) else v)
                    for k2, v in r.items()})

open(os.path.join(WORK, "step16g_output.txt"), "w", encoding="utf-8").write("\n".join(out))
print("DONE")
