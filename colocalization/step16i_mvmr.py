# step16i: MVMR (NT-proBNP + SBP -> AF) + defensive sensitivity
# Instruments: NPPA-locus cis SNPs, NT-proBNP p<5e-6 (SCALLOP), LD clump r2<0.2.
# MVMR-IVW: weighted least squares by ~ bx_NP + bx_SBP, w = 1/se_y^2 (no intercept).
import json, math, os, csv

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
LOG = os.path.join(WORK, "step16i_output.txt")
_lf = open(LOG, "w", encoding="utf-8")
def P(*a):
    s = " ".join(str(x) for x in a)
    _lf.write(s + "\n"); _lf.flush()

def load_efile(path):
    rows = {}
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 10: continue
            try:
                rows[p[1]] = dict(a1=p[4], a2=p[5], freq=float(p[6]),
                                  b=float(p[7]), se=float(p[8]), pval=float(p[9]))
            except ValueError: continue
    return rows

def load_assoc(path):
    d = json.load(open(path))
    out = {}
    for rs, e in d.items():
        out[rs] = dict(a1=e["ea"], a2=e["nea"], b=e["beta"], se=e["se"],
                       pval=e.get("p"), freq=e.get("eaf"))
    return out

def harmon(exp, med, out):
    """Align med & out betas to the exposure effect allele (efile a1)."""
    snps = []
    for rs, x in exp.items():
        if rs not in med or rs not in out: continue
        m, o = med[rs], out[rs]
        # mediator
        if m["a1"] == x["a1"] and m["a2"] == x["a2"]: bm = m["b"]
        elif m["a1"] == x["a2"] and m["a2"] == x["a1"]: bm = -m["b"]
        else: continue
        # outcome
        if o["a1"] == x["a1"] and o["a2"] == x["a2"]: by = o["b"]
        elif o["a1"] == x["a2"] and o["a2"] == x["a1"]: by = -o["b"]
        else: continue
        snps.append(dict(rsid=rs, bx=x["b"], sex=x["se"], pexp=x["pval"],
                         bm=bm, sem=m["se"], by=by, sey=o["se"]))
    return snps

def clump(snps, R, hdr, pthresh, r2thresh):
    idx = {rs: i for i, rs in enumerate(hdr)}
    cand = sorted([s for s in snps if s["pexp"] < pthresh and s["rsid"] in idx],
                  key=lambda s: s["pexp"])
    kept = []
    for s in cand:
        i = idx[s["rsid"]]
        ok = True
        for t in kept:
            j = idx[t["rsid"]]
            r = R[i][j]
            if r * r >= r2thresh: ok = False; break
        if ok: kept.append(s)
    return kept

def mvmr_ivw(snps):
    """WLS: by = b1*bx + b2*bm, weights 1/sey^2. Exact 2x2 solve."""
    k = len(snps)
    S11 = sum(s["bx"]**2 / s["sey"]**2 for s in snps)
    S12 = sum(s["bx"] * s["bm"] / s["sey"]**2 for s in snps)
    S22 = sum(s["bm"]**2 / s["sey"]**2 for s in snps)
    T1 = sum(s["bx"] * s["by"] / s["sey"]**2 for s in snps)
    T2 = sum(s["bm"] * s["by"] / s["sey"]**2 for s in snps)
    det = S11 * S22 - S12 * S12
    b1 = (T1 * S22 - T2 * S12) / det
    b2 = (T2 * S11 - T1 * S12) / det
    se1 = math.sqrt(S22 / det)
    se2 = math.sqrt(S11 / det)
    # residual Q (df = k-2)
    Q = sum((s["by"] - b1 * s["bx"] - b2 * s["bm"])**2 / s["sey"]**2 for s in snps)
    # conditional F (approximate): mean bx^2/sex^2 etc.
    F1 = sum(s["bx"]**2 / s["sex"]**2 for s in snps) / k
    F2 = sum(s["bm"]**2 / s["sem"]**2 for s in snps) / k
    return dict(k=k, b1=b1, se1=se1, b2=b2, se2=se2, Q=Q, df=k-2, F1=F1, F2=F2)

def p2(z):
    return math.erfc(abs(z) / math.sqrt(2))

def chi2_sf(x, k):
    if x <= 0: return 1.0
    a = k / 2.0
    z, t = x / 2.0, a
    if z < 1.0 or z < t:
        term = 1.0 / t; s = term
        for n in range(1, 500):
            term *= z / (t + n); s += term
            if abs(term) < 1e-14 * abs(s): break
        return max(0.0, min(1.0, 1.0 - s * math.exp(-z + t * math.log(z) - math.lgamma(t))))
    b_ = z + 1.0 - t; c = 1.0 / 1e-300; d = 1.0 / b_; h = d
    for i in range(1, 500):
        an = -i * (i - t); b_ += 2.0
        d = an * d + b_; d = 1e-300 if abs(d) < 1e-300 else d
        c = b_ + an / c; c = 1e-300 if abs(c) < 1e-300 else c
        d = 1.0 / d; de = d * c; h *= de
        if abs(de - 1.0) < 1e-14: break
    return max(0.0, min(1.0, math.exp(-z + t * math.log(z) - math.lgamma(t)) * h))

# ---- LD ----
hdr = open(os.path.join(WORK, "ld_nppa.csv")).readline().strip().split(",")
R = []
with open(os.path.join(WORK, "ld_nppa.csv")) as f:
    f.readline()
    for line in f:
        R.append([float(x) for x in line.strip().split(",")[1:]])

exp = load_efile(os.path.join(WORK, "SCALLOP.efile"))
med = load_assoc(os.path.join(WORK, "nppa__ieu-b-38_assoc.json"))
out = load_assoc(os.path.join(WORK, "nppa__ebi-a-GCST006061_assoc.json"))

snps = harmon(exp, med, out)
P("harmonized SNPs across SCALLOP/SBP/AF: %d" % len(snps))

for pt, rt in [(5e-6, 0.2), (5e-6, 0.1), (1e-4, 0.2)]:
    ivs = clump(snps, R, hdr, pt, rt)
    P("\n--- instruments p<%g r2<%g: k=%d [%s] ---" %
      (pt, rt, len(ivs), ", ".join(s["rsid"] for s in ivs)))
    if len(ivs) >= 3:
        r = mvmr_ivw(ivs)
        P("MVMR-IVW: k=%d df=%d" % (r["k"], r["df"]))
        P("  NT-proBNP direct (adj SBP): b=%.5f se=%.5f z=%.2f p=%.4g"
          % (r["b1"], r["se1"], r["b1"] / r["se1"], p2(r["b1"] / r["se1"])))
        P("  SBP direct (adj NT-proBNP): b=%.5f se=%.5f z=%.2f p=%.4g"
          % (r["b2"], r["se2"], r["b2"] / r["se2"], p2(r["b2"] / r["se2"])))
        P("  Q=%.2f (df=%d, p=%.3g); mean F: NP=%.1f SBP=%.1f"
          % (r["Q"], r["df"], chi2_sf(r["Q"], r["df"]), r["F1"], r["F2"]))
        for s in ivs:
            P("    %s bx=%.4f(se %.4f) bm=%.3f(se %.3f) by=%.4f(se %.4f)"
              % (s["rsid"], s["bx"], s["sex"], s["bm"], s["sem"], s["by"], s["sey"]))
    else:
        P("  insufficient instruments")

# univariable for comparison (same instrument set, p<5e-6 r2<0.2)
ivs = clump(snps, R, hdr, 5e-6, 0.2)
if len(ivs) >= 2:
    w = [1 / s["sey"]**2 for s in ivs]
    b_uv = sum(wi * s["by"] / s["bx"] for wi, s in zip(w, ivs)) / sum(wi / s["bx"]**2 * s["bx"]**2 for wi, s in zip(w, ivs))
    # simpler: ratio-of-sums IVW (by~bx through origin)
    b_uv = sum(s["bx"] * s["by"] / s["sey"]**2 for s in ivs) / sum(s["bx"]**2 / s["sey"]**2 for s in ivs)
    se_uv = math.sqrt(1 / sum(s["bx"]**2 / s["sey"]**2 for s in ivs))
    P("\nUnivariable IVW (same IVs, NT-proBNP only): b=%.5f se=%.5f z=%.2f p=%.4g"
      % (b_uv, se_uv, b_uv / se_uv, p2(b_uv / se_uv)))

_lf.close()
print("done")
