# -*- coding: utf-8 -*-
"""
step13q_smr_heidi.py
SMR & HEIDI test — faithful Python port of the OFFICIAL SMR v1.4.2 source code
(JianYang-Lab/SMR, GPL-3.0; Zhu et al. 2016 Nat Genet 48:481-487).

Rationale: the prebuilt Windows binary (1.3.1) crashes with 0xC0000005 access
violation on all analysis invocations in this environment; v1.4.2 has no
Windows build; WSL is blocked by security policy. The implementation below
reproduces, line-by-line, the algorithms in:
  - src_SMR_data.cpp :: SMR test            (lines 4095-4109)
  - src_SMR_data.cpp :: heidi_test_new      (lines 3790-3887)  <- v1.4.2 default
  - src_SMR_data.cpp :: bxy_hetero3         (lines 1292-1365)
  - src_SMR_data.cpp :: est_cov_bxy         (lines 1254-1264)
  - src_SMR_data.cpp :: rm_cor_sbat         (lines 3645-3683)
  - src_SMR_data.cpp :: update_snidx        (lines 3617-3624)
  - src_StatFunc.cpp :: pchisqsum/psadd/psatt/K/Kp/Kpp/Brent (lines 218-352)
Defaults confirmed from SMR.cpp: peqtl-smr 5e-8, peqtl-heidi 1.5654e-3,
ld-prune 0.9, ld-min 0.05, m-hetero 3, opt-hetero 20, MAX_NUM_LD 500,
diff-freq 0.2, new-heidi-mth TRUE.

No sample-overlap correction (default, sampleoverlap=false).

Inputs (all built in prior steps from OpenGWAS + 1000G EUR):
  {exposure}.efile  : ProbeID SNP Chr BP A1 A2 Freq b se p
  ma_{outcome}.ma   : SNP A1 A2 freq b se p n
  ld_{locus}.csv    : 1000G EUR LD correlation matrix (from step6)
"""
import os
import csv
import math

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"

# ----------------------------------------------------------------------------
# SMR defaults (SMR.cpp)
# ----------------------------------------------------------------------------
P_SMR_DEFAULT = 5e-8      # --peqtl-smr
P_SMR_RELAXED = 5e-6      # sensitivity (SMR supports --peqtl-smr)
P_HEIDI       = 1.5654e-3 # --peqtl-heidi
LD_PRUNE_R2   = 0.9       # --ld-prune (upper r^2 bound w.r.t. top SNP)
LD_MIN_R2     = 0.05      # --ld-min
M_HETERO      = 3         # --m-herero (min #SNPs for HEIDI)
OPT_HETERO    = 20        # --opt-herero (max #SNPs used in HEIDI)
MAX_NUM_LD    = 500       # hard cap in source
DIFF_FREQ     = 0.2       # --diff-freq

# ----------------------------------------------------------------------------
# Statistical functions (ported from src_StatFunc.cpp / Numerical Recipes)
# ----------------------------------------------------------------------------
def gammln(xx):
    cof = [76.18009172947146, -86.50532032941677, 24.01409824083091,
           -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5]
    x = xx
    y = xx
    tmp = x + 5.5
    tmp -= (x + 0.5) * math.log(tmp)
    ser = 1.000000000190015
    for j in range(6):
        y += 1
        ser += cof[j] / y
    return -tmp + math.log(2.5066282746310005 * ser / x)

def gser(a, x):
    """Series representation of P(a,x)."""
    ap = a
    summ = 1.0 / a
    delt = summ
    for _ in range(1000):
        ap += 1
        delt *= x / ap
        summ += delt
        if abs(delt) < abs(summ) * 1e-16:
            break
    return summ * math.exp(-x + a * math.log(x) - gammln(a))

def gcf(a, x):
    """Continued fraction for Q(a,x) (Lentz)."""
    tiny = 1e-300
    b = x + 1 - a
    c = 1 / tiny
    d = 1 / b
    h = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1 / d
        delt = d * c
        h *= delt
        if abs(delt - 1) < 1e-16:
            break
    return math.exp(-x + a * math.log(x) - gammln(a)) * h

def gammq(a, x):
    """Upper regularized incomplete gamma Q(a,x)."""
    if x < 0 or a <= 0:
        return float("nan")
    if x == 0:
        return 1.0
    if x < a + 1:
        return 1.0 - gser(a, x)
    return gcf(a, x)

def pchisq(x, df):
    """Chi-square SURVIVAL function (SMR naming: P(chi2_df > x)); real df OK."""
    if x <= 0:
        return 1.0
    return gammq(df / 2.0, x / 2.0)

def chi_val(df, p):
    """Critical value: x s.t. P(chi2_df > x) = p (inverse of pchisq; SMR chi_val).
    Bisection, as the cdfchi inversion in source does numerically."""
    lo, hi = 0.0, 1.0
    while pchisq(hi, df) > p:
        hi *= 2.0
        if hi > 1e12:
            break
    for _ in range(200):
        mid = (lo + hi) / 2
        if pchisq(mid, df) > p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

def pnorm_upper(x):
    """Upper-tail standard normal."""
    return 0.5 * math.erfc(x / math.sqrt(2.0))

# --- pchisqsum (Kuonen 1999 saddlepoint; src_StatFunc.cpp 218-352) -----------
def K(zeta, lam):
    return -0.5 * sum(math.log(1.0 - 2.0 * zeta * l) for l in lam)

def Kp(zeta, lam):
    return sum(l / (1.0 - 2.0 * zeta * l) for l in lam)

def Kpp(zeta, lam):
    return 2.0 * sum(l * l / (1.0 - 2.0 * zeta * l) ** 2 for l in lam)

def brents_kp_min_x(lam, x, lower, upper, tol=1e-8):
    """Brent's method root of Kp(zeta)-x. Faithful port incl. early exits."""
    a, b = lower, upper
    fa, fb = Kp(a, lam) - x, Kp(b, lam) - x
    if fa * fb >= 0:
        return a if fa < fb else b
    if abs(fa) < abs(fb):
        a, b = b, a
        fa, fb = fb, fa
    c, fc = a, fa
    d = 1.7976931348623157e308
    mflag = True
    i = 0
    while fb != 0 and abs(a - b) > tol:
        if fa != fc and fb != fc:
            s = (a * fb * fc / (fa - fb) / (fa - fc)
                 + b * fa * fc / (fb - fa) / (fb - fc)
                 + c * fa * fb / (fc - fa) / (fc - fb))
        else:
            s = b - fb * (b - a) / (fb - fa)
        tmp2 = (3 * a + b) / 4
        if (not ((tmp2 < s < b) or (b < s < tmp2))
                or (mflag and abs(s - b) >= abs(b - c) / 2)
                or (not mflag and abs(s - b) >= abs(c - d) / 2)):
            s = (a + b) / 2
            mflag = True
        else:
            if (mflag and abs(b - c) < tol) or (not mflag and abs(c - d) < tol):
                s = (a + b) / 2
                mflag = True
            else:
                mflag = False
        fs = Kp(s, lam) - x
        d = c
        c, fc = b, fb
        if fa * fs < 0:
            b, fb = s, fs
        else:
            a, fa = s, fs
        if abs(fa) < abs(fb):
            a, b = b, a
            fa, fb = fb, fa
        i += 1
        if i > 1000:
            return upper + 10
    return b

def psadd(x, lam):
    d = max(lam)
    if d <= 0:
        return 2.0
    lam = [l / d for l in lam]
    x = x / d
    m = min(lam)
    if m < 0:
        lmin = 0.499995 / m
    elif x > sum(lam):
        lmin = -0.01
    else:
        lmin = -0.5 * len(lam) / x
    lmax = 0.499995 / max(lam)
    hatzeta = brents_kp_min_x(lam, x, lmin, lmax, 1e-8)
    if hatzeta > lmax + 9:
        return 2.0
    sign = -1.0 if hatzeta < 0 else 1.0
    w = sign * math.sqrt(2 * (hatzeta * x - K(hatzeta, lam)))
    v = hatzeta * math.sqrt(Kpp(hatzeta, lam))
    if abs(hatzeta) < 1e-4:
        return 2.0
    return pnorm_upper(w + math.log(v / w) / w)

def psatt(x, lam):
    ssum = sum(lam)
    if ssum == 0:
        return 2.0
    sq_sum = sum(l * l for l in lam)
    a = sq_sum / ssum
    b = ssum * ssum / sq_sum
    return pchisq(x / a, b)

def pchisqsum(x, lam):
    pval = psadd(x, lam)
    if pval > 1.0:
        pval = psatt(x, lam)
    return pval

# ----------------------------------------------------------------------------
# SMR / HEIDI core (ported from src_SMR_data.cpp)
# ----------------------------------------------------------------------------
def smr_test(bx, sex, by, sey):
    """SMR test at one SNP (source 4095-4109), no sample overlap."""
    b_smr = by / bx
    se_smr = math.sqrt((sey * sey * bx * bx + sex * sex * by * by)
                       / (bx * bx * bx * bx))
    chisq = b_smr * b_smr / (se_smr * se_smr)
    return b_smr, se_smr, chisq, pchisq(chisq, 1)

def update_snidx(zxz, sn_ids, max_snp):
    """update_snidx (3617-3624): if > max_snp keep top-|z| (stable by |z| desc)."""
    if len(sn_ids) > max_snp:
        sn_ids = sorted(sn_ids, key=lambda i: -abs(zxz[i]))[:max_snp]
    return sn_ids

def rm_cor_sbat(R, cutoff, m):
    """rm_cor_sbat (3645-3683): greedy removal of higher-degree member of each
    pair with |R(i,j)| > cutoff. Returns sorted removed indices."""
    edges = []
    for i in range(1, m):
        for j in range(i):
            if abs(R[i][j]) > cutoff:
                edges.append((j, i))
    if not edges:
        return []
    degree = {}
    for a, b in edges:
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1
    removed = set()
    for a, b in edges:
        removed.add(a if degree[a] >= degree[b] else b)
    return sorted(removed)

def est_cov_bxy(zsxz, bxy, seyz, bxz, LD):
    """est_cov_bxy (1254-1264). Vectorized with plain lists->nested loops n<=21."""
    n = len(zsxz)
    cov = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            t1 = seyz[i] * seyz[j] / (bxz[i] * bxz[j])
            t2 = bxy[i] * bxy[j] / (zsxz[i] * zsxz[j])
            cov[i][j] = LD[i][j] * (t1 + t2) - bxy[i] * bxy[j] / (
                zsxz[i] * zsxz[j] * zsxz[i] * zsxz[j])
    return cov

def bxy_hetero3(byz, bxz, seyz, sexz, zsxz, LD):
    """bxy_hetero3 (1292-1365): returns (pdev, nsnp_eigenvalues)."""
    n = len(byz)
    if n <= 1:
        return -9.0, -9
    bxy = [byz[i] / bxz[i] for i in range(n)]
    dev = [0.0] * (n - 1)
    maxid = max(range(n), key=lambda i: abs(zsxz[i]))
    for j in range(maxid):
        dev[j] = bxy[maxid] - bxy[j]
    for j in range(maxid + 1, n):
        dev[j - 1] = bxy[maxid] - bxy[j]

    cov = est_cov_bxy(zsxz, bxy, seyz, bxz, LD)
    tmp1 = cov[maxid][maxid]
    # tmp3 = cov(maxid, non-top)
    tmp3 = [0.0] * (n - 1)
    for i in range(maxid):
        tmp3[i] = cov[maxid][i]
    for i in range(maxid + 1, n):
        tmp3[i - 1] = cov[maxid][i]

    # vdev(i,j) = tmp1 + cov(i,j) - tmp3[i] - tmp3[j]  (non-top submatrix)
    vdev = [[0.0] * (n - 1) for _ in range(n - 1)]
    for i in range(n - 1):
        for j in range(n - 1):
            ci = i if i < maxid else i + 1
            cj = j if j < maxid else j + 1
            vdev[i][j] = tmp1 + cov[ci][cj] - tmp3[i] - tmp3[j] + (
                1e-8 if i == j else 0.0)

    # vardev and chisq_dev
    sum_chisq = 0.0
    for i in range(n - 1):
        ci = i if i < maxid else i + 1
        vardev = tmp1 + cov[ci][ci] - 2 * tmp3[i] + 1e-8
        dev[i] = dev[i] * dev[i] / vardev
        sum_chisq += dev[i]

    # corr_dev = vdev / sqrt(diag x diag); symmetric eigenvalues via Jacobi
    m = n - 1
    diag = [vdev[i][i] for i in range(m)]
    corr = [[vdev[i][j] / math.sqrt(diag[i] * diag[j]) if diag[i] > 0 and diag[j] > 0
             else (1.0 if i == j else 0.0) for j in range(m)] for i in range(m)]
    # enforce exact symmetry & unit diagonal
    for i in range(m):
        corr[i][i] = 1.0
        for j in range(i + 1, m):
            v = 0.5 * (corr[i][j] + corr[j][i])
            corr[i][j] = corr[j][i] = v
    lam = jacobi_eigenvalues(corr)
    pdev = pchisqsum(sum_chisq, lam)
    return pdev, m

def jacobi_eigenvalues(A):
    """Symmetric eigenvalues (Jacobi rotations) — matches SelfAdjointEigenSolver."""
    n = len(A)
    a = [row[:] for row in A]
    for _ in range(100):
        off = 0.0
        for i in range(n - 1):
            for j in range(i + 1, n):
                off += a[i][j] * a[i][j]
        if off < 1e-14:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-300:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = (1 if theta >= 0 else -1) / (abs(theta) + math.sqrt(theta * theta + 1))
                c = 1 / math.sqrt(t * t + 1)
                s = t * c
                for k in range(n):
                    akp = a[k][p]
                    akq = a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk = a[p][k]
                    aqk = a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
    return sorted(a[i][i] for i in range(n))

def heidi_test_new(byz, bxz, seyz, sexz, zsxz, R_full, thresh):
    """heidi_test_new (3790-3887). R_full: LD r matrix over ALL merged SNPs.
    Returns (pdev, nsnp_used, note)."""
    n_all = len(zsxz)
    # 1) select z^2 > threshold + 1e-6
    sn_ids = [i for i in range(n_all) if zsxz[i] * zsxz[i] - thresh > 1e-6]
    if len(sn_ids) < M_HETERO:
        return -9.0, -9, "z2>10 SNPs < m_hetero(3)"
    # 2) cap MAX_NUM_LD by |z|
    sn_ids = update_snidx(zsxz, sn_ids, MAX_NUM_LD)
    maxid0 = max(sn_ids, key=lambda i: abs(zsxz[i]))
    # 3) r^2 filter w.r.t. top: keep 0.05 < r^2 < 0.9 (plus top itself)
    sn_ids2 = []
    for i in sn_ids:
        if i == maxid0:
            sn_ids2.append(i)
        else:
            r2 = R_full[i][maxid0] * R_full[i][maxid0]
            # source: ldr2tmp = r^2; keep (r2 < ld_prune) and (r2 > ld_min)
            if r2 < LD_PRUNE_R2 and r2 > LD_MIN_R2:
                sn_ids2.append(i)
    if len(sn_ids2) < M_HETERO:
        return -9.0, -9, "proxies (0.05<r2<0.9) < m_hetero(3)"
    maxid1 = max(sn_ids2, key=lambda i: abs(zsxz[i]))
    # 4) rm_cor_sbat with cutoff sqrt(ld_prune^2... note: source uses
    #    ld_top = sqrt(ldr2_top)=sqrt(0.9))
    m = len(sn_ids2)
    idxmap = sn_ids2
    subR = [[R_full[idxmap[i]][idxmap[j]] for j in range(m)] for i in range(m)]
    removed = rm_cor_sbat(subR, math.sqrt(LD_PRUNE_R2), m)
    if m - len(removed) < M_HETERO:
        return -9.0, -9, "after rm_cor_sbat < m_hetero(3)"
    keep = [idxmap[i] for i in range(m) if i not in set(removed)]
    # 5) cap opt_hetero by |z|
    keep = update_snidx(zsxz, keep, OPT_HETERO)
    # 6) bxy_hetero3 on final set with LD submatrix
    kk = len(keep)
    subLD = [[R_full[keep[i]][keep[j]] for j in range(kk)] for i in range(kk)]
    pdev, neig = bxy_hetero3([byz[i] for i in keep], [bxz[i] for i in keep],
                             [seyz[i] for i in keep], [sexz[i] for i in keep],
                             [zsxz[i] for i in keep], subLD)
    return pdev, kk, ""

# ----------------------------------------------------------------------------
# Data loading / harmonization
# ----------------------------------------------------------------------------
def load_efile(path):
    rows = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        assert header[:2] == ["ProbeID", "SNP"], header
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            rows[p[1]] = dict(a1=p[4], a2=p[5], freq=float(p[6]),
                              b=float(p[7]), se=float(p[8]), pval=float(p[9]))
    return rows

def load_ma(path):
    rows = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        assert header[0] == "SNP", header
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            try:
                rows[p[0]] = dict(a1=p[1], a2=p[2], freq=float(p[3]),
                                  b=float(p[4]), se=float(p[5]),
                                  pval=float(p[6]),
                                  n=(float(p[7]) if p[7] not in ("NA", "", "NaN")
                                     else None))
            except ValueError:
                continue
    return rows

def load_ld(path):
    with open(path) as f:
        cols = f.readline().rstrip("\n").split(",")
        names = cols[1:]
        idx = {s: i for i, s in enumerate(names)}
        mat = []
        for line in f:
            p = line.rstrip("\n").split(",")
            mat.append([float(x) for x in p[1:]])
    return names, idx, mat

def harmonize(efile, ma):
    """Align outcome (GWAS) beta/se/freq onto exposure A1; apply diff-freq.
    Returns (snps, rows) where rows[i] = dict(bx,sex,by,sey,freq)."""
    snps, rows = [], []
    ndiff_allele, nfreq_excl = 0, 0
    for rs, e in efile.items():
        if rs not in ma:
            continue
        g = ma[rs]
        if e["a1"] == g["a1"] and e["a2"] == g["a2"]:
            by, fy = g["b"], g["freq"]
        elif e["a1"] == g["a2"] and e["a2"] == g["a1"]:
            by, fy = -g["b"], 1.0 - g["freq"]
            ndiff_allele += 1
        else:
            continue  # allele-ambiguous / mismatch
        if abs(e["freq"] - fy) > DIFF_FREQ:
            nfreq_excl += 1
            continue
        if e["se"] <= 0 or g["se"] <= 0 or e["b"] == 0:
            continue
        snps.append(rs)
        rows.append(dict(bx=e["b"], sex=e["se"], by=by, sey=g["se"],
                         freq=e["freq"], p_eQTL=e["pval"], p_GWAS=g["pval"]))
    return snps, rows, ndiff_allele, nfreq_excl

# ----------------------------------------------------------------------------
# Run all 22 pairs
# ----------------------------------------------------------------------------
PAIRS = []
NPPA_OUT = ["AF", "CAD_cardio", "CAD_vdh", "CES", "HF", "stroke"]
for exp in ["SCALLOP", "INTERVAL", "NPPAeQTL"]:
    for out in NPPA_OUT:
        PAIRS.append((exp, out, "nppa"))
for out in ["CAD_cardio", "CAD_vdh", "MI_har", "MI_finn"]:
    PAIRS.append(("NPR3eQTL", out, "npr3"))

THRESH_HEIDI = chi_val(1, P_HEIDI)

def main():
    log = []
    def P(s=""):
        log.append(s)

    P("SMR & HEIDI — in-house Python port of official SMR v1.4.2 source code")
    P("Reference: Zhu et al. 2016 Nat Genet 48:481-487 (doi:10.1038/ng.3538)")
    P("Implementation: faithful port of JianYang-Lab/SMR v1.4.2 (GPL-3.0)")
    P("Reason: prebuilt Windows binary unavailable (1.3.1 access-violation crash;")
    P("        1.4.2 no Windows build; WSL blocked by security policy)")
    P("")
    P("Parameters (SMR 1.4.2 defaults):")
    P("  peqtl-smr=%g  peqtl-heidi=%g (thresh chi2=%g)  ld-prune=%g  ld-min=%g"
      % (P_SMR_DEFAULT, P_HEIDI, THRESH_HEIDI, LD_PRUNE_R2, LD_MIN_R2))
    P("  m-hetero=%d  opt-hetero=%d  MAX_NUM_LD=%d  diff-freq=%g  new-heidi-mth=TRUE"
      % (M_HETERO, OPT_HETERO, MAX_NUM_LD, DIFF_FREQ))
    P("  sample-overlap correction: OFF (default)")
    P("LD reference: 1000 Genomes EUR (phase3 v5b), +/-100kb cis region")
    P("")

    ld_cache = {}
    def get_ld(locus):
        if locus not in ld_cache:
            ld_cache[locus] = load_ld(os.path.join(WORK, "ld_%s.csv" % locus))
        return ld_cache[locus]

    results = []
    for exp, out, locus in PAIRS:
        efile = load_efile(os.path.join(WORK, exp + ".efile"))
        ma = load_ma(os.path.join(WORK, "ma_%s_%s.ma" % (out, locus)))
        names, lidx, R = get_ld(locus)

        snps, rows, nflip, nfreq = harmonize(efile, ma)
        # restrict to SNPs present in LD matrix
        keep = [(s, r) for s, r in zip(snps, rows) if s in lidx]
        n_ld_miss = len(snps) - len(keep)
        snps = [s for s, _ in keep]
        rows = [r for _, r in keep]

        P("=" * 100)
        P("Exposure: %-10s  Outcome: %-10s  Locus: %s" % (exp, out, locus.upper()))
        P("  eQTL/pQTL SNPs=%d  GWAS SNPs=%d  merged & allele-matched=%d "
          "(allele-flipped=%d, freq-excluded(diff>0.2)=%d, not-in-LD-matrix=%d)"
          % (len(efile), len(ma), len(snps), nflip, nfreq, n_ld_miss))
        if not snps:
            P("  -> no overlapping SNPs; skipped")
            results.append(dict(exposure=exp, outcome=out, locus=locus,
                                note="no overlapping SNPs"))
            continue

        bx = [r["bx"] for r in rows]
        sex = [r["sex"] for r in rows]
        by = [r["by"] for r in rows]
        sey = [r["sey"] for r in rows]
        zxz = [b / s for b, s in zip(bx, sex)]

        # R submatrix over merged SNPs
        ii = [lidx[s] for s in snps]
        Rsub = [[R[i][j] for j in ii] for i in ii]

        maxid = max(range(len(zxz)), key=lambda i: abs(zxz[i]))
        topsnp = snps[maxid]
        p_eQTL_top = pchisq(zxz[maxid] ** 2, 1)
        P("  Top cis instrument: %s  b_exp=%.5g  se_exp=%.5g  z=%.2f  p_eQTL=%.3g"
          % (topsnp, bx[maxid], sex[maxid], zxz[maxid], p_eQTL_top))

        b_smr, se_smr, chisq, psmr = smr_test(bx[maxid], sex[maxid],
                                              by[maxid], sey[maxid])
        P("  SMR test: b_SMR=%.5g  se_SMR=%.5g  chi2=%.3f  p_SMR=%.4g"
          % (b_smr, se_smr, chisq, psmr))

        note = ""
        if p_eQTL_top > P_SMR_DEFAULT:
            if p_eQTL_top <= P_SMR_RELAXED:
                note = ("top cis p_eQTL=%.3g > 5e-8 default; reported at "
                        "relaxed peqtl-smr=5e-6 sensitivity" % p_eQTL_top)
            else:
                note = ("top cis p_eQTL=%.3g exceeds even relaxed 5e-6 "
                        "threshold; instrument too weak" % p_eQTL_top)

        pdev, nsnp, hnote = heidi_test_new(by, bx, sey, sex, zxz, Rsub,
                                           THRESH_HEIDI)
        if pdev >= 0:
            P("  HEIDI (new method): p_HEIDI=%.4g  nsnp=%d%s%s"
              % (pdev, nsnp, "  [" + hnote + "]" if hnote else "",
                 "  NOTE: " + note if note else ""))
        else:
            P("  HEIDI (new method): NA (nsnp<3)%s"
              % ("  [" + hnote + "]" if hnote else ""))

        results.append(dict(
            exposure=exp, outcome=out, locus=locus, topSNP=topsnp,
            b_eQTL=bx[maxid], se_eQTL=sex[maxid], p_eQTL=p_eQTL_top,
            b_GWAS=by[maxid], se_GWAS=sey[maxid],
            b_SMR=b_smr, se_SMR=se_smr, p_SMR=psmr,
            p_HEIDI=(pdev if pdev >= 0 else "NA"),
            nsnp_HEIDI=(nsnp if pdev >= 0 else "NA"),
            note=note + ("; " + hnote if hnote else "")))

    P("")
    P("=" * 100)
    P("SUMMARY TABLE")
    P("%-10s %-11s %-6s %-12s %10s %10s %10s %8s %6s"
      % ("Exposure", "Outcome", "Locus", "topSNP", "p_eQTL", "p_SMR",
         "b_SMR", "p_HEIDI", "nsnp"))
    for r in results:
        P("%-10s %-11s %-6s %-12s %10.3g %10.3g %10.3g %8s %6s"
          % (r["exposure"], r["outcome"], r["locus"], r.get("topSNP", "-"),
             r.get("p_eQTL", float("nan")), r.get("p_SMR", float("nan")),
             r.get("b_SMR", float("nan")),
             ("%.3g" % r["p_HEIDI"]) if isinstance(r.get("p_HEIDI"), float) else "NA",
             r.get("nsnp_HEIDI", "NA")))

    with open(os.path.join(WORK, "step13q_smr_heidi_results.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "exposure", "outcome", "locus", "topSNP", "b_eQTL", "se_eQTL",
            "p_eQTL", "b_GWAS", "se_GWAS", "b_SMR", "se_SMR", "p_SMR",
            "p_HEIDI", "nsnp_HEIDI", "note"])
        w.writeheader()
        w.writerows(results)

    with open(os.path.join(WORK, "step13q_output.txt"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(log))
    print("done: %d pairs" % len(results))

if __name__ == "__main__":
    main()
