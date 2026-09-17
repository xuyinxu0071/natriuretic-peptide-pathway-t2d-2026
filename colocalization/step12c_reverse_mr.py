# -*- coding: utf-8 -*-
"""
Step 12c (P2): Reverse MR — AF/CAD/MI/HF → NT-proBNP pQTL
================================================================================
暴露(IV): step12_tophits_*.json (clumped p<5e-8, r2<0.001, 10Mb)
结局: NT-proBNP pQTL — SCALLOP (ebi-a-GCST90012082, n=21,758) 主
      INTERVAL SomaLogic (prot-a-2078, n=3,301) 敏感性
方法: per-SNP Wald; IVW 固定/随机效应; 加权中位数 (bootstrap CI, seed=42);
      MR-Egger (>=10 IV); 异质性 Q; mean F; 排除回文 IV 与 MHC(chr6:25-35Mb);
      排除 NPPA/NPPB/NPR3 位点 ±500kb IV (排除限制专项保护)。
"""
import json, math, os, time, urllib.request, urllib.error
import numpy as np

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
DIR = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"
API = "https://api.opengwas.io/api"
TOK = open(os.path.join(DIR, ".opengwas_token.txt")).read().strip()

EXPOSURES = {
    "ebi-a-GCST006061": "AF (Nielsen 2018)",
    "ieu-a-7": "CAD (CARDIoGRAMplusC4D)",
    "ebi-a-GCST005195": "CAD (van der Harst 2018)",
    "ebi-a-GCST011364": "MI (Hartiala 2021)",
    "ebi-a-GCST009541": "HF (HERMES)",
}
OUTCOMES = {
    "ebi-a-GCST90012082": "NT-proBNP (SCALLOP Olink, n=21,758)",
    "prot-a-2078": "NT-proBNP (INTERVAL SomaLogic, n=3,301)",
}
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
MHC = ("6", 25000000, 35000000)
LOCI_EXCL = [("1", 11805974), ("5", 32714270)]  # NPPA/NPPB, NPR3 lead 位置

def api_post(path, payload, retries=5):
    for i in range(retries):
        try:
            req = urllib.request.Request(API + path,
                data=json.dumps(payload).encode("utf-8"),
                headers={"Authorization": "Bearer " + TOK, "Content-Type": "application/json"},
                method="POST")
            return json.load(urllib.request.urlopen(req, timeout=300))
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504, 401, 403) and i < retries - 1:
                wait = 20 * (i + 1) if e.code == 429 else 10
                time.sleep(wait); continue
            return None
        except Exception:
            if i < retries - 1:
                time.sleep(15); continue
            return None

def pull_assoc(dsid, rsids):
    out, done, batch = {}, 0, 64
    while done < len(rsids):
        chunk = rsids[done:done + batch]
        res = api_post("/associations", {"variant": chunk, "id": [dsid]})
        if res is None:
            return None
        for a in res:
            if a.get("id") == dsid and a.get("beta") is not None and a.get("se") is not None:
                out[a["rsid"]] = a
        done += batch
        time.sleep(0.6)
    return out

def harmonize_iv(iv, orec):
    """把结局 beta 对齐到 IV 的 ea; 返回 (beta, se) 或 None"""
    ea_i, nea_i = iv["ea"].upper(), iv["nea"].upper()
    ea_o, nea_o = (orec.get("ea") or "").upper(), (orec.get("nea") or "").upper()
    if len(ea_o) != 1 or len(nea_o) != 1:
        return None
    if {ea_i, nea_i} in ({"A", "T"}, {"C", "G"}):
        return None  # 回文 IV: 无法确认链方向 -> 丢弃
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
    """卡方分布上尾概率 (正则化不完全 gamma Q(df/2, x/2), Numerical Recipes)"""
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

def weighted_median(b, se, nboot=1000, seed=42):
    w = 1.0 / np.asarray(se) ** 2
    b = np.asarray(b)
    def wm(bs, ws):
        o = np.argsort(bs)
        bs, ws = bs[o], ws[o]
        cw = np.cumsum(ws) / np.sum(ws)
        return bs[np.searchsorted(cw, 0.5)]
    est = wm(b, w)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(nboot):
        i = rng.integers(0, len(b), len(b))
        try:
            boots.append(wm(b[i], w[i]))
        except Exception:
            pass
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return est, lo, hi

def mr_egger(bx, by, sey):
    """Egger: by = a + slope*bx, 权重 1/sey^2"""
    bx, by, sey = np.asarray(bx), np.asarray(by), np.asarray(sey)
    w = 1.0 / sey ** 2
    X = np.column_stack([np.ones_like(bx), bx])
    W = np.diag(w)
    cov = np.linalg.inv(X.T @ W @ X)
    coef = cov @ X.T @ W @ by
    resid = by - X @ coef
    s2 = float(np.sum(w * resid ** 2) / max(len(bx) - 2, 1))
    cov *= s2
    se_a, se_sl = math.sqrt(cov[0, 0]), math.sqrt(cov[1, 1])
    z_a = coef[0] / se_a
    p_a = math.erfc(abs(z_a) / math.sqrt(2))
    return coef[0], se_a, p_a, coef[1], se_sl

def p2(z):
    return math.erfc(abs(z) / math.sqrt(2))

OUT = []
def p(s): OUT.append(str(s))

results = []
for odsid, olabel in OUTCOMES.items():
    for dsid, label in EXPOSURES.items():
        ivs = json.load(open(os.path.join(WORK, f"step12_tophits_{dsid}.json")))
        # 排除 MHC / 目标位点 ±500kb
        keep = []
        for iv in ivs:
            c = str(iv.get("chr"))
            pos = iv.get("position") or iv.get("pos")
            if c == MHC[0] and MHC[1] <= pos <= MHC[2]:
                continue
            if any(c == lc[0] and abs(pos - lc[1]) <= 500000 for lc in LOCI_EXCL):
                continue
            keep.append(iv)
        rsids = [iv["rsid"] for iv in keep]
        # 拉结局关联
        cache_f = os.path.join(WORK, f"step12_ivassoc_{dsid}__{odsid}.json")
        if os.path.exists(cache_f):
            orecs = json.load(open(cache_f))
        else:
            orecs = pull_assoc(odsid, rsids)
            if orecs is None:
                p(f"[{label} -> {olabel}] 拉取失败, 跳过"); continue
            json.dump(orecs, open(cache_f, "w"))
        # 协调
        rows, dropped = [], {"palindromic_or_allele": 0, "no_outcome_assoc": 0}
        for iv in keep:
            orec = orecs.get(iv["rsid"])
            if orec is None:
                dropped["no_outcome_assoc"] += 1; continue
            h = harmonize_iv(iv, orec)
            if h is None:
                dropped["palindromic_or_allele"] += 1; continue
            b_out, se_out = h
            b_exp, se_exp = iv["beta"], iv["se"]
            if b_exp == 0:
                dropped["palindromic_or_allele"] += 1; continue
            rows.append(dict(rsid=iv["rsid"], b_exp=b_exp, se_exp=se_exp,
                             b_out=b_out, se_out=se_out))
        k = len(rows)
        if k < 3:
            p(f"\n[{label} -> {olabel}] k={k} < 3, 不可做 IVW (丢弃: {dropped})"); continue
        bx = np.array([r["b_exp"] for r in rows])
        by = np.array([r["b_out"] for r in rows])
        sex = np.array([r["se_exp"] for r in rows])
        sey = np.array([r["se_out"] for r in rows])
        bxy = by / bx
        se_xy = np.abs(bxy) * np.sqrt((sey / by) ** 2 + (sex / bx) ** 2)
        F = (bx / sex) ** 2
        bf, se_f, br, se_r, Q, df, tau2 = ivw(bxy, se_xy)
        wm, wm_lo, wm_hi = weighted_median(bxy, se_xy)
        pQ = chi2_sf(Q, df) if df > 0 else float("nan")
        eg = mr_egger(bx, by, sey) if k >= 10 else None
        z_f, z_r = bf / se_f, br / se_r
        p(f"\n[{label} -> {olabel}]")
        p(f"  IV: 拉取 {len(ivs)} -> 排除 MHC/位点后 {len(keep)} -> 协调后 k={k} (丢弃 {dropped})")
        p(f"  mean F = {np.mean(F):.1f}, min F = {np.min(F):.1f}")
        p(f"  IVW 固定:  beta={bf:+.4f} (se {se_f:.4f}), z={z_f:+.2f}, p={p2(z_f):.2g}")
        p(f"  IVW 随机:  beta={br:+.4f} (se {se_r:.4f}), z={z_r:+.2f}, p={p2(z_r):.2g}, tau2={tau2:.5f}")
        p(f"  加权中位数: beta={wm:+.4f} [boot95CI {wm_lo:+.4f}, {wm_hi:+.4f}]")
        p(f"  异质性 Q={Q:.1f} (df={df}), p={pQ:.3g}")
        if eg:
            p(f"  MR-Egger: 截距={eg[0]:+.4f} (p={eg[2]:.3g}); 斜率={eg[3]:+.4f} (se {eg[4]:.4f})")
        # 方向性: 每个 SNP Wald 的符号一致性
        pos_frac = np.mean(np.sign(bxy) == np.sign(bxy[0]))
        p(f"  Wald 符号一致率: {pos_frac:.0%}")
        results.append(dict(exposure=dsid, exp_label=label, outcome=odsid,
            out_label=olabel, n_iv_pulled=len(ivs), n_iv=k, dropped=json.dumps(dropped),
            mean_F=f"{np.mean(F):.1f}", min_F=f"{np.min(F):.1f}",
            ivw_fixed_beta=f"{bf:+.5f}", ivw_fixed_se=f"{se_f:.5f}",
            ivw_fixed_z=f"{z_f:+.3f}", ivw_fixed_p=f"{p2(z_f):.3g}",
            ivw_re_beta=f"{br:+.5f}", ivw_re_se=f"{se_r:.5f}",
            ivw_re_z=f"{z_r:+.3f}", ivw_re_p=f"{p2(z_r):.3g}",
            wm_beta=f"{wm:+.5f}", wm_lo=f"{wm_lo:+.5f}", wm_hi=f"{wm_hi:+.5f}",
            Q=f"{Q:.2f}", Q_df=df, Q_p=f"{pQ:.3g}", tau2=f"{tau2:.5f}",
            egger_intercept=f"{eg[0]:+.5f}" if eg else "", egger_p=f"{eg[2]:.3g}" if eg else "",
            egger_slope=f"{eg[3]:+.5f}" if eg else "",
            wald_sign_consistency=f"{pos_frac:.2f}"))

import csv
with open(os.path.join(WORK, "step12_reverse_mr_results.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
    w.writeheader(); w.writerows(results)
with open(os.path.join(WORK, "step12c_output.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT))
print("DONE")
