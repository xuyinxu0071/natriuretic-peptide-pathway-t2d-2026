# -*- coding: utf-8 -*-
"""
Step 11 (P1): COJO 式条件分析 — rs5068/rs5066/rs198411/rs1421811 条件化
================================================================================
方法:
  z 尺度精确条件检验 (多元正态 + 1000G EUR LD 作为采样相关矩阵):
    z_t|S = (z_t - R_tS R_SS^{-1} z_S) / sqrt(1 - R_tS R_SS^{-1} R_St)
  beta 条件近似 (COJO 同款近似, 同数据集内同 N):
    b_t|S = b_t - sum_j (se_t/se_j) * gamma_j * b_j,  se_t|S = se_t * sqrt(1 - gamma R_St)
  其中 gamma = R_tS R_SS^{-1}。
  全部 beta 先按 step4 规则对齐到 dbSNP ALT 等位方向。
输出: step11_cojo_results.csv + step11_output.txt
"""
import csv, json, math, os
import numpy as np

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

def harmonize_np(variants, assoc):
    """step4 协调规则的无回文丢弃版: 仅当 ea/nea 与 dbSNP ref/alt 精确匹配(或互补)
    时链方向即已解决, 回文不再构成歧义 -> 保留。eaf QC 与 step4 相同。"""
    out = {}
    for v in variants:
        rsid = v["rsid"]
        a = assoc.get(rsid)
        if a is None or a.get("beta") is None or a.get("se") is None:
            continue
        ea, nea = (a.get("ea") or "").upper(), (a.get("nea") or "").upper()
        if len(ea) != 1 or len(nea) != 1:
            continue
        try:
            eaf_raw = float(a.get("eaf"))
        except (TypeError, ValueError):
            eaf_raw = None
        ref, alt = v["ref"], v["alt"]
        beta = a["beta"]
        if ea == ref and nea == alt:
            beta, eaf = -beta, (1 - eaf_raw) if eaf_raw is not None else None
        elif ea == alt and nea == ref:
            eaf = eaf_raw
        elif ea == COMPLEMENT.get(ref) and nea == COMPLEMENT.get(alt):
            beta, eaf = -beta, (1 - eaf_raw) if eaf_raw is not None else None
        elif ea == COMPLEMENT.get(alt) and nea == COMPLEMENT.get(ref):
            eaf = eaf_raw
        else:
            continue
        if eaf is None:
            eaf = v["af"]
        else:
            if abs(eaf - v["af"]) > 0.25:
                continue
        if eaf is None or not (0.0 < eaf < 1.0):
            continue
        out[rsid] = (v["pos"], beta, a["se"], eaf, a.get("n"), alt, ref)
    return out

LOCI = {
    "nppa": {
        "variants": "nppa_locus_variants.json",
        "ld": "ld_nppa.csv",
        "datasets": {
            "ebi-a-GCST90012082": "NT-proBNP pQTL (SCALLOP)",
            "prot-a-2078": "NT-proBNP pQTL (INTERVAL SomaLogic)",
            "eqtl-a-ENSG00000175206": "NPPA eQTL (eQTLGen whole blood)",
            "ieu-a-7": "CAD (CARDIoGRAMplusC4D)",
            "ebi-a-GCST005195": "CAD (van der Harst)",
            "ebi-a-GCST006061": "AF (Nielsen)",
            "ebi-a-GCST009541": "HF (HERMES)",
            "ebi-a-GCST005838": "Any stroke (MEGASTROKE)",
            "ebi-a-GCST006910": "Cardioembolic stroke (MEGASTROKE)",
        },
        "keys": ["rs5068", "rs5066", "rs198411"],
    },
    "npr3": {
        "variants": "npr3_locus_variants.json",
        "ld": "ld_npr3.csv",
        "datasets": {
            "eqtl-a-ENSG00000113389": "NPR3 eQTL (eQTLGen whole blood)",
            "ieu-b-38": "SBP (Evangelou)",
            "ieu-a-7": "CAD (CARDIoGRAMplusC4D)",
            "ebi-a-GCST005195": "CAD (van der Harst)",
            "ebi-a-GCST011364": "MI (Hartiala)",
            "finn-b-I9_MI": "MI (FinnGen R9)",
        },
        "keys": ["rs1421811"],
    },
}

def p_two_sided(z):
    return math.erfc(abs(z) / math.sqrt(2))

def load_ld(path):
    with open(path) as f:
        rd = list(csv.reader(f))
    ids = rd[0][1:]
    M = np.array([[float(x) for x in row[1:]] for row in rd[1:]], dtype=float)
    idx = {r: i for i, r in enumerate(ids)}
    # 对称化 + 对角=1
    M = (M + M.T) / 2
    np.fill_diagonal(M, 1.0)
    return ids, idx, M

def cond_test(M, idx, ds, target, condset):
    """ds: {rsid: (pos,beta,se,eaf,n,alt,ref)} 已对齐 ALT; 返回条件统计或 None"""
    if target not in ds or target not in idx:
        return None
    if any(c not in ds or c not in idx for c in condset):
        return None
    bt, st = ds[target][1], ds[target][2]
    zt = bt / st
    if not condset:
        return dict(r2red=0.0, z=zt, b=bt, se=st)
    S = [c for c in condset]
    R_SS = M[np.ix_([idx[c] for c in S], [idx[c] for c in S])]
    R_tS = M[[idx[target]], [idx[c] for c in S]]
    try:
        gam = R_tS @ np.linalg.inv(R_SS)   # 1 x k
    except np.linalg.LinAlgError:
        return None
    r2red = float(gam @ R_tS.T)            # R_tS R_SS^-1 R_St, 标量
    if r2red >= 0.999:
        return None
    zS = np.array([ds[c][1] / ds[c][2] for c in S])
    bS = np.array([ds[c][1] for c in S])
    seS = np.array([ds[c][2] for c in S])
    zc = (zt - float(gam @ zS)) / math.sqrt(1 - r2red)
    bc = bt - float(np.sum((st / seS) * gam * bS))
    sc = st * math.sqrt(1 - r2red)
    return dict(r2red=r2red, z=zc, b=bc, se=sc)

def main():
    out_rows = []
    log = []
    def P(s):
        log.append(s)

    for locus, cfg in LOCI.items():
        variants = json.load(open(os.path.join(WORK, cfg["variants"])))
        ids, idx, M = load_ld(os.path.join(WORK, cfg["ld"]))
        P(f"\n===== {locus.upper()} 位点: LD 面板 {len(ids)} 变异 =====")

        # 关键 SNP 两两 LD
        P("\n-- 关键变异两两 LD (r / r²):")
        kk = cfg["keys"]
        for i in range(len(kk)):
            for j in range(i + 1, len(kk)):
                if kk[i] in idx and kk[j] in idx:
                    r = M[idx[kk[i]], idx[kk[j]]]
                    P(f"   {kk[i]} x {kk[j]}: r={r:+.4f}  r²={r*r:.4f}")

        for dsid, label in cfg["datasets"].items():
            f = os.path.join(WORK, f"{locus}__{dsid}_assoc.json")
            if not os.path.exists(f):
                P(f"\n-- {label} [{dsid}]: assoc 文件缺失, 跳过")
                continue
            ds = harmonize_np(variants, json.load(open(f)))
            # 仅保留 LD 面板内
            ds = {k: v for k, v in ds.items() if k in idx}
            P(f"\n-- {label} [{dsid}] (LD 面板内 {len(ds)} 变异):")
            # 边际统计
            P("   边际: " + "; ".join(
                f"{k}: beta={ds[k][1]:+.4g}, z={ds[k][1]/ds[k][2]:+.2f}, p={p_two_sided(ds[k][1]/ds[k][2]):.2g}"
                for k in cfg["keys"] if k in ds))
            # 区域 lead (按 |z|)
            if ds:
                lead = max(ds, key=lambda r: abs(ds[r][1] / ds[r][2]))
                zl = ds[lead][1] / ds[lead][2]
                P(f"   区域 lead: {lead} (z={zl:+.2f}, p={p_two_sided(zl):.2g})")
                if lead in idx:
                    for k in cfg["keys"]:
                        if k in idx:
                            r = M[idx[lead], idx[k]]
                            P(f"      lead {lead} x {k}: r={r:+.4f}")
            else:
                lead = None

            # ---- 条件检验组合 ----
            tests = []
            if locus == "nppa":
                tests = [
                    ("rs5068|rs5066", "rs5068", ["rs5066"]),
                    ("rs5066|rs5068", "rs5066", ["rs5068"]),
                    ("rs5068|rs198411", "rs5068", ["rs198411"]),
                    ("rs198411|rs5068", "rs198411", ["rs5068"]),
                    ("rs5068|rs5066+rs198411", "rs5068", ["rs5066", "rs198411"]),
                ]
            else:
                if lead and lead != "rs1421811":
                    tests = [
                        ("rs1421811|lead(%s)" % lead, "rs1421811", [lead]),
                        ("lead(%s)|rs1421811" % lead, lead, ["rs1421811"]),
                    ]
            if lead and lead not in cfg["keys"] and locus == "nppa":
                tests += [(f"rs5068|lead({lead})", "rs5068", [lead]),
                          (f"lead({lead})|rs5068", lead, ["rs5068"])]

            for name, t, S in tests:
                res = cond_test(M, idx, ds, t, S)
                if res is None:
                    P(f"   {name}: 不可计算 (LD 奇异或缺失)")
                    continue
                z, b, se, r2r = res["z"], res["b"], res["se"], res["r2red"]
                pv = p_two_sided(z)
                out_rows.append(dict(
                    locus=locus, dataset=dsid, label=label, test=name, target=t,
                    conditioned_on="+".join(S) if S else "none",
                    r2_reduced_by_ld=f"{r2r:.4f}",
                    beta_marginal=f"{ds[t][1]:+.6g}", z_marginal=f"{ds[t][1]/ds[t][2]:+.4f}",
                    p_marginal=f"{p_two_sided(ds[t][1]/ds[t][2]):.3g}",
                    beta_cond=f"{b:+.6g}", se_cond=f"{se:.6g}", z_cond=f"{z:+.4f}",
                    p_cond=f"{pv:.3g}",
                    OR_cond=f"{math.exp(b):.4g}", OR_lo=f"{math.exp(b-1.96*se):.4g}",
                    OR_hi=f"{math.exp(b+1.96*se):.4g}"))
                P(f"   {name}: z {ds[t][1]/ds[t][2]:+.2f} -> {z:+.2f} (p {p_two_sided(ds[t][1]/ds[t][2]):.2g} -> {pv:.2g}); "
                  f"OR_cond={math.exp(b):.3f} [{math.exp(b-1.96*se):.3f},{math.exp(b+1.96*se):.3f}]")

            # ---- 单步二级信号扫描: 条件化于 lead ----
            if ds and lead and lead in idx:
                zs = {r: ds[r][1] / ds[r][2] for r in ds}
                zlead = zs[lead]
                il = idx[lead]
                second = []
                for r in ds:
                    if r == lead or r not in idx:
                        continue
                    rr = M[idx[r], il]
                    if abs(rr) >= 0.95:
                        continue  # 高 LD 处 LD 参考失配会被 1/sqrt(1-r²) 放大 -> 排除
                    zc = (zs[r] - rr * zlead) / math.sqrt(1 - rr * rr)
                    second.append((abs(zc), zc, r))
                second.sort(reverse=True)
                n_test = len(second)
                bonf = 0.05 / max(n_test, 1)
                top = second[:5]
                P(f"   二级信号扫描 (条件于 {lead}, {n_test} 变异, Bonferroni p<{bonf:.2g}): "
                  + "; ".join(f"{r}(z={zc:+.2f},p={p_two_sided(zc):.2g})" for _, zc, r in top))
                # 关键 SNP 条件于 lead 的残差信号
                for k in cfg["keys"]:
                    if k in ds and k in idx and k != lead:
                        rr = M[idx[k], idx[lead]]
                        if abs(rr) < 0.95:
                            zc = (zs[k] - rr * zlead) / math.sqrt(1 - rr * rr)
                            P(f"      {k} | {lead}: z={zc:+.2f} (p={p_two_sided(zc):.2g})")
                        else:
                            P(f"      {k} | {lead}: r={rr:+.3f} (同一单倍型, 条件化不适用)")
                # 记录 top 二级信号
                for rank, (az, zc, r) in enumerate(top[:3], 1):
                    b_c = ds[r][1]
                    out_rows.append(dict(
                        locus=locus, dataset=dsid, label=label,
                        test=f"secondary#{rank}|{lead}", target=r, conditioned_on=lead,
                        r2_reduced_by_ld="", beta_marginal=f"{ds[r][1]:+.6g}",
                        z_marginal=f"{zs[r]:+.4f}", p_marginal=f"{p_two_sided(zs[r]):.3g}",
                        beta_cond="", se_cond="", z_cond=f"{zc:+.4f}",
                        p_cond=f"{p_two_sided(zc):.3g}", OR_cond="", OR_lo="", OR_hi=""))

    with open(os.path.join(WORK, "step11_cojo_results.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    with open(os.path.join(WORK, "step11_output.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(log))
    print(f"DONE: {len(out_rows)} rows")

if __name__ == "__main__":
    main()
