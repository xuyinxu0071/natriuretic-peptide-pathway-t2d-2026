# -*- coding: utf-8 -*-
"""
step13r_refsmr.py
refSNP SMR & HEIDI at rs5068 (SMR --snp rs5068 functionality).

Faithful port of heidi_test_ref_new (src_SMR_data.cpp lines 4590-4685):
  - candidate proxies: {refid} + {i: z_eQTL^2 > 10}; < m_hetero(3) -> NA
  - cap MAX_NUM_LD(500) by |z|, then re-append refid if dropped
  - keep proxies with (r^2 - 0.9) < 1e-6 and r^2 > 0.05 w.r.t. refid
  - rm_cor_sbat pairwise |r| > sqrt(0.9); remaining < 3 -> NA
  - cap opt_hetero(20) by |z|; if refid not in top-20, replace last slot
  - bxy_hetero3 (internal top = argmax|z| within final set)
Per source line 4084, when refSNP is specified the p_smr filter on the
instrument is NOT applied (only the default top-SNP path is filtered).

Exposures: SCALLOP NT-proBNP, INTERVAL NT-proBNP, eQTLGen NPPA eQTL.
Outcomes: 6 NPPA-locus GWAS (AF, CAD_cardio, CAD_vdh, CES, HF, stroke).
"""
import os
import csv
import importlib.util

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"

spec = importlib.util.spec_from_file_location(
    "step13q", os.path.join(WORK, "step13q_smr_heidi.py"))
q = importlib.util.module_from_spec(spec)
spec.loader.exec_module(q)

REF = "rs5068"
PAIRS = [(e, o) for e in ["SCALLOP", "INTERVAL", "NPPAeQTL"]
         for o in ["AF", "CAD_cardio", "CAD_vdh", "CES", "HF", "stroke"]]


def heidi_test_ref_new(byz, bxz, seyz, sexz, zsxz, R_full, refid):
    """Port of heidi_test_ref_new. refid = index of the specified SNP."""
    n_all = len(zsxz)
    sn_ids = [i for i in range(n_all)
              if i == refid or zsxz[i] * zsxz[i] - q.THRESH_HEIDI > 1e-6]
    if len(sn_ids) < q.M_HETERO:
        return -9.0, -9, "candidates (ref+z2>10) < m_hetero(3)"
    sn_ids = q.update_snidx(zsxz, sn_ids, q.MAX_NUM_LD)
    if refid not in sn_ids:
        sn_ids.append(refid)
    # r^2 filter w.r.t. refid: (r2 - 0.9) < 1e-6 and r2 > 0.05
    sn_ids2 = []
    for i in sn_ids:
        if i == refid:
            sn_ids2.append(i)
        else:
            r2 = R_full[i][refid] * R_full[i][refid]
            if (r2 - q.LD_PRUNE_R2) < 1e-6 and r2 > q.LD_MIN_R2:
                sn_ids2.append(i)
    if len(sn_ids2) < q.M_HETERO:
        return -9.0, -9, "proxies (0.05<r2<=0.9) < m_hetero(3)"
    # rm_cor_sbat at |r| > sqrt(0.9)
    m = len(sn_ids2)
    idxmap = sn_ids2
    subR = [[R_full[idxmap[i]][idxmap[j]] for j in range(m)] for i in range(m)]
    removed = q.rm_cor_sbat(subR, math_sqrt(q.LD_PRUNE_R2), m)
    if m - len(removed) < q.M_HETERO:
        return -9.0, -9, "after rm_cor_sbat < m_hetero(3)"
    keep = [idxmap[i] for i in range(m) if i not in set(removed)]
    # cap opt_hetero by |z|; refid guaranteed in final set
    keep = q.update_snidx(zsxz, keep, q.OPT_HETERO)
    if refid not in keep and len(keep) > 0:
        keep[-1] = refid
    kk = len(keep)
    subLD = [[R_full[keep[i]][keep[j]] for j in range(kk)] for i in range(kk)]
    pdev, neig = q.bxy_hetero3([byz[i] for i in keep], [bxz[i] for i in keep],
                               [seyz[i] for i in keep], [sexz[i] for i in keep],
                               [zsxz[i] for i in keep], subLD)
    return pdev, kk, ""


def math_sqrt(x):
    import math
    return math.sqrt(x)


def main():
    log = []

    def P(s=""):
        log.append(s)

    P("refSNP SMR & HEIDI at %s — port of official SMR v1.4.2 "
      "(--snp rs5068 path; heidi_test_ref_new)" % REF)
    P("Note: with --snp specified, SMR does NOT apply the p_smr filter to the")
    P("      instrument (source line 4084: filter only when refSNP==nullptr).")
    P("")

    names, lidx, R = q.load_ld(os.path.join(WORK, "ld_nppa.csv"))
    results = []

    for exp, out in PAIRS:
        efile = q.load_efile(os.path.join(WORK, exp + ".efile"))
        ma = q.load_ma(os.path.join(WORK, "ma_%s_nppa.ma" % out))
        snps, rows, nflip, nfreq = q.harmonize(efile, ma)
        keep = [(s, r) for s, r in zip(snps, rows) if s in lidx]
        snps = [s for s, _ in keep]
        rows = [r for _, r in keep]

        P("=" * 100)
        P("Exposure: %-10s  Outcome: %-11s  refSNP: %s" % (exp, out, REF))
        if REF not in snps:
            P("  -> %s not present after harmonization; skipped" % REF)
            results.append(dict(exposure=exp, outcome=out, refSNP=REF,
                                note="refSNP not in merged data"))
            continue
        ri = snps.index(REF)
        bx = [r["bx"] for r in rows]
        sex = [r["sex"] for r in rows]
        by = [r["by"] for r in rows]
        sey = [r["sey"] for r in rows]
        zxz = [b / s for b, s in zip(bx, sex)]

        ii = [lidx[s] for s in snps]
        Rsub = [[R[i][j] for j in ii] for i in ii]

        P("  instrument: b_exp=%.5g  se_exp=%.5g  z=%.2f  p_eQTL=%.3g"
          % (bx[ri], sex[ri], zxz[ri], q.pchisq(zxz[ri] ** 2, 1)))
        P("  outcome:   b_GWAS=%.5g  se_GWAS=%.5g"
          % (by[ri], sey[ri]))
        b_smr, se_smr, chisq, psmr = q.smr_test(bx[ri], sex[ri], by[ri], sey[ri])
        P("  SMR test at %s: b_SMR=%.5g  se_SMR=%.5g  chi2=%.3f  p_SMR=%.4g"
          % (REF, b_smr, se_smr, chisq, psmr))

        pdev, nsnp, hnote = heidi_test_ref_new(by, bx, sey, sex, zxz, Rsub, ri)
        if pdev >= 0:
            P("  HEIDI (ref new): p_HEIDI=%.4g  nsnp=%d%s"
              % (pdev, nsnp, "  [" + hnote + "]" if hnote else ""))
        else:
            P("  HEIDI (ref new): NA%s" % ("  [" + hnote + "]" if hnote else ""))

        results.append(dict(
            exposure=exp, outcome=out, refSNP=REF,
            b_eQTL=bx[ri], se_eQTL=sex[ri],
            p_eQTL=q.pchisq(zxz[ri] ** 2, 1),
            b_GWAS=by[ri], se_GWAS=sey[ri],
            b_SMR=b_smr, se_SMR=se_smr, p_SMR=psmr,
            p_HEIDI=(pdev if pdev >= 0 else "NA"),
            nsnp_HEIDI=(nsnp if pdev >= 0 else "NA"),
            note=hnote))

    P("")
    P("=" * 100)
    P("SUMMARY (refSNP = %s)" % REF)
    P("%-10s %-11s %10s %10s %10s %8s %6s"
      % ("Exposure", "Outcome", "p_eQTL", "p_SMR", "b_SMR", "p_HEIDI", "nsnp"))
    for r in results:
        P("%-10s %-11s %10.3g %10.3g %10.3g %8s %6s"
          % (r["exposure"], r["outcome"], r.get("p_eQTL", float("nan")),
             r.get("p_SMR", float("nan")), r.get("b_SMR", float("nan")),
             ("%.3g" % r["p_HEIDI"]) if isinstance(r.get("p_HEIDI"), float) else "NA",
             r.get("nsnp_HEIDI", "NA")))

    with open(os.path.join(WORK, "step13r_refsmr_results.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "exposure", "outcome", "refSNP", "b_eQTL", "se_eQTL", "p_eQTL",
            "b_GWAS", "se_GWAS", "b_SMR", "se_SMR", "p_SMR", "p_HEIDI",
            "nsnp_HEIDI", "note"])
        w.writeheader()
        w.writerows(results)
    with open(os.path.join(WORK, "step13r_output.txt"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(log))
    print("done: %d pairs" % len(results))


if __name__ == "__main__":
    main()
