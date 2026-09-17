# -*- coding: utf-8 -*-
"""
step14f_gtex_harmonize.py
QC + harmonization + lead-eQTL Wald cis-MR for GTEx heart eQTL.

Inputs:  gtex_{atrial,ventricle}_{NPPA,NPR3}.tsv (step14c)
         ma_{outcome}_{locus}.ma (outcomes; NPPA locus 6, NPR3 locus 4)
Outputs: gtex_coloc_inputs/{pair}_{exposure,outcome}.csv   (step4 format)
         step14f_cismr_results.csv  (lead-eQTL Wald per pair)
         step14f_output.txt
"""
import csv
import math
import os
import importlib.util

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
spec = importlib.util.spec_from_file_location(
    "q", os.path.join(W, "step13q_smr_heidi.py"))
q = importlib.util.module_from_spec(spec)
spec.loader.exec_module(q)

TISSUES = ["atrial", "ventricle"]
GENES = {"NPPA": ("nppa", ["AF", "CAD_cardio", "CAD_vdh", "CES", "HF",
                            "stroke"]),
         "NPR3": ("npr3", ["CAD_cardio", "CAD_vdh", "MI_har", "MI_finn"])}
SAMPLE_N = {"atrial": 372, "ventricle": 382}
KEY_SNPS = {"NPPA": ["rs5068", "rs5066", "rs198411", "rs198367", "rs198389",
                     "rs12567119"],
            "NPR3": ["rs1421811", "rs112133948", "rs7271446"]}


def load_gtex(path):
    rows = []
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                continue
            try:
                rows.append(dict(
                    rsid=p[col["rsid"]], ref=p[col["ref"]], alt=p[col["alt"]],
                    pos=int(p[col["position"]]),
                    maf=float(p[col["maf"]]), beta=float(p[col["beta"]]),
                    se=float(p[col["se"]]), pval=float(p[col["pvalue"]]),
                    an=int(p[col["an"]]), ac=int(p[col["ac"]])))
            except ValueError:
                continue
    return rows


def main():
    log = []

    def P(s=""):
        log.append(s)

    os.makedirs(os.path.join(W, "gtex_coloc_inputs"), exist_ok=True)
    results = []

    for tissue in TISSUES:
        for gene, (locus, outcomes) in GENES.items():
            rows = load_gtex(os.path.join(W, "gtex_%s_%s.tsv" % (tissue, gene)))
            P("=" * 100)
            P("Tissue: %s (n=%d)  Gene: %s  variants: %d"
              % (tissue, SAMPLE_N[tissue], gene, len(rows)))

            # lead cis-eQTL
            lead = min(rows, key=lambda r: r["pval"])
            P("  Lead cis-eQTL: %s  chr%s:%d %s/%s  beta=%.4g  se=%.4g  "
              "z=%.2f  p=%.3g  maf=%.4g"
              % (lead["rsid"], "NA", lead["pos"], lead["ref"], lead["alt"],
                 lead["beta"], lead["se"], lead["beta"] / lead["se"],
                 lead["pval"], lead["maf"]))

            # key SNP lookups
            byrsid = {}
            for r in rows:
                byrsid.setdefault(r["rsid"], r)
            for s in KEY_SNPS[gene]:
                if s in byrsid:
                    r = byrsid[s]
                    P("  %-12s beta=%+.4g se=%.4g z=%+.2f p=%.3g maf=%.4g "
                      "an=%d" % (s, r["beta"], r["se"], r["beta"] / r["se"],
                                 r["pval"], r["maf"], r["an"]))
                else:
                    P("  %-12s (not present)" % s)

            for out in outcomes:
                ma = q.load_ma(os.path.join(W, "ma_%s_%s.ma" % (out, locus)))
                exp_rows = []   # harmonized exposure (GTEx, beta per ALT)
                out_rows = []
                ndropped = 0
                nfreq = 0
                for r in rows:
                    if r["rsid"] not in ma or r["se"] <= 0:
                        continue
                    g = ma[r["rsid"]]
                    # GTEx effect allele = ALT; align outcome to ALT
                    if r["alt"] == g["a1"] and r["ref"] == g["a2"]:
                        by, fy = g["b"], g["freq"]
                    elif r["alt"] == g["a2"] and r["ref"] == g["a1"]:
                        by, fy = -g["b"], 1.0 - g["freq"]
                    else:
                        ndropped += 1
                        continue
                    # exposure ALT-af = ac/an (verified empirically:
                    # maf column is minor-allele freq, ac is ALT count)
                    eaf_alt = r["ac"] / r["an"]
                    # freq QC vs outcome (|diff|>0.25 drop, as step4)
                    if abs(eaf_alt - fy) > 0.25:
                        nfreq += 1
                        continue
                    exp_rows.append(dict(rsid=r["rsid"], beta=r["beta"],
                                         se=r["se"], eaf=eaf_alt,
                                         n=SAMPLE_N[tissue],
                                         pval=r["pval"]))
                    out_rows.append(dict(rsid=r["rsid"], beta=by, se=g["se"],
                                         eaf=fy, n=g["n"], pval=g["pval"]))
                pair = "gtex_%s_%s_%s" % (tissue, gene.lower(), out)
                if len(exp_rows) < 50:
                    P("  %s: only %d overlapping SNPs; coloc skipped"
                      % (out, len(exp_rows)))
                    continue
                # write coloc inputs
                d = os.path.join(W, "gtex_coloc_inputs")
                with open(os.path.join(d, pair + "_exposure.csv"), "w",
                          newline="") as f:
                    wtr = csv.DictWriter(f, fieldnames=["rsid", "beta", "se",
                                                        "eaf", "n", "pval"])
                    wtr.writeheader()
                    wtr.writerows(exp_rows)
                with open(os.path.join(d, pair + "_outcome.csv"), "w",
                          newline="") as f:
                    wtr = csv.DictWriter(f, fieldnames=["rsid", "beta", "se",
                                                        "eaf", "n", "pval"])
                    wtr.writeheader()
                    wtr.writerows(out_rows)

                # lead-eQTL Wald: strongest eQTL among OVERLAPPING variants
                # (dataset lead may be a rare GTEx-specific variant absent
                # from the outcome GWAS)
                lead_h, lead_z = None, 0.0
                for r in exp_rows:
                    z = abs(r["beta"] / r["se"])
                    if z > lead_z:
                        lead_z = z
                        lead_h = r
                if lead_h and abs(lead_h["beta"]) > 0:
                    o = next(o for o in out_rows
                             if o["rsid"] == lead_h["rsid"])
                    r = lead_h
                    b = o["beta"] / r["beta"]
                    se = math.sqrt((o["se"] ** 2 * r["beta"] ** 2
                                    + r["se"] ** 2 * o["beta"] ** 2)
                                   / r["beta"] ** 4)
                    z = b / se
                    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
                    P("  %-11s lead-Wald at %s (z_eQTL=%.2f): b=%.4g "
                      "se=%.4g z=%.2f p=%.4g (n_overlap=%d, "
                      "allele-mismatch dropped=%d, freq-QC dropped=%d)"
                      % (out, r["rsid"], lead_z, b, se, z, p, len(exp_rows),
                         ndropped, nfreq))
                    results.append(dict(tissue=tissue, gene=gene,
                                        outcome=out, leadSNP=r["rsid"],
                                        lead_p_eQTL=r["pval"],
                                        b_Wald=b, se_Wald=se, p_Wald=p,
                                        nsnps=len(exp_rows)))
                else:
                    P("  %-11s no usable overlapping lead variant "
                      "(n_overlap=%d)" % (out, len(exp_rows)))
                # also rs5068-specific Wald if present (NPPA only)
                if gene == "NPPA" and "rs5068" in byrsid:
                    for r, o in zip(exp_rows, out_rows):
                        if r["rsid"] == "rs5068" and abs(r["beta"]) > 0:
                            b = o["beta"] / r["beta"]
                            se = math.sqrt(
                                (o["se"] ** 2 * r["beta"] ** 2
                                 + r["se"] ** 2 * o["beta"] ** 2)
                                / r["beta"] ** 4)
                            z = b / se
                            p = 2 * (1 - 0.5
                                     * (1 + math.erf(abs(z) / math.sqrt(2))))
                            P("  %-11s rs5068-Wald: b=%.4g se=%.4g z=%.2f "
                              "p=%.4g" % (out, b, se, z, p))

    with open(os.path.join(W, "step14f_cismr_results.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["tissue", "gene", "outcome",
                                          "leadSNP", "lead_p_eQTL", "b_Wald",
                                          "se_Wald", "p_Wald", "nsnps"])
        w.writeheader()
        w.writerows(results)
    with open(os.path.join(W, "step14f_output.txt"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(log))
    print("done: %d lead-Wald pairs" % len(results))


def _maf_of(eaf):
    return min(eaf, 1 - eaf)


if __name__ == "__main__":
    main()
