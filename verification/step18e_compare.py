# -*- coding: utf-8 -*-
"""
step18e: Cross-validation — official smr-1.3.1-win.exe results vs in-house
Python port (step13q top-SNP / step13r refSNP rs5068).
"""
import os
import csv

V = r"C:\Users\xuyin\WorkBuddy\2026-09-13-15-20-56\smr_verify"
SRC = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"

NPPA_OUT = ["AF", "CAD_cardio", "CAD_vdh", "CES", "HF", "stroke"]
NPR3_OUT = ["CAD_cardio", "CAD_vdh", "MI_har", "MI_finn"]


def read_smr(tag):
    p = os.path.join(V, tag + ".smr")
    if not os.path.exists(p):
        return None
    with open(p) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return rows[0] if rows else None


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def cmp_val(a, b, rel_tol=5e-3, label=""):
    """relative tolerance comparison; float32 container -> small drift ok"""
    if a is None or b is None:
        return ("MISS" if (a is None) != (b is None) else "bothNA")
    if a == 0 and b == 0:
        return "OK"
    denom = max(abs(a), abs(b), 1e-12)
    if abs(a - b) / denom <= rel_tol:
        return "OK"
    return "DIFF %.4g vs %.4g" % (a, b)


def main():
    lines = []

    def P(s=""):
        lines.append(s)
        print(s)

    # ---- load in-house port results ----
    port_top, port_ref = {}, {}
    with open(os.path.join(SRC, "step13q_smr_heidi_results.csv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            port_top[(r["exposure"], r["outcome"])] = r
    with open(os.path.join(SRC, "step13r_refsmr_results.csv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            port_ref[(r["exposure"], r["outcome"])] = r

    out_rows = []
    P("=" * 110)
    P("OFFICIAL smr-1.3.1-win.exe  vs  in-house Python port (v1.4.2 source)")
    P("tolerance: rel 5e-3 on b_SMR/se_SMR; p-values compared at 2 significant digits")
    P("=" * 110)

    for mode, port_map, pairs in (
            ("top", port_top, [(e, o) for e in ("SCALLOP", "INTERVAL", "NPPAeQTL")
                                for o in NPPA_OUT] + [("NPR3eQTL", o) for o in NPR3_OUT]),
            ("ref", port_ref, [(e, o) for e in ("SCALLOP", "INTERVAL", "NPPAeQTL")
                               for o in NPPA_OUT])):
        P("")
        P("### mode = %s (%s)" % (mode, "default top-SNP" if mode == "top"
                                  else "--target-snp rs5068"))
        for exp, out in pairs:
            tag = "%s_%s_x_%s" % (mode, exp, out)
            row = read_smr(tag)
            pr = port_map.get((exp, out))
            if row is None:
                P("%-10s x %-10s : NO OFFICIAL OUTPUT" % (exp, out))
                continue
            if mode == "top":
                inst = row.get("topSNP", "")
            else:
                inst = row.get("targetSNP", row.get("topSNP", ""))
            checks = [
                ("inst", inst, pr.get("topSNP") if mode == "top" else "rs5068"),
                ("b_SMR", f(row.get("b_SMR")), f(pr.get("b_SMR"))),
                ("se_SMR", f(row.get("se_SMR")), f(pr.get("se_SMR"))),
                ("p_SMR", f(row.get("p_SMR")), f(pr.get("p_SMR"))),
                ("p_HEIDI", f(row.get("p_HEIDI")), f(pr.get("p_HEIDI"))),
                ("nsnp", f(row.get("nsnp_HEIDI")), f(pr.get("nsnp_HEIDI"))),
            ]
            verdict = []
            for name, a, b in checks:
                if name == "inst":
                    verdict.append("%s=%s%s" % (name, a, "" if a == b else " (port %s)" % b))
                elif name in ("p_SMR", "p_HEIDI"):
                    if a is None or b is None:
                        verdict.append("%s=%s/%s" % (name, a, b))
                    else:
                        ok = (abs(a - b) / max(abs(a), abs(b), 1e-12)) <= 0.05
                        verdict.append("%s %s (%.4g vs %.4g)" % (name, "OK" if ok else "DIFF", a, b))
                else:
                    verdict.append("%s %s" % (name, cmp_val(a, b)))
            P("%-10s x %-10s | %s" % (exp, out, " | ".join(verdict)))
            out_rows.append(dict(mode=mode, exposure=exp, outcome=out,
                                 official_inst=inst, port_inst=pr.get("topSNP"),
                                 official_bSMR=row.get("b_SMR"),
                                 official_seSMR=row.get("se_SMR"),
                                 official_pSMR=row.get("p_SMR"),
                                 official_pHEIDI=row.get("p_HEIDI"),
                                 official_nsnp=row.get("nsnp_HEIDI"),
                                 port_bSMR=pr.get("b_SMR"), port_seSMR=pr.get("se_SMR"),
                                 port_pSMR=pr.get("p_SMR"), port_pHEIDI=pr.get("p_HEIDI"),
                                 port_nsnp=pr.get("nsnp_HEIDI")))

    keys = list(out_rows[0].keys())
    with open(os.path.join(V, "step18e_comparison.csv"), "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(out_rows)
    with open(os.path.join(V, "step18e_comparison_report.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print("\nSaved: step18e_comparison.csv / step18e_comparison_report.txt")


if __name__ == "__main__":
    main()
