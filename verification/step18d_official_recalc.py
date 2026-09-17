# -*- coding: utf-8 -*-
"""
step18d: Official SMR recalculation — build dense BESD from ESD, then run the
OFFICIAL smr-1.3.1-win.exe analysis path on all exposure x outcome pairs.

Background:
  - Official prebuilt Windows binary smr-1.3.1-win.exe runs fine, but its
    `--make-besd` (sparse) code path crashes with 0xC0000005 access violation
    (confirmed with a minimal 3-SNP ESD, pure-ASCII paths, single thread, and
    all OpenMP runtime variants; `--make-besd-dense` is a BESD->dense
    conversion route, not an ESD builder).
  - Solution: write the dense BESD container ourselves following the official
    write routine (make_full_besd, SMR_data.cpp v1.4.2/master: 10/16-int
    reserved header + per-probe [beta(snpNum), se(snpNum)] float32 blocks),
    using DENSE_FILE_TYPE_1 (gflag=0, no reserved units) which is read by both
    1.3.x and 1.4.x. All STATISTICS are then computed by the official binary.
  - Round-trip validation: the official analysis reproduces known values only
    if the container is byte-correct; additionally --query dumps container
    contents for direct comparison with the source ESD.

Data sources (all real, from coloc_analysis/):
  - {SCALLOP,INTERVAL,NPPAeQTL,NPR3eQTL}.efile  (cis-window eQTL/pQTL stats)
  - ma_{outcome}_{locus}.ma                     (outcome GWAS summary stats)
  - eur_{nppa,npr3}.bed/.bim/.fam               (1000G EUR LD reference)
Reference targets:
  - step13q_smr_heidi_results.csv (top-SNP SMR/HEIDI, in-house port)
  - step13r_refsmr_results.csv    (refSNP rs5068 SMR/HEIDI, in-house port)
"""
import os
import shutil
import struct
import subprocess
import csv

SRC = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
V = r"C:\Users\xuyin\WorkBuddy\2026-09-13-15-20-56\smr_verify"
EXE = os.path.join(V, "smr-1.3.1-win.exe")

NPPA_OUT = ["AF", "CAD_cardio", "CAD_vdh", "CES", "HF", "stroke"]
NPR3_OUT = ["CAD_cardio", "CAD_vdh", "MI_har", "MI_finn"]
GENE = {"SCALLOP": "NPPA", "INTERVAL": "NPPA", "NPPAeQTL": "NPPA",
        "NPR3eQTL": "NPR3"}

LOG = []


def P(s=""):
    LOG.append(str(s))
    print(s)


def dec(b):
    """Official exe prints UTF-16LE to a pipe."""
    for enc in ("utf-16", "utf-8", "gbk"):
        s = (b or b"").decode(enc, errors="ignore")
        if len(s) > 20:
            return s
    return ""


def parse_efile(path):
    rows = []
    with open(path, encoding="utf-8-sig") as f:
        header = f.readline().rstrip("\n\r").split("\t")
        assert header[:10] == ["ProbeID", "SNP", "Chr", "BP", "A1", "A2",
                               "Freq", "b", "se", "p"], header
        for line in f:
            if not line.strip():
                continue
            c = line.rstrip("\n\r").split("\t")
            rows.append(dict(probe=c[0], snp=c[1], chr=c[2], bp=int(c[3]),
                             a1=c[4], a2=c[5], freq=float(c[6]),
                             b=float(c[7]), se=float(c[8]), p=float(c[9])))
    return rows


def build_besd(exp):
    """Write exp.esi / exp.epi / exp.besd (dense, DENSE_FILE_TYPE_1=0)."""
    rows = parse_efile(os.path.join(SRC, exp + ".efile"))
    probes = sorted(set(r["probe"] for r in rows))
    assert len(probes) == 1, "expected single probe per efile"

    # esi order = order of first appearance
    seen, snps = set(), []
    for r in rows:
        if r["snp"] not in seen:
            seen.add(r["snp"])
            snps.append(r)
    by_snp = {r["snp"]: r for r in rows}

    with open(os.path.join(V, exp + ".esi"), "w") as f:
        for r in snps:
            f.write("\t".join([r["chr"], r["snp"], "0", str(r["bp"]),
                               r["a1"], r["a2"], "%.6f" % r["freq"]]) + "\n")

    probe_bp = sum(r["bp"] for r in snps) // len(snps)
    with open(os.path.join(V, exp + ".epi"), "w") as f:
        f.write("\t".join([snps[0]["chr"], probes[0], "0", str(probe_bp),
                           GENE[exp], "+"]) + "\n")

    with open(os.path.join(V, exp + ".besd"), "wb") as f:
        f.write(struct.pack("<I", 0))                    # DENSE_FILE_TYPE_1
        for r in snps:
            f.write(struct.pack("<f", by_snp[r["snp"]]["b"]))
        for r in snps:
            f.write(struct.pack("<f", by_snp[r["snp"]]["se"]))
    return len(snps)


def run_smr(args, tag):
    """Run official SMR; return (rc, parsed .smr rows or None)."""
    cmd = [EXE] + args
    r = subprocess.run(cmd, capture_output=True, timeout=1800, cwd=V)
    smrfile = os.path.join(V, tag + ".smr")
    parsed = None
    if os.path.exists(smrfile):
        with open(smrfile) as f:
            parsed = list(csv.DictReader(f, delimiter="\t"))
    tail = (dec(r.stdout) + " STDERR:" + dec(r.stderr)).replace("\r", "")
    P("  [%s] rc=%d smr_rows=%s" % (tag, r.returncode,
                                    len(parsed) if parsed else 0))
    if r.returncode != 0 or not parsed:
        P("    tail: " + tail.replace("\n", " | ")[-320:])
    return r.returncode, parsed


def main():
    P("=" * 100)
    P("OFFICIAL SMR RECALCULATION - smr-1.3.1-win.exe (Zhang/Zhu/Yang, UQ)")
    P("Workspace: " + V)
    P("")

    # 1) copy inputs (ASCII workspace)
    for e in GENE:
        shutil.copy(os.path.join(SRC, e + ".efile"), os.path.join(V, e + ".efile"))
    ma_needed = {"ma_%s_%s.ma" % (o, "nppa") for o in NPPA_OUT} | \
                {"ma_%s_%s.ma" % (o, "npr3") for o in NPR3_OUT}
    for m in sorted(ma_needed):
        shutil.copy(os.path.join(SRC, m), os.path.join(V, m))
    for loc in ("nppa", "npr3"):
        for ext in ("bed", "bim", "fam"):
            shutil.copy(os.path.join(SRC, "eur_%s.%s" % (loc, ext)),
                        os.path.join(V, "eur_%s.%s" % (loc, ext)))

    # 2) build BESD containers
    P("[Step 1] Building dense BESD containers (DENSE_FILE_TYPE_1) from ESD")
    for e in GENE:
        n = build_besd(e)
        P("  %s: %d SNPs, 1 probe -> %s.{esi,epi,besd}" % (e, n, e))
    P("")

    # 3) round-trip check via --query (optional; tolerate absence)
    P("[Step 2] Round-trip container check (--query)")
    for e in GENE:
        rc, _ = run_smr(["--beqtl-summary", e, "--query", "--out",
                         "q_" + e], "q_" + e)
    P("")

    # 4) official analyses
    P("[Step 3] Official SMR/HEIDI runs")
    results = []

    # 4a) NPPA pairs: default top-SNP + refSNP rs5068
    for exp in ("SCALLOP", "INTERVAL", "NPPAeQTL"):
        for out in NPPA_OUT:
            base = ["--bfile", "eur_nppa", "--gwas-summary",
                    "ma_%s_nppa.ma" % out, "--beqtl-summary", exp,
                    "--thread-num", "1"]
            tag = "top_%s_x_%s" % (exp, out)
            rc, rows = run_smr(base + ["--out", tag], tag)
            for row in (rows or []):
                results.append(dict(mode="top", exposure=exp, outcome=out,
                                    locus="nppa", **row))
            tag2 = "ref_%s_x_%s" % (exp, out)
            rc2, rows2 = run_smr(base + ["--target-snp", "rs5068",
                                         "--out", tag2], tag2)
            for row in (rows2 or []):
                results.append(dict(mode="ref", exposure=exp, outcome=out,
                                    locus="nppa", **row))

    # 4b) NPR3 pairs: relaxed peqtl-smr=1 (port computed despite weak IV)
    for out in NPR3_OUT:
        base = ["--bfile", "eur_npr3", "--gwas-summary",
                "ma_%s_npr3.ma" % out, "--beqtl-summary", "NPR3eQTL",
                "--peqtl-smr", "1", "--thread-num", "1"]
        tag = "top_NPR3eQTL_x_%s" % out
        rc, rows = run_smr(base + ["--out", tag], tag)
        for row in (rows or []):
            results.append(dict(mode="top-relaxed", exposure="NPR3eQTL",
                                outcome=out, locus="npr3", **row))

    # 5) save
    out_csv = os.path.join(V, "step18d_official_results.csv")
    if results:
        keys = list(results[0].keys())
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(results)
    P("")
    P("Official results written: %s (%d rows)" % (out_csv, len(results)))

    with open(os.path.join(V, "step18d_output.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(LOG))


if __name__ == "__main__":
    main()
