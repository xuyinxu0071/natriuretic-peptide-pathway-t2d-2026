# -*- coding: utf-8 -*-
"""
Step 13b (P4): 构建 SMR 官方分析所需输入
================================================================================
1) 1000G v5b EUR 区域 VCF (仅 dbSNP 匹配变异, ID=rsid) -> plink --vcf -> .bed/.bim/.fam
2) 结局 .ma 文件 (COJO/SMR 格式: SNP A1 A2 freq b se p n, ALT 方向)
3) 暴露 .efile 文件 (ProbeID SNP Chr BP A1 A2 Freq b se p) -> smr --make-besd
"""
import gzip, json, math, os, re, struct, time, urllib.request, zlib
from step4_harmonize import harmonize

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
BIN = os.path.join(WORK, "bin")
BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
PANEL = "integrated_call_samples_v3.20130502.ALL.panel"
VCF = {
    "nppa": "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    "npr3": "ALL.chr5.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
}
REGIONS = {
    "nppa": ("1", 11805974, 12005974),
    "npr3": ("5", 32614270, 32814270),
}

EXPOSURES = {
    "nppa": {
        "ebi-a-GCST90012082": ("NTproBNP_SCALLOP", "NT-proBNP pQTL (SCALLOP)"),
        "prot-a-2078": ("NTproBNP_INTERVAL", "NT-proBNP pQTL (INTERVAL)"),
        "eqtl-a-ENSG00000175206": ("NPPA_eQTL", "NPPA eQTL (eQTLGen)"),
    },
    "npr3": {
        "eqtl-a-ENSG00000113389": ("NPR3_eQTL", "NPR3 eQTL (eQTLGen)"),
    },
}
OUTCOMES = {
    "nppa": ["ieu-a-7", "ebi-a-GCST005195", "ebi-a-GCST006061",
             "ebi-a-GCST009541", "ebi-a-GCST005838", "ebi-a-GCST006910"],
    "npr3": ["ieu-a-7", "ebi-a-GCST005195", "ebi-a-GCST011364", "finn-b-I9_MI"],
}
DSNAME = {
    "ieu-a-7": "CAD_cardio", "ebi-a-GCST005195": "CAD_vdh",
    "ebi-a-GCST006061": "AF", "ebi-a-GCST009541": "HF",
    "ebi-a-GCST005838": "stroke", "ebi-a-GCST006910": "CES",
    "ebi-a-GCST011364": "MI_har", "finn-b-I9_MI": "MI_finn",
    "ebi-a-GCST90012082": "SCALLOP", "prot-a-2078": "INTERVAL",
    "eqtl-a-ENSG00000175206": "NPPAeQTL", "eqtl-a-ENSG00000113389": "NPR3eQTL",
}

def fetch(url, start=None, end=None, retries=5, timeout=1200):
    for i in range(retries):
        try:
            hdr = {"User-Agent": "coloc-pipeline/1.0"}
            if start is not None:
                hdr["Range"] = f"bytes={start}-{end}"
            req = urllib.request.Request(url, headers=hdr)
            return urllib.request.urlopen(req, timeout=timeout).read()
        except Exception:
            if i == retries - 1:
                raise
            time.sleep(8 * (i + 1))

def reg2bins(beg, end):
    lst = [0]; end -= 1
    for k in range(1 + (beg >> 26), 1 + (end >> 26) + 1): lst.append(k)
    for k in range(9 + (beg >> 23), 9 + (end >> 23) + 1): lst.append(k)
    for k in range(73 + (beg >> 20), 73 + (end >> 20) + 1): lst.append(k)
    for k in range(585 + (beg >> 17), 585 + (end >> 17) + 1): lst.append(k)
    for k in range(4681 + (beg >> 14), 4681 + (end >> 14) + 1): lst.append(k)
    return lst

def parse_tbi(path):
    raw = gzip.open(path, "rb").read()
    off = 4
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip = struct.unpack_from("<iiiiiii", raw, off); off += 28
    l_nm = struct.unpack_from("<i", raw, off)[0]; off += 4
    names = raw[off:off + l_nm].decode().rstrip("\0").split("\0"); off += l_nm
    refs = {}
    for r in range(n_ref):
        n_bin = struct.unpack_from("<i", raw, off)[0]; off += 4
        bins = []
        for b in range(n_bin):
            bin_no, n_chunk = struct.unpack_from("<Ii", raw, off); off += 8
            chunks = []
            for c in range(n_chunk):
                cb, ce = struct.unpack_from("<QQ", raw, off); off += 16
                chunks.append((cb, ce))
            bins.append((bin_no, chunks))
        n_intv = struct.unpack_from("<i", raw, off)[0]; off += 4
        off += 8 * n_intv
        refs[names[r]] = bins
    return refs

def decompress_bgzf(data):
    out, pos = [], 0
    while pos < len(data):
        d = zlib.decompressobj(47)
        out.append(d.decompress(data[pos:]))
        unused = d.unused_data
        if not unused:
            break
        pos = len(data) - len(unused)
    return b"".join(out)

import gzip as _gz

def get_eur_samples():
    p = os.path.join(WORK, "kg_panel.txt")
    if not os.path.exists(p):
        open(p, "wb").write(fetch(BASE + PANEL))
    return [l.split("\t")[0] for l in open(p).read().splitlines()
            if l and not l.startswith("sample") and len(l.split("\t")) >= 3 and l.split("\t")[2] == "EUR"]

def get_eur_indices(vcf_url, eur_set):
    data = fetch(vcf_url, 0, 400000)
    d = zlib.decompressobj(47)
    text = (d.decompress(data) + d.flush()).decode("ascii", "ignore")
    for line in text.split("\n"):
        if line.startswith("#CHROM"):
            samples = line.split("\t")[9:]
            idx = [i for i, s in enumerate(samples) if s in eur_set]
            return samples, idx
    raise RuntimeError("no #CHROM header")

OUT = []
def P(s):
    OUT.append(str(s))

eur = get_eur_samples(); eur_set = set(eur)
P(f"EUR samples: {len(eur)}")

for locus, (chrom, lo, hi) in REGIONS.items():
    ldvars = json.load(open(os.path.join(WORK, f"ld_variants_{locus}.json")))  # rsid -> pos/ref/alt
    keymap = {(v["pos"], v["ref"], v["alt"]): r for r, v in ldvars.items()}
    P(f"\n=== {locus}: chr{chrom}:{lo}-{hi}, {len(keymap)} 目标变异 ===")

    tbi = os.path.join(WORK, VCF[locus] + ".tbi")
    refs = parse_tbi(tbi)
    bins = refs[chrom]
    wbins = set(reg2bins(lo - 1, hi))
    chunks = []
    for bin_no, chs in bins:
        if bin_no in wbins:
            chunks.extend(chs)
    ranges = sorted(set((cb >> 16, (ce >> 16) + 70000) for cb, ce in chunks))
    merged = []
    for s, e in ranges:
        if merged and s <= merged[-1][1] + 100000:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    P(f"  {len(chunks)} chunks -> {len(merged)} ranges ~{sum(e-s for s,e in merged)/1e6:.1f} MB")

    samples, eur_idx = get_eur_indices(BASE + VCF[locus], eur_set)
    P(f"  VCF {len(samples)} 样本, EUR {len(eur_idx)}")
    eur_names = [samples[i] for i in eur_idx]

    got = {}  # rsid -> (pos, ref, alt, [gts])
    for s, e in merged:
        data = fetch(BASE + VCF[locus], s, e)
        text = decompress_bgzf(data).decode("ascii", "ignore")
        for line in text.split("\n"):
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 10 or parts[0] != chrom:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            if pos < lo or pos > hi:
                continue
            key = (pos, parts[3], parts[4])
            if key not in keymap:
                continue
            r = keymap[key]
            if r in got:
                continue
            gts = []
            ok = True
            for i in eur_idx:
                gt = parts[9 + i].split(":")[0]
                a1, a2 = (gt.replace("|", "/").split("/") + ["/"])[:2]
                if a1 in (".", "/") or a2 in (".", "/"):
                    ok = False; break
                gts.append(a1 + "|" + a2)
            if not ok:
                continue
            got[r] = (pos, parts[3], parts[4], gts)
        del text, data
    P(f"  匹配 {len(got)} 变异, 写 EUR VCF")

    vcf_path = os.path.join(WORK, f"eur_{locus}.vcf")
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID={chrom}>\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(eur_names) + "\n")
        for r in sorted(got, key=lambda x: got[x][0]):
            pos, ref, alt, gts = got[r]
            f.write(f"{chrom}\t{pos}\t{r}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")

    # ---- 暴露 .efile 与结局 .ma ----
    variants_full = json.load(open(os.path.join(WORK, f"{locus}_locus_variants.json")))
    for dsid, (probe, plabel) in EXPOSURES[locus].items():
        f = os.path.join(WORK, f"{locus}__{dsid}_assoc.json")
        if not os.path.exists(f):
            P(f"  {dsid}: 缺失"); continue
        h = harmonize(variants_full, json.load(open(f)))
        h = {r: v for r, v in h.items() if r in got}
        efile = os.path.join(WORK, f"{DSNAME[dsid]}.efile")
        with open(efile, "w") as fo:
            fo.write("ProbeID\tSNP\tChr\tBP\tA1\tA2\tFreq\tb\tse\tp\n")
            for r in sorted(h, key=lambda x: h[x][0]):
                pos, b, se, eaf, n, alt, ref = h[r]
                p = math.erfc(abs(b / se) / math.sqrt(2))
                fo.write(f"{probe}\t{r}\t{chrom}\t{pos}\t{alt}\t{ref}\t{eaf:.4f}\t{b:.6g}\t{se:.6g}\t{p:.6g}\n")
        P(f"  {plabel} [{dsid}]: {len(h)} 变异 -> {DSNAME[dsid]}.efile")

    for dsid in OUTCOMES[locus]:
        f = os.path.join(WORK, f"{locus}__{dsid}_assoc.json")
        if not os.path.exists(f):
            P(f"  {dsid}: 缺失"); continue
        h = harmonize(variants_full, json.load(open(f)))
        h = {r: v for r, v in h.items() if r in got}
        ma = os.path.join(WORK, f"ma_{DSNAME[dsid]}_{locus}.ma")
        with open(ma, "w") as fo:
            fo.write("SNP\tA1\tA2\tfreq\tb\tse\tp\tn\n")
            for r in sorted(h, key=lambda x: h[x][0]):
                pos, b, se, eaf, n, alt, ref = h[r]
                p = math.erfc(abs(b / se) / math.sqrt(2))
                fo.write(f"{r}\t{alt}\t{ref}\t{eaf:.4f}\t{b:.6g}\t{se:.6g}\t{p:.6g}\t{n if n else 'NA'}\n")
        P(f"  {DSNAME[dsid]} [{dsid}]: {len(h)} 变异 -> ma_{DSNAME[dsid]}_{locus}.ma")

with open(os.path.join(WORK, "step13b_output.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT))
print("DONE")
