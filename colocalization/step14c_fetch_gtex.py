# -*- coding: utf-8 -*-
"""
step14c_fetch_gtex.py
Fetch GTEx heart-tissue eQTL summary statistics (eQTL Catalogue r7,
tabix-indexed TSV over HTTPS) for the NPPA/NPPB and NPR3 cis regions.

Datasets (from dataset_metadata_r7.tsv):
  QTD000251  GTEx heart_atrial_appendage  ge  n=372
  QTD000256  GTEx heart_left_ventricle    ge  n=382

Regions: chr1:11805974-12005974 (NPPA ENSG00000175206 / NPPB ENSG00000183705)
         chr5:32614270-32814270 (NPR3 ENSG00000113389)

Output: gtex_{tissue}_{gene}.tsv (one row per variant, deduplicated)
"""
import gzip
import json
import os
import struct
import time
import urllib.request
import zlib

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/"

DATASETS = {
    "atrial": "QTD000251",
    "ventricle": "QTD000256",
}
REGIONS = [
    ("1", 11805974, 12005974, ["ENSG00000175206", "ENSG00000183705"]),
    ("5", 32614270, 32814270, ["ENSG00000113389"]),
]
GENE_NAMES = {
    "ENSG00000175206": "NPPA",
    "ENSG00000183705": "NPPB",
    "ENSG00000113389": "NPR3",
}


def fetch(url, start=None, end=None, retries=5, timeout=1800):
    for i in range(retries):
        try:
            hdr = {"User-Agent": "coloc-pipeline/1.0"}
            if start is not None:
                hdr["Range"] = "bytes=%d-%d" % (start, end)
            req = urllib.request.Request(url, headers=hdr)
            return urllib.request.urlopen(req, timeout=timeout).read()
        except Exception:
            if i == retries - 1:
                raise
            time.sleep(8 * (i + 1))


def reg2bins(beg, end):
    lst = [0]
    end -= 1
    for k in range(1 + (beg >> 26), 1 + (end >> 26) + 1):
        lst.append(k)
    for k in range(9 + (beg >> 23), 9 + (end >> 23) + 1):
        lst.append(k)
    for k in range(73 + (beg >> 20), 73 + (end >> 20) + 1):
        lst.append(k)
    for k in range(585 + (beg >> 17), 585 + (end >> 17) + 1):
        lst.append(k)
    for k in range(4681 + (beg >> 14), 4681 + (end >> 14) + 1):
        lst.append(k)
    return lst


def parse_tbi(path):
    raw = gzip.open(path, "rb").read()
    off = 4
    (n_ref, fmt, col_seq, col_beg, col_end, meta, skip) = struct.unpack_from(
        "<iiiiiii", raw, off)
    off += 28
    l_nm = struct.unpack_from("<i", raw, off)[0]
    off += 4
    names = raw[off:off + l_nm].decode().rstrip("\0").split("\0")
    off += l_nm
    refs = {}
    for r in range(n_ref):
        n_bin = struct.unpack_from("<i", raw, off)[0]
        off += 4
        bins = []
        for b in range(n_bin):
            bin_no, n_chunk = struct.unpack_from("<Ii", raw, off)
            off += 8
            chunks = []
            for c in range(n_chunk):
                cb, ce = struct.unpack_from("<QQ", raw, off)
                off += 16
                chunks.append((cb, ce))
            bins.append((bin_no, chunks))
        n_intv = struct.unpack_from("<i", raw, off)[0]
        off += 4
        off += 8 * n_intv
        refs[names[r]] = bins
    return refs, col_seq, col_beg, meta, skip


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


def main():
    log = open(os.path.join(W, "step14c_log.txt"), "w", encoding="utf-8")

    def P(s):
        print(s, flush=True)
        log.write(s + "\n")
        log.flush()

    for tissue, qtd in DATASETS.items():
        gz_url = BASE + qtd + "/" + qtd + ".all.tsv.gz"
        tbi_path = os.path.join(W, qtd + ".all.tsv.gz.tbi")
        if not os.path.exists(tbi_path) or os.path.getsize(tbi_path) < 100:
            P("downloading tbi for %s ..." % qtd)
            open(tbi_path, "wb").write(fetch(gz_url + ".tbi"))
        refs, col_seq, col_beg, meta, skip = parse_tbi(tbi_path)
        P("%s (%s): tbi refs=%s col_seq=%d col_beg=%d meta=%r skip=%d"
          % (tissue, qtd, list(refs.keys()), col_seq, col_beg, chr(meta), skip))

        for chrom, lo, hi, genes in REGIONS:
            need = set(GENE_NAMES[g] for g in genes)
            bins = refs.get(chrom, refs.get("chr" + chrom, []))
            if not bins:
                P("  chr%s not in tbi; skip" % chrom)
                continue
            wbins = set(reg2bins(lo - 1, hi))
            chunks = []
            for bin_no, chs in bins:
                if bin_no in wbins:
                    chunks.extend(chs)
            ranges = sorted(set((cb >> 16, (ce >> 16) + 70000)
                                for cb, ce in chunks))
            merged = []
            for s, e in ranges:
                if merged and s <= merged[-1][1] + 100000:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            total = sum(e - s for s, e in merged)
            P("  chr%s:%d-%d: %d chunks -> %d ranges, ~%.1f MB"
              % (chrom, lo, hi, len(chunks), len(merged), total / 1e6))

            rows = {}
            header = None
            for s, e in merged:
                data = fetch(gz_url, s, e)
                text = decompress_bgzf(data).decode("utf-8", "replace")
                for line in text.split("\n"):
                    if not line or line.startswith("#"):
                        if line.startswith("#") and header is None:
                            header = line.lstrip("#").split("\t")
                        continue
                    p = line.split("\t")
                    if len(p) < 12:
                        continue
                    # columns (0-based): 0=molecular_trait_id (gene),
                    # 1=chromosome, 2=position
                    try:
                        c = p[col_seq - 1]
                        pos = int(p[col_beg - 1])
                    except (ValueError, IndexError):
                        continue
                    if c != chrom and c != "chr" + chrom:
                        continue
                    if pos < lo or pos > hi:
                        continue
                    gid = p[0]
                    if gid not in genes:
                        continue
                    key = (gid, p[5])  # variant id dedup
                    if key not in rows:
                        rows[key] = p
                del text, data
            if header is None:
                header = ("molecular_trait_id chromosome position ref alt "
                          "variant ma_samples maf pvalue beta se type ac an "
                          "r2 molecular_trait_object_id gene_id median_tpm "
                          "rsid").split()
                P("  header not in chunks; using canonical eQTL Catalogue "
                  "column layout")
            P("  header: %s" % (header[:20] if header else "NA"))
            P("  rows kept (dedup): %d" % len(rows))

            for gene in genes:
                gid_rows = [p for (g, _), p in rows.items() if g == gene]
                if not gid_rows:
                    P("  %s: NO ROWS" % gene)
                    continue
                out = os.path.join(W, "gtex_%s_%s.tsv"
                                   % (tissue, GENE_NAMES[gene]))
                with open(out, "w", encoding="utf-8") as f:
                    if header:
                        f.write("\t".join(header) + "\n")
                    for p in sorted(gid_rows, key=lambda x: int(x[col_beg - 1])):
                        f.write("\t".join(p) + "\n")
                P("  saved %s: %d variants" % (out, len(gid_rows)))
    P("DONE")
    log.close()


if __name__ == "__main__":
    main()
