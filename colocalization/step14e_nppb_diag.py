# -*- coding: utf-8 -*-
"""step14e: diagnose NPPB absence in GTEx ge datasets."""
import gzip
import os
import struct
import urllib.request
import zlib

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/"
out = []


def fetch(url, start=None, end=None):
    hdr = {"User-Agent": "coloc-pipeline/1.0"}
    if start is not None:
        hdr["Range"] = "bytes=%d-%d" % (start, end)
    req = urllib.request.Request(url, headers=hdr)
    return urllib.request.urlopen(req, timeout=1800).read()


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
    return refs, col_seq, col_beg


def decompress_bgzf(data):
    res, pos = [], 0
    while pos < len(data):
        d = zlib.decompressobj(47)
        res.append(d.decompress(data[pos:]))
        unused = d.unused_data
        if not unused:
            break
        pos = len(data) - len(unused)
    return b"".join(res)


for qtd, tissue in [("QTD000251", "atrial"), ("QTD000256", "ventricle")]:
    gz = BASE + qtd + "/" + qtd + ".all.tsv.gz"
    refs, col_seq, col_beg = parse_tbi(os.path.join(W, qtd + ".all.tsv.gz.tbi"))
    bins = refs["1"]
    # wide region around both genes: 11,780,000-12,020,000
    lo, hi = 11780000, 12020000
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
    genes_seen = {}
    nppb_rows = []
    for s, e in merged:
        data = fetch(gz, s, e)
        text = decompress_bgzf(data).decode("utf-8", "replace")
        for line in text.split("\n"):
            if not line or line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 12 or p[1] != "1":
                continue
            try:
                pos = int(p[2])
            except ValueError:
                continue
            if pos < lo or pos > hi:
                continue
            genes_seen[p[0]] = genes_seen.get(p[0], 0) + 1
            if p[0] == "ENSG00000183705":
                nppb_rows.append(p)
        del text, data
    out.append("%s (%s): genes in region: %d" % (tissue, qtd, len(genes_seen)))
    out.append("  gene list: %s" % ", ".join(sorted(genes_seen)))
    out.append("  NPPB rows: %d" % len(nppb_rows))
    if nppb_rows:
        out.append("  first NPPB row: %s" % "\t".join(nppb_rows[0]))

open(os.path.join(W, "_nppb_diag.txt"), "w", encoding="utf-8").write(
    "\n".join(out))
print("ok")
