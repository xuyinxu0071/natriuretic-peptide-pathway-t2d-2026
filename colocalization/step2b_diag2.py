# -*- coding: utf-8 -*-
"""诊断 tbi 解析：打印区间偏移与解压首行"""
import gzip, os, struct, urllib.request, zlib
from step2_tabix import parse_tbi, reg2bins, merge_ranges, decompress_bgzf, BASE, VCF, WORK

refs = parse_tbi()
print("names:", sorted(refs.keys()))

chrom = "1"
lo, hi = 11700000, 12100000
bins = refs[chrom]
want = set(reg2bins(lo - 1, hi))
chunks = []
for bin_no, chs in bins:
    if bin_no in want:
        chunks.extend(chs)
print(f"chr{chrom} bins matched: {sum(1 for b,_ in bins if b in want)}, chunks: {len(chunks)}")
ranges = [(cb >> 16, (ce >> 16) + 70000) for cb, ce in chunks]
merged = merge_ranges(ranges)
print("merged ranges:", merged)

s, e = merged[0]
req = urllib.request.Request(BASE + VCF, headers={"Range": f"bytes={s}-{e}"})
data = urllib.request.urlopen(req, timeout=300).read()
text = decompress_bgzf(data).decode("utf-8", "ignore")
lines = [l for l in text.split("\n") if l and not l.startswith("#")]
print(f"range 0: got {len(data)} bytes -> {len(lines)} VCF lines")
for l in lines[:3]:
    p = l.split("\t")
    print(f"  chr{p[0]} pos={p[1]} id={p[2]} {p[3]}/{p[4]}")
# 尾部行
for l in lines[-2:]:
    p = l.split("\t")
    print(f"  ... chr{p[0]} pos={p[1]} id={p[2]}")
