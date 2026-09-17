# -*- coding: utf-8 -*-
"""诊断 2: 解压文本内容分析"""
import urllib.request
from step2_tabix import decompress_bgzf, BASE, VCF

s, e = 6330257, 6786166
req = urllib.request.Request(BASE + VCF, headers={"Range": f"bytes={s}-{e}"})
data = urllib.request.urlopen(req, timeout=300).read()
print("raw bytes:", len(data), "head:", data[:20].hex())
text = decompress_bgzf(data)
print("decompressed bytes:", len(text))
lines = text.split(b"\n")
print("total lines:", len(lines))
n_hash = sum(1 for l in lines if l.startswith(b"#"))
n_empty = sum(1 for l in lines if not l.strip())
n_short = sum(1 for l in lines if 0 < len(l.split(b"\t")) < 5)
print(f"# lines: {n_hash}, empty: {n_empty}, short(<5 fields): {n_short}")
data_lines = [l for l in lines if l and not l.startswith(b"#")]
print("data lines:", len(data_lines))
if data_lines:
    print("first 3 reprs:", [l[:80] for l in data_lines[:3]])
    print("last 2 reprs:", [l[:80] for l in data_lines[-2:]])
    import collections
    chroms = collections.Counter(l.split(b"\t")[0].decode() for l in data_lines if l.split(b"\t")[0])
    print("chrom distribution:", dict(list(chroms.items())[:5]))
