# -*- coding: utf-8 -*-
"""诊断 4: 完整解压后按位置分带统计"""
import urllib.request
from step2_tabix import decompress_bgzf, BASE, VCF, parse_tbi, reg2bins, merge_ranges

refs = parse_tbi()
chrom, lo, hi = "1", 11700000, 12100000
bins = refs[chrom]
want = set(reg2bins(lo - 1, hi))
chunks = []
for bin_no, chs in bins:
    if bin_no in want:
        chunks.extend(chs)
ranges = [(cb >> 16, (ce >> 16) + 70000) for cb, ce in chunks]
merged = merge_ranges(ranges)
print("ranges:", merged)

for s, e in merged:
    req = urllib.request.Request(BASE + VCF, headers={"Range": f"bytes={s}-{e}"})
    data = urllib.request.urlopen(req, timeout=300).read()
    text = decompress_bgzf(data)
    print(f"range {s}-{e}: raw {len(data)} -> decompressed {len(text)}")
    lines = [l for l in text.split(b"\n") if l and not l.startswith(b"#")]
    print(f"  data lines: {len(lines)}")
    if lines:
        def pos(l):
            try: return int(l.split(b"\t")[1])
            except Exception: return None
        ps = [p for p in (pos(l) for l in lines) if p]
        print(f"  pos range: {min(ps)} - {max(ps)}")
        inwin = [l for l in lines if pos(l) and lo <= pos(l) <= hi]
        print(f"  in window {lo}-{hi}: {len(inwin)}")
        # 检查 INFO 字段样例
        if inwin:
            sample = inwin[len(inwin)//2].decode("utf-8","ignore")
            print("  sample line:", sample[:200])
            import re
            print("  has EUR_AF:", "EUR_AF=" in sample)
            withrs = [l for l in inwin if l.split(b"\t")[2] != b"."]
            print(f"  in-window with rsID: {len(withrs)}")
