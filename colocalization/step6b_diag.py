# -*- coding: utf-8 -*-
"""Step 6b: 遍历 chr1 区域全部 chunk, 统计行覆盖/rsID 命中"""
import gzip, json, os, struct, time, urllib.request, zlib, glob

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
VCF = "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
TBI_LOCAL = os.path.join(WORK, VCF + ".tbi")

def fetch(url, start=None, end=None, timeout=1200):
    hdr = {"User-Agent": "coloc-pipeline/1.0"}
    if start is not None:
        hdr["Range"] = f"bytes={start}-{end}"
    req = urllib.request.Request(url, headers=hdr)
    return urllib.request.urlopen(req, timeout=timeout).read()

def reg2bins(beg, end):
    lst = [0]
    end -= 1
    for k in range(1 + (beg >> 26), 1 + (end >> 26) + 1): lst.append(k)
    for k in range(9 + (beg >> 23), 9 + (end >> 23) + 1): lst.append(k)
    for k in range(73 + (beg >> 20), 73 + (end >> 20) + 1): lst.append(k)
    for k in range(585 + (beg >> 17), 585 + (end >> 17) + 1): lst.append(k)
    for k in range(4681 + (beg >> 14), 4681 + (end >> 14) + 1): lst.append(k)
    return lst

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

raw = gzip.open(TBI_LOCAL, "rb").read()
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

lo, hi = 11805974, 12005974
wbins = set(reg2bins(lo - 1, hi))
chunks = []
for bin_no, chs in refs["1"]:
    if bin_no in wbins:
        chunks.extend(chs)
print("chunks:", len(chunks))

# want_ids
want = set()
for f in glob.glob(os.path.join(WORK, "nppa__*_assoc.json")):
    want.update(json.load(open(f)).keys())
print("want_ids:", len(want))

ranges = sorted(set((cb >> 16, (ce >> 16) + 70000) for cb, ce in chunks))
# 合并
merged = []
for s, e in ranges:
    if merged and s <= merged[-1][1] + 100000:
        merged[-1][1] = max(merged[-1][1], e)
    else:
        merged.append([s, e])
print("merged ranges:", merged, "total MB:", sum(e-s for s,e in merged)/1e6)

stats = dict(rows=0, in_region=0, with_rsid=0, in_want=0)
examples = []
for s, e in merged:
    data = fetch(BASE + VCF, s, e)
    text = decompress_bgzf(data).decode("ascii", "ignore")
    for line in text.split("\n"):
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 10 or parts[0] != "1":
            continue
        stats["rows"] += 1
        try:
            pos = int(parts[1])
        except ValueError:
            continue
        if lo <= pos <= hi:
            stats["in_region"] += 1
            rsid = parts[2]
            if rsid != ".":
                stats["with_rsid"] += 1
                if rsid in want:
                    stats["in_want"] += 1
                    if len(examples) < 5:
                        examples.append((rsid, pos, parts[3], parts[4], parts[7][:60]))
print(stats)
for ex in examples:
    print("EX:", ex)
