# -*- coding: utf-8 -*-
"""Step 6a: 诊断基因型 VCF 区域字节块"""
import gzip, os, struct, urllib.request, zlib

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

raw = gzip.open(TBI_LOCAL, "rb").read()
off = 4
n_ref, fmt, col_seq, col_beg, col_end, meta, skip = struct.unpack_from("<iiiiiii", raw, off); off += 28
print("tbi fmt/cols:", fmt, col_seq, col_beg, col_end, meta, skip)
l_nm = struct.unpack_from("<i", raw, off)[0]; off += 4
names = raw[off:off + l_nm].decode().rstrip("\0").split("\0"); off += l_nm
print("refs:", names[:5], "n_ref:", n_ref)
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
print("want bins:", len(wbins))
chunks = []
for bin_no, chs in refs["1"]:
    if bin_no in wbins:
        chunks.extend(chs)
print("chunks:", len(chunks))
for cb, ce in chunks[:6]:
    print(f"  chunk: file {cb>>16}..{ce>>16} (+{cb&0xffff} into block) span={((ce>>16)-(cb>>16))/1e3:.0f}KB")

# 取第一个 chunk, 从块边界开始拉一小段
cb, ce = sorted(chunks)[0]
s = cb >> 16
data = fetch(BASE + VCF, s, s + 200000)
print("fetched", len(data), "bytes from offset", s)
d = zlib.decompressobj(47)
text = (d.decompress(data) + d.flush()).decode("ascii", "ignore")
print("decompressed first member:", len(text))
lines = [l for l in text.split("\n") if l.strip()]
print("lines:", len(lines))
for l in lines[:3]:
    print("LINE[:250]:", l[:250])
# 位置分布
import collections
pos = []
for l in lines:
    p = l.split("\t")
    if len(p) > 2 and p[0] == "1":
        try:
            pos.append(int(p[1]))
        except ValueError:
            pass
print("chr1 rows:", len(pos))
if pos:
    print("pos range:", min(pos), "-", max(pos))
