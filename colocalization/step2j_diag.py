# -*- coding: utf-8 -*-
"""Step 2j: 诊断 dbSNP common VCF 区域字节块解压后的实际内容"""
import gzip, json, os, re, struct, time, urllib.request, zlib

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
BASE = "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/"
VCF = "common_all_20180423.vcf.gz"
TBI_LOCAL = os.path.join(WORK, "dbsnp_common.tbi")

def fetch_range(url, start, end, retries=5):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
            return urllib.request.urlopen(req, timeout=600).read()
        except Exception:
            if i == retries - 1:
                raise
            time.sleep(5 * (i + 1))

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

chrom, lo, hi = "1", 11405974, 12405974
want = set(reg2bins(lo - 1, hi))
chunks = []
for bin_no, chs in refs[chrom]:
    if bin_no in want:
        chunks.extend(chs)
print("chunks:", len(chunks))
# 取最小起始的一个 chunk, 只拉一小段看内容
chunks.sort()
cb, ce = chunks[0]
print(f"first chunk voff: {cb:#x} -> file {cb>>16} + {cb&0xffff}")
data = fetch_range(BASE + VCF, cb >> 16, (cb >> 16) + 70000)
print("fetched", len(data), "bytes")
try:
    text = decompress_bgzf(data).decode("utf-8", "ignore")
    print("decompressed", len(text), "chars")
    lines = [l for l in text.split("\n") if l.strip()]
    print("non-empty lines:", len(lines))
    for l in lines[:5]:
        print("LINE:", l[:200])
except Exception as e:
    print("decompress error:", e)
    print("raw head:", data[:30])
