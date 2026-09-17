# -*- coding: utf-8 -*-
"""
Step 2i: dbSNP b151 GRCh37 common VCF mini-tabix 提取两个 locus 的常见双等位 SNV
=================================================================================
数据源: https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz
       (全基因组常见变体 AF>=1%, 1.6GB, tbi 2.3MB — 仅下载索引+区域字节块)
区域(GRCh37, ±500kb):
  NPPA/NPPB: chr1:11,905,974 (rs5068) -> 11,405,974-12,405,974
  NPR3:      chr5:32,714,270 (rs1421811) -> 32,214,270-33,214,270
产出: nppa_locus_variants.json / npr3_locus_variants.json
     [ {rsid, chrom, pos, ref, alt, af}, ... ] 按 pos 排序
"""
import gzip, json, os, re, struct, time, urllib.request, zlib

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
BASE = "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/"
VCF = "common_all_20180423.vcf.gz"
TBI = VCF + ".tbi"
TBI_LOCAL = os.path.join(WORK, "dbsnp_common.tbi")

REGIONS = {
    "1": (11405974, 12405974),   # NPPA/NPPB ±500kb (rs5068 GRCh37=11,905,974)
    "5": (32214270, 33214270),   # NPR3 ±500kb (rs1421811 GRCh37=32,714,270)
}
LOCUS_NAME = {"1": "nppa", "5": "npr3"}

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

def parse_tbi():
    if not os.path.exists(TBI_LOCAL) or os.path.getsize(TBI_LOCAL) < 1e5:
        print("downloading dbSNP tbi (2.3MB) ...", flush=True)
        data = urllib.request.urlopen(urllib.request.Request(BASE + TBI), timeout=600).read()
        open(TBI_LOCAL, "wb").write(data)
    raw = gzip.open(TBI_LOCAL, "rb").read()
    off = 0
    magic = raw[:4]; off = 4
    assert magic == b"TBI\x01", magic
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip = struct.unpack_from("<iiiiiii", raw, off); off += 28
    l_nm = struct.unpack_from("<i", raw, off)[0]; off += 4
    names = raw[off:off + l_nm].decode().rstrip("\0").split("\0")
    off += l_nm
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
    return refs, names

def merge_ranges(ranges, pad=100000):
    ranges = sorted(ranges)
    merged = []
    for s, e in ranges:
        if merged and s <= merged[-1][1] + pad:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged

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

SNV_RE = re.compile(r"^[ACGT]$")

if __name__ == "__main__":
    refs, names = parse_tbi()
    print("tbi references:", len(refs), "| first 5:", names[:5], flush=True)

    for chrom, (lo, hi) in REGIONS.items():
        bins = refs.get(chrom)
        if bins is None:
            print(f"chr{chrom} not in tbi! available sample: {names[:10]}"); continue
        want = set(reg2bins(lo - 1, hi))
        chunks = []
        for bin_no, chs in bins:
            if bin_no in want:
                chunks.extend(chs)
        if not chunks:
            print(f"chr{chrom}: no chunks"); continue
        ranges = [(cb >> 16, (ce >> 16) + 70000) for cb, ce in chunks]
        merged = merge_ranges(ranges)
        total = sum(e - s for s, e in merged)
        print(f"chr{chrom} {lo}-{hi}: {len(chunks)} chunks -> {len(merged)} ranges, ~{total/1e3:.0f} KB", flush=True)

        variants, seen = [], set()
        for s, e in merged:
            data = fetch_range(BASE + VCF, s, e)
            text = decompress_bgzf(data).decode("utf-8", "ignore")
            for line in text.split("\n"):
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 8 or parts[0] != chrom:
                    continue
                try:
                    pos = int(parts[1])
                except ValueError:
                    continue
                if pos < lo or pos > hi:
                    continue
                rsid = parts[2]
                if not rsid.startswith("rs"):
                    continue
                ref, alt = parts[3].upper(), parts[4].upper()
                if "," in alt or alt == "." or not SNV_RE.match(ref) or not SNV_RE.match(alt):
                    continue
                m = re.search(r"(?:^|;)CAF=([0-9.,]+)(?:;|$)", parts[7])
                if not m:
                    continue
                vals = m.group(1).split(",")
                if len(vals) != 2 or vals[1] in (".", ""):
                    continue
                try:
                    af = float(vals[1])
                except ValueError:
                    continue
                if af > 1.0:
                    af /= 100.0
                maf = min(af, 1.0 - af)
                if maf < 0.01:
                    continue
                key = rsid
                if key in seen:
                    continue
                seen.add(key)
                variants.append(dict(rsid=rsid, chrom=chrom, pos=pos,
                                     ref=ref, alt=alt, af=round(af, 4)))
        variants.sort(key=lambda v: v["pos"])
        name = LOCUS_NAME[chrom]
        with open(os.path.join(WORK, f"{name}_locus_variants.json"), "w") as f:
            json.dump(variants, f)
        print(f"chr{chrom}: {len(variants)} common biallelic SNVs saved -> {name}_locus_variants.json", flush=True)
        for target in ("rs5068", "rs1421811"):
            hit = [v for v in variants if v["rsid"] == target]
            if hit:
                print(f"  {target} found: chr{chrom}:{hit[0]['pos']} {hit[0]['ref']}/{hit[0]['alt']} AF={hit[0]['af']}", flush=True)
    print("DONE", flush=True)
