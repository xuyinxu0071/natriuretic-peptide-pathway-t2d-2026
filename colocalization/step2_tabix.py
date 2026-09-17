# -*- coding: utf-8 -*-
"""
Step 2: mini-tabix 客户端——从 1000G WGS sites VCF 按 HTTP Range 提取区域变体（EUR MAF>=0.01）
=================================================================================
原理: .tbi 索引(bgzf压缩)给出每个 bin 的 chunk 虚拟偏移(高32位=压缩文件偏移,低32位=块内偏移)。
对目标区域求 reg2bins,合并 chunk 字节区间,Range 请求,解压 bgzf(级联 gzip 成员),过滤行。
产出: nppa_locus_variants.json / npr3_locus_variants.json
"""
import gzip, io, json, os, re, struct, urllib.request, zlib

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
VCF = "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz"
TBI = VCF + ".tbi"

REGIONS = {
    "1": (11700000, 12100000),   # NPPA/NPPB 宽窗(GRCh37)
    "5": (32500000, 33000000),   # NPR3 宽窗(GRCh37)
}

def fetch_range(url, start, end, retries=4):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
            return urllib.request.urlopen(req, timeout=300).read()
        except Exception as e:
            if i == retries - 1:
                raise
            import time; time.sleep(5 * (i + 1))

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
    """返回 {chrom_name: [(bin, [(voff_beg, voff_end), ...]), ...] 的字典}"""
    tbi_path = os.path.join(WORK, "sites.tbi")
    if not os.path.exists(tbi_path) or os.path.getsize(tbi_path) < 1e5:
        print("downloading tbi ...")
        data = urllib.request.urlopen(BASE + TBI, timeout=300).read()
        open(tbi_path, "wb").write(data)
    raw = gzip.open(tbi_path, "rb").read()
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
        off += 8 * n_intv  # 线性索引跳过
        refs[names[r]] = bins
    return refs

def merge_ranges(ranges, pad=200000):
    ranges = sorted(ranges)
    merged = []
    for s, e in ranges:
        if merged and s <= merged[-1][1] + pad:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged

def decompress_bgzf(data):
    """跨多个 bgzf(gzip) 成员循环解压"""
    out, pos = [], 0
    while pos < len(data):
        d = zlib.decompressobj(47)
        out.append(d.decompress(data[pos:]))
        unused = d.unused_data
        if not unused:
            break
        pos = len(data) - len(unused)
    return b"".join(out)

if __name__ == "__main__":
    refs = parse_tbi()
    print("tbi references:", len(refs))

    for chrom, (lo, hi) in REGIONS.items():
        bins = refs.get(chrom)
        if bins is None:
            print(f"chr{chrom} not in tbi!"); continue
        want = set(reg2bins(lo - 1, hi))
        chunks = []
        for bin_no, chs in bins:
            if bin_no in want:
                chunks.extend(chs)
        if not chunks:
            print(f"chr{chrom}: no chunks"); continue
        # 虚拟偏移 → 压缩字节区间
        ranges = [(cb >> 16, (ce >> 16) + 70000) for cb, ce in chunks]
        merged = merge_ranges(ranges)
        total = sum(e - s for s, e in merged)
        print(f"chr{chrom} {lo}-{hi}: {len(chunks)} chunks -> {len(merged)} byte-ranges, ~{total/1e3:.0f} KB")

        variants = []
        for s, e in merged:
            data = fetch_range(BASE + VCF, s, e)
            text = decompress_bgzf(data).decode("utf-8", "ignore")
            for line in text.split("\n"):
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if parts[0] != chrom:
                    continue
                pos = int(parts[1])
                if pos < lo or pos > hi:
                    continue
                rsid = parts[2]
                if rsid == ".":
                    continue
                ref, alt = parts[3], parts[4]
                if "," in alt or len(ref) != 1 or len(alt) != 1:
                    continue
                m = re.search(r"EUR_AF=([0-9.eE+-]+)", parts[7])
                if not m:
                    continue
                eaf = float(m.group(1))
                maf = min(eaf, 1.0 - eaf)
                if maf < 0.01:
                    continue
                variants.append(dict(rsid=rsid, chrom=chrom, pos=pos,
                                     ref=ref, alt=alt, eur_af=round(eaf, 4)))
        variants.sort(key=lambda v: v["pos"])
        name = "nppa" if chrom == "1" else "npr3"
        with open(os.path.join(WORK, f"{name}_locus_variants.json"), "w") as f:
            json.dump(variants, f)
        print(f"chr{chrom}: {len(variants)} EUR-common biallelic SNVs saved -> {name}_locus_variants.json")
        # 定位关键 SNP 真实 GRCh37 位置
        for target in ("rs5068", "rs1421811"):
            hit = [v for v in variants if v["rsid"] == target]
            if hit:
                print(f"  {target} at true GRCh37 chr{chrom}:{hit[0]['pos']}")
