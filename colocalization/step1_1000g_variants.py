# -*- coding: utf-8 -*-
"""Step 1: 1000G phase3 sites VCF 下载 + 区域变体列表提取（EUR MAF>=0.01）"""
import gzip, io, json, os, re, urllib.request

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"

FILES = {
    "1": "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz",
    "5": "ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz",
}

def download(chrnum):
    dest = os.path.join(WORK, f"kg_chr{chrnum}_sites.vcf.gz")
    if os.path.exists(dest) and os.path.getsize(dest) > 1e6:
        print(f"chr{chrnum}: already downloaded ({os.path.getsize(dest)/1e6:.1f} MB)")
        return dest
    url = BASE + FILES[chrnum]
    print(f"chr{chrnum}: downloading {url} ...")
    data = urllib.request.urlopen(url, timeout=600).read()
    with open(dest, "wb") as f:
        f.write(data)
    print(f"chr{chrnum}: {len(data)/1e6:.1f} MB saved")
    return dest

# 区域定义（GRCh37 = 1000G phase3 build）
# 先用宽窗口提取，rs5068/rs1421811 真实位置从文件内定位后再精确定义 coloc 区
REGIONS = {
    "1": (11700000, 12200000),   # NPPA/NPPB 簇宽窗（定位 rs5068）
    "5": (32400000, 33300000),   # NPR3 宽窗（定位 rs1421811）
}

if __name__ == "__main__":
    os.makedirs(WORK, exist_ok=True)
    for chrnum in FILES:
        path = download(chrnum)
        lo, hi = REGIONS[chrnum]
        out = []
        rs5068_pos = rs1421811_pos = None
        with gzip.open(path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                pos = int(parts[1])
                if pos < lo:
                    continue
                if pos > hi:
                    break
                rsid = parts[2]
                if rsid == "rs5068":
                    rs5068_pos = pos
                if rsid == "rs1421811":
                    rs1421811_pos = pos
                if rsid == ".":
                    continue
                ref, alt = parts[3], parts[4]
                if "," in alt or len(ref) != 1 or len(alt) != 1:
                    continue  # 仅双等位 SNV
                info = parts[7]
                m = re.search(r"EUR_AF=([0-9.eE+-]+)", info)
                if not m:
                    continue
                eaf = float(m.group(1))
                maf = min(eaf, 1 - eaf)
                if maf < 0.01:
                    continue
                out.append(dict(rsid=rsid, chrom=chrnum, pos=pos, ref=ref, alt=alt,
                                eur_af=eaf, eur_maf=maf))
        with open(os.path.join(WORK, f"variants_chr{chrnum}_wide.json"), "w") as f:
            json.dump(out, f)
        print(f"chr{chrnum}: {len(out)} common biallelic SNVs in {lo}-{hi}")
        if rs5068_pos: print(f"  rs5068 at chr1:{rs5068_pos} (true GRCh37)")
        if rs1421811_pos: print(f"  rs1421811 at chr5:{rs1421811_pos} (true GRCh37)")
