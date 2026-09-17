# -*- coding: utf-8 -*-
"""Step 6c: 对比 panel 样本名与 VCF 样本名格式"""
import os, re, struct, gzip, urllib.request, zlib

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"

p = os.path.join(WORK, "kg_panel.txt")
lines = open(p).read().splitlines()
print("panel total lines:", len(lines))
print("panel head 3:", [l[:60] for l in lines[:3]])
eur = [l.split("\t")[0] for l in lines if len(l.split("\t")) >= 3 and l.split("\t")[2] == "EUR"]
print("EUR count:", len(eur), "first 5:", eur[:5])

# VCF header: 需要读文件头部获取 #CHROM 行
# 但 header 在文件开头 (offset 0), 简单 fetch 0-200KB 并找 #CHROM 行
def fetch(url, start=None, end=None, timeout=600):
    hdr = {"User-Agent": "ld-diag/1.0"}
    if start is not None:
        hdr["Range"] = f"bytes={start}-{end}"
    req = urllib.request.Request(url, headers=hdr)
    return urllib.request.urlopen(req, timeout=timeout).read()

data = fetch(BASE + "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", 0, 300000)
d = zlib.decompressobj(47)
text = (d.decompress(data) + d.flush()).decode("ascii", "ignore")
for line in text.split("\n"):
    if line.startswith("#CHROM"):
        samples = line.split("\t")[9:]
        print("VCF n_samples:", len(samples))
        print("VCF first 5:", samples[:5])
        print("overlap with EUR set:", len(set(samples) & set(eur)))
        break
else:
    print("no #CHROM line in first 300KB; first lines:", text.split(chr(10))[:3][:1])
