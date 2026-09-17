# -*- coding: utf-8 -*-
"""
Step 2f: 裁决双build坐标 + 选定变体列表数据源
================================================================================
产出: step2f_output.txt (打印) — 坐标 + 1000G/dbSNP 目录清单 + 候选文件大小
"""
import json, re, urllib.request, urllib.error

OUT = open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis\step2f_output.txt", "w", encoding="utf-8")

def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.write(s + "\n"); OUT.flush()

def get(url, accept="application/json", timeout=90):
    req = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": "coloc-pipeline/1.0"})
    return urllib.request.urlopen(req, timeout=timeout)

# ---------- 1. Ensembl 主REST: 两个SNP的双build坐标 ----------
for rsid in ("rs5068", "rs1421811"):
    try:
        j = json.load(get(f"https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json"))
        P(f"[{rsid}] mappings:")
        for m in j.get("mappings", []):
            P(f"   {m.get('assembly_name')}: {m.get('seq_region_name')}:{m.get('start')}-{m.get('end')} allele={m.get('allele_string')}")
    except Exception as e:
        P(f"[{rsid}] FAILED:", repr(e)[:150])

# ---------- 2. 1000G release 目录完整清单 ----------
try:
    html = get("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/", accept="text/html").read().decode("utf-8", "ignore")
    files = sorted(set(re.findall(r'href="([^"]+)"', html)))
    files = [f for f in files if not f.startswith(("../", "?"))]
    P("\n[1000G 20130502] 全部文件:")
    for f in files:
        P("   ", f)
except Exception as e:
    P("[1000G listing] FAILED:", repr(e)[:150])

# ---------- 3. NCBI dbSNP GRCh37 目录 ----------
for path in (
    "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_grch37/p13/VCF/",
    "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_grch37/p13/VCF/common_notested/",
):
    try:
        html = get(path, accept="text/html").read().decode("utf-8", "ignore")
        files = sorted(set(re.findall(r'href="([^"]+\.vcf\.gz(?:\.tbi)?)"', html)))
        P(f"\n[dbSNP {path.split('/snp/')[-1]}] vcf.gz 文件 ({len(files)}):")
        for f in files[:60]:
            P("   ", f)
        if len(files) > 60:
            P(f"    ... 共 {len(files)} 个")
    except Exception as e:
        P(f"[dbSNP {path}] FAILED:", repr(e)[:150])

# ---------- 4. 候选文件 HEAD 大小 ----------
candidates = [
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz",
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz",
]
P("\n[HEAD 大小]")
for u in candidates:
    try:
        r = get(u, accept="application/octet-stream", timeout=60)
        # urllib GET 会开始下载 — 用 Range 只取1字节测可达性与 Content-Range
        pass
    except Exception:
        pass
# 改用 HEAD 方式
import urllib.request as ur
for u in candidates:
    try:
        req = ur.Request(u, method="HEAD")
        r = ur.urlopen(req, timeout=60)
        P(f"   {u.split('/')[-1]}: {int(r.headers['Content-Length'])/1e6:.1f} MB")
    except Exception as e:
        P(f"   {u.split('/')[-1]}: HEAD FAILED {repr(e)[:100]}")

OUT.close()
