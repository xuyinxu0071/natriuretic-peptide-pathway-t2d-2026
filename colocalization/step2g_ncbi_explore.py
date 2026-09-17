# -*- coding: utf-8 -*-
"""Step 2g: 探查 NCBI dbSNP 目录结构，找 GRCh37 common VCF"""
import re, urllib.request, urllib.error

OUT = open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis\step2g_output.txt", "w", encoding="utf-8")

def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.write(s + "\n"); OUT.flush()

def get(url, accept="application/json", timeout=90):
    req = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": "coloc-pipeline/1.0"})
    return urllib.request.urlopen(req, timeout=timeout)

def listdir(url, show=40):
    try:
        html = get(url, accept="text/html").read().decode("utf-8", "ignore")
        items = re.findall(r'href="([^"]+)"', html)
        items = [i for i in items if not i.startswith(("../", "http", "?"))]
        P(f"\n[{url}]  ({len(items)} 项)")
        for i in items[:show]:
            P("   ", i)
        if len(items) > show:
            P(f"    ... 共 {len(items)}")
        return items
    except Exception as e:
        P(f"\n[{url}] FAILED: {repr(e)[:130]}")
        return []

items = listdir("https://ftp.ncbi.nlm.nih.gov/snp/")
# 找 organisms 或 human 相关
for cand in [i for i in items if i.rstrip("/").endswith(("organisms", "archive", "Homo_sapiens", "human")) or "organism" in i.lower()]:
    pass

# 常见候选路径直接探
cands = [
    "https://ftp.ncbi.nlm.nih.gov/snp/organisms/",
    "https://ftp.ncbi.nlm.nih.gov/snp/archive/",
    "https://ftp.ncbi.nlm.nih.gov/snp/human/",
    "https://ftp.ncbi.nlm.nih.gov/snp/Homo_sapiens/",
    "https://ftp.ncbi.nlm.nih.gov/snp/latest_variance/",
]
for c in cands:
    listdir(c)

OUT.close()
