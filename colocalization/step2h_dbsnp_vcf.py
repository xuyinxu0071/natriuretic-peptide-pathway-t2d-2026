# -*- coding: utf-8 -*-
"""Step 2h: 定位 b151 GRCh37 的 common VCF"""
import re, urllib.request

OUT = open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis\step2h_output.txt", "w", encoding="utf-8")

def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

def get(url, accept="application/json", timeout=90):
    req = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": "coloc-pipeline/1.0"})
    return urllib.request.urlopen(req, timeout=timeout)

def listdir(url, show=50):
    try:
        html = get(url, accept="text/html").read().decode("utf-8", "ignore")
        items = [i for i in re.findall(r'href="([^"]+)"', html) if not i.startswith(("../", "http", "?"))]
        P(f"\n[{url}] ({len(items)} 项)")
        for i in items[:show]:
            P("   ", i)
        if len(items) > show:
            P(f"    ... 共 {len(items)}")
        return items
    except Exception as e:
        P(f"\n[{url}] FAILED: {repr(e)[:130]}")
        return []

base = "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/"
items = listdir(base)
for sub in [i for i in items if i.endswith("/") and "vcf" in i.lower()]:
    subitems = listdir(base + sub)
    # 探查 VCF 子目录里的 common / per-chr 文件
    for f in subitems:
        if "common" in f.lower() and (f.endswith(".gz") or f.endswith(".gz.tbi")):
            # HEAD 大小
            try:
                req = urllib.request.Request(base + sub + f, method="HEAD")
                r = urllib.request.urlopen(req, timeout=60)
                P(f"   SIZE {f}: {int(r.headers['Content-Length'])/1e6:.1f} MB")
            except Exception as e:
                P(f"   SIZE {f}: FAILED {repr(e)[:80]}")

OUT.close()
