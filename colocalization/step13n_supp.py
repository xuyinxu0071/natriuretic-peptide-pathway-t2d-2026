# -*- coding: utf-8 -*-
"""尝试多种 URL/头下载 Zhu 2016 补充材料 PDF"""
import urllib.request, os

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
OUT = []
cands = [
    ("supp1", "https://static-content.springer.com/esm/art%3A10.1038%2Fng.3538/MediaObjects/41588_2016_BFng3538_MOESM1_ESM.pdf",
     {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36", "Referer": "https://www.nature.com/"}),
    ("supp2", "https://www.nature.com/articles/ng.3538.pdf", {}),
]
hdr0 = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/120.0 Safari/537.36",
        "Accept": "application/pdf,*/*"}
for tag, url, extra in cands:
    h = dict(hdr0); h.update(extra)
    try:
        req = urllib.request.Request(url, headers=h)
        data = urllib.request.urlopen(req, timeout=300).read()
        dest = os.path.join(WORK, f"zhu2016_{tag}.pdf" if data[:4] == b"%PDF" else f"zhu2016_{tag}.html")
        open(dest, "wb").write(data)
        OUT.append(f"{tag}: {len(data):,} bytes, magic={data[:8]}")
    except Exception as e:
        OUT.append(f"{tag}: {e}")
open(os.path.join(WORK, "step13n_output.txt"), "w", encoding="utf-8").write("\n".join(OUT))
print("DONE")
