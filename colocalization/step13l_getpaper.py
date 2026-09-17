# -*- coding: utf-8 -*-
"""下载 Zhu 2016 SMR 论文正文与补充材料 PDF"""
import urllib.request, os

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
OUT = []
urls = {
    "zhu2016_smr_main.pdf": "https://www.nature.com/articles/ng.3538.pdf",
    "zhu2016_smr_supp.pdf": "https://static-content.springer.com/esm/art%3A10.1038%2Fng.3538/MediaObjects/41588_2016_BFng3538_MOESM1_ESM.pdf",
}
for name, url in urls.items():
    dest = os.path.join(WORK, name)
    if os.path.exists(dest) and os.path.getsize(dest) > 10000:
        OUT.append(f"{name}: 已存在"); continue
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        data = urllib.request.urlopen(req, timeout=300).read()
        open(dest, "wb").write(data)
        OUT.append(f"{name}: {len(data):,} bytes")
    except Exception as e:
        OUT.append(f"{name}: 失败 {e}")
open(os.path.join(WORK, "step13l_output.txt"), "w", encoding="utf-8").write("\n".join(OUT))
print("DONE")
