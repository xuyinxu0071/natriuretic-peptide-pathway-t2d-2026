# -*- coding: utf-8 -*-
"""抓 SMR 页面 HTML, 提取全部 .zip/.gz 下载链接"""
import re, urllib.request

OUT = []
for url in ("https://yanglab.westlake.edu.cn/software/smr/",
            "https://yanglab.westlake.edu.cn/software/smr/#Download"):
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        html = urllib.request.urlopen(req, timeout=120).read().decode("utf-8", "ignore")
        links = sorted(set(re.findall(r'href="([^"]+\.(?:zip|tar\.gz|gz))"', html)))
        OUT.append(f"{url}: {len(links)} links")
        for l in links:
            OUT.append("  " + l)
    except Exception as e:
        OUT.append(f"{url}: ERR {e}")

open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis\step13h_output.txt",
     "w", encoding="utf-8").write("\n".join(OUT))
print("DONE")
