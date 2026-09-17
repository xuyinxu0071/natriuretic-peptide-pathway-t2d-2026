# -*- coding: utf-8 -*-
"""下载 SMR 1.4.2 Windows 版"""
import os, zipfile, urllib.request

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
BIN = os.path.join(WORK, "bin")
OUT = []
for ver in ("1.4.2", "1.4.1", "1.3.2"):
    url = f"https://yanglab.westlake.edu.cn/software/smr/download/smr-{ver}-win-x86_64.zip"
    dest = os.path.join(BIN, f"smr-{ver}-win-x86_64.zip")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        data = urllib.request.urlopen(req, timeout=600).read()
        open(dest, "wb").write(data)
        OUT.append(f"{ver}: {len(data):,} bytes 下载成功")
        with zipfile.ZipFile(dest) as z:
            z.extractall(BIN)
            OUT.append(f"  解压: {', '.join(n for n in z.namelist() if not n.startswith('__'))[:300]}")
        break
    except Exception as e:
        OUT.append(f"{ver}: 失败 {e}")
open(os.path.join(WORK, "step13g_output.txt"), "w", encoding="utf-8").write("\n".join(OUT))
print("DONE")
