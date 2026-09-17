# -*- coding: utf-8 -*-
"""step13a: 下载 SMR 与 PLINK Windows 二进制"""
import os, zipfile, urllib.request

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
BIN = os.path.join(WORK, "bin")
os.makedirs(BIN, exist_ok=True)
OUT = []
def p(s): OUT.append(str(s))

URLS = {
    "smr.zip": "https://yanglab.westlake.edu.cn/software/smr/download/smr-1.3.1-win-x86_64.zip",
    "plink.zip": "http://s3.amazonaws.com/plink1-assets/plink_win64_20210606.zip",
}
for name, url in URLS.items():
    dest = os.path.join(BIN, name)
    if os.path.exists(dest):
        p(f"{name}: 已存在"); continue
    p(f"下载 {url} ...")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        data = urllib.request.urlopen(req, timeout=600).read()
        open(dest, "wb").write(data)
        p(f"  保存 {len(data):,} bytes")
        with zipfile.ZipFile(dest) as z:
            z.extractall(BIN)
            p("  解压: " + ", ".join(z.namelist()[:20]))
    except Exception as e:
        p(f"  失败: {e}")

with open(os.path.join(WORK, "step13a_output.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(OUT))
print("DONE")
