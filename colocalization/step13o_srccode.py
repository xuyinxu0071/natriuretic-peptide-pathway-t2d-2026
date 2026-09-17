# -*- coding: utf-8 -*-
"""下载 SMR 官方源码中实现 HEIDI 的核心文件"""
import os, urllib.request

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
SRC = os.path.join(WORK, "smr_src")
os.makedirs(SRC, exist_ok=True)
OUT = []
files = [
    "src/SMR_data.cpp", "src/SMR_data.hpp", "src/StatFunc.cpp", "src/StatFunc.hpp",
]
for f in files:
    url = f"https://raw.githubusercontent.com/JianYang-Lab/SMR/master/{f}"
    dest = os.path.join(SRC, f.replace("/", "_"))
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        data = urllib.request.urlopen(req, timeout=300).read()
        open(dest, "wb").write(data)
        OUT.append(f"{f}: {len(data):,} bytes")
    except Exception as e:
        OUT.append(f"{f}: {e}")
open(os.path.join(WORK, "step13o_output.txt"), "w", encoding="utf-8").write("\n".join(OUT))
print("DONE")
