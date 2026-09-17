# -*- coding: utf-8 -*-
"""下载 SMR.cpp 与 SMR.hpp 并抽取默认参数"""
import os, re, urllib.request

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
SRC = os.path.join(WORK, "smr_src")
OUT = []
for f in ("src/SMR.cpp", "src/SMR.hpp"):
    dest = os.path.join(SRC, f.replace("/", "_"))
    if os.path.exists(dest):
        OUT.append(f"{f}: 已存在"); continue
    url = f"https://raw.githubusercontent.com/JianYang-Lab/SMR/master/{f}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        data = urllib.request.urlopen(req, timeout=300).read()
        open(dest, "wb").write(data)
        OUT.append(f"{f}: {len(data):,} bytes")
    except Exception as e:
        OUT.append(f"{f}: {e}")

# 抽默认值
txt = open(os.path.join(SRC, "src_SMR.cpp"), "r", encoding="utf-8", errors="ignore").read()
pats = ["peqtl-heidi", "peqtl_heidi", "heidi-min-m", "heidi_max_m", "heidi-max-m",
        "ld-lower-limit", "ld_lower_limit", "ld-upper-limit", "ld_upper_limit",
        "ld_r2_top", "ldr2_top", "MAX_NUM_LD", "heidi-mtd", "diff-freq", "peqtl-smr",
        "rm_cor_sbat"]
for p in pats:
    for m in re.finditer(re.escape(p), txt):
        s = max(0, m.start() - 120); e = min(len(txt), m.end() + 120)
        OUT.append(f"--- {p} @{m.start()}: ...{txt[s:e].replace(chr(10), ' | ')}...")
        break  # 每个模式只取首处

# rm_cor_sbat 在 SMR_data.cpp 中
txt2 = open(os.path.join(SRC, "src_SMR_data.cpp"), "r", encoding="utf-8", errors="ignore").read()
m = re.search(r"void rm_cor_sbat[\s\S]{0,1500}", txt2)
if m:
    OUT.append("=== rm_cor_sbat ===\n" + m.group(0)[:1400])

open(os.path.join(WORK, "step13p_output.txt"), "w", encoding="utf-8").write("\n".join(OUT))
print("DONE")
