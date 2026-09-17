# -*- coding: utf-8 -*-
"""读 plink 日志并清洗输出"""
import re
p = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis\step13c_plink_test.txt"
raw = open(p, "rb").read()
text = raw.decode("utf-8", "ignore")
text = "".join(ch if (ch == "\n" or ch == "\r" or 32 <= ord(ch) < 127 or ord(ch) > 127) else " " for ch in text)
out = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis\step13c_plink_clean.txt"
open(out, "w", encoding="utf-8").write(text[-2000:])
print("ok")
