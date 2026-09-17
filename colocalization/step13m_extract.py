# -*- coding: utf-8 -*-
"""从 nature HTML 提取 SMR/HEIDI 方法的纯文本"""
import re, html as ihtml

SRC = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis\zhu2016_smr_main.pdf"
raw = open(SRC, "r", encoding="utf-8", errors="ignore").read()
# 去 script/style
raw = re.sub(r"<(script|style)[\s\S]*?</\1>", " ", raw)
# 块级标签换行
raw = re.sub(r"</(p|div|h\d|li|section)>", "\n", raw)
text = re.sub(r"<[^>]+>", " ", raw)
text = ihtml.unescape(text)
text = re.sub(r"[ \t]+", " ", text)
lines = [l.strip() for l in text.split("\n") if l.strip()]

OUT = []
# 找 Methods 段与含关键词的行
grab = False
for i, l in enumerate(lines):
    if re.search(r"HEIDI|heterogeneity in dependent|SMR test|b\s*xy|top associated SNP|pleiotropy", l, re.I):
        OUT.append(f"[{i}] {l}")

open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis\zhu2016_methods.txt",
     "w", encoding="utf-8").write("\n\n".join(OUT))
# 也存全文方法段
full = "\n".join(lines)
idx = full.find("Methods")
open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis\zhu2016_fulltext.txt",
     "w", encoding="utf-8").write(full[max(0,idx-200): idx+20000])
print("DONE")
