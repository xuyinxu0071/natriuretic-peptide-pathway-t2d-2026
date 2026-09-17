# -*- coding: utf-8 -*-
import re
BASE = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\Merged_Submission"
en = open(BASE + r"\Manuscript_EN.md", encoding="utf-8").read()
m = re.search(r"## Abstract\n(.*?)\n---", en, re.S)
abst = re.sub(r"\*\*(Background|Methods|Results|Conclusions)\.\*\*", "", m.group(1))
ws = len(abst.split())
alnum = len(re.findall(r"[A-Za-z0-9\-']+", abst))
open(BASE + r"\..\coloc_analysis\step9h_output.txt", "w", encoding="utf-8").write(
    f"whitespace-split: {ws}\nalnum-tokens: {alnum}\n")
