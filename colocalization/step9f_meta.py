# -*- coding: utf-8 -*-
"""取新数据集元信息 (作者/年份/PMID) 供引文使用"""
import json
WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
idx = json.load(open(WORK + r"\gwasinfo_full.json", encoding="utf-8"))
OUT = open(WORK + r"\step9f_output.txt", "w", encoding="utf-8")
for i in ["ebi-a-GCST90012082", "prot-a-2078", "prot-a-2076",
          "bbj-a-159", "bbj-a-109", "bbj-a-129", "bbj-a-71",
          "ebi-a-GCST90085762"]:
    e = idx.get(i)
    if e:
        OUT.write(f"{i}\n  trait={e.get('trait')}\n  author={e.get('author')}  year={e.get('year')}  pmid={e.get('pmid')}\n  n={e.get('sample_size')} ncase={e.get('ncase')} ncontrol={e.get('ncontrol')}\n  pop={e.get('population')}  group={e.get('group_name')}\n\n")
OUT.close()
