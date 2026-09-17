# -*- coding: utf-8 -*-
"""终验: EN 摘要词数 + 关键新数值在 EN/ZH 的对称性 + 陈旧表述扫描"""
import re
BASE = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\Merged_Submission"
OUT = open(BASE + r"\..\coloc_analysis\step9g_output.txt", "w", encoding="utf-8")
def P(*a):
    s = " ".join(str(x) for x in a); OUT.write(s + "\n"); OUT.flush(); print(s, flush=True)

en = open(BASE + r"\Manuscript_EN.md", encoding="utf-8").read()
zh = open(BASE + r"\Manuscript_ZH.md", encoding="utf-8").read()

# 1. 摘要词数 (Background..Conclusions, 去掉小节标签)
m = re.search(r"## Abstract\n(.*?)\n---", en, re.S)
abst = re.sub(r"\*\*(Background|Methods|Results|Conclusions)\.\*\*", "", m.group(1))
words = len(re.findall(r"[A-Za-z0-9\-']+", abst))
P(f"[1] EN abstract words: {words} (limit 350)")

# 2. 新数值对称性
checks = [
    ("PP4 ≤ 0.06 蛋白层默认", "PP4 ≤ 0.06", "PP4 ≤ 0.06"),
    ("PP1 0.78–0.95", "PP1 0.78–0.95", "PP1 0.78–0.95"),
    ("PP3 ≤ 0.22", "PP3 ≤ 0.22", "PP3 ≤ 0.22"),
    ("宽松 0.34–0.39", "0.34–0.39", "0.34–0.39"),
    ("协调变异 1,866–1,997", "1,866–1,997", "1,866–1,997"),
    ("SCALLOP z +4.30", "+4.30", "+4.30"),
    ("SomaLogic z +2.97", "+2.97", "+2.97"),
    ("ANP z +2.19", "+2.19", "+2.19"),
    ("CAD z −2.46", "−2.46", "−2.46"),
    ("MAF ≈ 0.01 / ≈1%", "MAF ≈ 0.01", "MAF ≈ 1%"),
    ("ANP PP0 0.62–0.85", "PP0 0.62–0.85", "PP0 0.62–0.85"),
    ("32 对 (表6)", "32 exposure–outcome pairs", "32 个暴露–结局对"),
    ("p12 = 1×10⁻⁴ 宽松", "p12 = 1×10⁻⁴", "p12 = 1×10⁻⁴"),
    ("三层表述", "expression, protein, and blood-pressure", "表达、蛋白与血压"),
    ("rs198389 pQTL lead", "rs198389", "rs198389"),
    ("Ishigaki/Low 引文", "Ishigaki K, Akiyama M", "Ishigaki K, Akiyama M"),
]
for label, a, b in checks:
    P(f"[{'OK' if (a in en and b in zh) else 'FAIL'}] {label}: EN={'Y' if a in en else 'N'} ZH={'Y' if b in zh else 'N'}")

# 3. 陈旧表述扫描 (不应再出现)
stale = ["no public pQTL", "formally untested", "not yet possible with public pQTL",
         "有待合适数据", "无公开 pQTL 资源提供 NPPA/NPPB", "尚不能经共定位裁决",
         "p12 = 1×10⁻⁵）敏感性", "permissive-prior (p12 = 1×10⁻⁵)", "permissive prior (p12 = 1×10⁻⁵)"]
for s in stale:
    hit = []
    if s in en: hit.append("EN")
    if s in zh: hit.append("ZH")
    for f in ["Supplement_EN.md", "Supplement_ZH.md"]:
        if s in open(BASE + "\\" + f, encoding="utf-8").read(): hit.append(f)
    P(f"[{'CLEAN' if not hit else 'STALE->'+','.join(hit)}] {s}")
OUT.close()
