# -*- coding: utf-8 -*-
"""step15d: 精确词边界检索 NP 级联 pQTL（覆盖 INTERVAL/UKB-PPP/SCALLOP/deCODE）"""
import json, os, re

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
info = json.load(open(os.path.join(W, "gwasinfo_cache.json"), encoding="utf-8"))
out = ["total: %d" % len(info)]

# 精确 trait 匹配模式（词边界）
PATTERNS = {
    "CORIN_exact":   r"\bcorin\b",
    "NEPRILYSIN":    r"neprilysin",
    "MME_exact":     r"\bmme\b",
    "DPP4":          r"dipeptidyl peptidase[- ]?(iv|4)|\bdpp4\b|\bcd26\b",
    "NPR1":          r"natriuretic peptide receptor[- ]?(1|a|one)|\bnpr1\b|guanylate cyclase a\b",
    "NPR2":          r"natriuretic peptide receptor[- ]?(2|b|two)|\bnpr2\b",
    "NPR3":          r"natriuretic peptide receptor[- ]?(3|c|three)|\bnpr3\b|clearance receptor",
    "NPPA_ANP":      r"atrial natriuretic|\bnppa\b|pro-?anp|nt-?proanp",
    "NPPB_BNP":      r"brain natriuretic|b-?type natriuretic|\bnppb\b|pro-?bnp|nt-?probnp|pro-?brain",
    "SCALLOP_any":   r"scallop",
    "UKBPPP_probe":  r"sun bb|elliott p|razq",  # 探针: INTERVAL 作者签名
}

for name, pat in PATTERNS.items():
    rx = re.compile(pat, re.I)
    hh = []
    for rid, rec in info.items():
        tr = rec.get("trait") or ""
        if rx.search(tr):
            hh.append((rid, tr, rec.get("sample_size", ""), rec.get("author", ""),
                       rec.get("year", ""), rec.get("population", "")))
    out.append("")
    out.append("=== %s (%d) ===" % (name, len(hh)))
    # 按样本量排序取前 20
    def n_of(x):
        try:
            return -int(x[2])
        except Exception:
            return 0
    for h in sorted(set(hh), key=n_of)[:20]:
        out.append("   %s | n=%s | %s | %s %s | %s" % (h[0], h[2], h[1][:80], h[3], h[4], h[5]))

open(os.path.join(W, "step15d_output.txt"), "w", encoding="utf-8").write("\n".join(out))
print("done")
