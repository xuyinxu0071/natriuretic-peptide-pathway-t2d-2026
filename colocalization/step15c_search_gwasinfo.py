# -*- coding: utf-8 -*-
"""step15c: 用 GET /gwasinfo 全库元数据检索 NP 加工酶级联 pQTL + T2D 结局数据集"""
import json, os, urllib.request

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
TOKEN = open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\.opengwas_token.txt").read().strip()
API = "https://api.opengwas.io/api"

CACHE = os.path.join(W, "gwasinfo_cache.json")
if os.path.exists(CACHE):
    info = json.load(open(CACHE, encoding="utf-8"))
else:
    req = urllib.request.Request(API + "/gwasinfo", headers={
        "Authorization": "Bearer " + TOKEN, "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=300) as r:
        info = json.loads(r.read().decode("utf-8"))
    json.dump(info, open(CACHE, "w", encoding="utf-8"))

out = ["total datasets in gwasinfo: %d" % len(info)]

# --- pQTL 检索：靶蛋白名 + pQTL 证据特征 ---
# SCALLOP: ebi-a-GCST90012* (Olink, trait 如 "NT-proBNP levels")
# INTERVAL: prot-a-* (SomaScan)
# UKB-PPP: ebi-a-GCST9024-9027* (Olink)
# deCODE: ebi-a-GCST9019-9020* (SomaScan)
TARGETS = {
    "CORIN":   ["corin", "corin levels", "atosin"],
    "FURIN":   ["furin", "pace4", "pcsks3"],
    "MME":     ["neprilysin", "mme", "enkephalinase", "cd10", "neprilysin levels"],
    "DPP4":    ["dpp4", "dipeptidyl peptidase", "cd26"],
    "NPR1":    ["npr1", "natriuretic peptide receptor 1", "guanylate cyclase a"],
    "NPR2":    ["npr2", "natriuretic peptide receptor 2"],
    "NPR3":    ["npr3", "natriuretic peptide receptor 3", "npr-c"],
    "NPPA":    ["natriuretic peptide a", "anp", "nppa", "nt-proanp", "proanp"],
    "NPPB":    ["natriuretic peptide b", "bnp", "nppb", "nt-probnp", "probnp"],
}

def is_pqtl_like(rid, rec):
    if rid.startswith("prot-a-"):
        return True
    tr = (rec.get("trait") or "").lower()
    if rid.startswith("ebi-a-GCST90") and any(
            k in tr for k in ["levels", "measurement", "protein"]):
        return True
    return False

hits = {}
for rid, rec in info.items():
    tr = (rec.get("trait") or "").lower()
    note = (rec.get("note") or "").lower()
    blob = tr + " " + note
    for gene, kws in TARGETS.items():
        for kw in kws:
            if kw in blob:
                if is_pqtl_like(rid, rec) or any(
                        a in blob for a in ["olink", "somalogic", "somascan",
                                            "plasma protein", "pqtl"]):
                    hits.setdefault(gene, []).append(
                        (rid, rec.get("trait", ""), rec.get("sample_size", ""),
                         rec.get("author", ""), rec.get("year", "")))
                break

out.append("")
out.append("=== pQTL datasets per NP-cascade target ===")
for gene in TARGETS:
    hh = hits.get(gene, [])
    # de-dup by (id)
    hh = sorted(set(hh), key=lambda x: -int(str(x[2]) or 0) if str(x[2]).isdigit() else 0)
    out.append("%s: %d" % (gene, len(hh)))
    for h in hh[:15]:
        out.append("   %s | n=%s | %s | %s %s" % (h[0], h[2], h[1][:70], h[3], h[4]))

# --- T2D / 血糖结局检索 ---
out.append("")
out.append("=== T2D outcome candidates ===")
t2d = []
for rid, rec in info.items():
    tr = (rec.get("trait") or "").lower()
    if ("type 2 diabetes" in tr or "type ii diabetes" in tr
            or tr.strip() in ["t2d"]) and (rec.get("ncase") or 0) > 2000:
        t2d.append((rid, rec.get("trait", ""), rec.get("ncase", ""),
                    rec.get("ncontrol", ""), rec.get("author", ""),
                    rec.get("year", ""), rec.get("population", "")))
t2d = sorted(set(t2d), key=lambda x: -int(str(x[2]) or 0) if str(x[2]).isdigit() else 0)
for h in t2d[:15]:
    out.append("   %s | cases=%s ctrl=%s | %s | %s %s | %s" %
               (h[0], h[2], h[3], h[1][:55], h[4], h[5], h[6]))

open(os.path.join(W, "step15c_output.txt"), "w", encoding="utf-8").write("\n".join(out))
print("done")
