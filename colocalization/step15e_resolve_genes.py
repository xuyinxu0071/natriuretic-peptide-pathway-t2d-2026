# -*- coding: utf-8 -*-
"""step15e: 解析 NP 级联基因 ENSG ID，确认 eqtl-a- 数据集存在，并探测 GTEx 心脏组织可用性"""
import json, os, urllib.request

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
info = json.load(open(os.path.join(W, "gwasinfo_cache.json"), encoding="utf-8"))
out = []

GENES = ["CORIN", "FURIN", "MME", "DPP4", "NPR1", "NPR2", "NPR3", "NPPA", "NPPB",
         "KLKB1", "MMP2"]  # KLKB1/MMP2 备选(ANP 降解相关)

# 1) Ensembl REST 解析 symbol -> ENSG + chr:pos
def ensembl_lookup(symbol):
    url = ("https://rest.ensembl.org/lookup/symbol/homo_sapiens/%s"
           "?expand=0" % symbol)
    req = urllib.request.Request(url, headers={
        "Content-Type": "application/json", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))

ensg = {}
for g in GENES:
    try:
        rec = ensembl_lookup(g)
        ensg[g] = rec
        out.append("%s -> %s chr%s:%d-%d" % (g, rec.get("id"),
                  rec.get("seq_region_name"), rec.get("start"), rec.get("end")))
    except Exception as e:
        out.append("%s -> FAIL %s" % (g, e))

json.dump(ensg, open(os.path.join(W, "np_cascade_genes.json"), "w"))

# 2) 检查 eqtl-a-ENSG... 在 gwasinfo 中是否存在
out.append("")
out.append("=== eqtl-a- availability (eQTLGen whole blood, via OpenGWAS) ===")
for g, rec in ensg.items():
    eid = "eqtl-a-" + rec.get("id", "")
    if eid in info:
        m = info[eid]
        out.append("%s: %s OK | trait=%s | n=%s | nsnp=%s" %
                   (g, eid, m.get("trait"), m.get("sample_size"), m.get("nsnp")))
    else:
        out.append("%s: %s MISSING" % (g, eid))

# 3) 顺带确认结局数据集元数据（AF/CAD/HF/T2D）
out.append("")
out.append("=== outcome datasets meta ===")
for oid in ["ebi-a-GCST006061", "ieu-a-7", "ebi-a-GCST005195",
            "ebi-a-GCST009541", "ebi-a-GCST007515", "finn-b-T2D",
            "ebi-a-GCST90018926", "prot-a-1150", "prot-a-2076"]:
    m = info.get(oid)
    if m:
        out.append("%s | %s | n=%s cases=%s ctrl=%s | %s %s" %
                   (oid, (m.get("trait") or "")[:60], m.get("sample_size"),
                    m.get("ncase"), m.get("ncontrol"),
                    m.get("author"), m.get("year")))
    else:
        out.append("%s | MISSING" % oid)

open(os.path.join(W, "step15e_output.txt"), "w", encoding="utf-8").write("\n".join(out))
print("done")
