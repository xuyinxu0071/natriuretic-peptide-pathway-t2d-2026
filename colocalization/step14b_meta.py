# -*- coding: utf-8 -*-
"""step14b: locate GTEx heart-tissue eQTL datasets in eQTL Catalogue
metadata (r7) and verify HTTPS accessibility of sumstats files."""
import urllib.request
import os
import csv
import io

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
META_URL = ("https://raw.githubusercontent.com/eQTL-Catalogue/"
            "eQTL-Catalogue-resources/master/data_tables/dataset_metadata_r7.tsv")


def fetch(url, binary=False):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        data = r.read()
    return data if binary else data.decode("utf-8", "replace")


out = []
meta_text = fetch(META_URL)
rows = list(csv.DictReader(io.StringIO(meta_text), delimiter="\t"))
out.append("metadata rows: %d; columns: %s" % (len(rows), list(rows[0].keys())))

heart = [r for r in rows
         if "heart" in r.get("tissue_label", "").lower()
         or "atrial" in r.get("tissue_label", "").lower()
         or "ventricle" in r.get("tissue_label", "").lower()]
out.append("heart-tissue datasets: %d" % len(heart))
for r in heart:
    out.append("  study=%s qtl_group=%s QTD=%s QTS=%s tissue=%s n=%s method=%s" % (
        r.get("study_label"), r.get("sample_group"), r.get("dataset_id"),
        r.get("study_id"), r.get("tissue_label"), r.get("sample_size"),
        r.get("quant_method")))

gtex = [r for r in rows if "gtex" in r.get("study_label", "").lower()]
out.append("GTEx datasets: %d" % len(gtex))
seen = {}
for r in gtex:
    seen.setdefault(r.get("sample_group"), r)
out.append("GTEx sample_groups (%d): %s" % (
    len(seen), ", ".join(sorted(seen.keys()))))
for qg in sorted(seen):
    if ("heart" in qg.lower() or "atrial" in qg.lower()
            or "ventri" in qg.lower()):
        r = seen[qg]
        out.append("  HEART group: %s -> QTD=%s QTS=%s n=%s method=%s" % (
            qg, r.get("dataset_id"), r.get("study_id"),
            r.get("sample_size"), r.get("quant_method")))

open(os.path.join(W, "_eqtlcat_meta.txt"), "w", encoding="utf-8").write(
    "\n".join(out))
print("ok")
