# -*- coding: utf-8 -*-
"""step14d: fetch header of eQTL Catalogue sumstats files to confirm columns."""
import os
import struct
import urllib.request
import zlib

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/"

out = []
for qtd in ["QTD000251", "QTD000256"]:
    url = BASE + qtd + "/" + qtd + ".all.tsv.gz"
    req = urllib.request.Request(url, headers={
        "User-Agent": "coloc-pipeline/1.0", "Range": "bytes=0-262143"})
    data = urllib.request.urlopen(req, timeout=600).read()
    d = zlib.decompressobj(47)
    text = (d.decompress(data) + d.flush()).decode("utf-8", "replace")
    lines = text.split("\n")
    out.append("%s: first lines:" % qtd)
    for l in lines[:4]:
        out.append("  " + l[:400])
open(os.path.join(W, "_gtex_header.txt"), "w", encoding="utf-8").write(
    "\n".join(out))
print("ok")
