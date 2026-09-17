# -*- coding: utf-8 -*-
import urllib.request
import os
import re

W = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
out = []


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace"), r.geturl()


for u in ["https://www.ebi.ac.uk/eqtl/data-access",
          "https://www.ebi.ac.uk/eqtl/api-info",
          "https://www.ebi.ac.uk/eqtl/release-notes"]:
    try:
        s, final = get(u)
        txt = re.sub(r"<[^>]+>", " ", s)
        txt = re.sub(r"\s+", " ", txt)
        out.append("OK %s (final=%s) len=%d" % (u, final, len(s)))
        out.append("text: " + txt[:3000])
        for m in re.finditer(r'href="([^"]+)"', s):
            h = m.group(1)
            if any(k in h.lower() for k in ["api", "download", "ftp",
                                             "sumstat", "rest", "json",
                                             "bucket", "s3", "opensearch"]):
                out.append("  link: " + h)
        out.append("=" * 80)
    except Exception as e:
        out.append("FAIL %s -> %s" % (u, e))

open(os.path.join(W, "_eqtlcat_probe4.txt"), "w", encoding="utf-8").write(
    "\n".join(out))
print("ok")
