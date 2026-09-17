# -*- coding: utf-8 -*-
"""Step 1b: Ensembl overlap rsID-SNV 密度测量（100kb 样区）"""
import json, urllib.request

def ensembl_get(path):
    req = urllib.request.Request("https://rest.ensembl.org" + path,
                                 headers={"Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=180))

if __name__ == "__main__":
    for name, region in [("NPPA簇北侧100kb", "1:11740000-11840000"),
                         ("NPPA簇本体120kb", "1:11840000-11960000")]:
        r = ensembl_get(f"/overlap/region/human/{region}?feature=variation")
        snv_rs = [v for v in r if str(v.get("id", "")).startswith("rs")
                  and len(v.get("alleles", [])) == 2
                  and all(len(a) == 1 for a in v["alleles"])]
        print(f"{name} {region}: total {len(r)}, rsID-双等位SNV {len(snv_rs)}")
