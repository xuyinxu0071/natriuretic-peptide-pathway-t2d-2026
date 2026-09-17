# -*- coding: utf-8 -*-
"""Step 0c: rs5068/rs1421811 双 build 坐标核实 + overlap 记录字段结构"""
import json, urllib.request

def ensembl_get(path):
    req = urllib.request.Request("https://rest.ensembl.org" + path,
                                 headers={"Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=120))

if __name__ == "__main__":
    v = ensembl_get("/variation/human/rs5068?populations=1000GENOMES:phase3:EUR")
    m = v["mappings"][0]
    print(f"rs5068: GRCh38 {m['seq_region_name']}:{m['start']}-{m['end']} allele={v.get('allele_string')}")
    for mm in v["mappings"]:
        print(f"  mapping: {mm['assembly']} {mm['seq_region_name']}:{mm['start']}")
    eur = v.get("populations", [])
    for p in eur:
        if "EUR" in str(p.get("population", "")):
            print(f"  EUR AF: {p.get('allele_frequency')} (allele {p.get('allele')})")
    print(f"  MAF: {v.get('MAF')} minor={v.get('minor_allele')}")

    v2 = ensembl_get("/variation/human/rs1421811?populations=1000GENOMES:phase3:EUR")
    m2 = v2["mappings"][0]
    print(f"rs1421811: GRCh38 {m2['seq_region_name']}:{m2['start']}-{m2['end']} allele={v2.get('allele_string')}")
    print(f"  MAF: {v2.get('MAF')} minor={v2.get('minor_allele')}")

    # overlap 记录字段结构
    r = ensembl_get("/overlap/region/human/1:11857000-11859000?feature=variation")
    rec = next(x for x in r if x.get("id", "").startswith("rs"))
    print("\noverlap record keys:", sorted(rec.keys()))
    print(json.dumps(rec, indent=1)[:800])
