# -*- coding: utf-8 -*-
"""Step 0d: overlap 记录结构 + MAF 分箱计数（评估 API 拉取规模）"""
import json, urllib.request

def ensembl_get(path):
    req = urllib.request.Request("https://rest.ensembl.org" + path,
                                 headers={"Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=120))

if __name__ == "__main__":
    v = ensembl_get("/variation/human/rs1421811")
    print("rs1421811 mappings:", [(m["seq_region_name"], m["start"], m.get("coord_system","")) for m in v["mappings"]])
    print("rs1421811 MAF:", v.get("MAF"), "minor:", v.get("minor_allele"))

    r = ensembl_get("/overlap/region/human/1:11857000-11859000?feature=variation")
    rec = next(x for x in r if str(x.get("id", "")).startswith("rs"))
    print("\noverlap record keys:", sorted(rec.keys()))
    print(json.dumps(rec, indent=1)[:900])

    rs_snv = [x for x in r if str(x.get("id", "")).startswith("rs")
              and x.get("var_class") == "SNV"
              and isinstance(x.get("MAF"), (int, float))]
    print(f"\n20kb 区域: 全部 {len(r)} | rs+SNV+MAF {len(rs_snv)} | MAF>=0.01: {sum(1 for x in rs_snv if x['MAF']>=0.01)} | MAF>=0.05: {sum(1 for x in rs_snv if x['MAF']>=0.05)}")
