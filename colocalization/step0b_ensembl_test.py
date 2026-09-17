# -*- coding: utf-8 -*-
"""Step 0b: 测试 Ensembl GRCh37 REST 区域变体列表端点"""
import json, urllib.request

def ensembl_get(path):
    req = urllib.request.Request("https://rest.ensembl.org" + path,
                                 headers={"Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=120))

if __name__ == "__main__":
    # 小区域测试: GRCh38 上 NPPA/NPPB 簇附近 20kb（rsID 匹配 OpenGWAS，build 无关）
    r = ensembl_get("/overlap/region/human/1:11857000-11877000?feature=variation")
    print(f"Total features: {len(r)}")
    snv_maf = [v for v in r if v.get("class") == "SNV" and isinstance(v.get("MAF"), (int, float))]
    print(f"SNV with MAF: {len(snv_maf)}")
    common = [v for v in snv_maf if v["MAF"] >= 0.01]
    print(f"SNV MAF>=0.01: {len(common)}")
    for v in common[:8]:
        print(f"  {v['id']} {v['seq_region_name']}:{v['start']} {v.get('allele_string')} MAF={v.get('MAF')}")
    # 基因坐标确认（主 REST = GRCh38）
    for sym in ("NPPA", "NPPB", "NPR3"):
        g = ensembl_get(f"/lookup/symbol/homo_sapiens/{sym}?expand=0")
        print(f"{sym}: chr{g['seq_region_name']}:{g['start']}-{g['end']} strand={g['strand']} (GRCh38)")
