# -*- coding: utf-8 -*-
"""Step 2e: 裁决 rs5068/rs1421811 双 build 坐标 + 检查 dbSNP FTP"""
import json, re, urllib.request

def get(url, timeout=120):
    req = urllib.request.Request(url, headers={"User-Agent": "python-urllib"})
    return urllib.request.urlopen(req, timeout=timeout).read()

# 1) NCBI refsnp API (双 build 坐标)
for rsid in ("5068", "1421811"):
    try:
        d = json.loads(get(f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{rsid}"))
        placements = d.get("primary_snapshot_data", {}).get("placements_with_allele", [])
        for p in placements:
            asm = p.get("placement_annot", {}).get("seq_id_traits_by_assembly", [])
            name = p.get("seq_id")
            pos = p.get("alleles", [{}])[0].get("allele", {}).get("spdi", {})
            for a in asm:
                an = a.get("assembly_name", "")
                if an in ("GRCh37.p13", "GRCh38.p13", "GRCh38", "GRCh37"):
                    print(f"rs{rsid}: {an} chr{name}:{pos.get('position', '?')} {pos.get('deleted_sequence')}>{pos.get('inserted_sequence')}")
    except Exception as e:
        print(f"rs{rsid} refsnp failed: {e}")

# 2) dbSNP GRCh37 VCF 目录
try:
    html = get("https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_grch37/VCF/").decode("utf-8", "ignore")
    names = sorted(set(re.findall(r'href="([^"]+\.vcf\.gz)"', html)))
    print("\ndbSNP GRCh37 VCF files:")
    for n in names[:20]:
        print(" ", n)
except Exception as e:
    print("dbSNP FTP failed:", e)
