# -*- coding: utf-8 -*-
"""诊断 Range 请求行为与解压结果"""
import urllib.request, zlib

BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
VCF = "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz"

req = urllib.request.Request(BASE + VCF, headers={"Range": "bytes=0-9999"})
resp = urllib.request.urlopen(req, timeout=120)
print("status:", resp.status)
print("Content-Range:", resp.headers.get("Content-Range"))
print("Content-Length:", resp.headers.get("Content-Length"))
data = resp.read()
print("bytes got:", len(data))
try:
    text = zlib.decompressobj(47).decompress(data).decode("utf-8", "ignore")
    print("decompressed chars:", len(text))
    print("first 500:", text[:500])
except Exception as e:
    print("decompress err:", e)
    print("raw head:", data[:50])
