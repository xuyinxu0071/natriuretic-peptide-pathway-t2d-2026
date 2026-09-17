# -*- coding: utf-8 -*-
"""列出 1000G release/20130502 目录内容"""
import re, urllib.request

url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
html = urllib.request.urlopen(url, timeout=120).read().decode()
names = sorted(set(re.findall(r'href="([^"]+)"', html)))
for n in names:
    if "chr1." in n or "chr5." in n or "sites" in n.lower():
        print(n)
print("---total entries:", len(names))
# 也检查 supported_variants 等子目录
sub = [n for n in names if n.endswith("/")]
print("subdirs:", sub[:30])
