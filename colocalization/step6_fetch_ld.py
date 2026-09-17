# -*- coding: utf-8 -*-
"""
Step 6 (v2): 1000G v5b 基因型 → EUR LD 矩阵 (±100kb, 按 pos+ref+alt 匹配)
================================================================================
v5b 文件无 rsID -> 用 (chrom,pos,ref,alt) 与 dbSNP b151 变体列表精确匹配。
INFO 的 EUR_AF 与 dbSNP CAF 做频率一致性 QC (|diff|<0.25)。
产出: ld_nppa.csv / ld_npr3.csv + ld_variants_*.json
"""
import gzip, json, math, os, struct, time, urllib.request, zlib, glob, re

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
PANEL = "integrated_call_samples_v3.20130502.ALL.panel"
VCF = {
    "nppa": "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    "npr3": "ALL.chr5.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
}
REGIONS = {
    "nppa": ("1", 11805974, 12005974),
    "npr3": ("5", 32614270, 32814270),
}

def fetch(url, start=None, end=None, retries=5, timeout=1200):
    for i in range(retries):
        try:
            hdr = {"User-Agent": "coloc-pipeline/1.0"}
            if start is not None:
                hdr["Range"] = f"bytes={start}-{end}"
            req = urllib.request.Request(url, headers=hdr)
            return urllib.request.urlopen(req, timeout=timeout).read()
        except Exception:
            if i == retries - 1:
                raise
            time.sleep(8 * (i + 1))

def reg2bins(beg, end):
    lst = [0]
    end -= 1
    for k in range(1 + (beg >> 26), 1 + (end >> 26) + 1): lst.append(k)
    for k in range(9 + (beg >> 23), 9 + (end >> 23) + 1): lst.append(k)
    for k in range(73 + (beg >> 20), 73 + (end >> 20) + 1): lst.append(k)
    for k in range(585 + (beg >> 17), 585 + (end >> 17) + 1): lst.append(k)
    for k in range(4681 + (beg >> 14), 4681 + (end >> 14) + 1): lst.append(k)
    return lst

def parse_tbi(path):
    raw = gzip.open(path, "rb").read()
    off = 4
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip = struct.unpack_from("<iiiiiii", raw, off); off += 28
    l_nm = struct.unpack_from("<i", raw, off)[0]; off += 4
    names = raw[off:off + l_nm].decode().rstrip("\0").split("\0"); off += l_nm
    refs = {}
    for r in range(n_ref):
        n_bin = struct.unpack_from("<i", raw, off)[0]; off += 4
        bins = []
        for b in range(n_bin):
            bin_no, n_chunk = struct.unpack_from("<Ii", raw, off); off += 8
            chunks = []
            for c in range(n_chunk):
                cb, ce = struct.unpack_from("<QQ", raw, off); off += 16
                chunks.append((cb, ce))
            bins.append((bin_no, chunks))
        n_intv = struct.unpack_from("<i", raw, off)[0]; off += 4
        off += 8 * n_intv
        refs[names[r]] = bins
    return refs

def decompress_bgzf(data):
    out, pos = [], 0
    while pos < len(data):
        d = zlib.decompressobj(47)
        out.append(d.decompress(data[pos:]))
        unused = d.unused_data
        if not unused:
            break
        pos = len(data) - len(unused)
    return b"".join(out)

def get_eur_samples():
    p = os.path.join(WORK, "kg_panel.txt")
    if not os.path.exists(p):
        open(p, "wb").write(fetch(BASE + PANEL))
    return [l.split("\t")[0] for l in open(p).read().splitlines()
            if l and not l.startswith("sample") and len(l.split("\t")) >= 3 and l.split("\t")[2] == "EUR"]

def get_eur_indices(vcf_url, eur_set):
    """从 VCF #CHROM 头行解析样本列序 -> EUR 列索引"""
    data = fetch(vcf_url, 0, 400000)
    d = zlib.decompressobj(47)
    text = (d.decompress(data) + d.flush()).decode("ascii", "ignore")
    for line in text.split("\n"):
        if line.startswith("#CHROM"):
            samples = line.split("\t")[9:]
            idx = [i for i, s in enumerate(samples) if s in eur_set]
            return len(samples), idx
    raise RuntimeError("no #CHROM header in first 400KB")

EUR_AF_RE = re.compile(r"EUR_AF=([0-9.eE+-]+)")

if __name__ == "__main__":
    log = open(os.path.join(WORK, "step6_log.txt"), "w", encoding="utf-8")
    def P(s):
        print(s, flush=True); log.write(s + "\n"); log.flush()

    eur = get_eur_samples()
    eur_set = set(eur)
    P(f"EUR samples: {len(eur)}")

    for locus, (chrom, lo, hi) in REGIONS.items():
        variants = json.load(open(os.path.join(WORK, f"{locus}_locus_variants.json")))
        # 窄窗内的变体, key=(pos,ref,alt)
        keymap = {}
        for v in variants:
            if lo <= v["pos"] <= hi:
                keymap[(v["pos"], v["ref"], v["alt"])] = v
        P(f"\n=== {locus}: chr{chrom}:{lo}-{hi}, {len(keymap)} target variants (narrow) ===")

        tbi = os.path.join(WORK, VCF[locus] + ".tbi")
        if not os.path.exists(tbi) or os.path.getsize(tbi) < 1e5:
            P("downloading tbi ...")
            open(tbi, "wb").write(fetch(BASE + VCF[locus] + ".tbi"))
        refs = parse_tbi(tbi)
        bins = refs[chrom]
        wbins = set(reg2bins(lo - 1, hi))
        chunks = []
        for bin_no, chs in bins:
            if bin_no in wbins:
                chunks.extend(chs)
        ranges = sorted(set((cb >> 16, (ce >> 16) + 70000) for cb, ce in chunks))
        merged = []
        for s, e in ranges:
            if merged and s <= merged[-1][1] + 100000:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        total = sum(e - s for s, e in merged)
        P(f"{len(chunks)} chunks -> {len(merged)} ranges, ~{total/1e6:.1f} MB")

        got = {}       # rsid -> [dosages]
        meta = {}      # rsid -> (pos, ref, alt)
        n_samples, eur_idx = get_eur_indices(BASE + VCF[locus], eur_set)
        P(f"  VCF samples: {n_samples}, EUR matched: {len(eur_idx)}")
        if len(eur_idx) < 400:
            P("  EUR match FAILED"); continue
        for s, e in merged:
            data = fetch(BASE + VCF[locus], s, e)
            text = decompress_bgzf(data).decode("ascii", "ignore")
            for line in text.split("\n"):
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 10 or parts[0] != chrom:
                    continue
                try:
                    pos = int(parts[1])
                except ValueError:
                    continue
                if pos < lo or pos > hi:
                    continue
                ref, alt = parts[3], parts[4]
                key = (pos, ref, alt)
                if key not in keymap:
                    continue
                v = keymap[key]
                if v["rsid"] in got:
                    continue
                # EUR_AF 频率 QC
                m = EUR_AF_RE.search(parts[7])
                if m:
                    try:
                        vaf = float(m.group(1))
                        if abs(vaf - v["af"]) > 0.25:
                            continue
                    except ValueError:
                        pass
                dose = []
                ok = True
                for i in eur_idx:
                    gt = parts[9 + i].split(":")[0]
                    if "|" in gt:
                        a1, a2 = gt.split("|")
                    elif "/" in gt:
                        a1, a2 = gt.split("/")
                    else:
                        a1 = a2 = gt
                    if a1 == "." or a2 == ".":
                        ok = False; break
                    dose.append(int(a1) + int(a2))
                if not ok:
                    continue
                got[v["rsid"]] = dose
                meta[v["rsid"]] = (pos, ref, alt)
            del text, data
        P(f"{locus}: matched {len(got)} variants with EUR dosages")

        rsids = sorted(got, key=lambda r: meta[r][0])
        n = len(eur_idx) if eur_idx else len(eur)
        def corr(x, y):
            mx = sum(x) / n; my = sum(y) / n
            sxy = sxx = syy = 0.0
            for a, b in zip(x, y):
                dx, dy = a - mx, b - my
                sxy += dx * dy; sxx += dx * dx; syy += dy * dy
            d = math.sqrt(sxx * syy)
            return sxy / d if d > 0 else 0.0
        m = len(rsids)
        P(f"computing LD {m}x{m} ...")
        t0 = time.time()
        R = [[1.0] * m for _ in range(m)]
        for i in range(m):
            xi = got[rsids[i]]
            for j in range(i + 1, m):
                c = round(corr(xi, got[rsids[j]]), 4)
                R[i][j] = c; R[j][i] = c
            if (i + 1) % 100 == 0:
                P(f"  row {i+1}/{m} ({time.time()-t0:.0f}s)")
        with open(os.path.join(WORK, f"ld_{locus}.csv"), "w") as f:
            f.write("," + ",".join(rsids) + "\n")
            for i in range(m):
                f.write(rsids[i] + "," + ",".join(f"{v:.4f}" for v in R[i]) + "\n")
        json.dump({r: dict(pos=meta[r][0], ref=meta[r][1], alt=meta[r][2]) for r in rsids},
                  open(os.path.join(WORK, f"ld_variants_{locus}.json"), "w"))
        # 关键 SNP 校验
        for t in ("rs5068", "rs1421811"):
            if t in got:
                P(f"  {t}: matched at chr{chrom}:{meta[t][0]} {meta[t][1]}/{meta[t][2]}")
        P(f"{locus}: LD matrix saved ({m} variants)")
    P("DONE")
    log.close()
