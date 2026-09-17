# -*- coding: utf-8 -*-
"""
Step 16a: 为 NP 加工酶级联 5 个新基因构建 locus 变体表 + 1000G EUR LD 矩阵
================================================================================
基因体(GRCh37) ±100kb 窗口，复用 step2i (dbSNP b151 common VCF) 与 step6 (1000G v5b) 机制。
产出: {gene}_locus_variants.json + ld_{gene}.csv + ld_variants_{gene}.json
"""
import gzip, json, math, os, re, struct, time, urllib.request, zlib

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
DBSNP_BASE = "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/"
DBSNP_VCF = "common_all_20180423.vcf.gz"
KG_BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
KG_PANEL = "integrated_call_samples_v3.20130502.ALL.panel"
TBI_DBSNP = os.path.join(WORK, "dbsnp_common.tbi")

GENES = ["corin", "furin", "mme", "dpp4", "npr2"]
FLANK = 100000

log = open(os.path.join(WORK, "step16a_log.txt"), "w", encoding="utf-8")
def P(s):
    print(s, flush=True); log.write(s + "\n"); log.flush()

# ---------- 通用 HTTP ----------
def fetch(url, start=None, end=None, retries=5, timeout=1800):
    for i in range(retries):
        try:
            hdr = {"User-Agent": "coloc-pipeline/1.0"}
            if start is not None:
                hdr["Range"] = "bytes=%d-%d" % (start, end)
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

def parse_tbi_file(path):
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
            for c in range(n_chunk):
                cb, ce = struct.unpack_from("<QQ", raw, off); off += 16
                bins.append((bin_no, [(cb, ce)] if False else None))
                # 重读: 手动展开
                bins[-1] = (bin_no, [])
            # 重新解析 chunks
            off -= 16 * n_chunk
            chunks = []
            for c in range(n_chunk):
                cb, ce = struct.unpack_from("<QQ", raw, off); off += 16
                chunks.append((cb, ce))
            bins[-1] = (bin_no, chunks) if isinstance(bins[-1][0], int) else bins[-1]
        n_intv = struct.unpack_from("<i", raw, off)[0]; off += 4
        off += 8 * n_intv
        refs[names[r]] = bins
    return refs, names

# 上面 parse_tbi_file 写得太绕, 直接照抄 step6 的干净实现
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

def merge_ranges(ranges, pad=100000):
    ranges = sorted(ranges)
    merged = []
    for s, e in ranges:
        if merged and s <= merged[-1][1] + pad:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged

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

# ---------- 1) GRCh37 基因坐标 ----------
def ensembl_grch37(symbol):
    url = "https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/%s" % symbol
    req = urllib.request.Request(url, headers={
        "Content-Type": "application/json", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))

REGIONS = {}
for g in GENES:
    rec = ensembl_grch37(g.upper())
    REGIONS[g] = (rec["seq_region_name"], rec["start"] - FLANK, rec["end"] + FLANK,
                  rec["id"])
    P("%s: %s chr%s:%d-%d (%s)" % (g.upper(), rec["id"], rec["seq_region_name"],
                                   rec["start"], rec["end"], rec.get("biotype")))
json.dump({g: dict(ensg=REGIONS[g][3], chrom=REGIONS[g][0], lo=REGIONS[g][1],
                   hi=REGIONS[g][2]) for g in GENES},
          open(os.path.join(WORK, "np_cascade_regions.json"), "w"))
P("regions: " + json.dumps(REGIONS))

# ---------- 2) dbSNP 变体表 ----------
if not os.path.exists(TBI_DBSNP) or os.path.getsize(TBI_DBSNP) < 1e5:
    P("downloading dbSNP tbi ...")
    open(TBI_DBSNP, "wb").write(fetch(DBSNP_BASE + DBSNP_VCF + ".tbi"))
dbsnp_refs = parse_tbi(TBI_DBSNP)
SNV_RE = re.compile(r"^[ACGT]$")

for gene, (chrom, lo, hi, ensg) in REGIONS.items():
    dest = os.path.join(WORK, "%s_locus_variants.json" % gene)
    if os.path.exists(dest):
        P("%s: variants exist, skip" % gene)
        continue
    bins = dbsnp_refs.get(chrom)
    want = set(reg2bins(lo - 1, hi))
    chunks = []
    for bin_no, chs in bins:
        if bin_no in want:
            chunks.extend(chs)
    ranges = [(cb >> 16, (ce >> 16) + 70000) for cb, ce in chunks]
    merged = merge_ranges(ranges)
    P("%s chr%s:%d-%d: %d chunks -> %d ranges, ~%d KB" %
      (gene, chrom, lo, hi, len(chunks), len(merged),
       sum(e - s for s, e in merged) // 1000))
    variants, seen = [], set()
    for s, e in merged:
        data = fetch(DBSNP_BASE + DBSNP_VCF, s, e)
        text = decompress_bgzf(data).decode("utf-8", "ignore")
        for line in text.split("\n"):
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8 or parts[0] != chrom:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            if pos < lo or pos > hi:
                continue
            rsid = parts[2]
            if not rsid.startswith("rs"):
                continue
            ref, alt = parts[3].upper(), parts[4].upper()
            if "," in alt or alt == "." or not SNV_RE.match(ref) or not SNV_RE.match(alt):
                continue
            m = re.search(r"(?:^|;)CAF=([0-9.,]+)(?:;|$)", parts[7])
            if not m:
                continue
            vals = m.group(1).split(",")
            if len(vals) != 2 or vals[1] in (".", ""):
                continue
            try:
                af = float(vals[1])
            except ValueError:
                continue
            if af > 1.0:
                af /= 100.0
            if min(af, 1.0 - af) < 0.01:
                continue
            if rsid in seen:
                continue
            seen.add(rsid)
            variants.append(dict(rsid=rsid, chrom=chrom, pos=pos,
                                 ref=ref, alt=alt, af=round(af, 4)))
        del text, data
    variants.sort(key=lambda v: v["pos"])
    json.dump(variants, open(dest, "w"))
    P("%s: %d common biallelic SNVs saved" % (gene, len(variants)))

# ---------- 3) 1000G EUR LD 矩阵 ----------
panel_path = os.path.join(WORK, "kg_panel.txt")
if not os.path.exists(panel_path):
    open(panel_path, "wb").write(fetch(KG_BASE + KG_PANEL))
eur = [l.split("\t")[0] for l in open(panel_path).read().splitlines()
       if l and not l.startswith("sample") and len(l.split("\t")) >= 3
       and l.split("\t")[2] == "EUR"]
eur_set = set(eur)
P("EUR samples: %d" % len(eur))
EUR_AF_RE = re.compile(r"EUR_AF=([0-9.eE+-]+)")

for gene, (chrom, lo, hi, ensg) in REGIONS.items():
    ld_csv = os.path.join(WORK, "ld_%s.csv" % gene)
    if os.path.exists(ld_csv):
        P("%s: LD exists, skip" % gene)
        continue
    variants = json.load(open(os.path.join(WORK, "%s_locus_variants.json" % gene)))
    keymap = {(v["pos"], v["ref"], v["alt"]): v for v in variants}
    P("%s: %d target variants in chr%s:%d-%d" % (gene, len(keymap), chrom, lo, hi))
    vcf = ("ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5b."
           "20130502.genotypes.vcf.gz" % chrom)
    tbi = os.path.join(WORK, vcf + ".tbi")
    if not os.path.exists(tbi) or os.path.getsize(tbi) < 1e5:
        P("downloading chr%s tbi ..." % chrom)
        open(tbi, "wb").write(fetch(KG_BASE + vcf + ".tbi"))
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
    P("  %d chunks -> %d ranges, ~%.1f MB" %
      (len(chunks), len(merged), sum(e - s for s, e in merged) / 1e6))

    got, meta = {}, {}
    data = fetch(KG_BASE + vcf, 0, 400000)
    d = zlib.decompressobj(47)
    text = (d.decompress(data) + d.flush()).decode("ascii", "ignore")
    samples = None
    for line in text.split("\n"):
        if line.startswith("#CHROM"):
            samples = line.split("\t")[9:]
            break
    del data, text
    if samples is None:
        P("  #CHROM header FAILED"); continue
    eur_idx = [i for i, s in enumerate(samples) if s in eur_set]
    P("  VCF samples: %d, EUR matched: %d" % (len(samples), len(eur_idx)))
    if len(eur_idx) < 400:
        P("  EUR match FAILED"); continue
    for s, e in merged:
        data = fetch(KG_BASE + vcf, s, e)
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
    P("  %s: matched %d variants with EUR dosages" % (gene, len(got)))

    rsids = sorted(got, key=lambda r: meta[r][0])
    n = len(eur_idx)
    def corr(x, y):
        mx = sum(x) / n; my = sum(y) / n
        sxy = sxx = syy = 0.0
        for a, b in zip(x, y):
            dx, dy = a - mx, b - my
            sxy += dx * dy; sxx += dx * dx; syy += dy * dy
        dd = math.sqrt(sxx * syy)
        return sxy / dd if dd > 0 else 0.0
    m = len(rsids)
    P("  computing LD %dx%d ..." % (m, m))
    t0 = time.time()
    R = [[1.0] * m for _ in range(m)]
    for i in range(m):
        xi = got[rsids[i]]
        for j in range(i + 1, m):
            c = round(corr(xi, got[rsids[j]]), 4)
            R[i][j] = c; R[j][i] = c
        if (i + 1) % 200 == 0:
            P("    row %d/%d (%.0fs)" % (i + 1, m, time.time() - t0))
    with open(ld_csv, "w") as f:
        f.write("," + ",".join(rsids) + "\n")
        for i in range(m):
            f.write(rsids[i] + "," + ",".join("%.4f" % v for v in R[i]) + "\n")
    json.dump({r: dict(pos=meta[r][0], ref=meta[r][1], alt=meta[r][2]) for r in rsids},
              open(os.path.join(WORK, "ld_variants_%s.json" % gene), "w"))
    P("  %s: LD matrix saved (%d variants)" % (gene, m))

P("DONE")
log.close()
