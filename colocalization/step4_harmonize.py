# -*- coding: utf-8 -*-
"""
Step 4: 等位基因协调 — 为每对 coloc 组装输入 CSV
================================================================================
规则:
  1. 按 rsid 匹配 dbSNP 变体列表与 OpenGWAS 关联记录;
  2. 等位对齐: OpenGWAS ea/nea vs dbSNP ref/alt — 反向则翻转 beta 与 eaf;
  3. 丢弃回文(A/T, C/G)与等位不匹配;
  4. eaf 缺失时用 dbSNP CAF(对齐后)兜底; eaf 与 CAF 偏差>0.25 时丢弃(QC);
  5. 产出 coloc_inputs/{pair}_exposure.csv 与 {pair}_outcome.csv
     列: rsid, pos, beta, se, eaf, n, ea, nea
"""
import json, os

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
INDIR = os.path.join(WORK, "coloc_inputs")
os.makedirs(INDIR, exist_ok=True)

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

# pair -> (locus, exposure_ds, outcome_ds, label)
PAIRS = [
    # NPPA eQTL × 结局
    ("nppa_eqtl_cad_vdh",    "nppa", "eqtl-a-ENSG00000175206", "ebi-a-GCST005195"),
    ("nppa_eqtl_cad_cardio", "nppa", "eqtl-a-ENSG00000175206", "ieu-a-7"),
    ("nppa_eqtl_af",         "nppa", "eqtl-a-ENSG00000175206", "ebi-a-GCST006061"),
    ("nppa_eqtl_hf",         "nppa", "eqtl-a-ENSG00000175206", "ebi-a-GCST009541"),
    ("nppa_eqtl_stroke",     "nppa", "eqtl-a-ENSG00000175206", "ebi-a-GCST005838"),
    ("nppa_eqtl_ces",        "nppa", "eqtl-a-ENSG00000175206", "ebi-a-GCST006910"),
    # NPR3 eQTL × 结局
    ("npr3_eqtl_cad_vdh",    "npr3", "eqtl-a-ENSG00000113389", "ebi-a-GCST005195"),
    ("npr3_eqtl_cad_cardio", "npr3", "eqtl-a-ENSG00000113389", "ieu-a-7"),
    ("npr3_eqtl_mi_har",     "npr3", "eqtl-a-ENSG00000113389", "ebi-a-GCST011364"),
    ("npr3_eqtl_mi_finn",    "npr3", "eqtl-a-ENSG00000113389", "finn-b-I9_MI"),
    # SBP × 结局 (NPR3 位点, 检验血压途径的共定位)
    ("npr3_sbp_cad_vdh",     "npr3", "ieu-b-38", "ebi-a-GCST005195"),
    ("npr3_sbp_cad_cardio",  "npr3", "ieu-b-38", "ieu-a-7"),
    ("npr3_sbp_mi_har",      "npr3", "ieu-b-38", "ebi-a-GCST011364"),
    ("npr3_sbp_mi_finn",     "npr3", "ieu-b-38", "finn-b-I9_MI"),
]

def load_json(p):
    return json.load(open(p))

def harmonize(variants, assoc):
    """variants: [{rsid,pos,ref,alt,af}], assoc: {rsid:{beta,se,p,ea,nea,eaf,n}}
    返回 {rsid: (pos, beta, se, eaf, n, ea, nea)} 对齐到 dbSNP ref/alt 方向"""
    out = {}
    for v in variants:
        rsid = v["rsid"]
        a = assoc.get(rsid)
        if a is None or a.get("beta") is None or a.get("se") is None:
            continue
        ea, nea = (a.get("ea") or "").upper(), (a.get("nea") or "").upper()
        if len(ea) != 1 or len(nea) != 1:
            continue
        try:
            eaf_raw = float(a.get("eaf"))
        except (TypeError, ValueError):
            eaf_raw = None
        ref, alt = v["ref"], v["alt"]
        beta, se = a["beta"], a["se"]
        if ea == ref and nea == alt:
            # OpenGWAS beta 以 ref 等位为效应方向 -> 翻转到 alt
            beta, eaf = -beta, (1 - eaf_raw) if eaf_raw is not None else None
        elif ea == alt and nea == ref:
            eaf = eaf_raw
        elif ea == COMPLEMENT.get(ref) and nea == COMPLEMENT.get(alt):
            beta, eaf = -beta, (1 - eaf_raw) if eaf_raw is not None else None
        elif ea == COMPLEMENT.get(alt) and nea == COMPLEMENT.get(ref):
            eaf = eaf_raw
        else:
            continue  # 等位不匹配
        # 回文丢弃
        if ref in "AT" and alt in "AT" and {ref, alt} == {"A", "T"}:
            continue
        if {ref, alt} == {"C", "G"}:
            continue
        # eaf 兜底与QC (eaf 定义为 alt 等位频率)
        if eaf is None:
            eaf = v["af"]
        else:
            if abs(eaf - v["af"]) > 0.25:
                continue  # 频率严重不一致 -> 疑似错配
        if eaf is None or not (0.0 < eaf < 1.0):
            continue
        n = a.get("n")
        out[rsid] = (v["pos"], beta, se, eaf, n, alt, ref)
    return out

if __name__ == "__main__":
    log = open(os.path.join(WORK, "step4_log.txt"), "w", encoding="utf-8")
    def P(s):
        print(s, flush=True); log.write(s + "\n"); log.flush()

    vcache = {}
    for pair, locus, exp_ds, out_ds in PAIRS:
        if locus not in vcache:
            vcache[locus] = load_json(os.path.join(WORK, f"{locus}_locus_variants.json"))
        exp_f = os.path.join(WORK, f"{locus}__{exp_ds}_assoc.json")
        out_f = os.path.join(WORK, f"{locus}__{out_ds}_assoc.json")
        if not (os.path.exists(exp_f) and os.path.exists(out_f)):
            P(f"{pair}: MISSING assoc file ({exp_ds} or {out_ds}) — skip")
            continue
        exp_h = harmonize(vcache[locus], load_json(exp_f))
        out_h = harmonize(vcache[locus], load_json(out_f))
        shared = sorted(set(exp_h) & set(out_h), key=lambda r: exp_h[r][0])
        # 仅保留两数据集都有的 rsid
        def rows(h, keys):
            return [",".join(str(x) for x in (r, h[r][0], h[r][1], h[r][2], h[r][3], h[r][4] or "", h[r][5], h[r][6])) for r in keys]
        with open(os.path.join(INDIR, f"{pair}_exposure.csv"), "w") as f:
            f.write("rsid,pos,beta,se,eaf,n,ea,nea\n" + "\n".join(rows(exp_h, shared)) + "\n")
        with open(os.path.join(INDIR, f"{pair}_outcome.csv"), "w") as f:
            f.write("rsid,pos,beta,se,eaf,n,ea,nea\n" + "\n".join(rows(out_h, shared)) + "\n")
        P(f"{pair}: exposure {len(exp_h)} / outcome {len(out_h)} -> shared {len(shared)} SNPs")
    P("DONE")
    log.close()
