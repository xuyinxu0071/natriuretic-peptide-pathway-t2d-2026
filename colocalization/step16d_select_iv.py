# -*- coding: utf-8 -*-
"""
Step 16d: M1 工具变量选择 + 结局关联拉取
================================================================================
暴露: corin/furin/mme/dpp4/npr2 eQTLGen cis-eQTL + furin INTERVAL protein (prot-a-1150)
选择: cis 窗口内 p<5e-8, LD clump r²<0.01 (1000G EUR), 保留 |z| 最强者
结局: AF, CAD(ieu-a-7), CAD(GCST005195), HF, T2D —— 仅拉取工具 SNP
"""
import json, math, os, time, urllib.error, urllib.request

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
TOKEN = open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\.opengwas_token.txt").read().strip()
API = "https://api.opengwas.io/api"

log = open(os.path.join(WORK, "step16d_log.txt"), "w", encoding="utf-8")
def P(s):
    print(s, flush=True); log.write(s + "\n"); log.flush()

def api_post(path, payload, retries=6):
    for i in range(retries):
        try:
            req = urllib.request.Request(
                API + path, data=json.dumps(payload).encode("utf-8"),
                headers={"Authorization": "Bearer " + TOKEN,
                         "Content-Type": "application/json"})
            return json.load(urllib.request.urlopen(req, timeout=300))
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and i < retries - 1:
                wait = 20 * (i + 1) if e.code == 429 else 10
                P("  [HTTP %d] retry in %ds" % (e.code, wait))
                time.sleep(wait); continue
            P("  [HTTP %d] %s" % (e.code, e.read().decode("utf-8", "ignore")[:200]))
            return None
        except Exception as e:
            if i < retries - 1:
                time.sleep(15); continue
            P("  [ERR] %s" % e); return None

# ---------- LD 矩阵载入 ----------
def load_ld(gene):
    path = os.path.join(WORK, "ld_%s.csv" % gene)
    with open(path) as f:
        hdr = f.readline().strip().split(",")
        idx = {r: i for i, r in enumerate(hdr)}
        R = {}
        for line in f:
            p = line.strip().split(",")
            R[p[0]] = [float(x) for x in p[1:]]
    return idx, R

def load_assoc(locus, dsid):
    return json.load(open(os.path.join(WORK, "%s__%s_assoc.json" % (locus, dsid))))

EXPOSURES = [
    ("corin", "eqtl-a-ENSG00000145244", "CORIN_eQTL"),
    ("furin", "eqtl-a-ENSG00000140564", "FURIN_eQTL"),
    ("furin", "prot-a-1150",            "FURIN_pQTL"),
    ("mme",   "eqtl-a-ENSG00000196549", "MME_eQTL"),
    ("dpp4",  "eqtl-a-ENSG00000197635", "DPP4_eQTL"),
    ("npr2",  "eqtl-a-ENSG00000159899", "NPR2_eQTL"),
]
OUTCOMES = {
    "AF":   "ebi-a-GCST006061",
    "CAD":  "ieu-a-7",
    "CAD2": "ebi-a-GCST005195",
    "HF":   "ebi-a-GCST009541",
    "T2D":  "ebi-a-GCST90018926",   # Sakaue 2021 EUR (GCST007515 覆盖率 ~1% 弃用)
    "T2D_FINN": "finn-b-T2D",       # FinnGen R9 复制
}
GWS = 5e-8
CLUMP_R2 = 0.01

if __name__ == "__main__":
    instruments_all = {}
    for locus, dsid, name in EXPOSURES:
        P("=== %s (%s x %s) ===" % (name, locus, dsid))
        assoc = load_assoc(locus, dsid)
        idx, R = load_ld(locus)
        # GWS 候选
        cands = [(rs, a) for rs, a in assoc.items()
                 if a.get("p") is not None and a["p"] < GWS
                 and rs in idx]
        cands.sort(key=lambda x: x[1]["p"])
        P("  GWS candidates (in LD panel): %d" % len(cands))
        # clump
        kept = []
        for rs, a in cands:
            ok = True
            for k in kept:
                r = R[rs][idx[k]]
                if r * r > CLUMP_R2:
                    ok = False; break
            if ok:
                kept.append(rs)
        P("  after clump r2<%.2f: %d -> %s" %
          (CLUMP_R2, len(kept), kept[:20]))
        instruments_all[name] = dict(
            locus=locus, dsid=dsid, ivs=[dict(
                rsid=rs, **assoc[rs]) for rs in kept])
        if not kept:
            # 无 GWS 时记录最强 cis 信号供参考
            best = sorted(assoc.items(),
                          key=lambda x: x[1].get("p") or 1)[:1]
            if best:
                P("  NO GWS; best cis signal: %s p=%.3g" %
                  (best[0][0], best[0][1].get("p", float("nan"))))
                instruments_all[name]["best_nongws"] = dict(
                    rsid=best[0][0], **best[0][1])
    json.dump(instruments_all,
              open(os.path.join(WORK, "step16d_instruments.json"), "w"))
    P("instruments saved")

    # ---------- 结局拉取(全部工具 SNP 并集) ----------
    all_iv = sorted({iv["rsid"] for x in instruments_all.values()
                     for iv in x["ivs"]})
    P("union of instrument SNPs: %d" % len(all_iv))
    for oname, oid in OUTCOMES.items():
        dest = os.path.join(WORK, "step16d_outcome__%s.json" % oname)
        if os.path.exists(dest):
            P("  %s: exists, skip" % oname); continue
        P("  pulling %s (%s) at %d SNPs ..." % (oname, oid, len(all_iv)))
        res = {}
        done = 0
        while done < len(all_iv):
            chunk = all_iv[done:done + 64]
            r = api_post("/associations", {"variant": chunk, "id": [oid]})
            if r is None:
                P("    FAILED at %d" % done); break
            for a in r:
                if a.get("id") == oid and a.get("beta") is not None:
                    res[a["rsid"]] = {k: a.get(k) for k in
                                      ("beta", "se", "p", "ea", "nea", "eaf", "n")}
            done += 64
            time.sleep(0.6)
        json.dump(res, open(dest, "w"))
        P("  %s: %d/%d hit" % (oname, len(res), len(all_iv)))
    P("ALL DONE")
    log.close()
