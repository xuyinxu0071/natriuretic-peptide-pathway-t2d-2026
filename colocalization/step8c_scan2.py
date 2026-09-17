# -*- coding: utf-8 -*-
import json, traceback
WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
OUT = open(WORK + r"\step8c_output.txt", "w", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a)
    OUT.write(s + "\n"); OUT.flush(); print(s, flush=True)
try:
    d = json.load(open(WORK + r"\gwasinfo_full.json", encoding="utf-8"))
    log("type:", type(d), "len:", len(d))
    x = d[0] if isinstance(d, list) else list(d.values())[0] if isinstance(d, dict) else None
    log("keys:", sorted(x.keys()))
    log("sample trait:", x.get("trait"))
    KEYS = ["natriuretic", "probnp", "nt-pro", "nppa", "nppb", "npr3", "corin",
            "neprilysin", "bnp", "guanylyl cyclase"]
    hits = [e for e in d.values() if any(k in str(e.get("trait", "")).lower() for k in KEYS)]
    log("trait-matched:", len(hits))
    for e in sorted(hits, key=lambda z: str(z.get("id"))):
        log(f"  {e.get('id')}  n={e.get('sample_size')}  {e.get('year')}  {str(e.get('population'))[:12]}  {str(e.get('trait'))[:75]}")
except Exception:
    log(traceback.format_exc())
OUT.close()
