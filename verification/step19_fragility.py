# -*- coding: utf-8 -*-
"""
step19_fragility.py — Fragility index for the three primary observational
hypotheses (Walsh 2014 concept, adapted to Cox models by event flipping).

Data frames and Cox specifications replicate main_analysis_v4.py exactly:
  H1: HRS DM-controlled stratum, per-IQR ln NT-proBNP -> all-cause death
      (12 covariates incl. ln cystatin C)                    n=510, 98 deaths
  H2: HRS DM-controlled & baseline CVD-free, per-IQR ln NT-proBNP -> incident
      CVD (11 covariates, base_cvd dropped by design)         n=330, 63 events
  H3: CHARLS DM-controlled & baseline CVD-free, per-IQR ln hsCRP -> incident
      CVD (7 covariates, cluster-robust by community)          n=361, 106 events

Fragility index = minimum number of events that must be reclassified as
non-events (censored at their observed time) to push the exposure's
two-sided Wald p >= 0.05, flipping at each step the single event whose
reclassification maximizes the resulting p-value (greedy exact search).
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import pyreadstat
from lifelines import CoxPHFitter

HRS_BASE = r"D:\2026\3. HRS  美国\HRS_美国"
CH_BASE = r"D:\2026\1. CHARLS  中国\CHARLS_中国"

LOG = []


def P(s=""):
    LOG.append(str(s))
    print(s, flush=True)


# ============================================================
# PART 1 — HRS frame (verbatim from main_analysis_v4.py)
# ============================================================
vbs, _ = pyreadstat.read_sav(HRS_BASE + r"\Raw_data\2016 静脉血研究VBS\HRS2016VBS\hrs2016vbs.sav")
vbs["hhidpn"] = vbs["HHID"].astype(int) * 1000 + vbs["PN"].astype(int)
vbs = vbs[["hhidpn", "PNTBNPE", "PCRP", "PCYSC", "PVBSWGTR"]].rename(
    columns={"PNTBNPE": "ntbnp", "PCRP": "crp", "PCYSC": "cysc", "PVBSWGTR": "vbswgt"})

dbs, _ = pyreadstat.read_dta(HRS_BASE + r"\Raw_data\2016 HRS\BIOMK16BL\BIOMK16BL_R.dta")
dbs["hhidpn"] = dbs["HHID"].astype(int) * 1000 + dbs["PN"].astype(int)
dbs = dbs[["hhidpn", "PA1C_ADJ"]].rename(columns={"PA1C_ADJ": "a1c"})

hrs, _ = pyreadstat.read_dta(HRS_BASE + r"\Working_data\hrs.dta")
w12 = hrs[hrs["wave"] == 12][["hhidpn", "ragender", "rabyear", "raracem", "raedyrs", "bmi",
                              "smoken", "drink", "hibpe", "diabe", "hearte", "stroke",
                              "cancre", "shlt"]].copy()
w12["age"] = 2016 - w12["rabyear"]

fu = hrs[hrs["wave"].isin([13, 14, 15])][["hhidpn", "wave", "hearte", "stroke"]]
fu = fu.sort_values(["hhidpn", "wave"])
first_ev = fu.groupby("hhidpn").apply(
    lambda g: g.loc[(g["hearte"] == 1) | (g["stroke"] == 1), "wave"].min(), include_groups=False)
last_ob = fu.groupby("hhidpn")["wave"].max()

trk, _ = pyreadstat.read_dta(HRS_BASE + r"\Raw_data\2022 HRS\2022 HRS Core\h22core\trk2022\trk2022tr_r.dta")
trk["hhidpn"] = trk["HHID"].astype(int) * 1000 + trk["PN"].astype(int)
trk = trk[["hhidpn", "EXDEATHYR", "EXDEATHMO"]]

d = (vbs.merge(dbs, on="hhidpn").merge(w12, on="hhidpn").merge(trk, on="hhidpn", how="left")
     .set_index("hhidpn"))
d = d.join(first_ev.rename("ev_wave")).join(last_ob.rename("last_wave"))


def glyrow(r):
    a = r["a1c"]
    if pd.isna(a):
        return np.nan
    if r["diabe"] == 1:
        return 1 if a < 7.0 else 2
    if a >= 6.5:
        return 3
    if a >= 5.7:
        return 4
    return 5


d["gly"] = d.apply(glyrow, axis=1)
d["death"] = ((d["EXDEATHYR"] >= 2016) & (d["EXDEATHYR"] <= 2023)).astype(int)
t_ev = np.where(d["EXDEATHMO"].notna(),
                d["EXDEATHYR"] + (d["EXDEATHMO"].astype(float) - 0.5) / 12.0,
                d["EXDEATHYR"] + 0.5)
d["t_death"] = np.where(d["death"] == 1, np.maximum(t_ev - 2016.5, 0.1), 7.5)
d["base_cvd"] = ((d["hearte"] == 1) | (d["stroke"] == 1)).astype(int)
d["new_cvd"] = 0
d.loc[(d["base_cvd"] == 0) & d["ev_wave"].notna(), "new_cvd"] = 1
d["t_cvd"] = np.where(d["new_cvd"] == 1,
                      (d["ev_wave"] - 12) * 2.0,
                      np.where(d["last_wave"].notna(), (d["last_wave"] - 12) * 2.0, np.nan))
d.loc[d["t_cvd"].isna() & (d["base_cvd"] == 0) & d["gly"].notna(), "t_cvd"] = 5.0
d["lnbnp"] = np.log(d["ntbnp"].clip(lower=1))
d["lncysc"] = np.log(d["cysc"].clip(lower=0.1))
d["female"] = (d["ragender"] == 0).astype(int)
d["smoke_now"] = (d["smoken"] == 1).astype(int)
d["drink_now"] = (d["drink"] == 1).astype(int)

mm = d[d["gly"].notna() & d["ntbnp"].notna()].copy()
goal = mm[mm["gly"] == 1].copy()
goal["bnp_iqr"] = goal["lnbnp"] / (goal["lnbnp"].quantile(0.75) - goal["lnbnp"].quantile(0.25))

COV_H = ["age", "female", "raracem", "raedyrs", "bmi", "smoke_now", "drink_now",
         "hibpe", "cancre", "shlt", "base_cvd", "lncysc"]

goal3 = goal[goal["base_cvd"] == 0].copy()
COV_H3 = [c for c in COV_H if c != "base_cvd"]

# ============================================================
# PART 2 — CHARLS frame (verbatim from main_analysis_v4.py)
# ============================================================
ch, _ = pyreadstat.read_dta(CH_BASE + r"\Working_data\charls.dta")
w3 = ch[ch["wave"] == 3][["ID", "ragender", "age", "bmi", "raeducl", "hibpe", "diabe",
                          "hearte", "stroke", "bl_crp", "bl_hbalc", "bl_cysc",
                          "bl_tg", "bl_hdl", "cancre",
                          "communityID", "bloodweight",
                          "smoken", "drinkl"]].copy()
fu4 = ch[ch["wave"].isin([4, 5])][["ID", "wave", "hearte", "stroke"]].sort_values(["ID", "wave"])
first_ev4 = fu4.groupby("ID").apply(
    lambda g: g.loc[(g["hearte"] == 1) | (g["stroke"] == 1), "wave"].min(), include_groups=False)
last_ob4 = fu4.groupby("ID")["wave"].max()
c = w3.set_index("ID").join(first_ev4.rename("ev_wave")).join(last_ob4.rename("last_wave"))


def glyrow_c(r):
    a = r["bl_hbalc"]
    if pd.isna(a):
        return np.nan
    if r["diabe"] == 1:
        return 1 if a < 7.0 else 2
    if a >= 6.5:
        return 3
    return 5


c["gly"] = c.apply(glyrow_c, axis=1)
c["base_cvd"] = ((c["hearte"] == 1) | (c["stroke"] == 1)).astype(int)
c["new_cvd"] = 0
c.loc[(c["base_cvd"] == 0) & c["ev_wave"].notna(), "new_cvd"] = 1
c["t_cvd"] = np.where(c["new_cvd"] == 1, (c["ev_wave"] - 3) * 2.5,
                      np.where(c["last_wave"].notna(), (c["last_wave"] - 3) * 2.5, np.nan))
c.loc[c["t_cvd"].isna() & (c["base_cvd"] == 0) & c["gly"].notna(), "t_cvd"] = 5.0

cm = c[c["gly"].notna() & c["bl_crp"].notna() & (c["base_cvd"] == 0)].copy()
cm["lncrp"] = np.log(cm["bl_crp"].clip(lower=0.1))
cm["crp_iqr"] = cm["lncrp"] / (cm["lncrp"].quantile(0.75) - cm["lncrp"].quantile(0.25))
cm["female"] = (cm["ragender"] == 0).astype(int)
cm["smoke_now"] = (cm["smoken"] == 1).astype(int)
cm["drink_now"] = (cm["drinkl"] == 1).astype(int)
goalc = cm[cm["gly"] == 1].copy()
COV_C = ["age", "female", "bmi", "raeducl", "hibpe", "smoke_now", "drink_now"]

# ============================================================
# PART 3 — fragility index (greedy exact event flipping)
# ============================================================
def fit_p(df, exposure, covs, event, time, cluster=None):
    cols = [exposure] + covs + [event, time] + ([cluster] if cluster else [])
    sub = df[cols].dropna()
    cph = CoxPHFitter()
    cph.fit(sub, duration_col=time, event_col=event, robust=True,
            cluster_col=cluster if cluster else None)
    return cph.summary.loc[exposure, "p"], len(sub), int(sub[event].sum())


def fragility(df, exposure, covs, event, time, cluster=None, label=""):
    cols = [exposure] + covs + [event, time] + ([cluster] if cluster else [])
    base = df[cols].dropna().copy()
    p0, n, ev = fit_p(df, exposure, covs, event, time, cluster)
    P("\n=== %s ===" % label)
    P("  n=%d, events=%d, baseline p=%.4g" % (n, ev, p0))
    if p0 >= 0.05:
        P("  FI = 0 (baseline not significant)")
        return 0, p0, p0
    work = base.copy()
    ev_idx = list(work.index[work[event] == 1])
    flipped = []
    p_cur = p0
    for step in range(len(ev_idx)):
        best_p, best_i = -1.0, None
        for i in ev_idx:
            work.loc[i, event] = 0
            try:
                p, _, _ = fit_p(work, exposure, covs, event, time, cluster)
            except Exception:
                p = 1.0
            work.loc[i, event] = 1
            if p > best_p:
                best_p, best_i = p, i
        work.loc[best_i, event] = 0
        ev_idx.remove(best_i)
        flipped.append(best_i)
        p_cur = best_p
        P("  flip %2d: %s -> p=%.4g" % (step + 1, best_i, p_cur))
        if p_cur >= 0.05:
            P("  FI = %d (p reached %.4g after %d flips of %d events)"
              % (step + 1, p_cur, step + 1, ev))
            return step + 1, p0, p_cur
    P("  FI > %d (all events flipped without reaching p>=0.05)" % ev)
    return None, p0, p_cur


res = []
fi, p0, p1 = fragility(goal, "bnp_iqr", COV_H, "death", "t_death",
                      label="H1 HRS controlled: lnNT-proBNP -> all-cause death")
res.append(("H1_HRS_death", fi, p0, p1))
fi, p0, p1 = fragility(goal3, "bnp_iqr", COV_H3, "new_cvd", "t_cvd",
                      label="H2 HRS CVD-free controlled: lnNT-proBNP -> incident CVD")
res.append(("H2_HRS_incidentCVD", fi, p0, p1))
fi, p0, p1 = fragility(goalc, "crp_iqr", COV_C, "new_cvd", "t_cvd",
                      cluster="communityID",
                      label="H3 CHARLS controlled CVD-free: lnhsCRP -> incident CVD")
res.append(("H3_CHARLS_incidentCVD", fi, p0, p1))

P("\n================ SUMMARY ================")
for name, fi, p0, p1 in res:
    P("%-24s FI=%-4s baseline p=%.4g -> final p=%.4g" % (name, fi, p0, p1))

with open(r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\step19_fragility_results.txt",
          "w", encoding="utf-8") as f:
    f.write("\n".join(LOG))
print("\nDONE -> step19_fragility_results.txt")
