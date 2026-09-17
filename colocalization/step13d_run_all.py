# -*- coding: utf-8 -*-
"""step13d: PLINK bfile + SMR BESD + SMR/HEIDI 全流程 (subprocess 驱动)"""
import os, subprocess

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
PLINK = os.path.join(WORK, "bin", "plink.exe")
SMR = os.path.join(WORK, "bin", "smr-1.3.1-win-x86_64", "smr-1.3.1-win.exe")
LOG = []
def P(s): LOG.append(str(s))

def run(cmd, logname):
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900,
                       encoding="utf-8", errors="ignore")
    with open(os.path.join(WORK, logname), "w", encoding="utf-8") as f:
        f.write((r.stdout or "") + "\n=== STDERR ===\n" + (r.stderr or ""))
    return r.returncode

# 1) PLINK bfile
for locus in ("nppa", "npr3"):
    if not os.path.exists(os.path.join(WORK, f"eur_{locus}.bed")):
        rc = run([PLINK, "--vcf", os.path.join(WORK, f"eur_{locus}.vcf"),
                  "--make-bed", "--out", os.path.join(WORK, f"eur_{locus}")],
                 f"step13d_plink_{locus}.log")
        P(f"plink {locus}: rc={rc}, bed exists={os.path.exists(os.path.join(WORK, f'eur_{locus}.bed'))}")
    else:
        P(f"plink {locus}: 已存在")
    n = sum(1 for _ in open(os.path.join(WORK, f"eur_{locus}.bim"))) if \
        os.path.exists(os.path.join(WORK, f"eur_{locus}.bim")) else 0
    P(f"  bim SNPs: {n}")

# 2) make-besd
for e in ("SCALLOP", "INTERVAL", "NPPAeQTL", "NPR3eQTL"):
    if not os.path.exists(os.path.join(WORK, f"{e}.besd")):
        rc = run([SMR, "--make-besd", os.path.join(WORK, f"{e}.efile"),
                  "--out", os.path.join(WORK, e)], f"step13d_besd_{e}.log")
        P(f"besd {e}: rc={rc}, exists={os.path.exists(os.path.join(WORK, f'{e}.besd'))}")
    else:
        P(f"besd {e}: 已存在")

# 3) SMR/HEIDI
PAIRS = [("nppa", "SCALLOP", o) for o in ("CAD_cardio","CAD_vdh","AF","HF","stroke","CES")] + \
        [("nppa", "INTERVAL", o) for o in ("CAD_cardio","CAD_vdh","AF","HF","stroke","CES")] + \
        [("nppa", "NPPAeQTL", o) for o in ("CAD_cardio","CAD_vdh","AF","HF","stroke","CES")] + \
        [("npr3", "NPR3eQTL", o) for o in ("CAD_cardio","CAD_vdh","MI_har","MI_finn")]

for loc, exp, out in PAIRS:
    tag = f"{exp}__x__{out}_{loc}"
    dest = os.path.join(WORK, f"smr_{tag}.smr")
    if os.path.exists(dest):
        P(f"smr {tag}: 已存在"); continue
    rc = run([SMR,
              "--bfile", os.path.join(WORK, f"eur_{loc}"),
              "--gwas-summary", os.path.join(WORK, f"ma_{out}_{loc}.ma"),
              "--beqtl-summary", os.path.join(WORK, exp),
              "--out", os.path.join(WORK, f"smr_{tag}"),
              "--thread-num", "2"], f"step13d_smr_{tag}.log")
    ok = os.path.exists(dest)
    P(f"smr {tag}: rc={rc}, result={'OK' if ok else 'NO OUTPUT'}")

with open(os.path.join(WORK, "step13d_output.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(LOG))
print("DONE")
