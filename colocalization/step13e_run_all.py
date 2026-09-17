# -*- coding: utf-8 -*-
"""step13e: SMR 全流程 v2 — cwd + ASCII 相对路径规避中文路径 argv 损坏"""
import os, subprocess

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
SMR = os.path.join(WORK, "bin", "smr-1.3.1-win-x86_64", "smr-1.3.1-win.exe")
LOG = []
def P(s): LOG.append(str(s))

def run_rel(args, logname):
    r = subprocess.run([SMR] + args, capture_output=True, text=True, timeout=900,
                       encoding="utf-8", errors="ignore", cwd=WORK)
    with open(os.path.join(WORK, logname), "w", encoding="utf-8") as f:
        f.write((r.stdout or "") + "\n=== STDERR ===\n" + (r.stderr or ""))
    return r.returncode

for e in ("SCALLOP", "INTERVAL", "NPPAeQTL", "NPR3eQTL"):
    if not os.path.exists(os.path.join(WORK, f"{e}.besd")):
        rc = run_rel(["--efile", f"{e}.efile", "--make-besd", "--out", e], f"step13e_besd_{e}.log")
        P(f"besd {e}: rc={rc}, exists={os.path.exists(os.path.join(WORK, f'{e}.besd'))}")
    else:
        P(f"besd {e}: 已存在")

PAIRS = [("nppa", "SCALLOP", o) for o in ("CAD_cardio","CAD_vdh","AF","HF","stroke","CES")] + \
        [("nppa", "INTERVAL", o) for o in ("CAD_cardio","CAD_vdh","AF","HF","stroke","CES")] + \
        [("nppa", "NPPAeQTL", o) for o in ("CAD_cardio","CAD_vdh","AF","HF","stroke","CES")] + \
        [("npr3", "NPR3eQTL", o) for o in ("CAD_cardio","CAD_vdh","MI_har","MI_finn")]

for loc, exp, out in PAIRS:
    tag = f"{exp}__x__{out}_{loc}"
    dest = os.path.join(WORK, f"smr_{tag}.smr")
    if os.path.exists(dest):
        P(f"smr {tag}: 已存在"); continue
    rc = run_rel(["--bfile", f"eur_{loc}",
                  "--gwas-summary", f"ma_{out}_{loc}.ma",
                  "--beqtl-summary", exp,
                  "--out", f"smr_{tag}",
                  "--thread-num", "2"], f"step13e_smr_{tag}.log")
    ok = os.path.exists(dest)
    P(f"smr {tag}: rc={rc}, result={'OK' if ok else 'NO OUTPUT'}")

with open(os.path.join(WORK, "step13e_output.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(LOG))
print("DONE")
