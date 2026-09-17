# -*- coding: utf-8 -*-
"""诊断 v2: 纯 ASCII cwd 测试 SMR"""
import os, shutil, subprocess

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
TMP = r"C:\smr_test"
os.makedirs(TMP, exist_ok=True)
OUT = []

smrdir = os.path.join(WORK, "bin", "smr-1.3.1-win-x86_64")
for fn in os.listdir(smrdir):
    if fn.lower().endswith((".exe", ".dll")) and not fn.startswith("."):
        shutil.copy2(os.path.join(smrdir, fn), os.path.join(TMP, fn))
shutil.copy2(os.path.join(WORK, "mini.efile"), os.path.join(TMP, "mini.efile"))

SMR = os.path.join(TMP, "smr-1.3.1-win.exe")
r = subprocess.run([SMR, "--efile", "mini.efile", "--make-besd", "--out", "minibesd",
                    "--thread-num", "1"], capture_output=True, text=True,
                   encoding="utf-8", errors="ignore", cwd=TMP, timeout=300)
OUT.append(f"ASCII cwd make-besd: rc={r.returncode}, besd={os.path.exists(os.path.join(TMP, 'minibesd.besd'))}")
OUT.append((r.stdout or "")[-500:])
OUT.append("STDERR: " + (r.stderr or "")[-300:])

open(os.path.join(WORK, "step13j_output.txt"), "w", encoding="utf-8").write("\n".join(OUT))
print("DONE")
