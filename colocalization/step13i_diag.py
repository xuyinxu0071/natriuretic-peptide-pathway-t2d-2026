# -*- coding: utf-8 -*-
"""诊断: 最小 efile 测试 + 下载 smr_Win_v1.03"""
import os, subprocess, zipfile, urllib.request

WORK = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估\coloc_analysis"
SMR = os.path.join(WORK, "bin", "smr-1.3.1-win-x86_64", "smr-1.3.1-win.exe")
OUT = []

# 1) 最小 efile
mini = os.path.join(WORK, "mini.efile")
with open(mini, "w") as f:
    f.write("ProbeID\tSNP\tChr\tBP\tA1\tA2\tFreq\tb\tse\tp\n")
    f.write("prb1\trs1001\t1\t1001\tG\tA\t0.45\t0.12\t0.02\t1e-6\n")
    f.write("prb1\trs1002\t1\t2002\tG\tA\t0.30\t0.08\t0.02\t6e-5\n")
    f.write("prb1\trs1003\t1\t3003\tT\tC\t0.25\t-0.05\t0.02\t0.012\n")
r = subprocess.run([SMR, "--efile", "mini.efile", "--make-besd", "--out", "minibesd",
                    "--thread-num", "1"], capture_output=True, text=True,
                   encoding="utf-8", errors="ignore", cwd=WORK, timeout=300)
OUT.append(f"mini efile make-besd: rc={r.returncode}, besd exists={os.path.exists(os.path.join(WORK, 'minibesd.besd'))}")
OUT.append((r.stdout or "")[-300:])
OUT.append("STDERR: " + (r.stderr or "")[-200:])

# 2) smr_Win_v1.03
url = "https://yanglab.westlake.edu.cn/software/smr/download/smr_Win_v1.03.zip"
dest = os.path.join(WORK, "bin", "smr_Win_v1.03.zip")
try:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    data = urllib.request.urlopen(req, timeout=600).read()
    open(dest, "wb").write(data)
    with zipfile.ZipFile(dest) as z:
        z.extractall(os.path.join(WORK, "bin", "smr103"))
        OUT.append(f"v1.03 下载 {len(data):,} bytes: {z.namelist()}")
except Exception as e:
    OUT.append(f"v1.03 失败: {e}")

open(os.path.join(WORK, "step13i_output.txt"), "w", encoding="utf-8").write("\n".join(OUT))
print("DONE")
