# -*- coding: utf-8 -*-
"""step18b: dense make-besd on real efiles; capture full stdout+stderr."""
import os, subprocess, glob

V = r"C:\Users\xuyin\WorkBuddy\2026-09-13-15-20-56\smr_verify"
EXE = os.path.join(V, "smr-1.3.1-win.exe")

def dec(b):
    for enc in ("utf-16", "utf-8", "gbk"):
        s = (b or b"").decode(enc, errors="ignore")
        if len(s) > 20:
            return s
    return ""

for name in ("SCALLOP",):
    r = subprocess.run([EXE, "--efile", name + ".efile", "--make-besd-dense",
                        "--out", name + "_d", "--thread-num", "1"],
                       capture_output=True, timeout=600, cwd=V)
    print(name, "rc =", r.returncode)
    out = dec(r.stdout)
    print("--- stdout ---")
    print(out[:2000])
    err = dec(r.stderr)
    print("--- stderr ---")
    print(err[:1000])
    print("--- new files ---")
    for f in sorted(glob.glob(os.path.join(V, name + "_d*"))):
        print("  ", os.path.basename(f), os.path.getsize(f))
