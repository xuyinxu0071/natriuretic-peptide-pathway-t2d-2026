# -*- coding: utf-8 -*-
"""step18a: official SMR make-besd crash triage in pure-ASCII workspace."""
import os, subprocess

V = r"C:\Users\xuyin\WorkBuddy\2026-09-13-15-20-56\smr_verify"
EXE = os.path.join(V, "smr-1.3.1-win.exe")

def dec(b):
    for enc in ("utf-16", "utf-8", "gbk"):
        s = (b or b"").decode(enc, errors="ignore")
        if len(s) > 20:
            return s
    return ""

def run(args, tag, outname):
    r = subprocess.run([EXE] + args, capture_output=True, timeout=300, cwd=V)
    ok = os.path.exists(os.path.join(V, outname + ".besd"))
    print(tag, "| rc =", r.returncode, "| besd =", ok)
    print("   tail:", dec(r.stdout).replace("\r", "").replace("\n", " | ")[-300:])
    return r.returncode, ok

# T3: dense format
run(["--efile", "mini.efile", "--make-besd-dense", "--out", "mini_t3", "--thread-num", "1"],
    "T3-mini-dense", "mini_t3")

# T5: LF-only line endings, dense
raw = open(os.path.join(V, "mini.efile"), "rb").read().replace(b"\r\n", b"\n")
open(os.path.join(V, "mini_lf.esd"), "wb").write(raw)
run(["--efile", "mini_lf.esd", "--make-besd", "--out", "mini_t5", "--thread-num", "1"],
    "T5-mini-LF-sparse", "mini_t5")
run(["--efile", "mini_lf.esd", "--make-besd-dense", "--out", "mini_t6", "--thread-num", "1"],
    "T6-mini-LF-dense", "mini_t6")
