# -*- coding: utf-8 -*-
"""step18c: sparse make-besd crash - OpenMP runtime variants."""
import os, subprocess, shutil

V = r"C:\Users\xuyin\WorkBuddy\2026-09-13-15-20-56\smr_verify"
EXE = os.path.join(V, "smr-1.3.1-win.exe")

def dec(b):
    for enc in ("utf-16", "utf-8", "gbk"):
        s = (b or b"").decode(enc, errors="ignore")
        if len(s) > 20:
            return s
    return ""

def trial(tag, env_extra=None, dll_state=None):
    # dll_state: None | ('hide','libomp.dll') | ('hide','libomp140.x86_64.dll')
    hidden = []
    if dll_state:
        src = os.path.join(V, dll_state[1])
        dst = src + ".hidden"
        shutil.move(src, dst)
        hidden.append((dst, src))
    env = dict(os.environ)
    if env_extra:
        env.update(env_extra)
    try:
        r = subprocess.run([EXE, "--efile", "mini.efile", "--make-besd",
                            "--out", "m_" + tag, "--thread-num", "1"],
                           capture_output=True, timeout=300, cwd=V, env=env)
        ok = os.path.exists(os.path.join(V, "m_" + tag + ".besd"))
        print(tag, "| rc =", r.returncode, "| besd =", ok)
        t = dec(r.stdout) + " STDERR:" + dec(r.stderr)
        print("   ", t.replace("\r", "").replace("\n", " | ")[-350:])
    finally:
        for dst, src in hidden:
            shutil.move(dst, src)

trial("A_default")
trial("B_omp1", env_extra={"OMP_NUM_THREADS": "1"})
trial("C_no140", dll_state=("hide", "libomp140.x86_64.dll"))
trial("D_nolibomp", dll_state=("hide", "libomp.dll"))
