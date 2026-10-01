#!/usr/bin/env python3
"""
Energy decomposition over the FIRST N particle dumps only.

WHY THIS EXISTS.  extract_paper.py's open_subset() picks n_keep dumps EVENLY
SPACED over the run, and its whole speed argument rests on n_keep being small:
nt2.Data() pays a cost proportional to how many .bp are in the series it opens.
Asking for ENERGY_FRAMES=500 on a 500-dump run keeps every dump, which is the
full series, which is the thing the symlink trick exists to avoid.  Job 3468041
did exactly that and hit the 2 h wall without finishing one run.

What is actually needed for the paper comparison is not uniform coverage -- it
is RESOLUTION BELOW wp t ~ 2000, where the paper's k != 0 burst rises five
decades between wp t = 130 and 600 and then plateaus.  The LZ runs dump every
40 time units, so the first 50 dumps cover 0..1960 at dt = 40: ~12 points inside
the burst, against the ONE point the 60-frame product has there (dt = 320).

Output: $OUT_DIR/<tag>_energy_early.npz, same keys as <tag>_energy.npz.
It is a SEPARATE product and does not touch the full-span 60-frame file.
Baselines (d_ran_rel_* etc.) are differences from the first dump, which is
t = 0 here, so they mean the same thing as in the full-span file.
"""
import os, sys, time, tempfile, shutil
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import extract_paper as E

NKEEP = int(os.environ.get("NKEEP", "50"))
TAGS  = [t for t in os.environ.get("TAGS", "lz1p0n,lz1p0_s").split(",") if t]
OUT   = E.out_dir


def first_n(tag, kind, n_keep):
    """open_subset(), but the first n_keep dumps instead of evenly spaced."""
    src   = os.path.join(E.rundir(tag), kind)
    dumps = sorted(d for d in os.listdir(src) if d.endswith(".bp"))
    keep  = dumps[:n_keep]
    tmp   = tempfile.mkdtemp(prefix=f"early_{tag}_{kind}_")
    leaf  = os.path.join(tmp, tag, kind)
    os.makedirs(leaf)
    for d in keep:
        os.symlink(os.path.join(src, d), os.path.join(leaf, d))
    csv = os.path.join(E.rundir(tag), f"{tag}_stats.csv")
    if os.path.exists(csv):
        os.symlink(csv, os.path.join(tmp, tag, f"{tag}_stats.csv"))
    import nt2
    return nt2.Data(os.path.join(tmp, tag)), tmp, keep


E.open_subset = first_n

os.makedirs(OUT, exist_ok=True)
for tag in TAGS:
    t0 = time.time()
    print(f"[{tag}] first {NKEEP} dumps ...", flush=True)
    en = E.read_energy(tag, NKEEP)
    path = os.path.join(OUT, f"{tag}_energy_early.npz")
    np.savez(path, **en)
    print(f"[{tag}] t = {en['t'][0]:.1f} .. {en['t'][-1]:.1f}  "
          f"n = {en['t'].size}  ->  {path}   ({time.time()-t0:.0f} s)", flush=True)
