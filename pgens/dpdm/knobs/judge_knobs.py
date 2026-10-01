#!/usr/bin/env python3
"""
Judge the fast-case knob scan against the HHS Fig. 2 target.

Reads ONLY <run>/fast/fast_stats.csv (no nt2 needed, runs on a login node).
The scan observable is dKE_e (total electron kinetic energy change from T00_1):
above wp t ~ 300 the coherent k=0 part is <~5% of it, so dKE_e stands in for
the paper's dE_e(k != 0) at the precision that matters here (a 2x deficit).

Targets (HHS Fig. 2 digitized, see $SCRATCH/dpdm_npy/data_2511/run_fast.txt):
    dE_e(k!=0):  t=300: 1.04e-2   t=500: 1.93e-2   t=1000: 2.60e-2  t=4800: 4.12e-2
Baseline fast_s (Sep 16 run):
    dKE_e:       t=300: 5.6e-3    t=500: 1.05e-2   t=1000: 1.49e-2  t=4800: 2.97e-2
"""
import os, glob
import numpy as np

SCRATCH = os.environ["SCRATCH"]
HHS = {300: 1.04e-2, 500: 1.93e-2, 1000: 2.60e-2, 2000: 2.77e-2, 4800: 4.12e-2}
TQ = sorted(HHS)

def read(path):
    raw = np.genfromtxt(path, delimiter=",", names=True, deletechars=" ",
                        autostrip=True)
    g = lambda k: np.asarray(raw[k], dtype=np.float64)
    t = g("time")
    return dict(t=t, eps_E=0.5*g("E1^2"), dKE_e=g("T00_1")-g("T00_1")[0],
                dKE_i=g("T00_2")-g("T00_2")[0])

def report(label, s):
    t = s["t"]
    # peak of the smoothed field energy envelope (window ~ 20 wp^-1)
    w = min(41, len(t)//2*2-1); ker = np.ones(w)/w
    eE = np.convolve(s["eps_E"], ker, mode="same")
    ipk = np.argmax(eE)
    row = [f"{label:8s} t_end={t[-1]:7.0f}  epsE_pk={eE[ipk]:.2e}@{t[ipk]:4.0f}"]
    for tq in TQ:
        if tq > t[-1]:
            row.append(f"t={tq}: --")
            continue
        v = np.interp(tq, t, s["dKE_e"])
        row.append(f"t={tq}: {v:.2e} ({v/HHS[tq]:.2f}x HHS)")
    print("  ".join(row))
    return {tq: np.interp(tq, t, s["dKE_e"]) for tq in TQ if tq <= t[-1]}

print("target   " + "  ".join(f"t={tq}: {HHS[tq]:.2e} (1.00x)" for tq in TQ))
base_csv = f"{SCRATCH}/dpdm/paper_shear/fast/fast/fast_stats.csv"
if os.path.exists(base_csv):
    report("fast_s", read(base_csv))
for d in sorted(glob.glob(f"{SCRATCH}/dpdm/knobs/*/fast/fast_stats.csv")):
    tag = d.split("/knobs/")[1].split("/")[0]
    try:
        report(tag, read(d))
    except Exception as e:
        print(f"{tag:8s} unreadable ({type(e).__name__}: {e})")
