#!/usr/bin/env python3
"""How far from Maxwellian is the post-burst distribution?

HHS arXiv:2510.13956 main.tex:203 shuts off resonant conversion with the clause
"...the damping rate on electrons would be polynomial suppressed IF THE ELECTRONS
WERE THERMALIZED WITH THEMSELVES". This script tests that clause against our
own runs. It cannot hold in their setup: the paper is 1+1D (main.tex:133,485),
and electron-electron collisions do not exist in 1V.

Input: the <tag>_phase.npz written by pgens/dpdm/paper/extract_paper.py
       (keys: t, u_edges_e, u_edges_i, H_e, H_i -- H is [frame, x, u]).

Usage:  python3 fdist_diag.py $SCRATCH/dpdm_npy/fast_s_phase.npz [n_frames]
"""
import sys
import numpy as np

tz = np.trapz if not hasattr(np, "trapezoid") else np.trapezoid


def shape_stats(counts, uc):
    """Moments and Maxwellian comparison for a 1D velocity histogram."""
    f = counts / tz(counts, uc)
    um = tz(f * uc, uc)
    var = tz(f * (uc - um) ** 2, uc)
    vth = np.sqrt(var)
    w = (uc - um) / vth
    g = np.exp(-w ** 2 / 2.0)
    g /= tz(g, uc)
    kurt = tz(f * (uc - um) ** 4, uc) / var ** 2 - 3.0     # 0 for a Maxwellian
    core = np.abs(w) < 1.0
    core_dev = tz(np.abs(f - g)[core], uc[core]) / tz(g[core], uc[core])
    return f, g, w, vth, kurt, core_dev


def report(label, H, u_edges, frames, extra=()):
    uc = 0.5 * (u_edges[1:] + u_edges[:-1])
    counts = H[frames].sum(axis=(0, 1)).astype(float)   # sum over frames and x
    f, g, w, vth, kurt, core_dev = shape_stats(counts, uc)
    print(f"\n=== {label} ===")
    print(f"  vth = {vth:.4e}   excess kurtosis = {kurt:+.3f} (Maxwellian: 0)")
    print(f"  integrated |f - Maxwellian| over |v| < vth : {core_dev * 100:.1f}%")
    print(f"  histogram spans |v| <= {abs(u_edges[-1]) / vth:.1f} vth")
    print(f"\n  {'v/vth':>7} {'counts':>9} {'stat err':>9} {'f/Maxwellian':>13}")
    for ww in (0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5):
        if ww > abs(w).max():
            continue
        j = int(np.argmin(np.abs(w - ww)))
        n = counts[j]
        err = 100.0 / np.sqrt(max(n, 1.0))
        ratio = f[j] / max(g[j], 1e-300)
        print(f"  {w[j]:7.2f} {n:9.0f} {err:8.1f}% {ratio:13.3e}")
    for name, v in extra:
        if abs(v) / vth > abs(w).max():
            print(f"\n  !! {name} sits at {v / vth:.2f} vth -- OUTSIDE the histogram "
                  f"range, cannot measure. Widen the extractor.")
        else:
            j = int(np.argmin(np.abs(w - v / vth)))
            print(f"\n  {name} at {v / vth:.2f} vth: f/Maxwellian = "
                  f"{f[j] / max(g[j], 1e-300):.3e} (counts {counts[j]:.0f})")


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "fast_s_phase.npz"
    nfr = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    p = np.load(path)
    t = p["t"]
    frames = slice(-nfr, None)
    print(f"{path}: averaging {nfr} frames, omega_p t = "
          f"{t[frames][0]:.0f} to {t[frames][-1]:.0f}")

    report("electrons", p["H_e"], p["u_edges_e"], frames)

    # ion-acoustic phase speed: electron Landau damping of the IAW samples
    # df_i/dv at v = c_s, and the exp(-T_e/2T_i) suppression HHS quote is a
    # MAXWELLIAN result -- so measuring f_i there is the point.
    ui = 0.5 * (p["u_edges_i"][1:] + p["u_edges_i"][:-1])
    ci = p["H_i"][frames].sum(axis=(0, 1)).astype(float)
    _, _, _, vthi, _, _ = shape_stats(ci, ui)
    Te_guess = 0.0744          # code units, fast case at omega_p t = 1e4
    cs = np.sqrt(Te_guess / 1836.0)
    report("ions", p["H_i"], p["u_edges_i"], frames, extra=[("c_s", cs)])
    print(f"\n  c_s / vth_i = {cs / vthi:.2f}  ->  Maxwellian ion-Landau "
          f"suppression exp(-(c_s/vth_i)^2/2) = {np.exp(-(cs / vthi) ** 2 / 2):.2e}")
