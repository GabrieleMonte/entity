#!/usr/bin/env python3
"""Validation gates for the quiet start.  module load python3/3.9.7; python3 check_quiet.py"""
import os, numpy as np, nt2

ROOT = os.path.expandvars("$SCRATCH/dpdm/qstest")
PPC, NX, DX = 64, 1000, 0.04
FLOOR = 2.06e-4 / PPC          # measured random-start floor, eps = 2.06e-4/ppc
ok = lambda b: "PASS" if b else "FAIL"

def load(tag):
    return nt2.Data(os.path.join(ROOT, tag, tag))

print("="*72)
print("GATE 1 -- t=0 loading.  Entity starts with E=0 identically and BOTH")
print("loadings co-locate the pair, so rho(0)=0 either way and Ex(0)=0 proves")
print("nothing. The discriminator is the TOTAL density n_e+n_i: that is the")
print("channel the quiet start exists to kill.")
print("="*72)
for tag in ("t0_quiet",):
    d = load(tag)
    n1 = np.asarray(d.fields["N_1"].isel(t=0).values, dtype=np.float64)
    n2 = np.asarray(d.fields["N_2"].isel(t=0).values, dtype=np.float64)
    # Rho is MASS density (1 + 1836 = 1837); Charge is the charge density
    rho = np.asarray(d.fields["Charge"].isel(t=0).values, dtype=np.float64)
    ntot = n1 + n2
    print(f"  {tag}: <n_e>={n1.mean():.6f}  dn_e/n_e = {n1.std()/n1.mean():.3e}")
    print(f"          <n_tot>={ntot.mean():.6f}  dn_tot/n_tot = {ntot.std()/ntot.mean():.3e}")
    print(f"          max|charge| = {np.abs(rho).max():.3e}   (co-location, exact)")
    print(f"    expected: random start would give dn/n ~ 1/sqrt(ppc) = {1/np.sqrt(PPC):.3f}")
    print(f"    {ok(n1.std()/n1.mean() < 1e-4)}  lattice uniformity")
    print(f"    {ok(np.abs(rho).max() < 1e-12)}  charge neutrality")

print()
print("="*72)
print("GATE 2 -- undriven noise floor, quiet vs random, same ppc0. This is the")
print("test that matters: does the ordering SURVIVE, or relax to the random")
print("floor within a few plasma periods (which is what happens if the velocity")
print("loading is left stochastic).")
print("="*72)
res = {}
for tag in ("undriven_q", "undriven_r"):
    csv = os.path.join(ROOT, tag, tag, f"{tag}_stats.csv")
    raw = np.genfromtxt(csv, delimiter=",", names=True, deletechars=" ", autostrip=True)
    t, e2 = np.asarray(raw["time"]), 0.5*np.asarray(raw["E1^2"])
    late = t > 100
    res[tag] = (t, e2, np.median(e2[late]))
    print(f"  {tag:12} plateau <eps_E> = {res[tag][2]:.3e}   "
          f"(measured random floor {FLOOR:.3e})")
r = res["undriven_q"][2] / res["undriven_r"][2]
print(f"\n  quiet/random = {r:.3e}")
print(f"  {ok(r < 0.1)}  quiet start suppresses the floor by >10x")
print(f"  {ok(abs(res['undriven_r'][2]/FLOOR - 1) < 1.0)}  random run reproduces "
      f"the known floor (sanity check on the comparison)")

print()
print("="*72)
print("GATE 3 -- no coherent seeded mode. A naive lattice (consecutive slots")
print("getting consecutive velocity quantiles) plants a v(x) ramp in every cell:")
print("a coherent mode, worse than the noise it removed. The stride permutation")
print("is what prevents it, and this is the check that it worked.")
print("="*72)
d = load("undriven_q")
Ex = np.asarray(d.fields["Ex"].isel(t=-1).values, dtype=np.float64)
k = 2*np.pi*np.fft.rfftfreq(NX, d=DX)
P = np.abs(np.fft.rfft(Ex)/NX)**2
j = int(np.argmax(P[1:])) + 1
print(f"  loudest mode: k = {k[j]:.2f} (k*lambda_D = {k[j]*np.sqrt(1e-3):.3f}),"
      f"  P = {P[j]:.3e}")
print(f"  median P over k>0 = {np.median(P[1:]):.3e}   peak/median = {P[j]/np.median(P[1:]):.1f}")
print(f"  {ok(P[j]/np.median(P[1:]) < 100)}  no single mode dominating the spectrum")

print()
print("="*72)
print("GATE 4 -- the drive still works. Quiet loading must not have broken the")
print("normalization: eps_E + dKE_e = (A0^2/8)(wp t)^2, measured C=0.12432.")
print("="*72)
csv = os.path.join(ROOT, "driven_q", "driven_q", "driven_q_stats.csv")
raw = np.genfromtxt(csv, delimiter=",", names=True, deletechars=" ", autostrip=True)
t = np.asarray(raw["time"]); A0 = 9.4868e-4
tot = 0.5*np.asarray(raw["E1^2"]) + (np.asarray(raw["T00_1"]) - raw["T00_1"][0])
m = (t > 10) & (t < 40)
C = np.median(tot[m]/(A0**2 * t[m]**2))
print(f"  C = {C:.5f}   (theory 1/8 = 0.125,  run34 measured 0.12432)")
print(f"  {ok(abs(C/0.125 - 1) < 0.05)}  drive normalization preserved")
