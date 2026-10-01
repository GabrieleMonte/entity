#!/usr/bin/env python3
"""Collision-rate hierarchy for the DPDM resonance, vs HHS arXiv:2510.13956.

Answers one question: how many electron-electron collision times elapse during
(a) the HHS simulations and (b) the physical resonance?

Calibration check: this reproduces the paper's own quoted tau_ei ~ 1e4 s at
recombination (main.tex:215) to within 25%, which fixes the normalization.

Usage:  python3 cosmo_rates.py
"""
import numpy as np

# --- cosmology -------------------------------------------------------------
h, Om, Orad = 0.674, 0.315, 9.2e-5
H0  = 2.13e-18 * h          # s^-1
nb0 = 2.503e-7              # cm^-3, Omega_b h^2 = 0.0224

def Hub(z):  return H0 * np.sqrt(Orad * (1 + z) ** 4 + Om * (1 + z) ** 3)
def Tgam(z): return 2.35e-4 * (1 + z)                       # eV

def Xe(z):
    """Crude ionization history: 1 before recombination, 2e-4 after."""
    if z > 1400: return 1.0
    if z > 700:  return 2e-4 + 0.9998 / (1 + np.exp((1100 - z) / 60.0))
    return 2e-4

# --- plasma / collision rates ---------------------------------------------
LNL = 20.0                  # Coulomb log, per paper main.tex:207
EPS = 1e-8                  # kinetic mixing; resonance lasts eps/H

def rates(z):
    ne = Xe(z) * nb0 * (1 + z) ** 3
    # electrons track photons until thermal decoupling (z ~ 150), then adiabatic
    Te = Tgam(z) if z > 150 else Tgam(150) * ((1 + z) / 151.0) ** 2
    wp   = 5.64e4 * np.sqrt(ne)                 # rad/s
    nuei = 2.91e-6 * ne * LNL * Te ** -1.5      # e-i MOMENTUM exchange, NRL
    nuee = nuei / np.sqrt(2)                    # like-particle, Z = 1
    tau_ei = 918.0 / nuei                       # e-i ENERGY exchange (HHS Eq. for tau_ei)
    return dict(ne=ne, Te=Te, wp=wp, nuee=nuee, nuei=nuei, tau_ei=tau_ei,
                H=Hub(z), dt_res=EPS / Hub(z))

HHS_LONGEST = 8e4           # omega_p t, longest HHS run (main.tex:203)

if __name__ == "__main__":
    hdr = (f"{'z':>8} {'n_e[cm^-3]':>11} {'T_e[eV]':>9} {'wp[1/s]':>10} "
           f"{'nu_ee/wp':>10} {'wp*tau_ee':>10} {'N_coll_res':>11} {'N_coll_sim':>11}")
    print(hdr); print("-" * len(hdr))
    for z in [1e6, 1e5, 3e4, 1e4, 3e3, 1400, 1100, 700, 300, 100, 30]:
        r = rates(z)
        n_res = r["nuee"] * r["dt_res"]                 # collision times per resonance
        n_sim = HHS_LONGEST * r["nuee"] / r["wp"]       # collision times per HHS run
        print(f"{z:8.0f} {r['ne']:11.3e} {r['Te']:9.3e} {r['wp']:10.3e} "
              f"{r['nuee']/r['wp']:10.3e} {r['wp']/r['nuee']:10.3e} "
              f"{n_res:11.3e} {n_sim:11.3e}")

    r = rates(1100)
    print(f"\ncalibration vs paper: tau_ei(z=1100) = {r['tau_ei']:.2e} s "
          f"(paper main.tex:215 quotes 'about 1e4 s')")
    print(f"tau_ei / tau_ee = {r['tau_ei'] * r['nuee']:.0f}")
    r = rates(1e4)
    print(f"\nnu_ee/wp is epoch-INDEPENDENT before recombination "
          f"(nu ~ n^1/2 T^-3/2, n ~ a^-3, T ~ a^-1):  {r['nuee']/r['wp']:.2e}")
    print(f"equivalently n_e * lambda_D^3 ~ 3e6 (weakly coupled)")
    print(f"\nIn the saturated state T_e/T_i ~ 50, the ion self-collision rate is")
    print(f"  nu_ii/nu_ee = (m_e/m_i)^1/2 (T_e/T_i)^3/2 = "
          f"{(1/1836.)**0.5 * 52.8**1.5:.1f}  -- ions are the MOST collisional species")
