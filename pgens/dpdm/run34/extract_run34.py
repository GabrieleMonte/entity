#!/usr/bin/env python3
"""
Reduce the run34 ADIOS2 output to .npz for plotting.  Run once; then every plot
script needs only numpy.

    pip install nt2py          # brings adios2 + xarray + dask
    python3 extract_run34.py

Reads  run34/<tag>/<tag>/{fields,particles}/*.bp  and  <tag>_stats.csv
Writes run34_npy/<tag>_{stats,fields,temps}.npz + summary.npz + meta.json
(~2.7 GB of .bp  ->  ~65 MB of npz)

ALL UNITS ARE ENTITY CODE UNITS
    omega_pe = 1, c = 1, n_e = n0, m_e = 1
    energies  in n_e m_e c^2       time in 1/omega_p       k in omega_p/c
    vth_e = sqrt(1e-3) = 0.0316    lambda_D = 0.0316 c/omega_p (k*lambda_D=1 at k=31.6)

TWO NORMALIZATION TRAPS, both measured and confirmed (DPDM_PIC_NOTES.md section 0):

  1. The linear-theory line is       eps_E + dKE_e = (A0^2/8)(omega_p t)^2
     NOT eps_E alone.  The resonantly driven Langmuir wave equipartitions, so
     the field carries only half the wave energy on cycle-average (and swings
     between 3% and 85% of it instantaneously at 2*omega_p).  Measured
     C = 0.12432 vs 1/8 = 0.125.  Plotting the paper's "purple line" against
     eps_E alone is wrong by a factor 2.

  2. The electron thermal reference is the 1V value  0.5*n_e*vth_e^2 = 5e-4,
     NOT T00_1 - 1 = 1.5e-3.  Entity carries a 3V isotropic Maxwellian and the
     paper is 1V.  dKE_e itself is safe: with no transverse forces in a 1D
     electrostatic run the perpendicular energy is constant and cancels in the
     difference.
"""
import json
import os

import numpy as np

# ----------------------------------------------------------------------------
run34_dir = "PATH/TO/run34"          # <-- the untarred directory
out_dir   = "run34_npy"
# ----------------------------------------------------------------------------

# Reading all 300 steps in one .values call makes dask open every .bp at once;
# each reserves a 16 MB ADIOS2 buffer and the reader dies with
# "FFS out of memory". Pull the time axis in slices instead.
CHUNK = 25

# v_q^D / vth_e for each run;  A0 = ratio * vth_e = ratio * sqrt(1e-3)
RATIO = {"r1e-3": 1e-3, "r3e-3": 3e-3, "r1e-2": 1e-2, "r3e-2": 3e-2,
         "r1e-1": 1e-1, "r3e-1": 3e-1, "conv500": 3e-3, "conv8000": 3e-3}

LADDER   = ["r1e-3", "r3e-3", "r1e-2", "r3e-2", "r1e-1", "r3e-1"]  # ppc0=2000, Plot 4
CONVERGE = ["conv500", "r3e-3", "conv8000"]                        # v_q/vth=3e-3, Plot 3

VTH_E          = np.sqrt(1e-3)   # 0.0316 c
EPS_THERMAL_1V = 0.5 * 1e-3      # 5e-4, the reference line for Plots 1 and 4
MI_ME          = 1836.0
T_ION          = np.sqrt(MI_ME)  # 1/omega_pi = 42.8, the fast-regime halt time


# ============================================================================
# PLOT 1  energy budget vs omega_p t          (paper Fig. 2)
# PLOT 3  Np convergence overlay
# PLOT 4  regime transition across the ladder
#     -- all three read only this, from the plain-text stats CSV. No nt2 needed.
# ============================================================================
def read_stats(tag):
    """
    -> dict of 1D arrays over the stats cadence (dt_stats = 1.0, or 0.5 for the
       two fastest runs), i.e. ~3000-6000 samples per run.

      t        omega_p t
      eps_E    0.5*<E1^2>   field energy. Only E1 is physical in 1D ES.
      dKE_e    T00_1 - T00_1[0]   electron kinetic energy change
      dKE_i    T00_2 - T00_2[0]   ion kinetic energy change
      eps_tot  eps_E + dKE_e      <-- THIS is what obeys the (A0^2/8)t^2 law
      JdotE    work done by the field on the particles
      e_ext    the drive itself, A0*cos(omega t); lets you check phase/amplitude
               against the response in the same file
    """
    path = os.path.join(run34_dir, tag, tag, f"{tag}_stats.csv")
    # deletechars=" " keeps the literal header names ("E1^2", "J.E", ...);
    # numpy's default deletechars would strip ^ and . and mangle them.
    raw = np.genfromtxt(path, delimiter=",", names=True, deletechars=" ",
                        autostrip=True)
    g = lambda k: np.asarray(raw[k], dtype=np.float64)

    T1, T2 = g("T00_1"), g("T00_2")
    eps_E  = 0.5 * g("E1^2")
    dKE_e  = T1 - T1[0]
    return dict(t=g("time"), eps_E=eps_E, dKE_e=dKE_e, dKE_i=T2 - T2[0],
                eps_tot=eps_E + dKE_e, T00_1=T1, T00_2=T2,
                JdotE=g("J.E"), e_ext=g("e_ext"))


# ============================================================================
# PLOT 2  E(k,t) and n(k,t) spectra           (paper Fig. 8, lower panels)
# PLOT 6  real-space snapshots                (paper Fig. 8, upper panels)
# ============================================================================
def read_fields(data):
    """
    -> dict with, on a (t=300, x=1000) grid and its rfft (k=501):

      t, x, k        axes.  k = 2*pi*rfftfreq(nx, dx), so k*lambda_D = 1 at k=31.6
      w_k            Parseval weights: sum(w_k * P) == mean(a^2) exactly.
                     1 at k=0 and at Nyquist, 2 elsewhere (rfft folds +/-k).
                     Use this to split k=0 from k!=0 without guessing factors:
                        P_k0   = w_k[0]  * PE[:, 0]
                        P_knz  = (w_k[1:] * PE[:, 1:]).sum(axis=1)
      PE, Pn1, Pn2   |rfft(.)/nx|^2 for Ex, N_1, N_2      -- PLOT 2
      Ex, N_1, N_2, Rho, Jx   real-space fields           -- PLOT 6

    PLOT 2 is the mechanism: the k!=0 power must start rising exactly where
    eps_tot departs the linear line in PLOT 1. Modulational daughters appear in
    PE; the ion-acoustic branch shows up in Pn2 (ions) at low k.

    float32 throughout -- these are for plotting, and it halves the download.
    """
    f  = data.fields
    t  = np.asarray(f["t"].values, dtype=np.float64)
    x  = np.asarray(f["x"].values, dtype=np.float64)
    nx = x.size
    dx = float(x[1] - x[0])
    k  = 2.0 * np.pi * np.fft.rfftfreq(nx, d=dx)

    w_k = np.full(k.size, 2.0)
    w_k[0] = 1.0
    if nx % 2 == 0:                 # Nyquist bin is unpaired for even nx
        w_k[-1] = 1.0

    out = {"t": t, "x": x, "k": k, "w_k": w_k}
    spec_of = {"PE": "Ex", "Pn1": "N_1", "Pn2": "N_2"}
    for var in ("Ex", "N_1", "N_2", "Rho", "Jx"):
        # chunked read: see CHUNK above
        parts = [np.asarray(f[var].isel(t=slice(i, i + CHUNK)).values,
                            dtype=np.float64)
                 for i in range(0, t.size, CHUNK)]
        a = np.concatenate(parts, axis=0)
        out[var] = a.astype(np.float32)
        for pname, src in spec_of.items():
            if src == var:
                out[pname] = (np.abs(np.fft.rfft(a, axis=-1) / nx) ** 2
                              ).astype(np.float32)
    return out


# ============================================================================
# PLOT 5  T_e/T_i vs time                     (paper Fig. 7)
# ============================================================================
def read_temperatures(data, n_t):
    """
    -> t, Te, Ti, u_e, u_i  (one value per particle dump)

    PARALLEL temperature only, from the ux variance, weighted by the
    macroparticle weight w. Entity carries an isotropic 3V Maxwellian while
    SHARP is 1V; the transverse degrees of freedom feel no force in a 1D
    electrostatic run, so including them would dilute T_e by 3x and destroy the
    comparison (DPDM_PIC_NOTES.md section 6.5).

    kT_par/(m_e c^2) = mass * <w (ux - <ux>_w)^2>_w.  Non-relativistic here
    (vth_e = 0.03 c) so u = gamma*v ~= v and no gamma correction is needed.

    The paper's claim: fast-growth runs reach T_e/T_i ~ O(30) with electrons
    O(10^2) times hotter than they started.
    """
    p = data.particles
    t, Te, Ti, u_e, u_i = [], [], [], [], []
    for i in range(n_t):
        snap = p.isel(t=i)
        mom = {}
        for sp in (1, 2):
            df = snap.sel(sp=sp).load(cols=["ux", "w"])
            u = df["ux"].to_numpy(dtype=np.float64)
            w = df["w"].to_numpy(dtype=np.float64)
            wsum = w.sum()
            mu = (w * u).sum() / wsum
            mom[sp] = (mu, (w * (u - mu) ** 2).sum() / wsum)
        u_e.append(mom[1][0]);      u_i.append(mom[2][0])
        Te.append(mom[1][1] * 1.0)  # m_e = 1
        Ti.append(mom[2][1] * MI_ME)
        t.append(float(np.atleast_1d(np.asarray(snap.t)).ravel()[0]))
    return dict(t=np.array(t), Te=np.array(Te), Ti=np.array(Ti),
                u_e=np.array(u_e), u_i=np.array(u_i))


def toml_value(tag, key):
    line = next(l for l in open(os.path.join(run34_dir, tag, f"{tag}.toml"))
                if l.strip().startswith(key))
    return float(line.split("=")[1])


def main():
    import nt2

    os.makedirs(out_dir, exist_ok=True)
    meta = {}

    for tag, ratio in RATIO.items():
        print(f"[{tag}]", end=" ", flush=True)

        st = read_stats(tag)
        np.savez_compressed(os.path.join(out_dir, f"{tag}_stats.npz"), **st)
        print(f"stats({st['t'].size})", end=" ", flush=True)

        data = nt2.Data(os.path.join(run34_dir, tag, tag))

        fl = read_fields(data)
        np.savez_compressed(os.path.join(out_dir, f"{tag}_fields.npz"), **fl)
        print(f"fields({fl['t'].size})", end=" ", flush=True)

        tp = read_temperatures(data, fl["t"].size)
        np.savez_compressed(os.path.join(out_dir, f"{tag}_temps.npz"), **tp)
        print(f"temps({tp['t'].size})", flush=True)

        A0   = ratio * VTH_E
        ppc0 = toml_value(tag, "ppc0")
        meta[tag] = dict(
            vq_over_vthe   = ratio,
            A0             = A0,
            ppc0           = ppc0,
            runtime        = toml_value(tag, "runtime"),
            # eps_driven reaches the 1V electron thermal energy at this time
            t_sat          = 0.0632 / A0,
            # measured shot-noise floor; the signal must sit well above this
            eps_noise      = 2.06e-4 / ppc0,
            eps_thermal_1V = EPS_THERMAL_1V,
            # peak field energy, normalized -- the PLOT 4 ordinate.
            # slow regime saturates below 1, fast regime overshoots above 1.
            eps_E_max_norm   = float(st["eps_E"].max()   / EPS_THERMAL_1V),
            eps_tot_max_norm = float(st["eps_tot"].max() / EPS_THERMAL_1V),
            eps_tot_end_norm = float(st["eps_tot"][-1]   / EPS_THERMAL_1V),
        )

    # ---- PLOT 4: one row per ladder amplitude, ready to scatter ------------
    np.savez_compressed(
        os.path.join(out_dir, "summary.npz"),
        tags             = np.array(LADDER),
        vq_over_vthe     = np.array([meta[t]["vq_over_vthe"]     for t in LADDER]),
        A0               = np.array([meta[t]["A0"]               for t in LADDER]),
        eps_E_max_norm   = np.array([meta[t]["eps_E_max_norm"]   for t in LADDER]),
        eps_tot_max_norm = np.array([meta[t]["eps_tot_max_norm"] for t in LADDER]),
        eps_tot_end_norm = np.array([meta[t]["eps_tot_end_norm"] for t in LADDER]),
        # predicted slow/fast boundary, (m_e/m_p)^(1/2)/2
        boundary         = np.array([0.5 / np.sqrt(MI_ME)]),
    )

    meta["_groups"] = {"ladder": LADDER, "converge": CONVERGE}
    meta["_constants"] = dict(
        vth_e=float(VTH_E), eps_thermal_1V=EPS_THERMAL_1V, mi_me=MI_ME,
        t_ion=float(T_ION), lambda_D=float(VTH_E), k_lambdaD_1_at=1.0 / VTH_E)
    meta["_linear_law"] = "eps_E + dKE_e = (A0^2/8)*(omega_p t)^2  -- not eps_E alone"
    with open(os.path.join(out_dir, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)

    print(f"\nwrote {out_dir}/")


if __name__ == "__main__":
    main()
