#!/usr/bin/env python3
"""
Reduce the paper-faithful production runs to .npz for plotting.

This is the companion of run34/extract_run34.py for the runs that run34 CANNOT
produce.  run34 is a random-start survey at a uniform ppc0 = 2000 out to
omega_p t = 3000.  The six runs here are quiet-start, at the paper's own particle
counts, and include the swept-frequency family.  Different initial conditions,
different N_pc, different physics.

    tag          family         v_q^D/vth   N_pc     omega_p t   what it is
    -----------  -------------  ----------  -------  ----------  --------------------
    fast         resonance      3e-2        2e2      1e4         paper Fig. 2, upper
    slow         resonance      1e-3        2e5      1e4         paper Fig. 2, lower
    lz1p0        landau_zener   1.0         2.7e3    2e4         paper Fig. 3
    heat         undriven       0           2.7e3    2e4         relaxation + heating gate
    undriven_q   undriven       0           64       400         quiet   } the A/B
    undriven_r   undriven       0           64       400         random  } pair

WHERE TO RUN THIS
    On Lonestar6, not on your laptop: it turns ~1.9 GB of .bp into ~50 MB of
    .npz, so reducing first means downloading 50 MB instead of 1.9 GB.  nt2py is
    already installed there (1.5.3).  Do NOT run it on a login node: opening the
    dump series is the memory-hungry step, not reading it.  Each .bp reserves a
    16 MB ADIOS2 buffer at open, and `fast` has 1000 field dumps against run34's
    300 -- an nt2.Data() call on `fast` from a login node dies with "FFS out of
    memory" before returning.  A development node has the memory for it:

        sbatch -p development -N 1 -n 1 -t 01:00:00 -A PHY23028 \
               --wrap "module load gcc/13.2.0; python3 extract_paper.py"

    Then copy paper_npy/ home.  It also runs unchanged on a laptop against a
    downloaded tree; set paper_dir below.

INPUT LAYOUT -- <paper_dir>/<tag>/<tag>/{fields,particles}/*.bp
                <paper_dir>/<tag>/<tag>/<tag>_stats.csv
                <paper_dir>/<tag>/<tag>.toml
    which is what you get by staging the science and leaving the checkpoints
    behind (paper/slow is 75 GB of .ckpt beside 483 MB of science):

        DEST=$SCRATCH/dpdm_export
        for src in paper/fast paper/slow lz/lz1p0 lz/heat \
                   qstest/undriven_q qstest/undriven_r; do
          tag=$(basename $src); mkdir -p $DEST/$tag
          cp -a $SCRATCH/dpdm/$src/$tag $SCRATCH/dpdm/$src/$tag.toml $DEST/$tag/
        done

OUTPUT -- paper_npy/
    <tag>_stats.npz    scalar time series          Plots A, B, C, D
    <tag>_fields.npz   real space + k space        Plots C, E
    <tag>_temps.npz    T_e, T_i per particle dump  Plots B, D
    <tag>_phase.npz    (x, ux) histograms          Plot F  -- the "display" one
    summary.npz        one row per run             Plots A, B
    meta.json          every constant, per run

ALL UNITS ARE ENTITY CODE UNITS
    omega_pe = 1, c = 1, n_e = n0, m_e = 1
    energies in n_e m_e c^2     time in 1/omega_p      k in omega_p/c
    vth_e = sqrt(1e-3) = 0.0316    lambda_D = 0.0316    k*lambda_D = 1 at k = 31.6

THE SAME TWO NORMALIZATION TRAPS AS run34 (dpdm_doc/main.tex, sec. "Deriving the
physical quantities"), restated because they change conclusions by factors:

  1. The linear-theory line is    eps_E + dKE_e = (A0^2/8)(omega_p t)^2
     NOT eps_E alone.  A resonantly driven Langmuir wave equipartitions, so the
     field holds half the wave energy on cycle average and swings between 3% and
     85% of it instantaneously at 2*omega_p.  Measured C = 0.1264 (fast),
     0.1237 (slow) against 1/8 = 0.125.

  2. The electron thermal reference is the 1V value 0.5*n_e*vth_e^2 = 5e-4,
     NOT T00_1 - 1 = 1.5e-3, which is the 3V rest-frame value.  dKE_e itself is
     safe: with no transverse forces the perpendicular energy is constant and
     cancels in the difference.  (These runs set one_v = true, so uy = uz = 0
     identically -- but keep the 5e-4 reference anyway, for comparability with
     run34 and with the paper.)
"""
import json
import os
import warnings

import numpy as np

# nt2 assigns into a DataFrame slice when loading particle columns, which
# raises one pandas SettingWithCopyWarning per species per dump -- thousands of
# lines of noise around the progress output. Harmless; silence it.
try:
    from pandas.errors import SettingWithCopyWarning
    warnings.simplefilter("ignore", SettingWithCopyWarning)
except ImportError:
    pass

# ----------------------------------------------------------------------------
# Both can be overridden from the environment, which is how extract_paper.sh
# drives this on the cluster without editing the file.
paper_dir = os.environ.get("PAPER_DIR", "PATH/TO/paper_runs")
out_dir   = os.environ.get("OUT_DIR",   "paper_npy")
# FORCE=1 re-extracts runs whose .npz already exist. Without it they are skipped,
# so a job that hit its walltime can be resubmitted and will only do what is left.
force     = os.environ.get("FORCE", "") == "1"
# ----------------------------------------------------------------------------

# Reading every dump in one .values call makes dask open all the .bp files at
# once; each reserves a 16 MB ADIOS2 buffer and the reader dies with "FFS out of
# memory".  Pull the time axis in slices.  (Same value as run34.)
CHUNK = 25

# Phase-space frames to keep per run.  A 240x240 float32 histogram is 230 kB, so
# 60 frames x 2 species is ~28 MB -- enough for a smooth movie, small enough to
# download.  Every particle dump is still used for the T_e/T_i time series.
PHASE_FRAMES = 60
PHASE_BINS   = 240

VTH_E          = np.sqrt(1e-3)   # 0.0316 c
EPS_THERMAL_1V = 0.5 * 1e-3      # 5e-4  -- the reference line, 1V not 3V
MI_ME          = 1836.0
T_ION          = np.sqrt(MI_ME)  # 1/omega_pi = 42.8

# Measured floors quoted in dpdm_doc/main.tex; used as reference lines, not
# recomputed here.
EPS_NOISE_COEFF   = 2.06e-4      # eps_noise = EPS_NOISE_COEFF / N_pc
QUIET_FLOOR_T0    = 1.31e-17     # undriven_q at t = 0
RANDOM_FLOOR_T0   = 3.03e-6      # undriven_r at t = 0
HEAT_PLATEAU      = 4.8e-8       # heat, reached by omega_p t ~ 5000
HEAT_RANDOM_FLOOR = 7.6e-8       # random start at the same N_pc

# Ordered cheapest first. The particle loop dominates and costs ~1 s per dump per
# species, so `fast` (1000 dumps) is most of the total: if a job runs out of time
# it will have finished the other five, and resubmitting picks up only `fast`.
RUNS = {
    "heat":       dict(family="undriven",     ratio=0.0,  loading="quiet"),
    "undriven_q": dict(family="undriven",     ratio=0.0,  loading="quiet"),
    "undriven_r": dict(family="undriven",     ratio=0.0,  loading="random"),
    "slow":       dict(family="resonance",    ratio=1e-3, loading="quiet"),
    "lz1p0":      dict(family="landau_zener", ratio=1.0,  loading="quiet"),
    "fast":       dict(family="resonance",    ratio=3e-2, loading="quiet"),
}

GROUPS = {
    "paper_fig2":      ["fast", "slow"],          # Plot A
    "landau_zener":    ["lz1p0"],                 # Plot B
    "quiet_vs_random": ["undriven_q", "undriven_r"],  # Plot C -- must be plotted together
    "relaxation":      ["heat"],                  # Plot D
}


def rundir(tag):
    return os.path.join(paper_dir, tag, tag)


# ============================================================================
# PLOT A  energy budget, paper Fig. 2          (fast, slow)
# PLOT B  Landau-Zener transition              (lz1p0)
# PLOT C  quiet vs random noise floor          (undriven_q, undriven_r)
# PLOT D  quiet-start relaxation               (heat)
#     -- all four read only this, from the plain-text stats CSV.  No nt2 needed,
#        so this half of the script runs anywhere in seconds.
# ============================================================================
def read_stats(tag):
    """
    -> dict of 1D arrays over the stats cadence.

      t        omega_p t
      eps_E    0.5*<E1^2>.  Only E1 is physical in a 1D electrostatic run;
               E2^2 and E3^2 are zero by construction.
      dKE_e    T00_1 - T00_1[0]      electron kinetic energy change
      dKE_i    T00_2 - T00_2[0]      ion kinetic energy change
      eps_tot  eps_E + dKE_e         <-- THIS obeys the (A0^2/8)t^2 law
      Te_rel   T00_1 / T00_1[0]      relative drift; the heating gate reads this
      JdotE    work done by the field on the particles
      e_ext    the drive itself, present only if the deck asked for it.  For the
               swept runs this is the cleanest way to see omega(t) advance --
               it is the actual field the pusher used, not a reconstruction.
    """
    path = os.path.join(rundir(tag), f"{tag}_stats.csv")
    # deletechars=" " keeps the literal header names ("E1^2", "J.E", ...);
    # numpy's default deletechars strips ^ and . and mangles them.
    raw = np.genfromtxt(path, delimiter=",", names=True, deletechars=" ",
                        autostrip=True)
    have = set(raw.dtype.names)
    g = lambda k: np.asarray(raw[k], dtype=np.float64)

    T1, T2 = g("T00_1"), g("T00_2")
    eps_E  = 0.5 * g("E1^2")
    dKE_e  = T1 - T1[0]
    out = dict(t=g("time"), eps_E=eps_E, dKE_e=dKE_e, dKE_i=T2 - T2[0],
               eps_tot=eps_E + dKE_e, T00_1=T1, T00_2=T2,
               Te_rel=T1 / T1[0], Ti_rel=T2 / T2[0])
    for opt in ("J.E", "e_ext"):
        if opt in have:
            out[{"J.E": "JdotE", "e_ext": "e_ext"}[opt]] = g(opt)
    return out


# ============================================================================
# PLOT C  density noise floor, real space      (undriven_q vs undriven_r)
# PLOT E  E(k,t) spectra at paper parameters   (fast, slow, lz1p0)
# ============================================================================
def read_fields(data):
    """
    -> dict on a (t, x=1000) grid plus its rfft (k=501):

      t, x, k        axes.  k = 2*pi*rfftfreq(nx, dx), so k*lambda_D = 1 at 31.6
      w_k            Parseval weights: sum(w_k * P) == mean(a^2) exactly.
                     1 at k=0 and at Nyquist, 2 elsewhere (rfft folds +/-k).
                     Splits k=0 from k!=0 without guessing factors:
                        P_k0  = w_k[0]  * PE[:, 0]
                        P_knz = (w_k[1:] * PE[:, 1:]).sum(axis=1)
      PE, Pn1, Pn2   |rfft(.)/nx|^2 of Ex, N_1, N_2
      dn_over_n      (N_1 + N_2)/<N_1 + N_2> - 1, the TOTAL density fluctuation.
                     This is the channel the quiet start exists to kill, and the
                     only honest discriminator between the two loaders: the
                     CHARGE density is zero for both, because both co-locate the
                     species, and Ex(t=0) = 0 for both because Entity
                     initialises E to zero identically.
      Ex, N_1, N_2, Rho, ...  whatever real-space fields the deck wrote.
                     NOTE Rho is MASS density (1 + 1836 = 1837 in code units);
                     charge density is the separate Charge field, which only the
                     validation decks requested.

    float32 for the arrays -- these are for plotting, and it halves the download.
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
    spec_of = {"Ex": "PE", "N_1": "Pn1", "N_2": "Pn2"}
    present = set(f.data_vars)
    dense   = {}

    # Field variables differ per deck: the production decks wrote E/J/Rho/N/T00,
    # the validation decks wrote Charge instead of J.  Take what is there.
    for var in ("Ex", "N_1", "N_2", "Rho", "Charge", "Jx"):
        if var not in present:
            continue
        parts = [np.asarray(f[var].isel(t=slice(i, i + CHUNK)).values,
                            dtype=np.float64)
                 for i in range(0, t.size, CHUNK)]           # see CHUNK above
        a = np.concatenate(parts, axis=0)
        out[var] = a.astype(np.float32)
        if var in spec_of:
            out[spec_of[var]] = (np.abs(np.fft.rfft(a, axis=-1) / nx) ** 2
                                 ).astype(np.float32)
        if var in ("N_1", "N_2"):
            dense[var] = a

    if "N_1" in dense and "N_2" in dense:
        ntot = dense["N_1"] + dense["N_2"]
        out["dn_over_n"] = (ntot / ntot.mean(axis=-1, keepdims=True)
                            - 1.0).astype(np.float32)
    return out


# ============================================================================
# PLOT B  T_e(t) through the sweep             (lz1p0)
# PLOT D  numerical heating over 1e6 steps     (heat -- no particle dumps, so
#                                               this falls back to the CSV)
# PLOT F  phase space (x, ux)                  (fast, slow, lz1p0)
# ============================================================================
def read_particles(data, x_max):
    """
    One pass over the particle dumps, producing two things.

    moments (every dump):  t, Te, Ti, u_e, u_i
        PARALLEL temperature only, from the ux variance weighted by the
        macroparticle weight w:
            kT_par/(m c^2) = mass * <w (ux - <ux>_w)^2>_w
        Non-relativistic here (vth_e = 0.03 c) so u = gamma*v ~= v and no gamma
        correction is needed.  These runs are 1V (uy = uz = 0), so this is the
        whole temperature, not a projection.

    phase (PHASE_FRAMES evenly spaced dumps):  t_phase, x_edges, u_edges_{e,i},
                                               H_e, H_i
        2D histograms of (x, ux), weighted by w.  The u range is fixed across
        time from a scan of the first/middle/last frames, so frames can be
        animated without the axes jumping.  Histogramming rather than storing
        particles keeps this at ~28 MB instead of gigabytes.
    """
    p = data.particles
    times = np.asarray(p.times, dtype=np.float64)
    n_t   = times.size

    frames = np.unique(np.linspace(0, n_t - 1, min(PHASE_FRAMES, n_t)).astype(int))

    # --- fix the velocity range once, from a scan, so the movie axes hold still
    urange = {1: 0.0, 2: 0.0}
    for i in (0, n_t // 2, n_t - 1):
        for sp in (1, 2):
            u = data.particles.isel(t=i).sel(sp=sp).load(cols=["ux"])["ux"]
            urange[sp] = max(urange[sp],
                             float(np.percentile(np.abs(np.asarray(u)), 99.9)))
    for sp in (1, 2):
        urange[sp] = urange[sp] * 1.15 or 1e-3

    x_edges = np.linspace(0.0, x_max, PHASE_BINS + 1)
    u_edges = {sp: np.linspace(-urange[sp], urange[sp], PHASE_BINS + 1)
               for sp in (1, 2)}

    mom = {1: ([], []), 2: ([], [])}
    H   = {1: [], 2: []}
    t_out, t_phase = [], []

    for i in range(n_t):
        snap = p.isel(t=i)
        want_frame = i in frames
        for sp in (1, 2):
            df = snap.sel(sp=sp).load(
                cols=["ux", "w", "x"] if want_frame else ["ux", "w"])
            u = df["ux"].to_numpy(dtype=np.float64)
            w = df["w"].to_numpy(dtype=np.float64)
            wsum = w.sum()
            mu   = (w * u).sum() / wsum
            var  = (w * (u - mu) ** 2).sum() / wsum
            mom[sp][0].append(mu)
            mom[sp][1].append(var)
            if want_frame:
                h, _, _ = np.histogram2d(
                    df["x"].to_numpy(dtype=np.float64), u,
                    bins=[x_edges, u_edges[sp]], weights=w)
                H[sp].append(h.astype(np.float32))
        t_out.append(float(times[i]))
        if want_frame:
            t_phase.append(float(times[i]))

    moments = dict(t=np.array(t_out),
                   Te=np.array(mom[1][1]) * 1.0,        # m_e = 1
                   Ti=np.array(mom[2][1]) * MI_ME,
                   u_e=np.array(mom[1][0]), u_i=np.array(mom[2][0]))
    phase = dict(t=np.array(t_phase), x_edges=x_edges,
                 u_edges_e=u_edges[1], u_edges_i=u_edges[2],
                 H_e=np.stack(H[1]), H_i=np.stack(H[2]))
    return moments, phase


def deck_value(tag, key, default=None):
    """Read a scalar out of <tag>.toml.  The decks are the record of what ran."""
    path = os.path.join(paper_dir, tag, f"{tag}.toml")
    for line in open(path):
        s = line.split("#")[0].strip()
        if s.startswith(key) and "=" in s and s.split("=")[0].strip() == key:
            v = s.split("=", 1)[1].strip().strip('"')
            try:
                return float(v)
            except ValueError:
                return v
    return default


def main():
    import nt2

    os.makedirs(out_dir, exist_ok=True)
    meta = {}

    for tag, spec in RUNS.items():
        if not os.path.isdir(rundir(tag)):
            print(f"[{tag}] not present, skipping")
            continue
        meta_only = (os.path.exists(os.path.join(out_dir, f"{tag}_fields.npz"))
                     and not force)
        if meta_only:
            # summary.npz still needs the stats, but they are cheap to re-read.
            st = read_stats(tag)
            print(f"[{tag}] already extracted, skipping (FORCE=1 to redo)",
                  flush=True)
        else:
            print(f"[{tag}]", end=" ", flush=True)

        if not meta_only:
            st = read_stats(tag)
            np.savez_compressed(os.path.join(out_dir, f"{tag}_stats.npz"), **st)
            print(f"stats({st['t'].size})", end=" ", flush=True)

        data = None if meta_only else nt2.Data(rundir(tag))

        if not meta_only:
            fl = read_fields(data)
            np.savez_compressed(os.path.join(out_dir, f"{tag}_fields.npz"), **fl)
            print(f"fields({fl['t'].size})", end=" ", flush=True)

            # heat wrote no particle dumps -- it is a stats-only run by design.
            try:
                mo, ph = read_particles(data, x_max=float(fl["x"][-1]))
                np.savez_compressed(os.path.join(out_dir, f"{tag}_temps.npz"), **mo)
                np.savez_compressed(os.path.join(out_dir, f"{tag}_phase.npz"), **ph)
                print(f"temps({mo['t'].size}) phase({ph['t'].size})", flush=True)
            except Exception as exc:
                print(f"no particles ({type(exc).__name__})", flush=True)

        A0   = deck_value(tag, "A0", 0.0)
        npc  = deck_value(tag, "ppc0")
        w0   = deck_value(tag, "omega0", 0.0)
        dwdt = deck_value(tag, "domega_dt", 0.0)
        meta[tag] = dict(
            family         = spec["family"],
            loading        = deck_value(tag, "loading", "random"),
            one_v          = deck_value(tag, "one_v", "true"),
            vq_over_vthe   = spec["ratio"],
            A0             = A0,
            N_pc           = npc,
            N_p_total      = npc * 1000.0,   # Nx = 1000; the shot-noise formula
                                             # wants the TOTAL, not the per-cell
            runtime        = deck_value(tag, "runtime"),
            omega0         = w0,
            domega_dt      = dwdt,
            # resonance is crossed when omega(t) = omega_p = 1
            t_cross        = ((1.0 - w0) / dwdt) if dwdt else None,
            # eps_driven reaches the 1V electron thermal energy at this time
            t_sat          = (0.0632 / A0) if A0 else None,
            eps_noise      = EPS_NOISE_COEFF / npc,
            eps_thermal_1V = EPS_THERMAL_1V,
            eps_E_max_norm   = float(st["eps_E"].max()   / EPS_THERMAL_1V),
            eps_tot_max_norm = float(st["eps_tot"].max() / EPS_THERMAL_1V),
            t_at_eps_E_max   = float(st["t"][int(np.argmax(st["eps_E"]))]),
            dKE_e_end_norm   = float(st["dKE_e"][-1] / EPS_THERMAL_1V),
            # the heating gate: fractional drift over the whole run
            dTe_frac         = float(st["Te_rel"][-1] - 1.0),
            dTi_frac         = float(st["Ti_rel"][-1] - 1.0),
        )

    # ---- one row per run, ready to tabulate or scatter ---------------------
    tags = [t for t in RUNS if t in meta]
    np.savez_compressed(
        os.path.join(out_dir, "summary.npz"),
        tags             = np.array(tags),
        family           = np.array([meta[t]["family"] for t in tags]),
        loading          = np.array([str(meta[t]["loading"]) for t in tags]),
        vq_over_vthe     = np.array([meta[t]["vq_over_vthe"] for t in tags]),
        A0               = np.array([meta[t]["A0"] for t in tags]),
        N_pc             = np.array([meta[t]["N_pc"] for t in tags]),
        eps_noise        = np.array([meta[t]["eps_noise"] for t in tags]),
        eps_E_max_norm   = np.array([meta[t]["eps_E_max_norm"] for t in tags]),
        eps_tot_max_norm = np.array([meta[t]["eps_tot_max_norm"] for t in tags]),
        t_at_eps_E_max   = np.array([meta[t]["t_at_eps_E_max"] for t in tags]),
        dKE_e_end_norm   = np.array([meta[t]["dKE_e_end_norm"] for t in tags]),
        dTe_frac         = np.array([meta[t]["dTe_frac"] for t in tags]),
        # the paper's slow/fast boundary, (m_e/m_p)^(1/2)/2
        boundary         = np.array([0.5 / np.sqrt(MI_ME)]),
    )

    meta["_groups"] = GROUPS
    meta["_constants"] = dict(
        vth_e=float(VTH_E), eps_thermal_1V=EPS_THERMAL_1V, mi_me=MI_ME,
        t_ion=float(T_ION), lambda_D=float(VTH_E),
        k_lambdaD_1_at=1.0 / VTH_E,
        eps_noise_coeff=EPS_NOISE_COEFF,
        quiet_floor_t0=QUIET_FLOOR_T0, random_floor_t0=RANDOM_FLOOR_T0,
        heat_plateau=HEAT_PLATEAU, heat_random_floor=HEAT_RANDOM_FLOOR)
    meta["_linear_law"] = "eps_E + dKE_e = (A0^2/8)*(omega_p t)^2  -- not eps_E alone"
    meta["_thermal_ref"] = "5e-4 is the 1V value; T00_1 - 1 = 1.5e-3 is 3V and wrong here"
    with open(os.path.join(out_dir, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)

    print(f"\nwrote {out_dir}/")


if __name__ == "__main__":
    main()
