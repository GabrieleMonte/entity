#!/bin/bash
#SBATCH -J dpdm_gate                         # Job name
#SBATCH -o sbatch_logs/out/gate.o%j          # Name of stdout output file
#SBATCH -e sbatch_logs/err/gate.e%j          # Name of stderr error file
#SBATCH -p development                       # ~20 min of work, well under the 2 h cap
#SBATCH -N 1                                 # Total # of nodes
#SBATCH -n 1                                 # Total # of tasks
#SBATCH -t 01:00:00                          # Run time
#SBATCH --mail-type=all                      # Send email at begin and end of job
#SBATCH -A PHY23028                          # Project/Allocation name
#SBATCH --mail-user=montefalcone@utexas.edu
#
# LOADER GATE. Undriven (A0 = 0), three loadings at identical settings.
#
# The paper's Fig. A00 shows undriven runs starting at eps_noise = 0 and rising
# to the analytic Poisson floor L^2/(12 Np Nx) within a few 1/omega_p, then
# sitting flat. Our old loader starts at 0 and stays ~12 orders below the floor
# for thousands of 1/omega_p, which seeds every driven run from round-off
# instead of shot noise.
#
# PASS for `quiet`:  eps_E(t=0) = 0 to round-off  AND  eps_E within ~2x of
#                    4.938e-8 by omega_p t ~ 10, flat thereafter.
# Expected `quiet_lattice`: starts at 0, stays far below the floor  (the bug).
# Expected `random`:        already at the floor at t=0, never zero (wrong IC).
#
#   sbatch run_gate_loader.sh

module load gcc/13.2.0     # binary is linked against GCC 13's libstdc++

export OMP_NUM_THREADS=112
export OMP_PROC_BIND=close
export OMP_PLACES=cores

ROOT=$SCRATCH/dpdm/gate
EXE=/work/09218/gab97/ls6/entity/build-dpdm/src/entity.xc

# N_pc = 2700, Nx = 1000, L = 40  ->  floor = 40^2/(12 * 2.7e6 * 1000) = 4.938e-8
for LOAD in quiet quiet_lattice random; do
  D=$ROOT/$LOAD
  mkdir -p $D/$LOAD
  cat > $D/$LOAD.toml <<EOF
[simulation]
  name = "$LOAD"
  engine = "srpic"
  runtime = 500.0
[grid]
  resolution = [1000]
  extent = [[0.0, 40.0]]
  [grid.metric]
    metric = "minkowski"
  [grid.boundaries]
    fields = [["PERIODIC"]]
    particles = [["PERIODIC"]]
[scales]
  skindepth0 = 1.0
  larmor0 = 1.0
[algorithms]
  current_filters = 0
  [algorithms.timestep]
    CFL = 0.5
[particles]
  ppc0 = 2700.0
  [[particles.species]]
    label = "e-"
    mass = 1.0
    charge = -1.0
    maxnpart = 3.2e6
    pusher = "Vay"
  [[particles.species]]
    label = "i+"
    mass = 1836.0
    charge = 1.0
    maxnpart = 3.2e6
    pusher = "Vay"
[setup]
  densities = [2.0]
  temperatures = [1.0e-3, 1.0e-3]
  drive = "resonance"
  loading = "$LOAD"
  one_v = true
  A0 = 0.0                     # UNDRIVEN: any field that appears is noise
  omega = 1.0
[diagnostics]
  interval = 2000
  log_level = "WARNING"
  colored_stdout = false
[output]
  interval_time = 50.0
  [output.fields]
    quantities = ["E", "N_1", "N_2", "Charge"]
  [output.stats]
    enable = true
    interval_time = 0.5        # fine: the rise happens in the first few 1/w_p
    quantities = ["E^2", "T00_1", "T00_2"]
  [output.particles]
    species = []
  [output.spectra]
    enable = false
EOF
  ( cd $D && $EXE -input $LOAD.toml > $LOAD.stdout 2>&1
    printf "%-14s exit=%s\n" "$LOAD" "$?" )
done

echo
echo "================= LOADER GATE ================="
python3 - <<'PY'
import numpy as np, os
L, Nx, NPC = 40.0, 1000, 2700
floor = L**2 / (12 * NPC*Nx * Nx)
print(f"analytic Poisson floor  L^2/(12 Np Nx) = {floor:.4e}\n")
print(f"{'loading':15s} {'eps_E(0)':>11s} {'wp t=1':>10s} {'10':>10s} {'50':>10s} {'200':>10s} {'500':>10s}")
root = os.environ["SCRATCH"] + "/dpdm/gate"
for tag in ("quiet", "quiet_lattice", "random"):
    f = f"{root}/{tag}/{tag}/{tag}_stats.csv"
    if not os.path.exists(f):
        print(f"{tag:15s} (missing)"); continue
    r = np.genfromtxt(f, delimiter=",", names=True, deletechars=" ", autostrip=True)
    t = np.asarray(r["time"], float); e = 0.5*np.asarray(r["E1^2"], float)
    vals = []
    for tt in (0, 1, 10, 50, 200, 500):
        i = min(int(np.searchsorted(t, tt)), len(e)-1)
        vals.append(e[i])
    print(f"{tag:15s} " + " ".join(f"{v:10.3e}" for v in vals))
    print(f"{'':15s} " + " ".join(f"{v/floor:10.2e}" for v in vals) + "   <- ratio to floor")
print("\nPASS for `quiet`: eps_E(0) ~ 0, and within ~2x of the floor by wp t ~ 10.")
PY
