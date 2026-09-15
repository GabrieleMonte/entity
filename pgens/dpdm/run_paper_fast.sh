#!/bin/bash
#SBATCH -J dpdm_fast                         # Job name
#SBATCH -o sbatch_logs/out/fast.o%j          # Name of stdout output file
#SBATCH -e sbatch_logs/err/fast.e%j          # Name of stderr error file
#SBATCH -p development                       # Queue: 41 min job, well under the 2 h dev cap
#SBATCH -N 1                                 # Total # of nodes
#SBATCH -n 1                                 # Total # of tasks
#SBATCH -t 01:30:00                          # Run time (measured 41 min)
#SBATCH --mail-type=all                      # Send email at begin and end of job
#SBATCH -A PHY23028                          # Project/Allocation name
#SBATCH --mail-user=montefalcone@utexas.edu
#
# Paper-faithful fast-growth run: v_q^D/vth_e = 3e-2, N_pc = 200, omega_p t = 1e4.
# Quiet start, 1V. 4e5 particles, 5e5 steps -> ~8 min on one node.
#
#   sbatch run_paper_fast.sh

module load gcc/13.2.0     # binary is linked against GCC 13's libstdc++

export OMP_NUM_THREADS=112
export OMP_PROC_BIND=close
export OMP_PLACES=cores

ROOT=$SCRATCH/dpdm/paper/fast
EXE=/work/09218/gab97/ls6/entity/build-dpdm/src/entity.xc
mkdir -p $ROOT

cat > $ROOT/fast.toml <<EOF
# Hook, Huang & Shalaby arXiv:2510.13956 -- fast-growth case, paper parameters.
#   v_q^D/vth_e = 3e-2  =>  A0 = 3e-2*sqrt(1e-3) = 9.48683e-4
#   t_sat = 0.0632/A0 = 67;  1/omega_pi = 42.8;  runtime = 1e4 = 150 t_sat
[simulation]
  name = "fast"
  engine = "srpic"
  runtime = 10000.0
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
  ppc0 = 200.0                 # N_pc, per cell per species (Np = N_pc*Nx = 2e5)
  [[particles.species]]
    label = "e-"
    mass = 1.0
    charge = -1.0
    maxnpart = 2.5e5
    pusher = "Vay"
  [[particles.species]]
    label = "i+"
    mass = 1836.0
    charge = 1.0
    maxnpart = 2.5e5
    pusher = "Vay"
[setup]
  densities = [2.0]
  temperatures = [1.0e-3, 1.0e-3]
  drive = "resonance"
  loading = "quiet"            # deterministic lattice, machine-precision neutral
  one_v = true                 # ux2 = ux3 = 0, matching SHARP
  A0 = 9.48683e-4
  omega = 1.0
[diagnostics]
  interval = 1000
  log_level = "WARNING"
  colored_stdout = false
[output]
  interval_time = 10.0         # 1000 dumps; t_sat = 67 is well resolved
  [output.fields]
    quantities = ["E", "N_1", "N_2", "Rho", "J", "T00_1", "T00_2"]
  [output.particles]
    species = [1, 2]
    stride = 20                # 1e4 particles/species/dump
  [output.stats]
    enable = true
    interval_time = 0.5
    quantities = ["E^2", "T00_1", "T00_2", "J.E"]
    custom = ["e_ext"]
  [output.spectra]
    enable = false
[checkpoint]
  interval_time = 2500.0
  keep = 2
EOF

cd $ROOT && $EXE -input fast.toml > fast.stdout 2>&1
echo "exit=$? $(date -Is)"
