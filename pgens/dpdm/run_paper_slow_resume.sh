#!/bin/bash
#SBATCH -J dpdm_slowr                         # Job name
#SBATCH -o sbatch_logs/out/slowr.o%j          # Name of stdout output file
#SBATCH -e sbatch_logs/err/slowr.e%j          # Name of stderr error file
#SBATCH -p gpu-a100-small                    # 1 A100 per node; 1/4 the billing weight of gpu-a100
#SBATCH -N 1                                 # Total # of nodes
#SBATCH -n 1                                 # Total # of tasks (single GPU, no MPI)
#SBATCH -t 24:00:00                          # 100k steps left at ~516 ms = 14.3-15.6 h
#SBATCH --mail-type=all                      # Send email at begin and end of job
#SBATCH -A PHY23028                          # Project/Allocation name
#SBATCH --mail-user=montefalcone@utexas.edu
#
# Paper-faithful SLOW-growth run: v_q^D/vth_e = 1e-3, N_pc = 2e5, omega_p t = 1e4.
# Quiet start, 1V. 4e8 particles, 5e5 steps on one A100.
#
# Measured on an A100 at this exact config: 265 ms/step, FLAT over 2000 steps
# (286 -> 265 ms, no growth) => 5e5 * 0.265 s = 36.8 h.
# The CPU showed +48% growth over the same range; the GPU does not, because
# Kokkos uses atomics into HBM there rather than per-thread grid duplication.
#
#   sbatch run_paper_slow.sh

module load gcc/13.2.0     # binary links GCC 13's libstdc++
module load cuda/12.8

ROOT=$SCRATCH/dpdm/paper/slow
EXE=/work/09218/gab97/ls6/entity/build-dpdm-cuda/src/entity.xc
# Pre-create Entity's <name>/ output dir (harmless here, required under MPI).
mkdir -p $ROOT/slow

nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader 2>&1 | head -3

cat > $ROOT/slow.toml <<EOF
# Hook, Huang & Shalaby arXiv:2510.13956 -- slow-growth case, paper parameters.
#   v_q^D/vth_e = 1e-3  =>  A0 = 1e-3*sqrt(1e-3) = 3.16228e-5
#   t_sat = 0.0632/A0 = 1999;  runtime = 1e4 = 5 t_sat
#   measured shot-noise floor 2.06e-4/N_pc = 1.03e-9, beaten by omega_p t ~ 3
[simulation]
  name = "slow"
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
  ppc0 = 200000.0              # N_pc, per cell per species (4e8 particles total)
  [[particles.species]]
    label = "e-"
    mass = 1.0
    charge = -1.0
    maxnpart = 2.05e8          # 23 GB of the A100's 40 GB
    pusher = "Vay"
  [[particles.species]]
    label = "i+"
    mass = 1836.0
    charge = 1.0
    maxnpart = 2.05e8
    pusher = "Vay"
[setup]
  densities = [2.0]
  temperatures = [1.0e-3, 1.0e-3]
  drive = "resonance"
  loading = "quiet"            # deterministic lattice; noise floor 1.3e-17 vs 3.0e-6 random
  one_v = true                 # ux2 = ux3 = 0, matching SHARP
  A0 = 3.16228e-5
  omega = 1.0
[diagnostics]
  interval = 2000
  log_level = "WARNING"
  colored_stdout = false
[output]
  interval_time = 20.0         # 500 dumps over the run; t_sat = 1999
  [output.fields]
    quantities = ["E", "N_1", "N_2", "Rho", "J", "T00_1", "T00_2"]
  [output.particles]
    species = [1, 2]
    stride = 20000             # 1e4 particles/species/dump
  [output.stats]
    enable = true
    interval_time = 1.0
    quantities = ["E^2", "T00_1", "T00_2", "J.E"]
    custom = ["e_ext"]
  [output.spectra]
    enable = false             # defaults ON; a full reduction we do not need
[checkpoint]
  interval_time = 1000.0       # 10 checkpoints; 36.8 h against a 48 h wall
  keep = 2
EOF

S=$(date +%s)
cd $ROOT && $EXE -input slow.toml -resume > slow_resume.stdout 2>&1
rc=$?
echo "SLOW rc=$rc wall=$(( $(date +%s) - S ))s  $(date -Is)"
