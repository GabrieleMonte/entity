#!/bin/bash
#SBATCH -J lz0p1
#SBATCH -o sbatch_logs/out/lz0p1.o%j
#SBATCH -e sbatch_logs/err/lz0p1.e%j
#SBATCH -p gpu-a100
#SBATCH -N 1
#SBATCH -n 3
#SBATCH -t 48:00:00
#SBATCH --mail-type=all
#SBATCH -A PHY23028
#SBATCH --mail-user=montefalcone@utexas.edu
#
# Landau-Zener, v_q^D/vth_e = 0.1  (paper Fig. LZcase)
#   omega(t) = 0.8 w_p + 0.1 w_p t/5000  ->  domega_dt = 2.0e-5   [paper line 551]
#   crossing at omega_p t = 10000;  nonlinearity kills transfer near 5000
#   N_pc = 270000.0 from the paper's LZ criterion t_noise <= 0.02
#     (driven energy 50x shot noise; t_noise = L/(A0 sqrt(3 Np Nx/2)), Np TOTAL)
#   Resume with:  $EXE -input lz0p1.toml -resume
#   DOMAIN COUNT IS LOCKED: a checkpoint written by N domains needs N to resume.

module load gcc/13.2.0
module load cuda/12.8

ROOT=$SCRATCH/dpdm/lz/lz0p1
EXE=/work/09218/gab97/ls6/entity/build-dpdm-cuda-mpi/src/entity.xc
# pre-create the <name>/ dir: under MPI every rank races to mkdir it
mkdir -p $ROOT/lz0p1

cat > $ROOT/lz0p1.toml <<EOF
[simulation]
  name    = "lz0p1"
  engine  = "srpic"
  runtime = 20000.0
[grid]
  resolution = [1000]
  extent     = [[0.0, 40.0]]
  [grid.metric]
    metric = "minkowski"
  [grid.boundaries]
    fields    = [["PERIODIC"]]
    particles = [["PERIODIC"]]
[scales]
  skindepth0 = 1.0
  larmor0    = 1.0
[algorithms]
  current_filters = 0
  [algorithms.timestep]
    CFL = 0.5
[particles]
  ppc0 = 270000.0
  [[particles.species]]
    label = "e-"
    mass = 1.0
    charge = -1.0
    maxnpart = 2.8e8
    pusher = "Vay"
  [[particles.species]]
    label = "i+"
    mass = 1836.0
    charge = 1.0
    maxnpart = 2.8e8
    pusher = "Vay"
[setup]
  densities    = [2.0]
  temperatures = [1.0e-3, 1.0e-3]
  drive        = "landau_zener"
  loading      = "quiet"
  one_v        = true
  A0           = 3.16228e-3
  omega0       = 0.8
  domega_dt    = 2.0e-5
[diagnostics]
  interval       = 2000
  log_level      = "WARNING"
  colored_stdout = false
[output]
  interval_time = 40.0
  [output.fields]
    quantities = ["E", "N_1", "N_2", "Rho", "J", "T00_1", "T00_2"]
  [output.particles]
    species = [1, 2]
    stride  = 27000
  [output.stats]
    enable        = true
    interval_time = 2.0
    quantities    = ["E^2", "T00_1", "T00_2", "J.E"]
    custom        = ["e_ext"]
  [output.spectra]
    enable = false
[checkpoint]
  interval_time = 1000.0
  keep          = 2
EOF

cd $ROOT && ibrun $EXE -input lz0p1.toml > lz0p1.stdout 2>&1
echo "exit=$? $(date -Is)"
