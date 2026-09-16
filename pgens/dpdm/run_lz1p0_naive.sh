#!/bin/bash
#SBATCH -J lz1p0n
#SBATCH -o sbatch_logs/out/lz1p0n.o%j
#SBATCH -e sbatch_logs/err/lz1p0n.e%j
#SBATCH -p gpu-a100-small
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 08:00:00
#SBATCH --mail-type=all
#SBATCH -A PHY23028
#SBATCH --mail-user=montefalcone@utexas.edu
#
# Landau-Zener, v_q^D/vth_e = 1.0  (paper Fig. LZcase)
#   omega(t) = 0.8 w_p + 0.1 w_p t/5000  ->  domega_dt = 2.0e-5   [paper line 551]
#   crossing at omega_p t = 10000;  nonlinearity kills transfer near 5000
#   N_pc = 2700.0 from the paper's LZ criterion t_noise <= 0.02
#     (driven energy 50x shot noise; t_noise = L/(A0 sqrt(3 Np Nx/2)), Np TOTAL)
#   Resume with:  $EXE -input lz1p0n.toml -resume
#   DOMAIN COUNT IS LOCKED: a checkpoint written by N domains needs N to resume.

module load gcc/13.2.0
module load cuda/12.8

ROOT=$SCRATCH/dpdm/lz/lz1p0n
EXE=/work/09218/gab97/ls6/entity/build-dpdm-cuda/src/entity.xc
# pre-create the <name>/ dir: under MPI every rank races to mkdir it
mkdir -p $ROOT/lz1p0n

cat > $ROOT/lz1p0n.toml <<EOF
[simulation]
  name    = "lz1p0n"
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
  ppc0 = 2700.0
  [[particles.species]]
    label = "e-"
    mass = 1.0
    charge = -1.0
    maxnpart = 3.0e6
    pusher = "Vay"
  [[particles.species]]
    label = "i+"
    mass = 1836.0
    charge = 1.0
    maxnpart = 3.0e6
    pusher = "Vay"
[setup]
  densities    = [2.0]
  temperatures = [1.0e-3, 1.0e-3]
  drive        = "landau_zener"
  loading      = "quiet"
  one_v        = true
  A0           = 3.16228e-2
  # THE ONLY DIFFERENCE FROM run_lz1p0.sh. "naive" uses phi = omega(t)*t, the
  # literal reading of the paper's "A_0 cos(omega t)" with
  # "omega(t) = 0.8 w_p + 0.1 w_p t/5000". Its instantaneous frequency is
  # 0.8 + 2*domega_dt*t, i.e. TWICE the nominal sweep, so resonance is crossed
  # at omega_p t = 5000 instead of 1e4.
  #
  # This is a SENSITIVITY TEST, not a correction: "integrated" is the correct
  # convention and remains the default. The question it answers is whether a 2x
  # change in sweep rate can move dKE_e by the ~500x that separates run lz1p0
  # from the paper's stated "not even changed by a factor of 2".
  sweep_phase  = "naive"
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
    stride  = 270
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

cd $ROOT && $EXE -input lz1p0n.toml > lz1p0n.stdout 2>&1
echo "exit=$? $(date -Is)"
