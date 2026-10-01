#!/bin/bash
#SBATCH -J k_nx4
#SBATCH -o ../sbatch_logs/out/k_nx4.o%j
#SBATCH -e ../sbatch_logs/err/k_nx4.e%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH --mail-type=all
#SBATCH -A PHY23028
#SBATCH --mail-user=montefalcone@utexas.edu
#
# Convergence ladder point 3: dx=0.01 (Nx=4000), run to t=5000 only -- the window where base and nx disagree: CFL=0.5, Nx=2000, ppc0=200.0. Estimated 2.7 h (dx = 0.02, dt = 0.01 via CFL).

module load gcc/13.2.0

export OMP_NUM_THREADS=112
export OMP_PROC_BIND=close
export OMP_PLACES=cores

ROOT=$SCRATCH/dpdm/knobs/nx4
EXE=/work/09218/gab97/ls6/entity/build-dpdm-res/src/entity.xc
mkdir -p $ROOT

cat > $ROOT/fast.toml <<TOML
# fast resonance knob run 'nx' -- baseline is run_paper_fast_shear.sh
[simulation]
  name = "fast"
  engine = "srpic"
  runtime = 5000.0
[grid]
  resolution = [4000]
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
  ppc0 = 200.0
  [[particles.species]]
    label = "e-"
    mass = 1.0
    charge = -1.0
    maxnpart = 9.0e5
    pusher = "Vay"
  [[particles.species]]
    label = "i+"
    mass = 1836.0
    charge = 1.0
    maxnpart = 9.0e5
    pusher = "Vay"
[setup]
  densities = [2.0]
  temperatures = [1.0e-3, 1.0e-3]
  drive = "resonance"
  loading = "quiet"
  one_v = true
  A0 = 9.48683e-4
  omega = 1.0
[diagnostics]
  interval = 1000
  log_level = "WARNING"
  colored_stdout = false
[output]
  interval_time = 10.0
  [output.fields]
    quantities = ["E", "N_1", "N_2", "Rho", "J", "T00_1", "T00_2"]
  [output.particles]
    species = [1, 2]
    stride = 80
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
TOML

cd $ROOT && $EXE -input fast.toml > fast.stdout 2>&1
echo "exit=$? $(date -Is)"
