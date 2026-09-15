#!/bin/bash
#SBATCH -J dpdm_diag
#SBATCH -o sbatch_logs/out/diag.o%j
#SBATCH -e sbatch_logs/err/diag.e%j
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -A PHY23028
# Which kernel eats the time at ppc0 = 2e5?  blocking_timers gives a per-kernel
# breakdown. 10 steps is enough -- we are looking at the per-step split.
set -uo pipefail
module load gcc/13.2.0
EXE=/work/09218/gab97/ls6/entity/build-dpdm/src/entity.xc
ROOT=$SCRATCH/dpdm/bench/diag
rm -rf $ROOT; mkdir -p $ROOT
export OMP_NUM_THREADS=112 OMP_PROC_BIND=close OMP_PLACES=cores
cat > $ROOT/diag.toml <<EOF
[simulation]
  name = "diag"
  engine = "srpic"
  runtime = 0.2
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
  ppc0 = 200000.0
  [[particles.species]]
    label = "e-"
    mass = 1.0
    charge = -1.0
    maxnpart = 2.05e8
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
  loading = "quiet"
  one_v = true
  A0 = 3.16228e-5
  omega = 1.0
[diagnostics]
  interval = 1
  blocking_timers = true
  log_level = "WARNING"
  colored_stdout = false
[output]
  interval_time = 1000.0
  [output.fields]
    quantities = ["E"]
  [output.particles]
    enable = false
  [output.stats]
    enable = false
  [output.spectra]
    enable = false
[checkpoint]
  interval_time = 0.0
  keep = 0
EOF
S=$(date +%s)
( cd $ROOT && $EXE -input diag.toml > diag.stdout 2>&1 ); rc=$?
echo "DIAG rc=$rc wall=$(( $(date +%s) - S ))s for 10 steps"
echo "DIAG_DONE"
