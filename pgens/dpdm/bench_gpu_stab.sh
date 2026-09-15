#!/bin/bash
#SBATCH -J dpdm_gstab
#SBATCH -o sbatch_logs/out/gstab.o%j
#SBATCH -e sbatch_logs/err/gstab.e%j
#SBATCH -p gpu-a100-dev
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:45:00
#SBATCH -A PHY23028
#
# Is a 1000-cell grid with 2e5 particles/cell viable on an A100?
# Huge particle count favours the GPU; a tiny grid punishes it, because Kokkos
# ScatterView defaults to ScatterAtomic on CUDA (duplication on CPU), so every
# deposit is an atomic onto one of only 1000 locations.
# Two ppc0 values: contention per cell scales with ppc0, so the pair separates
# "GPU is slow here" from "GPU is slow because of atomics".
set -uo pipefail
module load gcc/13.2.0
module load cuda/12.8

echo "=== GPU hardware ==="; nvidia-smi --query-gpu=index,name,memory.total --format=csv 2>&1 | head -5

EXE=/work/09218/gab97/ls6/entity/build-dpdm-cuda/src/entity.xc
ROOT=$SCRATCH/dpdm/bench/gstab
rm -rf $ROOT; mkdir -p $ROOT

deck () {   # $1 tag  $2 ppc0  $3 maxnpart  $4 runtime
  mkdir -p $ROOT/$1
  cat > $ROOT/$1/$1.toml <<EOF
[simulation]
  name = "$1"
  engine = "srpic"
  runtime = $4
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
  ppc0 = $2
  [[particles.species]]
    label = "e-"
    mass = 1.0
    charge = -1.0
    maxnpart = $3
    pusher = "Vay"
  [[particles.species]]
    label = "i+"
    mass = 1836.0
    charge = 1.0
    maxnpart = $3
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
  interval = 200
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
    enable = true
    interval_time = 1000.0
    quantities = ["E^2"]
[checkpoint]
  interval_time = 0.0
  keep = 0
EOF
}

run () {   # $1 label  $2 ppc0  $3 maxnpart
  for R in 40.0; do
    T="$1_$(echo $R | tr -d '.')"
    deck $T $2 $3 $R
    S=$(date +%s.%N)
    ( cd $ROOT/$T && $EXE -input $T.toml > $T.stdout 2>&1 ); rc=$?
    E=$(date +%s.%N)
    echo "RESULT gpu label=$1 ppc0=$2 runtime=$R rc=$rc wall=$(echo "$E - $S" | bc)"
    [ $rc -ne 0 ] && { echo "  --- tail ---"; tail -4 $ROOT/$T/$T.stdout; }
  done
}

run long 200000.0 2.05e8

echo "GPUBENCH_DONE"
