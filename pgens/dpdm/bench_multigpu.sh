#!/bin/bash
#SBATCH -J dpdm_mgpu
#SBATCH -o sbatch_logs/out/mgpu.o%j
#SBATCH -e sbatch_logs/err/mgpu.e%j
#SBATCH -p gpu-a100-dev
#SBATCH -t 01:45:00
#SBATCH -A PHY23028
#
# Multi-GPU scaling at the real slow-run config (ppc0=2e5, 4e8 particles).
#   sbatch -N 1 -n 3 bench_multigpu.sh     # 3 GPUs, one node
#   sbatch -N 2 -n 6 bench_multigpu.sh     # 6 GPUs, two nodes
#
# 1000 steps per point, not 200: the CPU showed a real +48% rise in per-step
# cost that only plateaued around step 300, so a short run reads too fast.
# blocking_timers every 100 steps gives the curve, so this measures the
# steady-state rate AND settles whether the GPU has the same rise -- one queue
# wait instead of two.
set -uo pipefail
module load gcc/13.2.0
module load cuda/12.8

NG=${SLURM_NTASKS:-1}
EXE=/work/09218/gab97/ls6/entity/build-dpdm-cuda-mpi/src/entity.xc
ROOT=$SCRATCH/dpdm/bench/mgpu/g$NG
rm -rf $ROOT
# Pre-create Entity's <name>/ output dir: every rank races to mkdir it and the
# losers abort with "cannot create directory: File exists".
mkdir -p $ROOT/run

echo "=== $NG GPU(s), $SLURM_NNODES node(s) ==="
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader 2>&1 | head -3

cat > $ROOT/run.toml <<EOF
[simulation]
  name = "run"
  engine = "srpic"
  runtime = 20.0
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
  interval = 100
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

S=$(date +%s.%N)
( cd $ROOT && ibrun $EXE -input run.toml > run.stdout 2>&1 ); rc=$?
E=$(date +%s.%N)
echo "RESULT gpus=$NG nodes=$SLURM_NNODES rc=$rc wall=$(echo "$E - $S" | bc) for 1000 steps"
[ $rc -ne 0 ] && { echo "--- stdout tail ---"; tail -6 $ROOT/run.stdout; }
echo "--- per-step curve ---"
grep -aoE "Step: [0-9]+|Timestep duration: [0-9]+s? ?[0-9]*m?s?" $ROOT/run/run.out 2>/dev/null | paste - - | tail -10
echo "MGPU_DONE gpus=$NG"
