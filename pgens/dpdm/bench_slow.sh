#!/bin/bash
#SBATCH -J dpdm_bench
#SBATCH -o sbatch_logs/out/bench.o%j
#SBATCH -e sbatch_logs/err/bench.e%j
#SBATCH -p development
#SBATCH -t 00:50:00
#SBATCH -A PHY23028
#
# Measure the REAL rate at the slow run's ppc0 = 2e5 (4e8 particles).
# Two step counts, differenced, so initialization is excluded.
#   sbatch -N 1 -n 1 -c 128 bench_slow.sh
#   sbatch -N 8 -n 8 -c 128 bench_slow.sh
#
# -c 128 IS MANDATORY. Without it SLURM gives each task ONE cpu, ibrun binds the
# task to that single core, and all OMP_NUM_THREADS threads pile onto it: 112x
# slower, which is how the first attempt hit its 40-minute walltime. Invoking
# the binary directly (no ibrun) hides the bug -- it is not confined.
set -uo pipefail
module load gcc/13.2.0

NN=${SLURM_NNODES:-1}
ROOT=$SCRATCH/dpdm/bench/n$NN
MPI=/work/09218/gab97/ls6/entity/build-dpdm-mpi/src/entity.xc
rm -rf $ROOT; mkdir -p $ROOT
export OMP_NUM_THREADS=112 OMP_PROC_BIND=close OMP_PLACES=cores

deck () {   # $1 dir  $2 tag  $3 runtime
  # Pre-create BOTH the run dir and Entity's own <name>/ output subdir: under
  # MPI every rank races to mkdir the latter and the losers abort with
  # "filesystem error: cannot create directory: File exists". Not rank-guarded
  # upstream; only shows up in MPI mode.
  mkdir -p $1 $1/$2
  cat > $1/$2.toml <<EOF
[simulation]
  name = "$2"
  engine = "srpic"
  runtime = $3
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
    maxnpart = 2.6e8
    pusher = "Vay"
  [[particles.species]]
    label = "i+"
    mass = 1836.0
    charge = 1.0
    maxnpart = 2.6e8
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
  [output.spectra]
    enable = false
[checkpoint]
  interval_time = 0.0
  keep = 0
EOF
}

for R in 2.0 6.0; do
  T=s$(echo $R | tr -d '.')
  deck $ROOT/$T $T $R
  S=$(date +%s.%N)
  ( cd $ROOT/$T && ibrun $MPI -input $T.toml > $T.stdout 2>&1 )
  E=$(date +%s.%N)
  echo "RESULT nodes=$NN runtime=$R wall=$(echo "$E - $S" | bc)"
  # Verify the backend AT RUNTIME. A Serial-backend build looks identical from
  # the outside but runs ~112x slower; that cost a 40-minute walltime once.
  grep -aoE "Kokkos: v[0-9.]+|OpenMP|Serial" $ROOT/$T/$T/$T.out 2>/dev/null | sort -u | tr "\n" " " | sed "s/^/  backend: /"; echo
  grep -aE "Timestep duration" $ROOT/$T/$T/$T.out 2>/dev/null | tail -1 | sed "s/^/  /"
done
echo "BENCH_DONE nodes=$NN"
