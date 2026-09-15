#!/bin/bash
#SBATCH -J dpdm_mpival
#SBATCH -o sbatch_logs/out/mpival.o%j
#SBATCH -e sbatch_logs/err/mpival.e%j
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 00:30:00
#SBATCH -A PHY23028
#
# Does the quiet start survive MPI domain decomposition?
# The loading is keyed on velocity CLASS, never on cell index, so 8 subdomains
# of 125 cells must reproduce the single-domain run. If it doesn't, the 134
# node-hour slow run is not trustworthy on 8 nodes.
set -uo pipefail
module load gcc/13.2.0

ROOT=$SCRATCH/dpdm/mpival
SER=/work/09218/gab97/ls6/entity/build-dpdm/src/entity.xc
MPI=/work/09218/gab97/ls6/entity/build-dpdm-mpi/src/entity.xc
rm -rf $ROOT; mkdir -p $ROOT/serial $ROOT/mpi8

deck () {
  cat > $1/$2.toml <<EOF
[simulation]
  name = "$2"
  engine = "srpic"
  runtime = 200.0
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
  ppc0 = 64.0
  [[particles.species]]
    label = "e-"
    mass = 1.0
    charge = -1.0
    maxnpart = 1.0e5
    pusher = "Vay"
  [[particles.species]]
    label = "i+"
    mass = 1836.0
    charge = 1.0
    maxnpart = 1.0e5
    pusher = "Vay"
[setup]
  densities = [2.0]
  temperatures = [1.0e-3, 1.0e-3]
  drive = "resonance"
  loading = "quiet"
  one_v = true
  A0 = 9.4868e-4
  omega = 1.0
[diagnostics]
  interval = 1000
  log_level = "WARNING"
  colored_stdout = false
[output]
  interval_time = 200.0
  [output.fields]
    quantities = ["E", "N_1", "N_2", "Charge"]
  [output.particles]
    enable = false
  [output.stats]
    enable = true
    interval_time = 0.5
    quantities = ["E^2", "T00_1", "T00_2"]
[checkpoint]
  interval_time = 0.0
  keep = 0
EOF
}
deck $ROOT/serial serial
deck $ROOT/mpi8   mpi8

echo "--- serial (1 domain, 112 threads) ---"
( cd $ROOT/serial && OMP_NUM_THREADS=112 $SER -input serial.toml > serial.stdout 2>&1; echo "exit=$?" )

echo "--- mpi8 (8 domains x 125 cells, 14 threads each) ---"
( cd $ROOT/mpi8 && OMP_NUM_THREADS=14 ibrun -n 8 $MPI -input mpi8.toml > mpi8.stdout 2>&1; echo "exit=$?" )
echo "MPIVAL_DONE"
