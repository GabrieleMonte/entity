#!/bin/bash
# Staged validation of the quiet start. All four runs are seconds long.
#   bash validate_quiet.sh
set -uo pipefail
module load gcc/13.2.0 2>/dev/null
EXE=/work/09218/gab97/ls6/entity/build-dpdm/src/entity.xc
ROOT=$SCRATCH/dpdm/qstest
rm -rf $ROOT; mkdir -p $ROOT

deck () {   # $1 tag  $2 loading  $3 A0  $4 runtime  $5 ppc0  $6 out_int
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
  ppc0 = $5
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
  loading = "$2"
  one_v = true
  A0 = $3
  omega = 1.0
[diagnostics]
  interval = 500
  log_level = "WARNING"
  colored_stdout = false
[output]
  interval_time = $6
  [output.fields]
    quantities = ["E", "N_1", "N_2", "Rho", "Charge"]
  [output.particles]
    species = [1, 2]
    stride = 1
  [output.stats]
    enable = true
    interval_time = 0.5
    quantities = ["E^2", "T00_1", "T00_2"]
[checkpoint]
  interval_time = 0.0
  keep = 0
EOF
}

#     tag          loading  A0         runtime ppc0 out_int
deck  t0_quiet quiet 0.0 0.04 64.0   0.02
deck  undriven_q quiet 0.0 400.0 64.0   20.0
deck  undriven_r random 0.0 400.0 64.0   20.0
deck  driven_q quiet 9.4868e-4 100.0 64.0   20.0

for t in t0_quiet undriven_q undriven_r driven_q; do
  ( cd $ROOT/$t && OMP_NUM_THREADS=8 taskset -c 0-7 $EXE -input $t.toml > $t.stdout 2>&1
    printf "%-12s exit=%s\n" "$t" "$?" )
done
echo; echo "runs complete -> $ROOT"
