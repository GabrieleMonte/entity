#!/bin/bash
#SBATCH -J dpdm_run34                        # Job name
#SBATCH -o sbatch_logs/out/launcher.o%j      # Name of stdout output file
#SBATCH -e sbatch_logs/err/launcher.e%j      # Name of stderr error file
#SBATCH -p normal                            # Queue (partition) name
#SBATCH -N 1                                 # Total # of nodes
#SBATCH -n 1                                 # Total # of tasks
#SBATCH -t 06:00:00                          # Run time (hh:mm:ss)
#SBATCH --mail-type=all                      # Send email at begin and end of job
#SBATCH -A PHY23028                          # Project/Allocation name
#SBATCH --mail-user=montefalcone@utexas.edu
#
# DPDM resonant conversion (arXiv:2510.13956) -- Runs 3 & 4.
# Writes its own input decks, then works through them two at a time: the node
# has two 64-core EPYC sockets, so one run per socket. ~3.7 h.
#
#   sbatch run34.sh

# entity.xc is linked against GCC 13's libstdc++; the default LS6 environment
# has gcc/9.4.0 and the binary then dies at load time (GLIBCXX_3.4.31 not found)
# before it ever reads the input file.
module load gcc/13.2.0

export OMP_NUM_THREADS=64
export OMP_PROC_BIND=close
export OMP_PLACES=cores

ROOT=$SCRATCH/dpdm/run34
EXE=/work/09218/gab97/ls6/entity/build-dpdm/src/entity.xc

# A0 = (v_q^D/vth_e) * sqrt(1e-3).  Longest run first: the dynamic scheduler
# must not strand conv8000 behind short jobs.
#      tag       A0           ppc0  maxnpart  runtime  fld_int  stat_int
RUNS="
conv8000  9.48683e-5   8000.0  1.0e7   3000.0   10.0   1.0
r1e-3     3.16228e-5   2000.0  2.5e6   6000.0   20.0   1.0
r3e-2     9.48683e-4   2000.0  2.5e6   3000.0   10.0   1.0
r1e-1     3.16228e-3   2000.0  2.5e6   3000.0   10.0   0.5
r3e-3     9.48683e-5   2000.0  2.5e6   3000.0   10.0   1.0
r1e-2     3.16228e-4   2000.0  2.5e6   3000.0   10.0   1.0
r3e-1     9.48683e-3   2000.0  2.5e6   3000.0   10.0   0.5
conv500   9.48683e-5    500.0  6.5e5   3000.0   10.0   1.0
"

mkdir -p $ROOT && : > $ROOT/queue
while read -r tag A0 ppc0 maxnpart runtime fld_int stat_int; do
  [ -z "$tag" ] && continue
  mkdir -p $ROOT/$tag
  cat > $ROOT/$tag/$tag.toml <<EOF
[simulation]
  name    = "$tag"
  engine  = "srpic"
  runtime = $runtime

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
  current_filters = 0          # keep 0: saturation physics lives at high k
  [algorithms.timestep]
    CFL = 0.5                  # dt = 0.02 / omega_p

[particles]
  ppc0 = $ppc0
  [[particles.species]]
    label = "e-"
    mass = 1.0
    charge = -1.0
    maxnpart = $maxnpart
    pusher = "Vay"
  [[particles.species]]
    label = "i+"
    mass = 1836.0
    charge = 1.0
    maxnpart = $maxnpart
    pusher = "Vay"

[setup]
  densities    = [2.0]         # n_e = n0 => omega_pe = 1 exactly
  temperatures = [1.0e-3, 1.0e-3]
  drive        = "resonance"
  A0           = $A0
  omega        = 1.0

[diagnostics]
  interval       = 200         # .out progress bar: default 1 = 28 MB / 2e4 steps
  log_level      = "WARNING"   # .log VERB trace: default = 4 kB PER STEP
  colored_stdout = false

[output]
  interval_time = $fld_int
  [output.fields]
    quantities = ["E", "N_1", "N_2", "Rho", "J", "T00_1", "T00_2"]
  [output.particles]
    species = [1, 2]
    stride  = 200
  [output.stats]
    enable        = true
    interval_time = $stat_int
    quantities    = ["E^2", "T00_1", "T00_2", "J.E"]
    custom        = ["e_ext"]
  [output.spectra]
    enable = false

[checkpoint]
  interval_time = 1000.0
  keep          = 2
EOF
  echo "$tag" >> $ROOT/queue
done <<< "$RUNS"

# Two workers, one per socket, both pulling from the same open file descriptor:
# whichever finishes first takes the next run. taskset is not optional --
# OMP_PROC_BIND is per-process, so without it both workers would bind their
# threads to the same cores (measured 8x slowdown).
exec 3< $ROOT/queue
for sock in 0 1; do
  (
    lo=$(( sock * 64 )); hi=$(( lo + 63 ))
    while read -r tag <&3; do
      echo "socket $sock -> $tag  $(date +%T)"
      ( cd $ROOT/$tag && taskset -c $lo-$hi $EXE -input $tag.toml > $tag.stdout 2>&1 )
      echo "socket $sock    $tag done (exit $?)  $(date +%T)"
    done
  ) &
done
wait
echo "all runs finished $(date +%T)"
