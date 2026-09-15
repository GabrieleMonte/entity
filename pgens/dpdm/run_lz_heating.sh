#!/bin/bash
#SBATCH -J lz_heat                           # Job name
#SBATCH -o sbatch_logs/out/lzheat.o%j        # stdout
#SBATCH -e sbatch_logs/err/lzheat.e%j        # stderr
#SBATCH -p gpu-a100-small                    # 1 A100
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 04:00:00                          # ~1.5 h expected at 5.4e6 particles
#SBATCH --mail-type=all
#SBATCH -A PHY23028
#SBATCH --mail-user=montefalcone@utexas.edu
#
# Step 1 of the Landau-Zener plan: undriven (A0 = 0) heating check at the LZ
# particle count and the FULL LZ duration, omega_p t = 2e4 = 1e6 steps.
#
# Why: our numerical-heating number (-4.5e-10 per 1/omega_p, scaling 1/Np) was
# measured with RANDOM loading over 2e4 steps. An ordered lattice can behave
# differently, and 1e6 steps is 50x longer. For LZ this is the dangerous case:
# a slow spurious drift in T_e shifts omega_p and would detune the sweep -- it
# would look exactly like the physics the run exists to measure.
#
# Also re-checks that the quiet start's noise floor holds over 1e6 steps; the
# slow run showed the per-step COST changes after saturation, so "it was flat
# for 2000 steps" is not evidence about 1e6.

module load gcc/13.2.0
module load cuda/12.8

ROOT=$SCRATCH/dpdm/lz/heat
EXE=/work/09218/gab97/ls6/entity/build-dpdm-cuda/src/entity.xc
mkdir -p $ROOT/heat

cat > $ROOT/heat.toml <<EOF
# Undriven quiet-start heating check at LZ duration.
#   N_pc = 2.7e3 (the v_q/vth=1 LZ requirement, t_noise <= 0.02)
#   A0 = 0: nothing should happen. Any drift in T00_1 is numerical.
[simulation]
  name    = "heat"
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
  drive        = "resonance"
  loading      = "quiet"
  one_v        = true
  A0           = 0.0                 # undriven
  omega        = 1.0
[diagnostics]
  interval       = 2000
  log_level      = "WARNING"
  colored_stdout = false
[output]
  interval_time = 2000.0             # 10 field dumps; this run is about the stats
  [output.fields]
    quantities = ["E", "N_1"]
  [output.particles]
    enable = false
  [output.stats]
    enable        = true
    interval_time = 10.0             # 2000 rows: T_e drift needs fine sampling
    quantities    = ["E^2", "T00_1", "T00_2"]
  [output.spectra]
    enable = false
[checkpoint]
  interval_time = 5000.0
  keep          = 2
EOF

cd $ROOT && $EXE -input heat.toml > heat.stdout 2>&1
echo "exit=$? $(date -Is)"
