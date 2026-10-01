#!/bin/bash
# Knob-scan on the paper's fast resonance case (see dpdm-hhs-burst-deficit):
# one numerics knob per run, judged on dE_e(k!=0) at wp t = 500 moving from
# our 1.0e-2 toward HHS's 1.9e-2. Baseline control re-runs fast_s unchanged
# on the current binary (build hash differs from the Sep 16 run).
#
# Generates sbatch scripts; run this, then submit the ones you want.
set -e
cd "$(dirname "$0")"

make_job () {
  local tag=$1 queue=$2 walltime=$3 cfl=$4 nx=$5 ppc=$6 maxn=$7 stride=$8 est=$9
  cat > run_knob_${tag}.sh <<EOF
#!/bin/bash
#SBATCH -J k_${tag}
#SBATCH -o ../sbatch_logs/out/k_${tag}.o%j
#SBATCH -e ../sbatch_logs/err/k_${tag}.e%j
#SBATCH -p ${queue}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t ${walltime}
#SBATCH --mail-type=all
#SBATCH -A PHY23028
#SBATCH --mail-user=montefalcone@utexas.edu
#
# Knob run '${tag}': CFL=${cfl}, Nx=${nx}, ppc0=${ppc}. Estimated ${est}.

module load gcc/13.2.0

export OMP_NUM_THREADS=112
export OMP_PROC_BIND=close
export OMP_PLACES=cores

ROOT=\$SCRATCH/dpdm/knobs/${tag}
EXE=/work/09218/gab97/ls6/entity/build-dpdm/src/entity.xc
mkdir -p \$ROOT

cat > \$ROOT/fast.toml <<TOML
# fast resonance knob run '${tag}' -- baseline is run_paper_fast_shear.sh
[simulation]
  name = "fast"
  engine = "srpic"
  runtime = 10000.0
[grid]
  resolution = [${nx}]
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
    CFL = ${cfl}
[particles]
  ppc0 = ${ppc}
  [[particles.species]]
    label = "e-"
    mass = 1.0
    charge = -1.0
    maxnpart = ${maxn}
    pusher = "Vay"
  [[particles.species]]
    label = "i+"
    mass = 1836.0
    charge = 1.0
    maxnpart = ${maxn}
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
    stride = ${stride}
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

cd \$ROOT && \$EXE -input fast.toml > fast.stdout 2>&1
echo "exit=\$? \$(date -Is)"
EOF
  echo "wrote run_knob_${tag}.sh"
}

#        tag   queue        walltime  CFL   Nx    ppc     maxn   stride est
make_job base  development  01:30:00  0.5   1000  200.0   2.5e5  20     "41 min (control, = fast_s)"
make_job dt    development  02:00:00  0.25  1000  200.0   2.5e5  20     "82 min (dt = 0.01)"
make_job nx    normal       05:00:00  0.5   2000  200.0   4.5e5  40     "2.7 h (dx = 0.02, dt = 0.01 via CFL)"
make_job npc   normal       10:00:00  0.5   1000  2000.0  2.2e6  200    "7 h (eps_noise /10)"
