#!/bin/bash
#SBATCH -J k_nx4b
#SBATCH -o ../sbatch_logs/out/k_nx4b.o%j
#SBATCH -e ../sbatch_logs/err/k_nx4b.e%j
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH --mail-type=all
#SBATCH -A PHY23028
#SBATCH --mail-user=montefalcone@utexas.edu
#
# nx4 leg B: resume from the last checkpoint leg A wrote before its 2 h wall
# (stats CSV is appended on resume, header only written on fresh start).

module load gcc/13.2.0
export OMP_NUM_THREADS=112
export OMP_PROC_BIND=close
export OMP_PLACES=cores

ROOT=$SCRATCH/dpdm/knobs/nx4
cd $ROOT && /work/09218/gab97/ls6/entity/build-dpdm-res/src/entity.xc -input fast.toml -resume > fast_resume.stdout 2>&1
echo "exit=$? $(date -Is)"
