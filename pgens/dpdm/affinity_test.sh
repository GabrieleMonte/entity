#!/bin/bash
#SBATCH -J dpdm_aff
#SBATCH -o sbatch_logs/out/aff.o%j
#SBATCH -e sbatch_logs/err/aff.e%j
#SBATCH -p development
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -t 00:10:00
#SBATCH -A PHY23028
# What CPU set does each MPI rank actually get? Entity only sees the cores its
# process is bound to, so this decides whether hybrid MPI+OpenMP works at all.
set -uo pipefail
module load gcc/13.2.0
export OMP_NUM_THREADS=112

probe='echo "rank=${PMI_RANK:-${MV2_COMM_WORLD_RANK:-?}} host=$(hostname -s) nproc=$(nproc) affinity=$(taskset -cp $$ 2>/dev/null | sed "s/.*: //" | cut -c1-40)"'

echo "### A: plain ibrun"
ibrun bash -c "$probe" 2>&1 | grep -E "^rank=" | sort

echo "### B: ibrun with OMP_PROC_BIND=spread OMP_PLACES=cores"
OMP_PROC_BIND=spread OMP_PLACES=cores ibrun bash -c "$probe" 2>&1 | grep -E "^rank=" | sort

echo "### C: srun --cpus-per-task=128 --cpu-bind=none"
srun -N 2 -n 2 --cpus-per-task=128 --cpu-bind=none bash -c "$probe" 2>&1 | grep -E "^rank=" | sort

echo "### D: srun --cpu-bind=none (no -c)"
srun -N 2 -n 2 --cpu-bind=none bash -c "$probe" 2>&1 | grep -E "^rank=" | sort

echo "### E: what does OpenMP actually spawn under ibrun"
ibrun bash -c 'export OMP_NUM_THREADS=112; echo "rank=${PMI_RANK:-?} omp_threads=$(OMP_NUM_THREADS=112 python3 -c "import os;print(len(os.sched_getaffinity(0)))")"' 2>&1 | grep -E "^rank=" | sort
echo "AFF_DONE"
