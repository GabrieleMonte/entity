#!/bin/bash
#SBATCH -J dpdm_early
#SBATCH -o sbatch_logs/out/early.o%j
#SBATCH -e sbatch_logs/err/early.e%j
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:45:00
#SBATCH -A PHY23028
#SBATCH --mail-type=all
#SBATCH --mail-user=montefalcone@utexas.edu
set -u
export DPDM_ROOT=$SCRATCH/dpdm
export OUT_DIR=$SCRATCH/dpdm_npy
python3 extract_energy_early.py
echo "early exit=$? $(date -Is)"
