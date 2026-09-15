#!/bin/bash
#SBATCH -J dpdm_extract                      # Job name
#SBATCH -o sbatch_logs/out/extract.o%j       # Name of stdout output file
#SBATCH -e sbatch_logs/err/extract.e%j       # Name of stderr error file
#SBATCH -p development                       # CPU only; reading .bp needs no GPU
#SBATCH -N 1                                 # Total # of nodes
#SBATCH -n 1                                 # Total # of tasks
#SBATCH -t 02:00:00                          # development cap; see note below
#SBATCH --mail-type=all                      # Send email at begin and end of job
#SBATCH -A PHY23028                          # Project/Allocation name
#SBATCH --mail-user=montefalcone@utexas.edu
#
# Reduce the six paper-faithful runs to .npz for local plotting.
#
#   cd $WORK/ls6/entity/pgens/dpdm/paper
#   mkdir -p sbatch_logs/out sbatch_logs/err
#   sbatch extract_paper.sh
#
# Why a job and not a login node: nt2.Data() reserves a 16 MB ADIOS2 buffer per
# .bp at open, and `fast` has 1000 field dumps. On a login node that dies with
# "FFS out of memory" before returning. A compute node has 256 GB.
#
# Walltime: the particle loop dominates at ~1 s per dump per species, and there
# are 3000 dump-species in total, so this is expected to take 60-90 min against
# the 2 h development cap. The runs are ordered cheapest-first and already
# extracted runs are skipped, so if it does hit the wall just resubmit -- it
# picks up where it stopped. FORCE=1 redoes everything.
#
# Output: $SCRATCH/dpdm_npy/  (~50 MB), plus a tarball to scp home.

set -u

DPDM=$SCRATCH/dpdm
STAGE=$DPDM/_export_stage                    # symlinks only; copies nothing
export PAPER_DIR=$STAGE
export OUT_DIR=$SCRATCH/dpdm_npy

# The six runs live under three different parents, and extract_paper.py wants
# them side by side as <dir>/<tag>/<tag>/. Symlink rather than copy: the script
# only reads, and copying would move 1.9 GB for no reason.
mkdir -p $STAGE $OUT_DIR
for src in paper/fast paper/slow lz/lz1p0 lz/heat \
           qstest/undriven_q qstest/undriven_r; do
  tag=$(basename $src)
  if [ -d "$DPDM/$src/$tag" ]; then
    ln -sfn $DPDM/$src $STAGE/$tag
  else
    echo "WARNING: $DPDM/$src/$tag missing, will be skipped"
  fi
done

echo "staged:"; ls -l $STAGE | sed 's/^/  /'
echo

python3 extract_paper.py
rc=$?
echo "extract exit=$rc $(date -Is)"
[ $rc -ne 0 ] && exit $rc

echo; echo "=== $OUT_DIR ==="
ls -lh $OUT_DIR
du -sh $OUT_DIR

# One tarball to fetch. This is ~50 MB, not the 1.9 GB of .bp behind it.
tar -czf $SCRATCH/dpdm_npy.tar.gz -C $SCRATCH dpdm_npy
echo
echo "download with:"
echo "  scp $USER@ls6.tacc.utexas.edu:$SCRATCH/dpdm_npy.tar.gz ."
