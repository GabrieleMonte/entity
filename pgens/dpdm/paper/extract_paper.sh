#!/bin/bash
#SBATCH -J dpdm_extract                      # Job name
#SBATCH -o sbatch_logs/out/extract.o%j       # Name of stdout output file
#SBATCH -e sbatch_logs/err/extract.e%j       # Name of stderr error file
#SBATCH -p vm-small                          # CPU only; reading .bp needs no GPU.
#                                            # vm-small over development: development has
#                                            # MaxJobsPU = 1, so this queues behind any other
#                                            # dev job, and its 2 h cap is tight once a SCHEMA
#                                            # bump forces every product to be rebuilt.
#SBATCH -N 1                                 # Total # of nodes
#SBATCH -n 1                                 # Total # of tasks
#SBATCH -t 04:00:00                          # headroom for a full rebuild
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
# Walltime: two jobs (3443843, 3445746) hit the 2 h wall before the cause was
# measured rather than guessed at. It is not the volume of data. nt2.Data() pays
# a cost proportional to how many .bp are IN the series, and every isel(t=i)
# inherits it. Measured on `slow` (540 field + 500 particle dumps):
#
#     full series             open 1049 s,  ~43 s per particle load
#     8-dump symlink subset   open  9.2 s,   1.2 s per particle load
#
# So the script now symlinks just the dumps it wants into a scratch directory and
# opens that. `slow`'s phase stage went from 90 min unfinished to 25 s for 8
# frames. FIELD_FRAMES (250) and PHASE_FRAMES (60) set how many are kept; both
# are finer than a plot can show and both are tunable from the environment.
#
# Resume is per PRODUCT and schema-versioned: a file written by an older script
# is rebuilt rather than skipped, so a partial job can always be resubmitted.
#
# Each stage prints its own timing -- open/read per stage -- so the log says
# where the time went.
# Output: $SCRATCH/dpdm_npy/  (~50 MB), plus a tarball to scp home.
# PHASE_FRAMES=30 sbatch extract_paper.sh  halves the particle work if needed.

set -u

export DPDM_ROOT=$SCRATCH/dpdm            # every run lives under here
export OUT_DIR=$SCRATCH/dpdm_npy
mkdir -p $OUT_DIR

# No staging: extract_paper.py addresses each run by (path, name) from its RUNS
# table. The old symlink tree keyed on basename could not represent the re-runs,
# since paper/fast and paper_shear/fast are both named "fast".

# SKIP a run that is still writing:  SKIP=slow_s sbatch extract_paper.sh
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
