#!/bin/bash
# Optimal policies (scripts/run_optimal.jl): for the rows whose system has
# NMIN <= n <= NMAX, the best r-out-of-n and the optimal count-based policies,
# and the optimal state-based policies if n <= 27. Output: output/optimal/row_<id>.csv
# (one directory for all sizes). Finished rows are skipped, so re-submitting the same
# command is safe; a task killed by --time loses only the row in progress.
#
# Memory and time depend strongly on n, so give them on the command line together
# with NMIN, NMAX and NCHUNKS (= size of --array, kept fixed between submissions for
# the same NMIN/NMAX, since it defines which rows each task owns), e.g.
#   sbatch --job-name=opt-n22 --array=1-50 --mem=2G --time=02:00:00 \
#          --export=ALL,NMIN=22,NMAX=22,NCHUNKS=50 templates/optimal.sbatch
# Pilot on one row (prints the time of every state-based Dinkelbach run):
#   sbatch --job-name=pilot-n27 --array=1 --mem=28G --time=08:00:00 \
#          --export=ALL,NMIN=27,NMAX=27,NCHUNKS=1,EXTRA="--rows=8009:8009 --verbose" templates/optimal.sbatch
# See README, "Optimal policies", for the commands of the production run.
#SBATCH --job-name=optimal
#SBATCH --partition=compute-slim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out

: "${NMIN:?set NMIN with --export}"
: "${NMAX:?set NMAX with --export}"
: "${NCHUNKS:?set NCHUNKS with --export (equal to the size of --array)}"

module purge
module load softwares
module load Julia/1.10.4-linux-x86_64

set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"

echo "host: $(hostname)  job: ${SLURM_JOB_ID}  task: ${SLURM_ARRAY_TASK_ID}/${NCHUNKS}  n: ${NMIN}-${NMAX}  start: $(date -Is)"
julia --project=. scripts/run_optimal.jl "${SLURM_ARRAY_TASK_ID}" "${NCHUNKS}" \
    --min-n="${NMIN}" --max-n="${NMAX}" ${EXTRA:-}
echo "end: $(date -Is)"
