#!/bin/bash
# Optimal policies for the largest systems: state-based policies for the rows whose
# system has NMIN <= n <= NMAX (24: system 20; 26: system 14; 27: system 17;
# 30: systems 13 and 19), 496 rows per system. Each row also recomputes the best
# r-out-of-n and count-based results. Output: output/optimal_large/row_<id>.csv,
# which scripts/summarize_optimal.jl prefers. Finished rows are skipped, so
# re-submitting is safe.
#
# Memory and time depend strongly on n, so they are NOT fixed here: give them on
# the command line, together with NMIN, NMAX and NCHUNKS (= size of --array), e.g.
#   sbatch --array=1-50 --mem=4G --time=12:00:00 \
#          --export=ALL,NMIN=24,NMAX=24,NCHUNKS=50 templates/optimal_xl.sbatch
# A pilot on one row prints the time of every Dinkelbach run:
#   sbatch --array=1 --mem=16G --time=06:00:00 \
#          --export=ALL,NMIN=26,NMAX=26,NCHUNKS=1,EXTRA="--rows=6449:6449 --verbose" templates/optimal_xl.sbatch
# Keep NCHUNKS fixed between submissions for the same NMIN/NMAX (it defines which
# rows each array task owns). See README, "Largest systems".
#SBATCH --job-name=optimal-xl
#SBATCH --partition=compute-slim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
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
    --min-n="${NMIN}" --max-n="${NMAX}" --state-max-n="${NMAX}" --out=large ${EXTRA:-}
echo "end: $(date -Is)"
