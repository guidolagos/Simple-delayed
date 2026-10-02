#!/bin/bash
# Optimal policies, main run: every row of data/grid.csv, best r-out-of-n and
# optimal count-based policies for all systems, state-based for n <= 14.
# Output: output/optimal/row_<id>.csv. Finished rows are skipped, so
# re-submitting is safe. Keep NCHUNKS fixed; --array must be 1-NCHUNKS.
#
# Submit from the repository root, after precompiling once (see README):
#   sbatch templates/optimal_array.sbatch
# Pilot (first 40 rows, i.e. system 1, in chunk 1 only):
#   sbatch --array=1 --export=ALL,EXTRA="--rows=1:40" templates/optimal_array.sbatch
#SBATCH --job-name=optimal
#SBATCH --partition=compute-slim
#SBATCH --array=1-20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out

NCHUNKS=20

module purge
module load softwares
module load Julia/1.10.4-linux-x86_64

set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"

echo "host: $(hostname)  job: ${SLURM_JOB_ID}  task: ${SLURM_ARRAY_TASK_ID}/${NCHUNKS}  start: $(date -Is)"
julia --project=. scripts/run_optimal.jl "${SLURM_ARRAY_TASK_ID}" "${NCHUNKS}" --state-max-n=14 ${EXTRA:-}
echo "end: $(date -Is)"
