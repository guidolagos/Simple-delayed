#!/bin/bash
# Optimal policies, large run: state-based policies for the systems with
# 15 <= n <= 22 (systems 7 and 16 with n = 21, system 18 with n = 22: 1488 rows).
# Each row also recomputes the best r-out-of-n and count-based results.
# Output: output/optimal_large/row_<id>.csv (scripts/summarize_optimal.jl prefers
# these files). Finished rows are skipped, so re-submitting is safe.
# Keep NCHUNKS fixed; --array must be 1-NCHUNKS.
#
# Cost: one state-based pass visits the 2^n states with O(n^2) work each, and each
# row runs Dinkelbach for 36 cost pairs. Run the pilot of the README first and
# adjust --time (and NCHUNKS) to its timing.
#
# Submit from the repository root:
#   sbatch templates/optimal_large.sbatch
#SBATCH --job-name=optimal-large
#SBATCH --partition=compute-slim
#SBATCH --array=1-150
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.out

NCHUNKS=150

module purge
module load softwares
module load Julia/1.10.4-linux-x86_64

set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"

echo "host: $(hostname)  job: ${SLURM_JOB_ID}  task: ${SLURM_ARRAY_TASK_ID}/${NCHUNKS}  start: $(date -Is)"
julia --project=. scripts/run_optimal.jl "${SLURM_ARRAY_TASK_ID}" "${NCHUNKS}" \
    --min-n=15 --max-n=22 --state-max-n=22 --out=large ${EXTRA:-}
echo "end: $(date -Is)"
