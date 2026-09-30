#!/bin/bash
# LTMC experiments: every row of data/grid.csv (systems x LFMO models x betas x delays).
# Each array task processes one chunk of rows (row_id mod NCHUNKS); rows whose
# output/ltmc/row_<id>.csv exists are skipped, so re-submitting is safe.
#
# Submit from the repository root, after precompiling once (see README):
#   sbatch templates/run_array.sbatch
# Pilot (a few rows of chunk 1 only):
#   sbatch --array=1 --export=ALL,EXTRA="--rows=1:400" templates/run_array.sbatch
# Re-run only some chunks:
#   sbatch --array=3,7 templates/run_array.sbatch
#
# Keep NCHUNKS fixed for a given grid (it defines which rows each chunk owns);
# the --array range must be 1-NCHUNKS.
#SBATCH --job-name=ltmc
#SBATCH --partition=compute-slim
#SBATCH --array=1-20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%A_%a.out

NCHUNKS=20

module purge
module load softwares
module load Julia/1.10.4-linux-x86_64

set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"

echo "host: $(hostname)  job: ${SLURM_JOB_ID}  task: ${SLURM_ARRAY_TASK_ID}/${NCHUNKS}  start: $(date -Is)"
julia --project=. scripts/run_experiment.jl "${SLURM_ARRAY_TASK_ID}" "${NCHUNKS}" ${EXTRA:-}
echo "end: $(date -Is)"
