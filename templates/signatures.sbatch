#!/bin/bash
# Structural signatures of the 22 systems of Yeh et al. (2024).
# One array task per system; task id = system id (1..22).
#
# Submit from the repository root:
#   sbatch templates/signatures.sbatch
# Re-submitting is safe: systems whose data/signatures/<name>.csv exists are skipped.
# To run only some systems:  sbatch --array=13,19,21 templates/signatures.sbatch
#
# Cost: 2^n connectivity checks. n <= 24 takes seconds; the three n = 30
# systems (13, 19, 21) take a few minutes on one core.
#SBATCH --job-name=signatures
#SBATCH --partition=compute-slim
#SBATCH --array=1-22
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%A_%a.out

module purge
module load softwares
module load Julia/1.10.4-linux-x86_64

set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"

echo "host: $(hostname)  job: ${SLURM_JOB_ID}  task: ${SLURM_ARRAY_TASK_ID}  start: $(date -Is)"
julia --project=. --threads="${SLURM_CPUS_PER_TASK:-1}" \
    scripts/compute_signatures.jl "${SLURM_ARRAY_TASK_ID}"
echo "end: $(date -Is)"
