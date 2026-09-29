#!/bin/bash
#SBATCH --job-name=simple-delayed
#SBATCH --partition=compute-slim
#SBATCH --array=1-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%A_%a.out

module purge
module load softwares
module load Julia/1.10.4-linux-x86_64

julia --project=. scripts/run_experiment.jl $SLURM_ARRAY_TASK_ID
