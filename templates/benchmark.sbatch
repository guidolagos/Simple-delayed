#!/bin/bash
# Runtime benchmark for the paper (scripts/benchmark.jl): one job, one core, timing
# every step of the r-out-of-n, count-based and state-based computations on the
# rows LFMO-C beta = 0.2 with delays none, exp 0.1, det 0.1 and insp 0.1 of every
# system (state-based for n <= 22 by default). Writes output/benchmark/benchmark.csv.
# Takes about one hour. Run it when the cluster is quiet, and keep the whole node
# type fixed (compute-slim) so that the times are comparable.
#
# Submit from the repository root:
#   sbatch templates/benchmark.sbatch
# Larger state spaces (one cost pair, e.g. system 17 with n = 27; see README):
#   sbatch --mem=40G --time=08:00:00 --export=ALL,EXTRA="--systems=17 --delays=det:0.1 --pairs=1.0:0.1 --state-max-n=27 --out=output/benchmark/benchmark_n27.csv" templates/benchmark.sbatch
#SBATCH --job-name=benchmark
#SBATCH --partition=compute-slim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out

module purge
module load softwares
module load Julia/1.10.4-linux-x86_64

set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"

echo "host: $(hostname)  job: ${SLURM_JOB_ID}  start: $(date -Is)"
lscpu | grep -E "Model name|Socket|Core\(s\)|Thread\(s\)|MHz" || true
julia --project=. scripts/benchmark.jl ${EXTRA:-}
echo "end: $(date -Is)"
