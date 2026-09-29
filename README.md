# Simple-delayed

Computational experiments for the paper on repair policies with delays for
semi-coherent systems under LFMO models.

## Setup
    julia --project=. -e 'using Pkg; Pkg.instantiate()'

## Run tests
    julia --project=. -e 'using Pkg; Pkg.test()'

## Run on the cluster (Slurm)
    cp templates/run_array.sbatch job.sbatch   # edit array range, time, memory
    sbatch job.sbatch

## Layout
- `src/` reusable code (signatures, LTMC formulas, simulation)
- `scripts/` entry points
- `templates/` Slurm job templates
- `data/` small inputs
- `output/` results (git-ignored)
