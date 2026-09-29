# Instructions for coding agents

- Language: Julia 1.10. Activate the project with `--project=.`.
- Put reusable logic in `src/`; keep `scripts/` as thin wrappers.
- Run tests with `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Never commit anything in `output/` or `logs/`.
- Never run heavy computations on the cluster login node; use `sbatch`.
- Randomness: use `StableRNGs`, seeded from the `seed` column of the `data/grid.csv` row being processed. Never seed from the Slurm task ID or array index, because one task processes a chunk of rows and chunking can change between runs.
- Each grid row writes one output file atomically (write to a temp file, then rename) and is skipped if that file already exists, so jobs are resumable.
