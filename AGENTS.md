# Instructions for coding agents

- Language: Julia 1.10. Activate the project with `--project=.`.
- Put reusable logic in `src/`; keep `scripts/` as thin wrappers.
- Run tests with `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Never commit anything in `output/` or `logs/`.
- Never run heavy computations on the cluster login node; use `sbatch`.
- Seed the RNG with the Slurm task ID for reproducibility.
