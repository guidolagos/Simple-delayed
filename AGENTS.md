# Instructions for coding agents

- Language: Julia 1.10. Activate the project with `--project=.`.
- Put reusable logic in `src/`; keep `scripts/` as thin wrappers.
- Run tests with `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Never commit anything in `output/` or `logs/`.
- Never run heavy computations on the cluster login node; use `sbatch`.
- Randomness: use `StableRNGs`, seeded from the `seed` column of the `data/grid.csv` row being processed. Never seed from the Slurm task ID or array index, because one task processes a chunk of rows and chunking can change between runs.
- Each grid row writes one output file atomically (write to a temp file, then rename) and is skipped if that file already exists, so jobs are resumable.
- Systems are two-terminal networks whose components are the edges (`data/systems/<family>_edges.csv`); component `e` is bit `e-1` of a state, and bit = 1 means working. Do not reorder edges: signature files and later experiments rely on this numbering.
- Signatures live in `data/signatures/<name>.csv` and are committed. They are computed by `scripts/compute_signatures.jl` (cluster: `templates/signatures.sbatch`) and validated by `scripts/check_signatures.jl` against `test/reference/`.
- Engines: the production values come from the closed forms of the paper in BigFloat (`src/closed_forms.jl`); the survivor-count chain engine (`*_chain` in `src/ltmc.jl` and `src/optimal.jl`) is an independent Float64 check whose largest difference is recorded in every row file. Do not evaluate the closed forms in Float64: they lose about n log10(3) digits.
- LTMC row files carry `# engine_version`; bump `ENGINE_VERSION` in `src/grid.jl` whenever the content of `output/ltmc` changes, and `OPTIMAL_VERSION` in `src/optimal.jl` for `output/optimal*`.
- `test/reference/*.py` are independent Python implementations (mpmath) used only to generate the reference CSV files compared by `Pkg.test()`; they are not run on the cluster. `optimal_port.py` and `ctmc_port.py` are Float64 ports of the Julia algorithms used to validate designs.
