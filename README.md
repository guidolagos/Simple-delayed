# Simple-delayed

Computational experiments for the paper on repair policies with delays for
semi-coherent systems under LFMO models.

## Setup
    julia --project=. -e 'using Pkg; Pkg.instantiate()'

On the cluster, run the setup (and the tests) through `srun`, never on the
login node:

    module purge; module load softwares; module load Julia/1.10.4-linux-x86_64
    srun --partition=compute-slim --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:30:00 \
        julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); Pkg.test()'

Precompiling once before submitting an array job avoids every array task
compiling the package at the same time.

## Run tests
    julia --project=. -e 'using Pkg; Pkg.test()'

## Structural signatures of the test systems

The systems of Yeh et al. (2024) (`\cite{yeh2024new}`; Figure `fig:yeh` of
the subsystems paper) are two-terminal networks: the **components are the
edges**, and the system works iff there is a path of working edges from
`Start` to `End`. Systems 6 and 7 (labels `6D`, `7D`) are directed; the
others are undirected.

Data files (`data/systems/`):

- `yeh2024_raw.txt`: the original edge lists (`SistemasYeh_corrected.txt`),
  kept for provenance.
- `yeh2024_systems.csv`: one row per system:
  `system, name, label, directed, source, target, n_components, n_nodes`.
- `yeh2024_edges.csv`: one row per component: `system, component, from, to`.
  Component `e` of a system is its `e`-th edge in the raw file; in a state,
  bit `e-1` is 1 iff component `e` works.

The test suite checks that the CSV files coincide with the raw file.

Compute the signatures on the cluster (one array task per system; existing
outputs are skipped, so re-submitting is safe):

    sbatch templates/signatures.sbatch
    squeue -u $USER                                   # wait until done
    srun --partition=compute-slim --ntasks=1 --mem=1G --time=00:10:00 \
        julia --project=. scripts/check_signatures.jl

(On your own machine, `julia --project=. scripts/check_signatures.jl` directly.)

Output: `data/signatures/<name>.csv` (e.g. `yeh2024_07.csv`), with header
comments (system, directedness, run time, host) and one row per
`k = 0..n`:

| column | meaning |
|---|---|
| `W` | number of working states with exactly `k` failed components |
| `binomial` | `binomial(n, k)` |
| `tail_signature` | `W / binomial(n, k)`, exact (`num/den`) |
| `signature` | `s_k = tail_{k-1} - tail_k`, exact; `s_0 := 0` |
| `signature_float` | `s_k` as Float64 |

`scripts/check_signatures.jl` compares every `W` column with an independent
C++ computation (`test/reference/yeh2024_W_reference.csv`) and exits with
code 1 on any discrepancy. The signature files are small and are committed.

## Run on the cluster (Slurm)
    cp templates/run_array.sbatch job.sbatch   # edit array range, time, memory
    sbatch job.sbatch

## Layout
- `src/` reusable code (systems, signatures, LTMC formulas, simulation)
- `scripts/` entry points
- `templates/` Slurm job templates
- `data/systems/` system definitions; `data/signatures/` computed signatures (committed)
- `test/reference/` independent reference values used by the tests
- `output/` results (git-ignored)
