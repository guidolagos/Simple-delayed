# Results of Simple-delayed v1.0.0

This file is the `README.md` of the results archive
`Simple-delayed-v1.0.0-results.tar.gz`, deposited on Zenodo
(https://doi.org/10.5281/zenodo.23841153) together with the source code of release v1.0.0 of
https://github.com/guidolagos/Simple-delayed.
These are the results that underlie the numbers of the paper
G. Lagos and H. Olivero, "Simple repair policies with delays for systems with
simultaneous failures: exact costs and the value of knowing which components
failed". They are released under CC BY 4.0.

## Contents

The archive extracts to a directory `Simple-delayed-v1.0.0-results/` with the
layout of the repository, so that its `output/` and `logs/` directories can be
copied into a clone of the repository and read by the scripts there.

| path | content |
|---|---|
| `output/ltmc/row_00001.csv` ... `row_09920.csv` | cost-free quantities of all r-out-of-n policies, one file per row of `data/grid.csv` (`# engine_version = 3`) |
| `output/optimal/row_00001.csv` ... `row_09920.csv` | best r-out-of-n, optimal count-based and (for n ≤ 27) optimal state-based policies for the 36 cost pairs, one file per row (`# optimal_version = 3`) |
| `output/optimal_summary.csv`, `output/optimal_rows.csv` | written from `output/optimal` by `scripts/summarize_optimal.jl` |
| `output/benchmark/` | timings of every step written by `scripts/benchmark.jl`; the files with `n30` in their name hold the state-based timings and memory for n = 30 (one cost pair, on a node with 512 GB) |
| `logs/` | Slurm logs of the runs that produced these files |
| `SHA256SUMS` | SHA-256 checksums of all the files above |

The columns of the row files are described in the repository's `README.md`
(sections "LTMC experiments" and "Optimal policies"). The lines starting with
`#` record the configuration, the engine version, the largest difference
between the closed forms and the chain engine (`engine_check_*`), the
computing times and the compute node.

`output/summary.csv` is not included because it is large; it is rebuilt from
`output/ltmc` with `julia --project=. scripts/summarize.jl`.

## Provenance

All files were computed in October 2026 on the HPC cluster of Universidad
Adolfo Ibáñez (nodes with two Intel Xeon E5-2650 processors and 128 GB of RAM;
the n = 30 benchmark on a node with four Intel Xeon E5-4650 processors and
512 GB), with Julia 1.10.4, one configuration per core. The code that computes
them is that of release v1.0.0; the parts that write `output/ltmc` have not
changed since commit 83b5a10, those that write `output/optimal` since commit
0af339e, and the n = 30 benchmark was computed with commit 85c1299, whose
state-based pass is the same as in v1.0.0.

The optimal state-based policies of the two networks with n = 30 (systems 13
and 19) were not computed in all configurations and are not part of the paper;
in their row files `ltmc_state` is `NaN` and `iters_state` is 0.
