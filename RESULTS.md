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
| `output/benchmark/benchmark_n30_s13.csv`, `benchmark_n30_s19.csv` | state-based timings and memory for n = 30, systems 13 and 19 (one cost pair, on a node with 512 GB), written by `scripts/benchmark.jl` |
| `output/benchmark/benchmark_engines.csv` | timings of every step of the r-out-of-n and count-based computations for all systems, written by `scripts/benchmark.jl` |
| `output/paper_numbers.txt` | Tables 3, 4 and 5 and the statistics of Sections 7.5 to 7.7 of the paper, written from these files by `scripts/paper_numbers.jl` |
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

The row files and the benchmarks were computed in October 2026 on the HPC
cluster of Universidad Adolfo Ibáñez (nodes with two Intel Xeon E5-2650
processors and 128 GB of RAM; the n = 30 benchmarks on a node with four Intel
Xeon E5-4650 processors and 512 GB), with Julia 1.10.4, one configuration per
core. The code that computes them is that of release v1.0.0: the parts that
write `output/ltmc` have not changed since commit 83b5a10, and those that write
`output/optimal` since commit 0af339e; `benchmark_engines.csv` was computed with
commit 83b5a10, and the n = 30 benchmarks with commit 85c1299, whose state-based
pass is the same as in v1.0.0. The summaries and `paper_numbers.txt` were
written from the row files by the scripts of v1.0.0. `logs/` holds the logs of
the Slurm jobs of these runs, and of the computation of the signatures in
`data/signatures/` (September 2026).

The optimal state-based policies of the two networks with n = 30 (systems 13
and 19) were not computed in all configurations and are not part of the paper;
in their row files `ltmc_state` is `NaN` and `iters_state` is 0.
