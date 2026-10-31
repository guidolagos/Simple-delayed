# Simple-delayed

Code and results of the computational experiments of the paper

> G. Lagos and H. Olivero, "Simple repair policies with delays for systems with
> simultaneous failures: exact costs and the value of knowing which components
> failed", submitted to *Reliability Engineering & System Safety* (2026).

The paper gives closed-form long-run mean costs of r-out-of-n repair policies
with a delay between the call to the repair crew and the repair (an exponential
or a constant delay, or a fixed inspection schedule), for semi-coherent systems
whose components fail according to a Lévy-frailty Marshall-Olkin (LFMO) model,
which allows simultaneous failures. It also computes the optimal policies that
observe only the number of failed components (count-based) and those that
observe which components have failed (state-based), to measure how much an
r-out-of-n policy loses by ignoring which components have failed.

This repository computes, for the 20 benchmark networks of Yeh et al. (2024),
31 failure models and 16 delays (9,920 configurations):

- the structural signatures of the networks (`data/signatures/`);
- the cost-free quantities of all r-out-of-n policies, from which the long-run
  mean cost follows for any cost pair (`output/ltmc/`);
- the best r-out-of-n policy and the optimal count-based policy for 36 cost
  pairs in all configurations, and the optimal state-based policy in the 8,928
  configurations of the 18 networks with at most 27 components (`output/optimal/`);
- the computing times of every step (`output/benchmark/` and the row files).

## Archived version and citation

Version v1.0.0 of this repository and the results that underlie the numbers of
the paper are archived together on Zenodo (see `CITATION.cff`). If you use the
code or the results, please cite the archive and the paper. The results are
not in this repository (`output/` and `logs/` are git-ignored); `RESULTS.md`
describes the archived results file.

## License

The code is released under the MIT License (`LICENSE`). The data files
(`data/`) and the archived results are released under the Creative Commons
Attribution 4.0 International License (CC BY 4.0). The networks in
`data/systems/` are taken from Figure 7 of W.-C. Yeh, G. Yang, C.-L. Huang,
"A New Node-Based Algorithm for Identifying the Complete Minimal Cut Set",
*Electronics* 13 (2024) 603; see `data/systems/README.md`.

## Requirements

- Julia 1.10. The experiments of the paper were run with Julia 1.10.4.
  The package `SimpleDelayed` uses only Julia's standard library (`BigFloat`,
  through the MPFR library bundled with Julia), so `Manifest.toml` lists no
  external packages.
- Python 3 with `mpmath`, only to regenerate the reference values in
  `test/reference/` (not needed to run the experiments or the tests).
- For the full experiments, a Slurm cluster (`templates/`); the full state-based
  computations need up to about 14 GB of memory per configuration for 27
  components, and about 1,800 CPU hours in total (see the paper).

## Quick start on a laptop

    git clone https://github.com/guidolagos/Simple-delayed.git
    cd Simple-delayed
    julia --project=. -e 'using Pkg; Pkg.instantiate()'
    julia --project=. -e 'using Pkg; Pkg.test()'

The tests compare both engines (below) with independent reference values and
with each other on networks with up to 30 components, so they take a while.

The first 16 rows of `data/grid.csv` are network 1 (5 components) under
independent exponential lifetimes ([Exp]) with the 16 delays. To compute them:

    julia --project=. scripts/run_experiment.jl 1 1 --rows=1:16
    julia --project=. scripts/run_optimal.jl 1 1 --rows=1:16

This writes `output/ltmc/row_00001.csv`, ..., `row_00016.csv` and
`output/optimal/row_00001.csv`, ..., `row_00016.csv`, which can be compared
with the same files in the archived results (the lines starting with `#`
record times, hosts and engine checks, and differ between runs). Any other set
of rows can be selected with `--rows=a:b`; network k owns rows
(k-1)·496+1 to k·496.

## Reproducing the experiments of the paper

The sections below describe each step: signatures, r-out-of-n quantities,
optimal policies and runtime benchmark. The Slurm templates in `templates/`
are written for the HPC cluster of Universidad Adolfo Ibáñez (Rocky Linux,
Slurm, Julia loaded with `module load softwares; module load
Julia/1.10.4-linux-x86_64`, partitions `compute-slim` and `compute-fat`);
on other clusters, change the `module` lines and the partitions.

Set up the package once through `srun`, never on the login node:

    module purge; module load softwares; module load Julia/1.10.4-linux-x86_64
    srun --partition=compute-slim --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:30:00 \
        julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); Pkg.test()'

Precompiling once before submitting an array job avoids every array task
compiling the package at the same time.

## From the results to the numbers of the paper

| results | script | used in the paper |
|---|---|---|
| `data/signatures/*.csv` | `scripts/compute_signatures.jl`, `scripts/check_signatures.jl` | Section 7.1 (systems) |
| `output/ltmc/row_*.csv` | `scripts/summarize.jl` → `output/summary.csv` (best trigger and its saving over corrective maintenance, per row and cost pair) | Section 7.5, Table 3 |
| `output/optimal/row_*.csv` | `scripts/summarize_optimal.jl` → `output/optimal_summary.csv` (per row and cost pair: best r-out-of-n, count- and state-based costs, and the savings `gain_*`) | Section 7.6, Table 4 |
| `output/optimal/row_*.csv` | `scripts/summarize_optimal.jl` → `output/optimal_rows.csv` (per row: times of every step, passes, peak memory) | Section 7.7, Table 5; Supplementary Material, computational setup and computational effort |
| `output/benchmark/*n30*.csv` | `scripts/benchmark.jl` | state-based times and memory for n = 30 (Section 7.7, Table 5; Supplementary Material) |
| all row files | `scripts/engine_check.jl` | agreement between the closed forms and the chain engine |

The statistics of the paper (medians, percentiles, maxima and shares) are taken
over the rows and cost pairs stated in each table caption, e.g. Table 4 gives the
median and the 90th percentile of `gain_state_vs_r` over the 18 networks with
n ≤ 27, the models with a given β, the delays of a given type and the 36 cost
pairs with `setup` > 0, and the maximum of `gain_count_vs_r` over all networks.
`output/summary.csv` is large and is not archived: `scripts/summarize.jl`
rebuilds it from `output/ltmc`.

## Paper references

The comments in the code and the text below refer to the paper by the labels of
its LaTeX source. In the submitted manuscript (v1.0.0 of this repository):

| label | in the paper | label | in the paper |
|---|---|---|---|
| `thm:repairdelay` | Theorem 7 | `sec:LFMO` | Section 5 |
| `rmk:call0` | Remark 4 | `def:lambda` | Definition 23 |
| `def:policyindep`, `def:policydelta` | Definitions 5, 6 | `sec:optimal` | Section 6 |
| `thm:LFMOrepairdelay` | Theorem 14 | `alg:dinkelbach` | Algorithm 1 |
| `lem:pbbs and exps` | Lemma 16 | `rmk:count vs rout` | Remark 26 |
| `thm:LFMO`, `thm:determinstic`, `thm:LFMOsch` | Theorems 18, 19, 20 | `sec:computational` | Section 7 |

## Structural signatures of the test systems

The systems of Yeh et al. (2024) (Figure 7 of their paper) are two-terminal networks: the **components are the
edges**, and the system works iff there is a path of working edges from
`Start` to `End`. We use the 20 benchmark networks of Yeh et al. (2024),
all undirected, as in their computations. `data/systems/README.md` documents
the differences with the raw file (systems 21 and 22 dropped as duplicates;
systems 6 and 7 read as undirected) and the checks against Table 9 of Yeh et al.

Data files (`data/systems/`):

- `yeh2024_raw.txt`: the original edge lists (`SistemasYeh_corrected.txt`),
  kept for provenance.
- `yeh2024_systems.csv`: one row per system:
  `system, name, label, directed, source, target, n_components, n_nodes`.
- `yeh2024_edges.csv`: one row per component: `system, component, from, to`.
  Component `e` of a system is its `e`-th edge in the raw file; in a state,
  bit `e-1` is 1 iff component `e` works.

The test suite checks that the CSV files coincide with the raw file (except
for the documented differences) and that the numbers of minimal cuts agree
with Table 9 of Yeh et al. (2024).

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

## LTMC experiments

The experiment grid `data/grid.csv` (committed; regenerate with
`julia --project=. scripts/make_grid.jl`) has one row per

- system: the 20 Yeh et al. (2024) systems;
- failure model: `Exp` (β = 0) and `LFMO-C`, `LFMO-G`, `LFMO-IG`, `LFMO-S`,
  `BFR`, `LSh`, each calibrated to Ψ(1) = 1 and pairwise beta factor
  β = 2 − Ψ(2) ∈ {0.01, 0.05, 0.1, 0.2, 0.5} (`src/models.jl`);
- delay between the trigger T_{r:n} and the repair start: `none`, or `exp`
  (exponential with that mean), `det` (constant) or `insp` (inspections every
  that many time units), each with value ∈ {0.01, 0.03, 0.1, 0.3, 1}
  (time unit = mean component lifetime);

that is 20 × 31 × 16 = 9920 rows. Each row writes
`output/ltmc/row_<row_id>.csv` with, for every trigger r = 1..n,
`p_fail = P(T_rep = T_fail)`, `e_nfail = E[N(T_rep)]` and `e_trep = E[T_rep]`
(Theorem thm:repairdelay of the paper), plus r = 0 (the crew is called right
after each repair; Remark rmk:call0) for the `exp` and `det` delays. Row files
carry `# engine_version = 3` (closed-form engine, below); the scripts refuse
files of versions 1 (no r = 0) and 2 (chain engine).
In the summary, r = 0 is admissible only when setup > 0. The cost enters only afterwards: with
c_sys = 1 and c_cmp(j) = c_setup + j c_cmp,

    LTMC(r) = (c_setup + c_cmp e_nfail + c_sys p_fail) / e_trep,   c_cmp = ratio c_sys / n,

for every ratio = n c_cmp/c_sys in `COST_RATIOS` (1/10, 1/5, 1/2, 2/3, 1,
4/3, 2, 4, 10), computed by `scripts/summarize.jl`. The summary also adds a
fixed cost per repair visit, c_cmp(j) = c_setup + j c_cmp with
setup = c_setup/c_sys in `COST_SETUPS` (0, 0.05, 0.1, 0.2, 0.5): without it,
repairing at every single failure (r = 1) is almost always optimal. The cost
grid costs nothing to extend.

Method (`src/closed_forms.jl`): the closed forms of the paper (Lemma
lem:pbbs and exps, Theorems thm:LFMO, thm:determinstic, thm:LFMOsch, and the
signature decomposition thm:repairdelay), evaluated in BigFloat with
`CLOSED_PRECISION` = 256 bits and rounded to Float64. They are alternating
binomial sums whose terms add up, in absolute value, to about 3^n times the
result, so Float64 is not enough (errors up to ~10% at n = 30); 256 bits leave
more than 200 bits for n ≤ 30 (`Pkg.test()` checks 256 against 512 bits).

Check (`src/ltmc.jl`, `*_chain`): an independent Float64 engine on the Markov
chain of the number of failed components, frozen at k failures for each
k-out-of-n system; all solves are triangular substitutions with nonnegative
terms (resolvent for exponential delays, uniformization for constant delays and
inspection periods). Every row file records the largest difference between both
engines (`engine_check_prob`: absolute, probabilities; `engine_check_rel`:
relative, expectations), and a row fails if either exceeds 1e-9;
`scripts/engine_check.jl` collects them.

Validation: `test/reference/ltmc_reference.py` evaluates the paper's closed forms
in 80-digit arithmetic and checks them against a Monte Carlo simulation of the
actual networks; `test/reference/ltmc_reference.csv` (43 cases, with r = 0 for
the exp and det ones) is what `Pkg.test()` compares both Julia engines with
(relative tolerance 1e-10). `Pkg.test()` also compares both engines with each
other for all delays on systems with up to n = 30 components (1e-12).

### Run on the cluster

    git pull
    module purge; module load softwares; module load Julia/1.10.4-linux-x86_64
    srun --partition=compute-slim --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:30:00 \
        julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); Pkg.test()'
    # pilot: 20 rows of chunk 1, check the log and the seconds in the output headers
    sbatch --array=1 --export=ALL,EXTRA="--rows=1:400" templates/run_array.sbatch
    # everything (20 chunks); rows already done are skipped. After an
    # ENGINE_VERSION change, move the old outputs away first:
    #   mv output/ltmc output/ltmc_v2
    sbatch templates/run_array.sbatch
    squeue -u $USER
    grep -l "failed [1-9]" logs/ltmc_*.out     # any failed rows?
    ls output/ltmc | wc -l                     # 9920 when complete

Then, on your own machine:

    rsync -av uai:Simple-delayed/output/ltmc/ output/ltmc/
    julia --project=. scripts/summarize.jl     # writes output/summary.csv

## Optimal policies (Section sec:optimal of the paper)

For every grid row and every cost pair (ratio in `OPT_RATIOS` = `COST_RATIOS`,
setup in `OPT_SETUPS` = 0.05, 0.1, 0.2, 0.5; setup > 0 because calling the crew
before any failure must cost something, Remark rmk:call0), `src/optimal.jl`
computes

- the best r-out-of-n policy, r ∈ {0..n} for `exp` and `det`, {1..n} otherwise;
- the optimal **count-based** policy (Dinkelbach, alg:dinkelbach, with the passes
  alg:count exp, its deterministic version and alg:count fixed), for all systems;
- the optimal **state-based** policy (alg:state exp, alg:state det,
  alg:state fixed), which visits the 2^n states with O(n²) work each, for the 18
  systems with n ≤ `OPT_STATE_MAX_N` = 27; not for n = 30 (systems 13, 19).
  The 32 [Exp] rows of the two n = 30 systems (5953-5968, 8929-8944) were run
  separately on compute-fat (~1-4 days and ~107 GB per row) with
  `--state-max-n=30 --out=output/optimal_n30`, which keeps them out of
  `output/optimal`.

The `none` delay is included as continuous inspection with D = 0 (calling the
crew is repairing at once, and it is not allowed with no failed component).
Each row writes `output/optimal/row_<id>.csv`, one line per cost pair with
`best_r, ltmc_best_r, ltmc_corrective, ltmc_count, iters_count, seconds_count,
count_policy` (`C` call, `.` continue, `-` unreachable, for f = 0..n-1),
`count_threshold` (whether the count policy is an r-out-of-n one), `ltmc_state`,
`iters_state`, `state_mixed_levels` (levels f where the state policy calls at
some but not all working states), `state_calls` (calls/working states per f) and
`seconds_state`. Every row checks state ≤ count ≤ best r-out-of-n.

Engine. Everything the policies use comes from the closed forms in BigFloat: the
r-out-of-n quantities (`rout_quantities`), and the after-call and one-period
tables (`call_tables`: prop:after call det for the constant delay, eq:EminDT,
eq:Ddet_Tkn and eq:PNcDdetT for the exponential delay, lem:period for
inspections), which feed both the count-based and the state-based passes; the
transition rates of def:lambda are computed in 512 bits. The chain engine
(`*_chain`: resolvent and uniformization) is run afterwards in every row only as
a check (`engine_check_*`, `seconds_chain`); none of its values is used.

Timings. Every row records, in full precision, `seconds_rout` (closed forms,
computed afresh in every row, delay forms and r-out-of-n quantities),
`seconds_tables` (after-call tables and count model), `seconds_states` (state
space; 0 when cached from the previous row of the same system, see
`states_cached`), `seconds_state_pre`, per cost pair `seconds_count` and
`seconds_state` (Dinkelbach's method, `iters_*` passes), and `maxrss_mb`, the
peak memory of the process when the row is finished. Each array task first runs
four rows of system 1 in a temporary directory (`optimal_warmup`), so that the
times exclude compilation.

Versions. `# optimal_version = 3`. Version 1 (commit 85c1299) computed the tables
with the chain; version 2 (commit 83b5a10) with the closed forms, but kept the
state-based policies for n ≥ 21 of version 1 in a separate directory. Version 3
computes all of them in one run, in `output/optimal/`. Readers refuse other
versions, except `scripts/compare_optimal.jl`, which compares new and old outputs.

Validation: `test/reference/optimal_reference.py` is an independent mpmath
implementation that follows the paper literally (recursion eq:ccall exp, closed
forms of prop:after call det and lem:period, sub-state sums enumerated
directly), checks the count-based optimum by brute force over all stopping sets,
and checks the computed policies by Monte Carlo simulation of the networks;
`optimal_reference.csv` (64 cases) is compared by `Pkg.test()`. `Pkg.test()` also
checks the closed forms against the chain (all models for n = 5 and 9, all delays
up to n = 30), that the passes with stopping set {r, ..., n-1} reproduce the
r-out-of-n costs (count- and state-based, every r), that [Exp] with continuous
inspection has an optimum equivalent to a threshold policy, with the cost of the best r (all 36 cost pairs), that the exponential-delay
tables agree with the recursion eq:ccall exp, and that the row files contain
exactly the closed-form values. `optimal_port.py` is a Python port of
`src/optimal.jl` used to validate the design.

### Run on the cluster

Sizes: n ≤ 14: 12 systems, 5,952 rows; n = 21: systems 7 and 16 (rows 2977-3472,
7441-7936); n = 22: system 18 (8433-8928); n = 24: system 20 (9425-9920);
n = 26: system 14 (6449-6944); n = 27: system 17 (7937-8432); n = 30: systems 13
and 19 (5953-6448, 8929-9424), count-based only. Measured on compute-slim
(version 2/1 runs): per row, ~2 min for n = 21-22, ~20 min (n = 24), ~30 min
(26), ~3 h (27); peak memory 0.6, 0.9, 2.6, 6.8 and 13.5 GB.

    git pull
    module purge; module load softwares; module load Julia/1.10.4-linux-x86_64
    srun --partition=compute-slim --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:45:00 \
        julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); Pkg.test()'
    # old outputs out of the way (kept for scripts/compare_optimal.jl)
    mv output/optimal output/optimal_v2; mv output/optimal_large output/optimal_large_v1
    # pilots: one row per size (LFMO-C, β = 0.2, constant delay 0.1), then sacct
    sbatch --job-name=pilot-n27 --array=1 --mem=28G --time=08:00:00 --export=ALL,NMIN=27,NMAX=27,NCHUNKS=1,EXTRA="--rows=8009:8009 --verbose" templates/optimal.sbatch
    # production: see templates/optimal.sbatch; one array per size
    sbatch --job-name=opt-n27 --array=1-248%40 --mem=28G --time=16:00:00 --export=ALL,NMIN=27,NMAX=27,NCHUNKS=248 templates/optimal.sbatch
    # checks
    grep -l "failed [1-9]" logs/opt-*.out
    ls output/optimal | wc -l                                  # 9920 when complete
    julia --project=. scripts/compare_optimal.jl output/optimal_large_v1 output/optimal_v2

(the last command through `srun`). Then, on your own machine:

    rsync -av uai:Simple-delayed/output/optimal/ output/optimal/
    julia --project=. scripts/summarize_optimal.jl    # output/optimal_summary.csv, output/optimal_rows.csv

The state-based timings for n = 30 come from the benchmark below (one cost pair,
on compute-fat); the full runs would take ~1.5-2 days per row and ~106 GB.

### Runtime benchmark

`scripts/benchmark.jl` times every step separately (transition rates, closed
forms independent of the delay and depending on it, r-out-of-n quantities and
after-call tables with both engines, evaluation of the LTMC of all r for one cost
pair, one count-based pass, count- and state-based Dinkelbach, state space and
per-state precomputation) on representative rows, after a warm-up, and writes
`output/benchmark/benchmark.csv`:

    sbatch templates/benchmark.sbatch                       # all systems, state-based n <= 22, ~1 h
    # one cost pair with a larger state space (n = 30 needs compute-fat, ~115 GB)
    sbatch --partition=compute-fat --mem=250G --time=12:00:00 \
        --export=ALL,EXTRA="--systems=13 --delays=det:0.1 --pairs=1.0:0.1 --state-max-n=30 --out=output/benchmark/benchmark_n30.csv" \
        templates/benchmark.sbatch

## Layout
- `src/` reusable code (systems, signatures, LFMO models, LTMC engine, experiment grid, optimal policies)
- `scripts/` entry points
- `templates/` Slurm job templates
- `data/systems/` system definitions; `data/signatures/` computed signatures; `data/grid.csv` experiment grid (all committed)
- `test/reference/` independent reference values used by the tests
- `output/` results (git-ignored)
