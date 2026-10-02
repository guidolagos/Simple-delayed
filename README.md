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
carry `# engine_version = 2`; the scripts refuse files of version 1 (no r = 0).
In the summary, r = 0 is admissible only when setup > 0. The cost enters only afterwards: with
c_sys = 1 and c_cmp(j) = c_setup + j c_cmp,

    LTMC(r) = (c_setup + c_cmp e_nfail + c_sys p_fail) / e_trep,   c_cmp = ratio c_sys / n,

for every ratio = n c_cmp/c_sys in `COST_RATIOS` (1/10, 1/5, 1/2, 2/3, 1,
4/3, 2, 4, 10), computed by `scripts/summarize.jl`. The summary also adds a
fixed cost per repair visit, c_cmp(j) = c_setup + j c_cmp with
setup = c_setup/c_sys in `COST_SETUPS` (0, 0.05, 0.1, 0.2, 0.5): without it,
repairing at every single failure (r = 1) is almost always optimal. The cost
grid costs nothing to extend.

Method and validation (`src/ltmc.jl`): Markov chain of the number of failed
components, frozen at k failures for each k-out-of-n system of the signature
decomposition; all solves are triangular substitutions with nonnegative terms
(no cancellation), and the transition rates are computed in 512-bit BigFloat.
`test/reference/ltmc_reference.py` evaluates the paper's closed forms
(thm:LFMO, thm:determinstic, thm:LFMOsch) in 80-digit arithmetic and checks
them against a Monte Carlo simulation of the actual networks;
`test/reference/ltmc_reference.csv` (43 cases, with r = 0 for the exp and det
ones) is what `Pkg.test()` compares the Julia engine with (relative tolerance 1e-10).

### Run on the cluster

    git pull
    module purge; module load softwares; module load Julia/1.10.4-linux-x86_64
    srun --partition=compute-slim --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:30:00 \
        julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); Pkg.test()'
    # pilot: 20 rows of chunk 1, check the log and the seconds in the output headers
    sbatch --array=1 --export=ALL,EXTRA="--rows=1:400" templates/run_array.sbatch
    # everything (20 chunks); rows already done are skipped
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
  alg:state fixed), which visits the 2^n states with O(n²) work each: for
  n ≤ 14 in the main run, for 15 ≤ n ≤ 22 (systems 7, 16, 18) in a separate
  large run; not computed for n ≥ 24 (systems 13, 14, 17, 19, 20).

The `none` delay is included as continuous inspection with D = 0 (calling the
crew is repairing at once, and it is not allowed with no failed component).
Each row writes `output/optimal/row_<id>.csv` (or `output/optimal_large/`), one
line per cost pair with `best_r, ltmc_best_r, ltmc_corrective, ltmc_count,
count_policy` (`C` call, `.` continue, `-` unreachable, for f = 0..n-1),
`count_threshold` (whether the count policy is an r-out-of-n one), `ltmc_state`,
`state_mixed_levels` (levels f where the state policy calls at some but not all
working states) and `state_calls` (calls/working states per f). Every row checks
state ≤ count ≤ best r-out-of-n.

Validation: `test/reference/optimal_reference.py` is an independent mpmath
implementation that follows the paper literally (recursion eq:ccall exp, closed
forms of prop:after call det and lem:period, sub-state sums enumerated
directly), checks the count-based optimum by brute force over all stopping sets,
and checks the computed policies by Monte Carlo simulation of the networks;
`optimal_reference.csv` (64 cases) is compared by `Pkg.test()`, which also checks
that forcing r-out-of-n actions in the passes reproduces the LTMC engine, that
[Exp] with continuous inspection has a threshold optimum, and the exponential
delay recursion. `optimal_port.py` is a Python port of `src/optimal.jl` used to
validate the design (agrees with the reference to 1e-15, also at n = 30 against
60-digit closed forms).

### Run on the cluster

    git pull
    module purge; module load softwares; module load Julia/1.10.4-linux-x86_64
    srun --partition=compute-slim --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:30:00 \
        julia --project=. -e 'using Pkg; Pkg.precompile(); Pkg.test()'
    # 1. LTMC grid again (adds r = 0); old outputs of engine version 1 must go
    mv output/ltmc output/ltmc_v1
    sbatch templates/run_array.sbatch
    # 2. optimal policies, main run (state-based for n <= 14)
    sbatch --array=1 --export=ALL,EXTRA="--rows=1:40" templates/optimal_array.sbatch   # pilot
    sbatch templates/optimal_array.sbatch
    # 3. state-based policies for n = 21, 22: time one row first (system 18 is rows 8433-8928)
    srun --partition=compute-slim --ntasks=1 --cpus-per-task=1 --mem=2G --time=02:00:00 \
        julia --project=. scripts/run_optimal.jl 1 1 --rows=8433:8433 --min-n=15 --max-n=22 \
        --state-max-n=22 --out=large --verbose
    sbatch templates/optimal_large.sbatch      # adjust --time to the pilot first
    # checks
    grep -l "failed [1-9]" logs/ltmc_*.out logs/optimal*.out
    ls output/ltmc | wc -l; ls output/optimal | wc -l; ls output/optimal_large | wc -l   # 9920, 9920, 1488

Then, on your own machine:

    rsync -av uai:Simple-delayed/output/ltmc/ output/ltmc/
    rsync -av uai:Simple-delayed/output/optimal/ output/optimal/
    rsync -av uai:Simple-delayed/output/optimal_large/ output/optimal_large/
    julia --project=. scripts/summarize.jl            # output/summary.csv
    julia --project=. scripts/summarize_optimal.jl    # output/optimal_summary.csv

## Layout
- `src/` reusable code (systems, signatures, LFMO models, LTMC engine, experiment grid, optimal policies)
- `scripts/` entry points
- `templates/` Slurm job templates
- `data/systems/` system definitions; `data/signatures/` computed signatures; `data/grid.csv` experiment grid (all committed)
- `test/reference/` independent reference values used by the tests
- `output/` results (git-ignored)
