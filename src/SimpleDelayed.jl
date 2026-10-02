module SimpleDelayed

# Paths of the data shipped with the repository.
const REPO_DIR = normpath(joinpath(@__DIR__, ".."))
const DATA_DIR = joinpath(REPO_DIR, "data")
const SYSTEMS_DIR = joinpath(DATA_DIR, "systems")
const SIGNATURES_DIR = joinpath(DATA_DIR, "signatures")

include("systems.jl")
include("signatures.jl")
include("models.jl")
include("ltmc.jl")
include("grid.jl")
include("optimal.jl")

# systems
export TwoTerminalSystem, ncomponents, nnodes, load_systems, get_system, parse_yeh_raw
# signatures
export CompiledSystem, compile_system, works, state_from_working,
       count_working_states, tail_signature, signature, check_counts,
       signature_path, write_signature_file, read_signature_file

# LFMO models, LTMC of the delayed r-out-of-n policy, experiment grid
export MODELS, LFMOModel, lfmo_model, psi, transition_rates, mean_batch,
       Delay, policy_quantities, policy_quantities_r0, ltmc, COST_RATIOS, COST_SETUPS,
       make_grid, write_grid, read_grid, run_row, row_output_path, read_row_output,
       GRID_PATH, OUTPUT_DIR, ENGINE_VERSION
# optimal count-based and state-based policies (sec:optimal of the paper)
export Costs, repair_cost, failure_cost, call_allowed, CallTables, call_tables,
       CountModel, count_pass, dinkelbach, is_threshold, policy_string,
       StateSpace, working_counts, StatePre, state_precompute, state_pass,
       state_calls_string, mixed_levels, rout_quantities, run_optimal_row, read_optimal_output,
       OPT_RATIOS, OPT_SETUPS, OPTIMAL_DIR, OPTIMAL_LARGE_DIR, OPTIMAL_HEADER, OPTIMAL_VERSION

end # module
