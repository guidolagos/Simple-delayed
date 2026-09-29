module SimpleDelayed

# Paths of the data shipped with the repository.
const REPO_DIR = normpath(joinpath(@__DIR__, ".."))
const DATA_DIR = joinpath(REPO_DIR, "data")
const SYSTEMS_DIR = joinpath(DATA_DIR, "systems")
const SIGNATURES_DIR = joinpath(DATA_DIR, "signatures")

include("systems.jl")
include("signatures.jl")

# systems
export TwoTerminalSystem, ncomponents, nnodes, load_systems, get_system, parse_yeh_raw
# signatures
export CompiledSystem, compile_system, works, state_from_working,
       count_working_states, tail_signature, signature, check_counts,
       signature_path, write_signature_file, read_signature_file

# TODO: LTMC formulas for delayed r-out-of-n policies, simulation

end # module
