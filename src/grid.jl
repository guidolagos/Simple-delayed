# Experiment grid and per-row processing.
#
# data/grid.csv has one row per (system, model, beta, delay). Each row writes
# output/ltmc/row_<row_id>.csv atomically, with one line per trigger r = 1..n,
# plus r = 0 (call the crew right after each repair) when the delay is
# exponential or deterministic (continuous inspection; Remark rmk:call0 of the
# paper); rows whose output file exists are skipped, so array jobs are resumable.
# The cost ratio n c_cmp/c_sys is not a grid dimension: the LTMC for any cost
# is computed afterwards from (p_fail, e_nfail, e_trep) by scripts/summarize.jl.

const GRID_PATH = joinpath(DATA_DIR, "grid.csv")
const OUTPUT_DIR = joinpath(REPO_DIR, "output", "ltmc")
# Version of the row files written by `run_row`; scripts refuse older outputs.
# 1: r = 1..n (commit dbd5263). 2: adds r = 0 for exponential and deterministic delays.
const ENGINE_VERSION = 2

# Experimental design (see the paper, Section "Computational experiments").
const GRID_BETAS = ["0.01", "0.05", "0.1", "0.2", "0.5"]
const GRID_DELAY_VALUES = ["0.01", "0.03", "0.1", "0.3", "1"]
const GRID_DELAY_KINDS = ["exp", "det", "insp"]

"""
    make_grid(; family = "yeh2024") -> Vector{NamedTuple}

All rows of the experiment, in a fixed order: systems × models (Exp with
β = 0, every other model with every β in `GRID_BETAS`) × delays (none, plus
every kind in `GRID_DELAY_KINDS` with every value in `GRID_DELAY_VALUES`).
`seed` equals `row_id`; the computations are deterministic and do not use it.
"""
function make_grid(; family::AbstractString = "yeh2024")
    systems = load_systems(family)
    models = [("Exp", "0")]
    for name in MODELS[2:end], b in GRID_BETAS
        push!(models, (name, b))
    end
    delays = [("none", "0")]
    for kind in GRID_DELAY_KINDS, v in GRID_DELAY_VALUES
        push!(delays, (kind, v))
    end
    rows = NamedTuple[]
    id = 0
    for sys in systems, (name, b) in models, (kind, v) in delays
        id += 1
        push!(rows, (row_id = id, system = sys.id, model = name, beta = b,
                     delay_type = kind, delay = v, seed = id))
    end
    return rows
end

const GRID_HEADER = ["row_id", "system", "model", "beta", "delay_type", "delay", "seed"]

function write_grid(path::AbstractString, rows)
    io = IOBuffer()
    println(io, join(GRID_HEADER, ","))
    for r in rows
        println(io, join((r.row_id, r.system, r.model, r.beta, r.delay_type, r.delay, r.seed), ","))
    end
    return atomic_write(path, take!(io))
end

function read_grid(path::AbstractString = GRID_PATH)
    header, rows = read_simple_csv(path)
    header == GRID_HEADER || error("$(path): unexpected header $(header)")
    return [(row_id = parse(Int, r[1]), system = parse(Int, r[2]), model = r[3], beta = r[4],
             delay_type = r[5], delay = r[6], seed = parse(Int, r[7])) for r in rows]
end

row_output_path(row_id::Integer; dir::AbstractString = OUTPUT_DIR) =
    joinpath(dir, "row_" * lpad(row_id, 5, '0') * ".csv")

"""
    run_row(row; systems, dir = OUTPUT_DIR, cache = Dict()) -> Symbol

Compute one grid row and write its output file; returns `:done`, or `:skipped`
if the file already exists. `cache` memoizes calibrated models and rates.
"""
function run_row(row; systems::Vector{TwoTerminalSystem}, dir::AbstractString = OUTPUT_DIR,
                 cache::Dict = Dict())
    path = row_output_path(row.row_id; dir = dir)
    isfile(path) && return :skipped
    t0 = time()
    sys = get_system(systems, row.system)
    n = ncomponents(sys)
    sigfile = read_signature_file(signature_path(sys))
    s = Float64.(sigfile.s)
    model = get!(cache, (:model, row.model, row.beta)) do
        lfmo_model(row.model, row.beta)
    end
    R, exits = get!(cache, (:rates, row.model, row.beta, n)) do
        transition_rates(model, n)
    end
    delay = Delay(row.delay_type, parse(Float64, row.delay))
    pf, en, et = policy_quantities(s, R, exits, delay)
    rs = collect(1:n)
    if delay.kind in (:exp, :det)
        pf0, en0, et0 = policy_quantities_r0(s, R, exits, delay)
        pushfirst!(rs, 0); pushfirst!(pf, pf0); pushfirst!(en, en0); pushfirst!(et, et0)
    end
    _check_quantities(rs, pf, en, et, n, row)
    seconds = time() - t0

    io = IOBuffer()
    println(io, "# SimpleDelayed LTMC quantities of the r-out-of-n policy with delay")
    for (key, val) in pairs(row)
        println(io, "# ", key, " = ", val)
    end
    println(io, "# name = ", sys.name)
    println(io, "# n_components = ", n)
    println(io, "# param = ", string(Float64(model.param)))
    println(io, "# mean_batch_n = ", mean_batch(model, n))
    println(io, "# seconds = ", round(seconds; digits = 4))
    println(io, "# engine_version = ", ENGINE_VERSION)
    println(io, "# julia_version = ", VERSION)
    println(io, "# host = ", _hostname())
    println(io, "r,p_fail,e_nfail,e_trep")
    for (i, r) in enumerate(rs)
        println(io, r, ",", repr(pf[i]), ",", repr(en[i]), ",", repr(et[i]))
    end
    atomic_write(path, take!(io))
    return :done
end

# Sanity checks that hold for every system, model and delay. The vectors are
# indexed like `rs` (triggers in increasing order, possibly starting at r = 0).
function _check_quantities(rs, pf, en, et, n, row)
    tol = 1e-9
    for (i, r) in enumerate(rs)
        -tol <= pf[i] <= 1 + tol || error("row $(row.row_id): p_fail[r=$(r)] = $(pf[i]) outside [0,1]")
        (r == 0 ? -tol : 1 - tol) <= en[i] <= n + tol ||
            error("row $(row.row_id): e_nfail[r=$(r)] = $(en[i]) out of range")
        et[i] > 0 && isfinite(et[i]) || error("row $(row.row_id): e_trep[r=$(r)] = $(et[i])")
    end
    # the last trigger (r = n) never fires before the system fails
    abs(pf[end] - 1) <= tol || error("row $(row.row_id): p_fail[r=n] = $(pf[end]) ≠ 1")
    # a later trigger never repairs earlier: T_rep, N(T_rep) and the event of a
    # system failure are nondecreasing in r for each sample path, for every
    # delay (also from r = 0 to r = 1, and for inspections, whose times are fixed)
    for v in (pf, en, et), i in 1:(length(rs) - 1)
        v[i] <= v[i + 1] + 1e-10 * max(1.0, abs(v[i + 1])) ||
            error("row $(row.row_id): quantities not nondecreasing in r at r = $(rs[i])")
    end
    return true
end

"""
    read_row_output(path) -> (meta, r, p_fail, e_nfail, e_trep)

Errors if the file was written by an engine older than `ENGINE_VERSION`.
"""
function read_row_output(path::AbstractString)
    meta = Dict{String,String}()
    for line in eachline(path)
        startswith(line, "#") || continue
        m = match(r"^#\s*([A-Za-z_]+)\s*=\s*(.*)$", line)
        m === nothing || (meta[m.captures[1]] = strip(m.captures[2]))
    end
    get(meta, "engine_version", "1") == string(ENGINE_VERSION) ||
        error("$(path) was written by engine version $(get(meta, "engine_version", "1")), expected ",
              "$(ENGINE_VERSION); move the old output/ltmc away and re-run templates/run_array.sbatch")
    header, rows = read_simple_csv(path)
    header == ["r", "p_fail", "e_nfail", "e_trep"] || error("$(path): unexpected header")
    col(i, T) = [parse(T, x[i]) for x in rows]
    return (meta = meta, r = col(1, Int), p_fail = col(2, Float64),
            e_nfail = col(3, Float64), e_trep = col(4, Float64))
end
