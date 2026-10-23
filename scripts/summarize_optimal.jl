# Consolidate the optimal-policy results of output/optimal (version OPTIMAL_VERSION).
#
# Usage (from the repository root; cheap, run it on your own machine after
# rsync-ing output/optimal from the cluster):
#   julia --project=. scripts/summarize_optimal.jl
#
# Writes two files:
#
# output/optimal_summary.csv: one line per grid row and cost pair (ratio = n c_comp/c_sys,
# setup = c_set/c_sys), with the row metadata, the columns of the row files (except
# state_calls), and
#   gain_count_vs_r     1 - ltmc_count / ltmc_best_r     value of departing from a threshold
#   gain_state_vs_count 1 - ltmc_state / ltmc_count      value of knowing which components failed
#   gain_state_vs_r     1 - ltmc_state / ltmc_best_r
#   gain_r_vs_corr      1 - ltmc_best_r / ltmc_corrective
# (state columns are NaN where state-based policies were not computed, i.e. n = 30).
#
# output/optimal_rows.csv: one line per grid row with the timings and memory of the
# row file (seconds_rout, seconds_tables, seconds_states, states_cached,
# seconds_state_pre, maxrss_mb, seconds_chain, seconds) and, summed over the 36 cost
# pairs, seconds_count, passes_count, seconds_state and passes_state.
#
# Every row must have its file, of version OPTIMAL_VERSION (read_optimal_output refuses
# others); missing rows are listed and the script exits with code 1.

using SimpleDelayed

function main()
    grid = read_grid()
    systems = load_systems()
    outdir = dirname(OPTIMAL_DIR)
    io = IOBuffer(); iorows = IOBuffer()
    keep = ["ratio", "setup", "best_r", "ltmc_best_r", "ltmc_corrective", "ltmc_count", "iters_count",
            "seconds_count", "count_policy", "count_threshold", "ltmc_state", "iters_state",
            "state_mixed_levels", "seconds_state"]
    println(io, "row_id,system,n,model,beta,param,mean_batch_n,delay_type,delay,state_based,",
            join(keep, ","), ",gain_count_vs_r,gain_state_vs_count,gain_state_vs_r,gain_r_vs_corr")
    rowmeta = ["seconds_rout", "seconds_tables", "seconds_states", "states_cached", "seconds_state_pre",
               "maxrss_mb", "seconds_chain", "seconds", "n_working_states", "host"]
    println(iorows, "row_id,system,n,model,beta,delay_type,delay,state_based,", join(rowmeta, ","),
            ",seconds_count,passes_count,seconds_state,passes_state")
    missing_rows = Int[]; nstate = 0
    for row in grid
        path = row_output_path(row.row_id; dir = OPTIMAL_DIR)
        if !isfile(path)
            push!(missing_rows, row.row_id)
            continue
        end
        res = read_optimal_output(path)
        n = ncomponents(get_system(systems, row.system))
        parse(Int, res.meta["n_components"]) == n || error("$(path): wrong n")
        length(res.rows) == length(OPT_RATIOS) * length(OPT_SETUPS) || error("$(path): wrong number of cost pairs")
        col(name) = findfirst(==(name), res.header)
        state_based = res.meta["state_based"]
        state_based == string(n <= OPT_STATE_MAX_N) || error("$(path): state_based = $(state_based) for n = $(n)")
        state_based == "true" && (nstate += 1)
        sc = 0.0; pc = 0; ss = 0.0; ps = 0
        for r in res.rows
            fl(name) = parse(Float64, r[col(name)])
            gc, gs, gr, gcorr = fl("ltmc_count"), fl("ltmc_state"), fl("ltmc_best_r"), fl("ltmc_corrective")
            sc += fl("seconds_count"); pc += parse(Int, r[col("iters_count")])
            ss += fl("seconds_state"); ps += parse(Int, r[col("iters_state")])
            println(io, join((row.row_id, row.system, n, row.model, row.beta, res.meta["param"],
                              res.meta["mean_batch_n"], row.delay_type, row.delay, state_based,
                              (r[col(k)] for k in keep)...,
                              repr(1 - gc / gr), repr(1 - gs / gc), repr(1 - gs / gr), repr(1 - gr / gcorr)), ","))
        end
        println(iorows, join((row.row_id, row.system, n, row.model, row.beta, row.delay_type, row.delay,
                              state_based, (res.meta[k] for k in rowmeta)..., repr(sc), pc, repr(ss), ps), ","))
    end
    out = joinpath(outdir, "optimal_summary.csv"); outrows = joinpath(outdir, "optimal_rows.csv")
    mkpath(outdir)
    write(out, take!(io)); write(outrows, take!(iorows))
    println("wrote $(out) and $(outrows): $(length(grid) - length(missing_rows)) of $(length(grid)) rows ",
            "($(nstate) with state-based policies) × $(length(OPT_RATIOS) * length(OPT_SETUPS)) cost pairs")
    if !isempty(missing_rows)
        println("missing $(length(missing_rows)) rows, e.g. ", join(first(missing_rows, 10), ", "))
        exit(1)
    end
end

main()
