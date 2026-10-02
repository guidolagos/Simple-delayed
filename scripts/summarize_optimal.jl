# Consolidate the optimal-policy results into output/optimal_summary.csv.
#
# Usage (from the repository root; cheap, run it on your own machine after
# rsync-ing output/optimal and output/optimal_large from the cluster):
#   julia --project=. scripts/summarize_optimal.jl
#
# One line per grid row and cost pair (ratio = n c_cmp/c_sys, setup = c_set/c_sys).
# For each row the file in output/optimal_large/ is used if present (it has the
# state-based results for the larger systems), otherwise the one in output/optimal/.
# Columns, besides the row metadata and those of the row files:
#   gain_count_vs_r     1 - ltmc_count / ltmc_best_r     value of departing from a threshold
#   gain_state_vs_count 1 - ltmc_state / ltmc_count      value of knowing which components failed
#   gain_state_vs_r     1 - ltmc_state / ltmc_best_r
#   gain_r_vs_corr      1 - ltmc_best_r / ltmc_corrective
# (state columns are NaN where state-based policies were not computed).

using SimpleDelayed

function main()
    grid = read_grid()
    systems = load_systems()
    out = joinpath(dirname(OPTIMAL_DIR), "optimal_summary.csv")
    io = IOBuffer()
    keep = ["ratio", "setup", "best_r", "ltmc_best_r", "ltmc_corrective", "ltmc_count", "iters_count",
            "count_policy", "count_threshold", "ltmc_state", "iters_state", "state_mixed_levels"]
    println(io, "row_id,system,n,model,beta,param,mean_batch_n,delay_type,delay,state_based,",
            join(keep, ","), ",gain_count_vs_r,gain_state_vs_count,gain_state_vs_r,gain_r_vs_corr")
    missing_rows = Int[]; nstate = 0
    for row in grid
        path_large = row_output_path(row.row_id; dir = OPTIMAL_LARGE_DIR)
        path = isfile(path_large) ? path_large : row_output_path(row.row_id; dir = OPTIMAL_DIR)
        if !isfile(path)
            push!(missing_rows, row.row_id)
            continue
        end
        res = read_optimal_output(path)
        n = ncomponents(get_system(systems, row.system))
        parse(Int, res.meta["n_components"]) == n || error("$(path): wrong n")
        res.meta["state_based"] == "true" && (nstate += 1)
        col(name) = findfirst(==(name), res.header)
        for r in res.rows
            v(name) = r[col(name)]
            fl(name) = parse(Float64, v(name))
            gc, gs, gr, gcorr = fl("ltmc_count"), fl("ltmc_state"), fl("ltmc_best_r"), fl("ltmc_corrective")
            println(io, join((row.row_id, row.system, n, row.model, row.beta, res.meta["param"],
                              res.meta["mean_batch_n"], row.delay_type, row.delay, res.meta["state_based"],
                              (v(k) for k in keep)...,
                              repr(1 - gc / gr), repr(1 - gs / gc), repr(1 - gs / gr), repr(1 - gr / gcorr)), ","))
        end
    end
    mkpath(dirname(out))
    write(out, take!(io))
    println("wrote $(out): $(length(grid) - length(missing_rows)) of $(length(grid)) rows ",
            "($(nstate) with state-based policies) × $(length(OPT_RATIOS) * length(OPT_SETUPS)) cost pairs")
    isempty(missing_rows) || println("missing $(length(missing_rows)) rows, e.g. ",
                                     join(first(missing_rows, 10), ", "))
end

main()
