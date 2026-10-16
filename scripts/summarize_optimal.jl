# Consolidate the optimal-policy results into output/optimal_summary.csv.
#
# Usage (from the repository root; cheap, run it on your own machine after
# rsync-ing output/optimal and output/optimal_large from the cluster):
#   julia --project=. scripts/summarize_optimal.jl
#
# One line per grid row and cost pair (ratio = n c_comp/c_sys, setup = c_set/c_sys).
# Every row needs its file in output/optimal/ (version OPTIMAL_VERSION: closed-form
# engine), which gives the r-out-of-n and count-based columns, and the state-based
# columns for n <= 14. For the larger systems the state-based columns come from
# output/optimal_large/, whose files may be of version 1 (tables from the
# survivor-count chain; computing them again would take ~1,900 CPU-hours). This is
# the only place where version 1 is accepted, and only after checking that the
# r-out-of-n and count-based columns of both files agree (relative tolerance 1e-9);
# the largest difference is printed. The tables of both engines agree to ~1e-15
# (engine_check_* in the row metadata; scripts/engine_check.jl).
# Columns, besides the row metadata and those of the row files:
#   state_source        main (output/optimal), large (output/optimal_large) or none
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
    state_cols = ("ltmc_state", "iters_state", "state_mixed_levels")
    println(io, "row_id,system,n,model,beta,param,mean_batch_n,delay_type,delay,state_based,state_source,",
            join(keep, ","), ",gain_count_vs_r,gain_state_vs_count,gain_state_vs_r,gain_r_vs_corr")
    missing_rows = Int[]; nstate = 0; nlarge = 0; nlarge_v1 = 0; maxdiff_large = 0.0
    for row in grid
        path = row_output_path(row.row_id; dir = OPTIMAL_DIR)
        if !isfile(path)
            push!(missing_rows, row.row_id)
            continue
        end
        res = read_optimal_output(path)
        n = ncomponents(get_system(systems, row.system))
        parse(Int, res.meta["n_components"]) == n || error("$(path): wrong n")
        col(name) = findfirst(==(name), res.header)
        path_large = row_output_path(row.row_id; dir = OPTIMAL_LARGE_DIR)
        large = isfile(path_large) ? read_optimal_output(path_large; versions = (1, OPTIMAL_VERSION)) : nothing
        if large !== nothing
            nlarge += 1
            large.meta["optimal_version"] == "1" && (nlarge_v1 += 1)
            parse(Int, large.meta["n_components"]) == n || error("$(path_large): wrong n")
            large.header == res.header && length(large.rows) == length(res.rows) ||
                error("$(path_large): not comparable with $(path)")
            for (a, b) in zip(res.rows, large.rows)
                a[col("ratio")] == b[col("ratio")] && a[col("setup")] == b[col("setup")] ||
                    error("$(path_large): cost pairs differ from $(path)")
                for name in ("ltmc_best_r", "ltmc_corrective", "ltmc_count")
                    x, y = parse(Float64, a[col(name)]), parse(Float64, b[col(name)])
                    d = abs(x - y) / abs(x)
                    d <= 1e-9 || error("row $(row.row_id): $(name) differs between $(path) and $(path_large) ($(d))")
                    maxdiff_large = max(maxdiff_large, d)
                end
            end
        end
        state_based = large !== nothing ? large.meta["state_based"] : res.meta["state_based"]
        source = large !== nothing && large.meta["state_based"] == "true" ? "large" :
                 res.meta["state_based"] == "true" ? "main" : "none"
        state_based == "true" && (nstate += 1)
        for (i, r) in enumerate(res.rows)
            rl = large === nothing ? r : large.rows[i]
            v(name) = name in state_cols ? rl[col(name)] : r[col(name)]
            fl(name) = parse(Float64, v(name))
            gc, gs, gr, gcorr = fl("ltmc_count"), fl("ltmc_state"), fl("ltmc_best_r"), fl("ltmc_corrective")
            println(io, join((row.row_id, row.system, n, row.model, row.beta, res.meta["param"],
                              res.meta["mean_batch_n"], row.delay_type, row.delay, state_based, source,
                              (v(k) for k in keep)...,
                              repr(1 - gc / gr), repr(1 - gs / gc), repr(1 - gs / gr), repr(1 - gr / gcorr)), ","))
        end
    end
    mkpath(dirname(out))
    write(out, take!(io))
    println("wrote $(out): $(length(grid) - length(missing_rows)) of $(length(grid)) rows ",
            "($(nstate) with state-based policies) × $(length(OPT_RATIOS) * length(OPT_SETUPS)) cost pairs")
    println("output/optimal_large: $(nlarge) rows ($(nlarge_v1) of version 1); largest relative difference ",
            "with output/optimal in the r-out-of-n and count-based columns: $(maxdiff_large)")
    isempty(missing_rows) || println("missing $(length(missing_rows)) rows, e.g. ",
                                     join(first(missing_rows, 10), ", "))
end

main()
