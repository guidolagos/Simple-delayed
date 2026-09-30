# Summarize the LTMC results: optimal trigger r* for every grid row and cost ratio.
#
# Usage (from the repository root; cheap, run it on your own machine after
# rsync-ing output/ltmc from the cluster):
#   julia --project=. scripts/summarize.jl
#
# For every grid row, every ratio ρ = n c_cmp/c_sys in COST_RATIOS and every
# setup σ = c_setup/c_sys in COST_SETUPS, with c_sys = 1 and repair cost
# c_cmp(j) = c_setup + j c_cmp, writes one line of output/summary.csv:
#   r_opt, ltmc_opt           optimal trigger (smallest minimizer) and its long-term mean cost
#   ltmc_r1                   repair at the first failure (r = 1)
#   ltmc_corrective           corrective only (r = n: repair when the system fails)
#   saving_vs_corrective      1 - ltmc_opt / ltmc_corrective
#   corrective_equivalent     true if r_opt is as good as corrective (no preventive gain)
# Missing row outputs are reported and skipped.

using SimpleDelayed

function main()
    grid = read_grid()
    systems = load_systems()
    out = joinpath(dirname(OUTPUT_DIR), "summary.csv")
    io = IOBuffer()
    println(io, "row_id,system,n,model,beta,param,mean_batch_n,delay_type,delay,ratio,setup,",
            "r_opt,ltmc_opt,ltmc_r1,ltmc_corrective,saving_vs_corrective,corrective_equivalent")
    missing_rows = Int[]
    for row in grid
        path = row_output_path(row.row_id)
        if !isfile(path)
            push!(missing_rows, row.row_id)
            continue
        end
        res = read_row_output(path)
        n = length(res.r)
        n == ncomponents(get_system(systems, row.system)) || error("$(path): wrong number of rows")
        for ρ in COST_RATIOS, σ in COST_SETUPS
            L = ltmc(res.p_fail, res.e_nfail, res.e_trep, n; ratio = ρ, setup = σ)
            ropt = argmin(L)
            println(io, join((row.row_id, row.system, n, row.model, row.beta, res.meta["param"],
                              res.meta["mean_batch_n"], row.delay_type, row.delay, ρ, σ,
                              ropt, repr(L[ropt]), repr(L[1]), repr(L[n]),
                              repr(1 - L[ropt] / L[n]), L[ropt] >= L[n] * (1 - 1e-12)), ","))
        end
    end
    mkpath(dirname(out))
    write(out, take!(io))
    println("wrote $(out): $(length(grid) - length(missing_rows)) of $(length(grid)) rows × $(length(COST_RATIOS)) ratios × $(length(COST_SETUPS)) setups")
    isempty(missing_rows) || println("missing $(length(missing_rows)) rows, e.g. ",
                                     join(first(missing_rows, 10), ", "))
end

main()
