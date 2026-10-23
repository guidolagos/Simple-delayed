# Compare the row files of output/optimal (version OPTIMAL_VERSION) with older outputs
# of the same rows, e.g. those of version 1 (state-based policies for n >= 21 with the
# chain tables) and version 2 (closed forms; state-based for n <= 14 only).
#
# Usage (from the repository root; cheap):
#   julia --project=. scripts/compare_optimal.jl OLD_DIR [OLD_DIR ...] [--rows=a:b]
# e.g., on the cluster after the pilots:
#   julia --project=. scripts/compare_optimal.jl output/optimal_large_v1 output/optimal_v2 --rows=2977:9920
#
# This is the only script that reads row files of other versions (versions 1 to
# OPTIMAL_VERSION are accepted in OLD_DIR). For every OLD_DIR and every row with a file
# in both directories, it compares, per cost pair:
#   relative differences of ltmc_best_r, ltmc_corrective, ltmc_count, and ltmc_state
#   (when both files have state-based policies);
#   equality of best_r, count_policy and state_calls (calls/working states per number
#   of failed components: the row files do not record the states themselves).
# It prints the largest differences by system size, and every pair whose best trigger
# or policy differs, with the relative difference of the corresponding cost: a
# different policy with costs equal to ~1e-13 is a near-tie.

using SimpleDelayed

function main(args)
    dirs = filter(a -> !startswith(a, "--"), args)
    isempty(dirs) && error("usage: compare_optimal.jl OLD_DIR [OLD_DIR ...] [--rows=a:b]")
    rowrange = 1:typemax(Int)
    for a in filter(a -> startswith(a, "--"), args)
        if startswith(a, "--rows=")
            lo, hi = parse.(Int, split(a[8:end], ':'))
            rowrange = lo:hi
        else
            error("unknown option $(a)")
        end
    end
    grid = [r for r in read_grid() if r.row_id in rowrange]
    systems = load_systems()
    nof = Dict(s.id => ncomponents(s) for s in systems)
    costcols = ("ltmc_best_r", "ltmc_corrective", "ltmc_count", "ltmc_state")
    for dir in dirs
        isdir(dir) || error("$(dir) is not a directory")
        # per n: rows compared, rows with state in both, max relative difference per cost column
        stats = Dict{Int,Vector{Float64}}()
        diffs = String[]
        ndiff = Dict("best_r" => 0, "count_policy" => 0, "state_calls" => 0)
        versions = Set{String}()
        for row in grid
            pnew = row_output_path(row.row_id; dir = OPTIMAL_DIR)
            pold = row_output_path(row.row_id; dir = dir)
            (isfile(pnew) && isfile(pold)) || continue
            new = read_optimal_output(pnew)
            old = read_optimal_output(pold; versions = Tuple(1:OPTIMAL_VERSION))
            push!(versions, old.meta["optimal_version"])
            n = nof[row.system]
            length(new.rows) == length(old.rows) || error("row $(row.row_id): different number of cost pairs")
            cn(name) = findfirst(==(name), new.header)
            co(name) = findfirst(==(name), old.header)
            both_state = new.meta["state_based"] == "true" && old.meta["state_based"] == "true"
            st = get!(stats, n, zeros(2 + length(costcols)))
            st[1] += 1; st[2] += both_state
            for (a, b) in zip(new.rows, old.rows)
                a[cn("ratio")] == b[co("ratio")] && a[cn("setup")] == b[co("setup")] ||
                    error("row $(row.row_id): cost pairs differ")
                rel = Dict{String,Float64}()
                for (k, name) in enumerate(costcols)
                    name == "ltmc_state" && !both_state && continue
                    x, y = parse(Float64, a[cn(name)]), parse(Float64, b[co(name)])
                    rel[name] = abs(x - y) / abs(y)
                    st[2 + k] = max(st[2 + k], rel[name])
                end
                checks = [("best_r", "ltmc_best_r"), ("count_policy", "ltmc_count")]
                both_state && push!(checks, ("state_calls", "ltmc_state"))
                for (name, cost) in checks
                    a[cn(name)] == b[co(name)] && continue
                    ndiff[name] += 1
                    push!(diffs, "row $(row.row_id) (n = $(n), $(row.model) β=$(row.beta), $(row.delay_type) " *
                                 "$(row.delay)), (ρ, σ) = ($(a[cn("ratio")]), $(a[cn("setup")])): $(name) differs " *
                                 "(new $(a[cn(name)]), old $(b[co(name)])); relative difference of $(cost): $(rel[cost])")
                end
            end
        end
        println("== $(dir) (versions $(join(sort(collect(versions)), ", "))) vs $(OPTIMAL_DIR) (version $(OPTIMAL_VERSION))")
        println("    n   rows  state  max rel. diff: best_r  corrective   count   state")
        for n in sort(collect(keys(stats)))
            st = stats[n]
            println(lpad(n, 5), lpad(Int(st[1]), 7), lpad(Int(st[2]), 7), "  ",
                    join((rpad(string(round(v; sigdigits = 2)), 10) for v in st[3:end]), " "))
        end
        println("pairs with a different best trigger: $(ndiff["best_r"]), count-based policy: ",
                "$(ndiff["count_policy"]), state-based calls per level: $(ndiff["state_calls"])")
        for d in first(diffs, 50)
            println("  ", d)
        end
        length(diffs) > 50 && println("  ... and $(length(diffs) - 50) more")
    end
end

main(ARGS)
