# Largest differences between the closed-form engine (BigFloat) and the
# survivor-count chain engine (Float64), from the metadata of the row files.
#
# Usage (from the repository root; cheap, after rsync-ing output/ltmc and output/optimal):
#   julia --project=. scripts/engine_check.jl
#
# Every row file of engine version >= 3 (output/ltmc) or optimal version >= 2
# (output/optimal) records
#   engine_check_prob  largest absolute difference of the probabilities: P(T_rep = T_fail)
#                      of every trigger r and, in output/optimal, the after-call tables
#                      pf, pN and enF (E[N(T_{i:w}); T_{i:w} < D])
#   engine_check_rel   largest relative difference of E[N(T_rep)], E[T_rep] and, in
#                      output/optimal, the tables enA and ell
# Writes output/engine_check.csv (one line per row file) and prints the maxima
# by directory and system size.

using SimpleDelayed

function scan(dir, kind, systems, grid, io)
    nof = Dict(s.id => ncomponents(s) for s in systems)
    worst = Dict{Int,Tuple{Float64,Float64,Int}}()
    nfiles = 0; missing_rows = 0
    for row in grid
        path = row_output_path(row.row_id; dir = dir)
        if !isfile(path)
            missing_rows += 1
            continue
        end
        meta = kind == "ltmc" ? read_row_output(path).meta : read_optimal_output(path).meta
        p = parse(Float64, meta["engine_check_prob"]); r = parse(Float64, meta["engine_check_rel"])
        n = nof[row.system]
        println(io, join((kind, row.row_id, row.system, n, row.model, row.beta, row.delay_type, row.delay,
                          repr(p), repr(r), meta["seconds"], meta["seconds_chain"]), ","))
        a, b, c = get(worst, n, (0.0, 0.0, 0))
        worst[n] = (max(a, p), max(b, r), c + 1)
        nfiles += 1
    end
    println("$(dir): $(nfiles) row files ($(missing_rows) missing)")
    for n in sort(collect(keys(worst)))
        a, b, c = worst[n]
        println("  n = $(lpad(n, 2)): $(lpad(c, 4)) rows, largest difference: probabilities $(a), relative $(b)")
    end
    isempty(worst) || println("  all: probabilities $(maximum(v[1] for v in values(worst))), ",
                              "relative $(maximum(v[2] for v in values(worst)))")
end

function main()
    grid = read_grid()
    systems = load_systems()
    io = IOBuffer()
    println(io, "dir,row_id,system,n,model,beta,delay_type,delay,engine_check_prob,engine_check_rel,seconds,seconds_chain")
    scan(OUTPUT_DIR, "ltmc", systems, grid, io)
    scan(OPTIMAL_DIR, "optimal", systems, grid, io)
    out = joinpath(dirname(OUTPUT_DIR), "engine_check.csv")
    write(out, take!(io))
    println("wrote $(out)")
end

main()
