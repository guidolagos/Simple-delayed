# Optimal count-based and state-based policies for one chunk of the rows of data/grid.csv.
#
# Usage (from the repository root):
#   julia --project=. scripts/run_optimal.jl CHUNK NCHUNKS [options]
#
# Options:
#   --state-max-n=K  compute state-based policies for systems with n <= K (default 14)
#   --min-n=K        only rows whose system has n >= K
#   --max-n=K        only rows whose system has n <= K
#   --out=large      write to output/optimal_large/ instead of output/optimal/
#   --rows=a:b       only row ids a..b (pilot runs)
#   --verbose        print the time of every state-based Dinkelbach run
#
# Chunk i (1-based) owns the selected rows with (position - 1) mod NCHUNKS == i - 1,
# where position is the index among the selected rows in grid order. Keep NCHUNKS and
# the selection options fixed between submissions; finished rows are skipped.
#
# Two runs cover the experiment (see README):
#   main:  all rows, state-based for n <= 14            -> output/optimal/
#   large: rows with 15 <= n <= 22, state-based for all -> output/optimal_large/

using SimpleDelayed

function main(args)
    positional = filter(a -> !startswith(a, "--"), args)
    length(positional) == 2 || error("usage: run_optimal.jl CHUNK NCHUNKS [options]")
    chunk, nchunks = parse.(Int, positional)
    1 <= chunk <= nchunks || error("chunk $(chunk) is outside 1..$(nchunks)")
    state_max_n = 14; min_n = 0; max_n = typemax(Int); dir = OPTIMAL_DIR
    rowrange = 1:typemax(Int); verbose = false
    for a in args
        startswith(a, "--") || continue
        key, val = occursin('=', a) ? split(a[3:end], '='; limit = 2) : (a[3:end], "")
        if key == "state-max-n"
            state_max_n = parse(Int, val)
        elseif key == "min-n"
            min_n = parse(Int, val)
        elseif key == "max-n"
            max_n = parse(Int, val)
        elseif key == "out"
            val == "large" || error("--out accepts only 'large'")
            dir = OPTIMAL_LARGE_DIR
        elseif key == "rows"
            lo, hi = parse.(Int, split(val, ':'))
            rowrange = lo:hi
        elseif key == "verbose"
            verbose = true
        else
            error("unknown option $(a)")
        end
    end

    systems = load_systems()
    nof = Dict(s.id => ncomponents(s) for s in systems)
    selected = [row for row in read_grid() if min_n <= nof[row.system] <= max_n && row.row_id in rowrange]
    mine = [row for (pos, row) in enumerate(selected) if (pos - 1) % nchunks == chunk - 1]
    mkpath(dir)
    println("chunk $(chunk)/$(nchunks): $(length(mine)) of $(length(selected)) selected rows; ",
            "state-based for n <= $(state_max_n); output in $(dir); threads = $(Threads.nthreads())")
    flush(stdout)

    cache = Dict()
    done = skipped = failed = 0
    t0 = time()
    for row in mine
        t1 = time()
        try
            status = run_optimal_row(row; systems = systems, dir = dir, state_max_n = state_max_n,
                                     cache = cache, verbose = verbose)
            if status == :done
                done += 1
                println("[done] row $(row.row_id) (system $(row.system), n = $(nof[row.system]), $(row.model) ",
                        "β=$(row.beta), $(row.delay_type) $(row.delay)): $(round(time() - t1; digits = 1)) s")
            else
                skipped += 1
            end
        catch err
            failed += 1
            println("[fail] row $(row.row_id) $(row): ", sprint(showerror, err))
        end
        flush(stdout)
    end
    println("chunk $(chunk)/$(nchunks): done $(done), skipped $(skipped), failed $(failed), ",
            round(time() - t0; digits = 1), " s")
    failed == 0 || exit(1)
end

main(ARGS)
