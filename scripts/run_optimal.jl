# Optimal count-based and state-based policies for one chunk of the rows of data/grid.csv.
#
# Usage (from the repository root):
#   julia --project=. scripts/run_optimal.jl CHUNK NCHUNKS [options]
#
# Options:
#   --state-max-n=K  compute state-based policies for systems with n <= K
#                    (default OPT_STATE_MAX_N = 27: all systems but the two with n = 30)
#   --min-n=K        only rows whose system has n >= K
#   --max-n=K        only rows whose system has n <= K
#   --rows=a:b       only row ids a..b (pilot runs)
#   --verbose        print the time of every state-based Dinkelbach run
#   --no-warmup      skip the warm-up (see below)
#
# Chunk i (1-based) owns the selected rows with (position - 1) mod NCHUNKS == i - 1,
# where position is the index among the selected rows in grid order. Keep NCHUNKS and
# the selection options fixed between submissions; finished rows are skipped.
#
# Before the first row, the task runs four rows of the smallest system in a temporary
# directory (optimal_warmup), so that the times recorded in the row files exclude
# compilation. All rows are written to output/optimal/; the production run is split
# into arrays by system size only to request the right memory (templates/optimal.sbatch).

using SimpleDelayed

function main(args)
    positional = filter(a -> !startswith(a, "--"), args)
    length(positional) == 2 || error("usage: run_optimal.jl CHUNK NCHUNKS [options]")
    chunk, nchunks = parse.(Int, positional)
    1 <= chunk <= nchunks || error("chunk $(chunk) is outside 1..$(nchunks)")
    state_max_n = OPT_STATE_MAX_N; min_n = 0; max_n = typemax(Int); dir = OPTIMAL_DIR
    rowrange = 1:typemax(Int); verbose = false; warmup = true
    for a in args
        startswith(a, "--") || continue
        key, val = occursin('=', a) ? split(a[3:end], '='; limit = 2) : (a[3:end], "")
        if key == "state-max-n"
            state_max_n = parse(Int, val)
        elseif key == "min-n"
            min_n = parse(Int, val)
        elseif key == "max-n"
            max_n = parse(Int, val)
        elseif key == "rows"
            lo, hi = parse.(Int, split(val, ':'))
            rowrange = lo:hi
        elseif key == "verbose"
            verbose = true
        elseif key == "no-warmup"
            warmup = false
        else
            error("unknown option $(a)")
        end
    end

    systems = load_systems()
    nof = Dict(s.id => ncomponents(s) for s in systems)
    selected = [row for row in read_grid() if min_n <= nof[row.system] <= max_n && row.row_id in rowrange]
    mine = [row for (pos, row) in enumerate(selected) if (pos - 1) % nchunks == chunk - 1]
    todo = count(row -> !isfile(row_output_path(row.row_id; dir = dir)), mine)
    mkpath(dir)
    println("chunk $(chunk)/$(nchunks): $(length(mine)) of $(length(selected)) selected rows ($(todo) to do); ",
            "state-based for n <= $(state_max_n); output in $(dir); optimal_version $(OPTIMAL_VERSION); ",
            "threads = $(Threads.nthreads())")
    flush(stdout)
    if warmup && todo > 0
        secs = optimal_warmup(systems; state_max_n = state_max_n)
        println("warm-up: $(round(secs; digits = 1)) s")
        flush(stdout)
    end

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
                        "β=$(row.beta), $(row.delay_type) $(row.delay)): $(round(time() - t1; digits = 1)) s, ",
                        "maxrss $(round(Sys.maxrss() / 2^30; digits = 2)) GB")
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
            round(time() - t0; digits = 1), " s, maxrss $(round(Sys.maxrss() / 2^30; digits = 2)) GB")
    failed == 0 || exit(1)
end

main(ARGS)
