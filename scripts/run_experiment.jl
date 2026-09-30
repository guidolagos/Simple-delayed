# Compute the LTMC quantities for one chunk of the rows of data/grid.csv.
#
# Usage (from the repository root):
#   julia --project=. scripts/run_experiment.jl CHUNK NCHUNKS [--rows=a:b]
#
# Chunk i (1-based) processes the rows with (row_id - 1) mod NCHUNKS == i - 1,
# so that every chunk gets a similar mix of small and large systems. NCHUNKS
# must stay the same between submissions of the same grid; rows whose output
# file (output/ltmc/row_<id>.csv) exists are skipped, so re-submitting is safe.
# --rows=a:b restricts to row ids a..b (useful for a pilot run).

using SimpleDelayed

function main(args)
    positional = filter(a -> !startswith(a, "--"), args)
    length(positional) == 2 || error("usage: run_experiment.jl CHUNK NCHUNKS [--rows=a:b]")
    chunk, nchunks = parse.(Int, positional)
    1 <= chunk <= nchunks || error("chunk $(chunk) is outside 1..$(nchunks)")
    rowrange = 1:typemax(Int)
    for a in args
        if startswith(a, "--rows=")
            lo, hi = parse.(Int, split(a[length("--rows=")+1:end], ':'))
            rowrange = lo:hi
        elseif startswith(a, "--")
            error("unknown option $(a)")
        end
    end

    grid = read_grid()
    mine = [row for row in grid if (row.row_id - 1) % nchunks == chunk - 1 && row.row_id in rowrange]
    systems = load_systems()
    mkpath(OUTPUT_DIR)
    println("chunk $(chunk)/$(nchunks): $(length(mine)) rows; output in $(OUTPUT_DIR)")
    flush(stdout)

    cache = Dict()
    done = skipped = failed = 0
    t0 = time()
    for row in mine
        try
            status = run_row(row; systems = systems, cache = cache)
            status == :done ? (done += 1) : (skipped += 1)
        catch err
            failed += 1
            println("[fail] row $(row.row_id) $(row): ", sprint(showerror, err))
            flush(stdout)
        end
    end
    println("chunk $(chunk)/$(nchunks): done $(done), skipped $(skipped), failed $(failed), ",
            round(time() - t0; digits = 1), " s")
    failed == 0 || exit(1)
end

main(ARGS)
