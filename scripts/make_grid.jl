# Write data/grid.csv, the list of experiments (one row per system, model, beta, delay).
#
# Usage (from the repository root, on your own machine; the file is committed):
#   julia --project=. scripts/make_grid.jl            # refuses to overwrite a different grid
#   julia --project=. scripts/make_grid.jl --force    # overwrite
#
# The design (betas, delay kinds and values) is set in src/grid.jl. Changing the
# grid after results exist changes the meaning of row ids: move output/ltmc away first.

using SimpleDelayed

function main(args)
    rows = make_grid()
    if isfile(GRID_PATH) && !("--force" in args)
        old = read_grid(GRID_PATH)
        same = length(old) == length(rows) &&
               all(Tuple(a) == Tuple(b) for (a, b) in zip(old, rows))
        if same
            println("$(GRID_PATH) is up to date ($(length(rows)) rows)")
            return
        end
        error("$(GRID_PATH) exists and differs from the design in src/grid.jl; use --force to overwrite")
    end
    write_grid(GRID_PATH, rows)
    println("wrote $(GRID_PATH): $(length(rows)) rows")
end

main(ARGS)
