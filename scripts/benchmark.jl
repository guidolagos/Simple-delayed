# Timing of each step of the r-out-of-n, count-based and state-based computations,
# for the runtime analysis of the paper. Not part of the production runs: it
# writes nothing in output/optimal*, only one CSV with timings.
#
# Usage (from the repository root, on a compute node, e.g. templates/benchmark.sbatch):
#   julia --project=. scripts/benchmark.jl [options]
#
# Options:
#   --systems=a,b,...   systems to time (default: all 20)
#   --model=NAME        LFMO model of the timed rows (default LFMO-C)
#   --beta=B            its beta, as written in data/grid.csv (default 0.2)
#   --delays=k:v,...    delays, as delay_type:delay of data/grid.csv
#                       (default none:0,exp:0.1,det:0.1,insp:0.1)
#   --pairs=all | r:s,...  cost pairs (ratio:setup) for the optimal policies
#                       (default all 36 of OPT_RATIOS x OPT_SETUPS)
#   --state-max-n=K     time state-based policies for n <= K (default 22)
#   --reps=K            repetitions of the cheap steps; the minimum is kept (default 5)
#   --out=PATH          output CSV (default output/benchmark/benchmark.csv)
#
# Every step is timed after a warm-up on the smallest system, so the times
# exclude compilation. Columns of the output:
#   system,n,model,beta,delay_type,delay,step,ratio,setup,seconds,iterations,maxrss_mb
# Steps (seconds are per call):
#   rates            transition rates of the survivor-count chain (512-bit), once per (model, n)
#   closed           delay-independent closed forms (closed_forms, includes `rates`), once per (model, n)
#   delay_forms      delay-dependent closed forms (BigFloat), once per (model, n, delay)
#   rout             p_fail, E N, E T of all r-out-of-n policies (cost-free) from the closed
#                    forms, given `delay_forms`, once per row
#   rout_chain       the same with the survivor-count chain engine (Float64 check)
#   ltmc_eval        long-term mean cost of all r for one cost pair, and its minimum
#   tables           after-call / one-period tables from the closed forms, given
#                    `delay_forms`, and the CountModel, once per row
#   tables_chain     the same with the survivor-count chain engine (Float64 check)
#   count_pass       one count-based pass (alg:count *), at the best r-out-of-n cost
#   count_opt        count-based Dinkelbach, warm-started at the best r-out-of-n (iterations = passes)
#   states           state space (2^n structure-function evaluations), once per system
#   state_pre        per-state signatures and after-call quantities, once per row
#   state_opt        state-based Dinkelbach, warm-started at the count-based optimum
#                    (iterations = passes; seconds / iterations = time of one pass)
# maxrss_mb is the peak resident memory of the process so far (Sys.maxrss).

using SimpleDelayed

function tmin(f, reps::Integer)
    best = Inf
    out = nothing
    for _ in 1:reps
        t = time_ns()
        out = f()
        best = min(best, (time_ns() - t) / 1e9)
    end
    return best, out
end

maxrss_mb() = round(Sys.maxrss() / 2^20; digits = 1)

"Long-term mean cost of every r-out-of-n policy for one cost pair, and its minimum."
function best_r(pf, en, et, costs)
    Lb = Inf
    ib = 0
    for i in eachindex(pf)
        v = (costs.c_set + costs.c_cmp * en[i] + costs.c_sys * pf[i]) / et[i]
        if v < Lb
            Lb = v
            ib = i
        end
    end
    return Lb, ib
end

function parse_args(args)
    o = Dict{String,String}()
    for a in args
        startswith(a, "--") || error("unexpected argument $(a)")
        key, val = occursin('=', a) ? split(a[3:end], '='; limit = 2) : (a[3:end], "")
        o[String(key)] = String(val)
    end
    known = ("systems", "model", "beta", "delays", "pairs", "state-max-n", "reps", "out")
    for k in keys(o)
        k in known || error("unknown option --$(k)")
    end
    return o
end

function main(args)
    o = parse_args(args)
    systems = load_systems()
    sysids = haskey(o, "systems") ? parse.(Int, split(o["systems"], ',')) : [s.id for s in systems]
    model = get(o, "model", "LFMO-C")
    beta = get(o, "beta", "0.2")
    delays = [Tuple(split(d, ':')) for d in split(get(o, "delays", "none:0,exp:0.1,det:0.1,insp:0.1"), ',')]
    pairs = get(o, "pairs", "all") == "all" ? [(ρ, σ) for ρ in OPT_RATIOS for σ in OPT_SETUPS] :
            [Tuple(parse.(Float64, split(p, ':'))) for p in split(o["pairs"], ',')]
    state_max_n = parse(Int, get(o, "state-max-n", "22"))
    reps = parse(Int, get(o, "reps", "5"))
    out = get(o, "out", joinpath(dirname(OPTIMAL_DIR), "benchmark", "benchmark.csv"))

    grid = read_grid()
    rows = [r for r in grid if r.system in sysids && r.model == model && r.beta == beta &&
            any(d -> r.delay_type == d[1] && r.delay == d[2], delays)]
    isempty(rows) && error("no grid row matches the selection")
    nof = Dict(s.id => ncomponents(s) for s in systems)
    sort!(rows; by = r -> (nof[r.system], r.row_id))
    println("benchmark: $(length(rows)) rows, model $(model) β=$(beta), $(length(pairs)) cost pairs, ",
            "state-based for n <= $(state_max_n), threads = $(Threads.nthreads()), host $(Libc.gethostname())")
    flush(stdout)

    lines = String["system,n,model,beta,delay_type,delay,step,ratio,setup,seconds,iterations,maxrss_mb"]
    function record(row, n, step, ρ, σ, secs, its)
        push!(lines, join((row.system, n, row.model, row.beta, row.delay_type, row.delay, step,
                           ρ, σ, secs, its, maxrss_mb()), ","))
    end

    # one row: all steps; `rec = false` for the warm-up
    function run_one(row, rec::Bool)
        sys = get_system(systems, row.system)
        n = ncomponents(sys)
        sigfile = read_signature_file(signature_path(sys))
        s = Float64.(sigfile.s)
        m = lfmo_model(row.model, row.beta)
        delay = Delay(row.delay_type, parse(Float64, row.delay))
        t, (R, exits) = tmin(() -> transition_rates(m, n), 1)
        rec && record(row, n, "rates", "", "", t, "")
        t, cf = tmin(() -> closed_forms(m, n), 1)
        rec && record(row, n, "closed", "", "", t, "")
        df = nothing
        if delay.kind != :none
            t, df = tmin(() -> delay_forms(cf, delay), reps)
            rec && record(row, n, "delay_forms", "", "", t, "")
        end
        t, (rs, pf, en, et) = tmin(() -> rout_quantities(sigfile.s, cf, delay; df = df), reps)
        rec && record(row, n, "rout", "", "", t, "")
        t, _ = tmin(() -> rout_quantities_chain(s, R, exits, delay), reps)
        rec && record(row, n, "rout_chain", "", "", t, "")
        t, cm = tmin(() -> CountModel(sigfile.W, R, exits, call_tables(cf, delay; df = df)), reps)
        rec && record(row, n, "tables", "", "", t, "")
        t, _ = tmin(() -> CountModel(sigfile.W, R, exits, delay), reps)
        rec && record(row, n, "tables_chain", "", "", t, "")
        do_state = n <= state_max_n
        S = pre = nothing
        if do_state
            GC.gc()
            t, S = tmin(() -> StateSpace(sys), 1)
            rec && record(row, n, "states", "", "", t, "")
            t, pre = tmin(() -> state_precompute(S, cm), 1)
            rec && record(row, n, "state_pre", "", "", t, "")
        end
        for (ρ, σ) in pairs
            costs = Costs(n; ratio = ρ, setup = σ)
            t, (Lb, ib) = tmin(() -> best_r(pf, en, et, costs), 100 * reps)
            rec && record(row, n, "ltmc_eval", ρ, σ, t, "")
            t, _ = tmin(() -> count_pass(cm, costs, Lb), reps)
            rec && record(row, n, "count_pass", ρ, σ, t, 1)
            t, (gc, _, itc) = tmin(() -> dinkelbach(g -> count_pass(cm, costs, g), Lb), reps)
            rec && record(row, n, "count_opt", ρ, σ, t, itc)
            if do_state
                t, (gs, _, its) = tmin(() -> dinkelbach(g -> state_pass(S, cm, pre, costs, g), gc), 1)
                rec && record(row, n, "state_opt", ρ, σ, t, its)
                gs <= gc * (1 + 1e-9) || error("row $(row.row_id): state-based $(gs) > count-based $(gc)")
            end
            rec || break                      # warm-up: one cost pair is enough
        end
        return n
    end

    # compile everything on a row of the smallest system (same model and delay kind),
    # so that the warm-up is cheap even when only large systems are selected
    small = argmin(s -> ncomponents(s), systems).id
    warm = only([r for r in grid if r.system == small && r.model == rows[1].model &&
                 r.beta == rows[1].beta && r.delay_type == rows[1].delay_type && r.delay == rows[1].delay])
    run_one(warm, false)
    t0 = time()
    for row in rows
        t1 = time()
        n = run_one(row, true)
        println("[done] row $(row.row_id) (system $(row.system), n = $(n), $(row.delay_type) $(row.delay)): ",
                round(time() - t1; digits = 1), " s, maxrss $(maxrss_mb()) MB")
        flush(stdout)
    end
    mkpath(dirname(out))
    tmp = out * ".tmp"
    write(tmp, join(lines, "\n") * "\n")
    mv(tmp, out; force = true)
    println("wrote $(out): $(length(lines) - 1) timings in ", round(time() - t0; digits = 1), " s")
end

main(ARGS)
