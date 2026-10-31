# Rebuild the numbers of the computational section of the paper (Tables 3, 4 and 5, and the
# statistics in the text of Sections 7.5 to 7.7) from the raw result files.
#
# Usage (from the repository root, with the results in output/, e.g. those of the Zenodo
# archive; cheap, about a minute):
#   julia --project=. scripts/paper_numbers.jl
#
# Reads output/ltmc/row_*.csv (version ENGINE_VERSION), output/optimal/row_*.csv (version
# OPTIMAL_VERSION), data/signatures/ and output/benchmark/benchmark_n30_s13.csv,
# benchmark_n30_s19.csv, and writes everything it prints to output/paper_numbers.txt.
# It does not use output/summary.csv nor output/optimal_summary.csv.
#
# Conventions, as in the paper:
# - long-term mean cost with c_sys = 1, LTMC(r) = (c_set + c_comp E N + P(fail)) / E T, for the
#   cost pairs (n c_comp/c_sys, c_set/c_sys) of COST_RATIOS x COST_SETUPS; r = 0 is admissible
#   only when c_set > 0; the best trigger is the smallest minimizer;
# - savings: 1 - g*_r / g_corr (Table 3), 1 - g*_count / g*_r and 1 - g*_state / g*_r (Table 4);
# - medians and percentiles with linear interpolation (type 7, as Statistics.quantile);
# - Table 5: medians over the configurations of each n of seconds_rout (r-out-of-n),
#   seconds_tables + the 36 seconds_count (count-based), seconds_state / iters_state over the
#   cost pairs (one pass), mean iters_state per cost pair (passes), the 36 seconds_state
#   (36 pairs) and maxrss_mb of the first configuration of a system in each task
#   (states_cached = false; memory, GiB); for n = 30, the state-based values come from the
#   one-cost-pair benchmarks, and the 36 pairs are extrapolated with the smallest and largest
#   mean number of passes for 24 <= n <= 27.

using SimpleDelayed
using Printf

const OUTPUT_ROOT = dirname(OPTIMAL_DIR)
const OUT = joinpath(OUTPUT_ROOT, "paper_numbers.txt")
const BENCH_N30 = [joinpath(OUTPUT_ROOT, "benchmark", "benchmark_n30_s$(s).csv") for s in (13, 19)]
const BETAS = [0.0, 0.01, 0.05, 0.1, 0.2, 0.5]
const DELAY_TYPES = ["none", "exp", "det", "insp"]
const DELAY_NAMES = Dict("none" => "None", "exp" => "Exponential", "det" => "Constant",
                         "insp" => "Inspection")
const LOG = IOBuffer()

say(xs...) = (println(xs...); println(LOG, xs...))

"Quantile with linear interpolation between order statistics (type 7)."
function quantile7(x::AbstractVector{<:Real}, p::Real)
    isempty(x) && return NaN
    y = sort(x)
    h = (length(y) - 1) * p + 1
    lo = floor(Int, h)
    hi = min(lo + 1, length(y))
    return y[lo] + (h - lo) * (y[hi] - y[lo])
end
med(x) = quantile7(x, 0.5)
share(f, x) = count(f, x) / length(x)
pct(x) = @sprintf("%.1f", 100 * x)
sig3(x) = x >= 100 ? @sprintf("%.0f", x) : @sprintf("%.3g", x)
cell(x) = "\$" * pct(med(x)) * "\$ \$(" * pct(quantile7(x, 0.9)) * ")\$"

"Ranks, with the average rank for ties."
function ranks(x::AbstractVector{<:Real})
    p = sortperm(x)
    r = zeros(length(x))
    i = 1
    while i <= length(x)
        j = i
        while j < length(x) && x[p[j + 1]] == x[p[i]]
            j += 1
        end
        for k in i:j
            r[p[k]] = (i + j) / 2
        end
        i = j + 1
    end
    return r
end

function spearman(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    rx, ry = ranks(x), ranks(y)
    dx, dy = rx .- sum(rx) / length(rx), ry .- sum(ry) / length(ry)
    return sum(dx .* dy) / sqrt(sum(dx .^ 2) * sum(dy .^ 2))
end

"Number of local minima of c (strictly below the left neighbour, not above the right one)."
function local_minima(c::AbstractVector{<:Real})
    m = 0
    for i in eachindex(c)
        left = i > firstindex(c) ? c[i - 1] : Inf
        right = i < lastindex(c) ? c[i + 1] : Inf
        (c[i] < left && c[i] <= right) && (m += 1)
    end
    return m
end

# ---------------------------------------------------------------------------------------------
# r-out-of-n policies (output/ltmc): best trigger for every row and cost pair

struct RoutRec
    system::Int
    n::Int
    model::String
    beta::Float64
    delay_type::String
    delay::Float64
    ratio::Float64
    setup::Float64
    r_opt::Int
    g_opt::Float64         # g*_r, cost of the best trigger
    g_corr::Float64        # g_corr, cost of corrective maintenance (r = n)
    saving::Float64        # 1 - g*_r / g_corr
    nonunimodal::Bool      # more than one local minimum in r
end

function rout_records(grid, systems)
    recs = RoutRec[]
    for row in grid
        res = read_row_output(row_output_path(row.row_id))
        n = ncomponents(get_system(systems, row.system))
        res.r == (res.r[1] == 0 ? collect(0:n) : collect(1:n)) ||
            error("row $(row.row_id): unexpected triggers $(res.r)")
        for ρ in COST_RATIOS, σ in COST_SETUPS
            L = ltmc(res.p_fail, res.e_nfail, res.e_trep, n; ratio = ρ, setup = σ)
            adm = [i for i in eachindex(res.r) if res.r[i] >= 1 || σ > 0]
            iopt = adm[argmin(L[adm])]
            push!(recs, RoutRec(row.system, n, row.model, parse(Float64, row.beta), row.delay_type,
                                parse(Float64, row.delay), ρ, σ, res.r[iopt], L[iopt], L[end], 1 - L[iopt] / L[end],
                                local_minima(L[adm]) > 1))
        end
    end
    return recs
end

# ---------------------------------------------------------------------------------------------
# optimal policies (output/optimal)

struct OptRec
    system::Int
    n::Int
    model::String
    beta::Float64
    delay_type::String
    delay::Float64
    state::Bool
    ratio::Float64
    setup::Float64
    g_r::Float64
    g_corr::Float64
    g_count::Float64
    g_state::Float64
    iters_count::Int
    seconds_count::Float64
    iters_state::Int
    seconds_state::Float64
    mixed::Int
    count_policy::String
end

struct OptRow
    n::Int
    state::Bool
    seconds_rout::Float64
    seconds_count36::Float64   # seconds_tables + the 36 seconds_count
    seconds_state36::Float64
    iters_state36::Int
    cached::Bool
    maxrss_mb::Float64
end

function optimal_records(grid, systems)
    recs = OptRec[]
    rows = OptRow[]
    for row in grid
        out = read_optimal_output(joinpath(OPTIMAL_DIR, "row_" * lpad(row.row_id, 5, '0') * ".csv"))
        n = ncomponents(get_system(systems, row.system))
        col(name) = findfirst(==(name), out.header)
        c = Dict(name => col(name) for name in OPTIMAL_HEADER)
        state = out.meta["state_based"] == "true"
        length(out.rows) == length(OPT_RATIOS) * length(OPT_SETUPS) ||
            error("row $(row.row_id): $(length(out.rows)) cost pairs")
        sc = 0.0; ss = 0.0; its = 0
        for f in out.rows
            F(name) = parse(Float64, f[c[name]])
            N(name) = isempty(f[c[name]]) ? 0 : parse(Int, f[c[name]])   # empty without state-based policies
            push!(recs, OptRec(row.system, n, row.model, parse(Float64, row.beta), row.delay_type,
                               parse(Float64, row.delay), state, F("ratio"), F("setup"),
                               F("ltmc_best_r"), F("ltmc_corrective"), F("ltmc_count"), F("ltmc_state"),
                               N("iters_count"), F("seconds_count"), N("iters_state"), F("seconds_state"),
                               N("state_mixed_levels"), f[c["count_policy"]]))
            sc += F("seconds_count"); ss += F("seconds_state"); its += N("iters_state")
        end
        M(name) = parse(Float64, out.meta[name])
        push!(rows, OptRow(n, state, M("seconds_rout"), M("seconds_tables") + sc, ss, its,
                           out.meta["states_cached"] == "true", M("maxrss_mb")))
    end
    return recs, rows
end

"Rows of a benchmark CSV with step state_opt: (one pass in seconds, peak memory in MiB)."
function bench_state(path::AbstractString)
    lines = [split(l, ',') for l in eachline(path) if !isempty(strip(l))]
    h = lines[1]
    i(name) = findfirst(==(name), h)
    for f in lines[2:end]
        f[i("step")] == "state_opt" || continue
        return (pass = parse(Float64, f[i("seconds")]) / parse(Int, f[i("iterations")]),
                mem = parse(Float64, f[i("maxrss_mb")]))
    end
    error("$(path): no state_opt step")
end

# ---------------------------------------------------------------------------------------------

function main()
    grid = read_grid()
    systems = load_systems()
    t0 = time()
    R = rout_records(grid, systems)
    O, rows = optimal_records(grid, systems)
    say("Read $(length(grid)) rows of output/ltmc and output/optimal in ",
        round(time() - t0; digits = 1), " s.")

    # ---------------- Section 7.5 and Table 3
    say("\n== Section 7.5: the best r-out-of-n policy")
    z = filter(x -> x.setup == 0, R)
    say("c_set = 0: best trigger r = 1 or corrective (r = n): ",
        pct(share(x -> x.r_opt == 1 || x.r_opt == x.n || x.saving <= 0, z)), "%")
    P = filter(x -> x.setup > 0, R)
    s = [x.saving for x in P]
    say("c_set > 0 ($(length(P)) combinations): saving over corrective, median ", pct(med(s)),
        "%, p90 ", pct(quantile7(s, 0.9)), "%, max ", pct(maximum(s)), "%, > 1% in ",
        pct(share(>(0.01), s)), "%")
    say("  no saving at n c_comp/c_sys = 10: ",
        pct(share(x -> x.saving <= 0, filter(x -> x.ratio == 10, P))), "%")
    for σ in OPT_SETUPS
        g = filter(x -> x.setup == σ, P)
        say("  c_set/c_sys = $(σ): median r*/n ", @sprintf("%.2f", med([x.r_opt / x.n for x in g])),
            ", r* = 1 in ", pct(share(x -> x.r_opt == 1, g)), "%")
    end
    for b in BETAS
        g = [x.saving for x in P if x.beta == b]
        say("  beta = $(b): median saving ", pct(med(g)), "%, no saving in ", pct(share(<=(0), g)), "%")
    end
    say("  no delay: median saving ", pct(med([x.saving for x in P if x.delay_type == "none"])), "%")
    for dt in ("det", "exp", "insp"), d in (0.01, 1.0)
        say("  $(dt), d = $(d): median saving ",
            pct(med([x.saving for x in P if x.delay_type == dt && x.delay == d])), "%")
    end
    say("  cost not unimodal in r (more than one local minimum): ", pct(share(x -> x.nonunimodal, P)), "%")
    r0 = filter(x -> x.r_opt == 0, P)
    say("  r = 0 optimal in $(length(r0)) combinations, delays: ",
        join(sort(unique(["$(x.delay_type) $(x.delay)" for x in r0])), ", "))

    say("\n== Table 3 (tab:rout): saving (%) of the best r-out-of-n policy over corrective maintenance,")
    say("   median (p90) by c_set/c_sys (rows) and beta = ", join(BETAS, ", "), " (columns)")
    for σ in OPT_SETUPS
        say(rpad("\$$(σ)\$", 7), " & ",
            join([cell([x.saving for x in R if x.setup == σ && x.beta == b]) for b in BETAS], " & "), " \\\\")
    end

    # ---------------- Section 7.6 and Table 4
    say("\n== Section 7.6: comparison with the optimal policies")
    gc = [1 - x.g_count / x.g_r for x in O]
    S = filter(x -> x.state, O)
    gs = [1 - x.g_state / x.g_r for x in S]
    say("ordering, largest relative violations: count > best r ", @sprintf("%.1e", maximum(x.g_count / x.g_r - 1 for x in O)),
        ", state > count ", @sprintf("%.1e", maximum(x.g_state / x.g_count - 1 for x in S)))
    say("best r-out-of-n cost recomputed from output/ltmc vs output/optimal: largest relative difference ",
        @sprintf("%.1e", rout_check(R, O)))
    say("count-based ($(length(O)) combinations): saving > 1e-9 in ", pct(share(>(1e-9), gc)),
        "%, > 0.1% in ", @sprintf("%.2f", 100 * share(>(0.001), gc)), "%, > 1% in ",
        @sprintf("%.3f", 100 * share(>(0.01), gc)), "%, max ", @sprintf("%.2f", 100 * maximum(gc)), "%")
    k = argmax(gc)
    xm = O[k]
    say("  largest: system $(xm.system), $(xm.model), beta = $(xm.beta), $(xm.delay_type) $(xm.delay), ",
        "cost pair ($(xm.ratio), $(xm.setup)), policy $(xm.count_policy) (f = 0, 1, ...), ",
        "saving ", @sprintf("%.2f", 100 * gc[k]), "% over best r, best r = corrective: ", xm.g_r == xm.g_corr)
    E = filter(x -> x.model == "Exp" && x.delay_type != "insp", O)
    say("  [Exp] with continuous inspection: largest |g*_count / g*_r - 1| = ",
        @sprintf("%.1e", maximum(abs(x.g_count / x.g_r - 1) for x in E)))
    say("state-based ($(length(S)) combinations, n <= 27): saving > 1e-9 in ", pct(share(>(1e-9), gs)),
        "%, median ", @sprintf("%.2f", 100 * med(gs)), "%, p90 ", @sprintf("%.2f", 100 * quantile7(gs, 0.9)),
        "%, p99 ", @sprintf("%.2f", 100 * quantile7(gs, 0.99)), "%, max ", @sprintf("%.2f", 100 * maximum(gs)),
        "%, > 5% in ", pct(share(>(0.05), gs)), "%")
    say("  constant delay d = 1: median saving ",
        pct(med([1 - x.g_state / x.g_r for x in S if x.delay_type == "det" && x.delay == 1.0])), "%")
    say("  median saving by n c_comp/c_sys = ", join(OPT_RATIOS, ", "), ": ",
        join([@sprintf("%.2f", 100 * med([1 - x.g_state / x.g_r for x in S if x.ratio == ρ])) for ρ in OPT_RATIOS], ", "))
    say("  at least one mixed level: ", pct(share(x -> x.mixed >= 1, S)), "%, median number of mixed levels ",
        med([x.mixed for x in S]))
    sysids = sort(unique([x.system for x in S]))
    medsys = [med([1 - x.g_state / x.g_r for x in S if x.system == i]) for i in sysids]
    EK = Float64[]
    nn = Float64[]
    for i in sysids
        sys = get_system(systems, i)
        sg = read_signature_file(signature_path(sys)).s
        push!(EK, Float64(sum(j * sg[j] for j in eachindex(sg))))
        push!(nn, ncomponents(sys))
    end
    say("  $(length(sysids)) systems, Spearman of the median saving with the signature mean sum k s_k: ",
        @sprintf("%.3f", spearman(medsys, EK)), ", with n: ", @sprintf("%.3f", spearman(medsys, nn)))

    say("\n== Table 4 (tab:optimal): saving (%) of the optimal state-based policy over the best r-out-of-n")
    say("   policy, median (p90) by delay type (rows) and beta = ", join(BETAS, ", "),
        "; last column: max saving (%) of the count-based policy (all 20 systems)")
    for dt in vcat(DELAY_TYPES, "all")
        sel(x) = dt == "all" || x.delay_type == dt
        cells = [cell([1 - x.g_state / x.g_r for x in S if sel(x) && x.beta == b]) for b in BETAS]
        cmax = maximum(1 - x.g_count / x.g_r for x in O if sel(x))
        say(rpad(get(DELAY_NAMES, dt, "All"), 12), "& ", join(cells, " & "), " & \$", pct(cmax), "\$ \\\\")
    end

    # ---------------- Section 7.7 and Table 5
    say("\n== Table 5 (tab:runtimes): medians per configuration of each n")
    say(rpad("n", 4), rpad("configs", 9), rpad("r-out (ms)", 12), rpad("count36 (ms)", 14),
        rpad("pass (s)", 11), rpad("passes", 8), rpad("36 pairs (s)", 14), "memory (GB)")
    T = Dict{Int,Any}()
    for n in sort(unique([r.n for r in rows]))
        g = filter(r -> r.n == n, rows)
        gs_ = filter(r -> r.state, g)
        st = [x.seconds_state / x.iters_state for x in S if x.n == n && x.iters_state > 0]
        rout_ms = 1e3 * med([r.seconds_rout for r in g])
        cnt_ms = 1e3 * med([r.seconds_count36 for r in g])
        mem = med([r.maxrss_mb for r in g if !r.cached]) / 1024
        if isempty(gs_)
            T[n] = (rout = rout_ms,)
            say(rpad(n, 4), rpad(length(g), 9), rpad(@sprintf("%.3g", rout_ms), 12),
                rpad(@sprintf("%.3g", cnt_ms), 14), "(state-based: see n = 30 below)")
        else
            passes = sum(r.iters_state36 for r in gs_) / (36 * length(gs_))
            s36 = med([r.seconds_state36 for r in gs_])
            T[n] = (rout = rout_ms, passes = passes, s36 = s36, mem = mem)
            say(rpad(n, 4), rpad(length(g), 9), rpad(@sprintf("%.3g", rout_ms), 12),
                rpad(@sprintf("%.3g", cnt_ms), 14), rpad(sig3(med(st)), 11),
                rpad(@sprintf("%.2f", passes), 8), rpad(sig3(s36), 14), @sprintf("%.2f", mem))
        end
    end
    b = [bench_state(p) for p in BENCH_N30]
    pr = [T[n].passes for n in (24, 26, 27)]
    lo = 36 * minimum(x.pass for x in b) * minimum(pr)
    hi = 36 * maximum(x.pass for x in b) * maximum(pr)
    say("n = 30 (benchmarks, systems 13 and 19): one pass ", join([@sprintf("%.0f", x.pass) for x in b], " and "),
        " s, memory ", join([@sprintf("%.1f", x.mem / 1024) for x in b], " and "), " GB; 36 pairs with ",
        @sprintf("%.2f", minimum(pr)), " to ", @sprintf("%.2f", maximum(pr)), " passes: ",
        @sprintf("%.2g", lo), " to ", @sprintf("%.2g", hi), " s (", @sprintf("%.2f", lo / 86400), " to ",
        @sprintf("%.2f", hi / 86400), " days)")
    say("state-based 36 pairs vs r-out-of-n: ", @sprintf("%.2g", T[22].s36 / (T[22].rout / 1e3)), " times for n = 22, ",
        @sprintf("%.2g", T[27].s36 / (T[27].rout / 1e3)), " times for n = 27")
    say("from n = 22 to n = 27: time x ", @sprintf("%.0f", T[27].s36 / T[22].s36), ", memory x ",
        @sprintf("%.0f", T[27].mem / T[22].mem))

    mkpath(dirname(OUT))
    write(OUT, take!(LOG))
    println("\nwrote $(OUT)")
end

"Largest relative difference between the best r-out-of-n cost of output/ltmc and output/optimal."
function rout_check(R::Vector{RoutRec}, O::Vector{OptRec})
    best = Dict((x.system, x.model, x.beta, x.delay_type, x.delay, x.ratio, x.setup) => x.g_opt for x in R)
    return maximum(abs(x.g_r / best[(x.system, x.model, x.beta, x.delay_type, x.delay, x.ratio, x.setup)] - 1)
                   for x in O)
end

main()
