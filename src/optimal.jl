# Optimal count-based and state-based repair policies (Section sec:optimal of
# the paper): Dinkelbach's method (alg:dinkelbach) with the passes
#   count-based: alg:count exp, its deterministic-delay version
#                (prop:after call det) and alg:count fixed (lem:period);
#   state-based: alg:state exp, alg:state det and alg:state fixed;
# plus the no-delay case (continuous inspection, D = 0: calling is repairing).
#
# Implementation choices (validated against an independent mpmath
# implementation, brute force over stopping sets and Monte Carlo simulation of
# the networks; see test/reference/optimal_reference.py and optimal_port.py):
#
# * What happens after a call, and during one inspection period, is written
#   with the conditional signature s^(f) (count-based) or the per-state
#   signature s^x (state-based) and the numbers a^(w)_i, b^(w)_i of
#   prop:after call det, for every delay (the paper notes the closed forms
#   also hold for the exponential delay). These numbers (`CallTables`) are
#   computed with the closed forms of lem:pbbs and exps in BigFloat
#   (`call_tables`, closed_forms.jl). The Float64 engine on the survivor-count
#   chain of a residual system with w components, frozen at i failures
#   (resolvent for an exponential delay, uniformization for a constant delay
#   or an inspection period; nonnegative terms only) is kept as a check
#   (`call_tables_chain`); each row records the largest difference.
# * The cost-free parts (probability of a system failure, expected number of
#   failed components, expected time) are stored separately: they are computed
#   once per (system, model, delay) and reused for every cost pair and pass.
# * State-based passes sweep the states by number of working components with
#   lem:level sums, keeping two levels in memory; at failed states the level
#   sums have closed forms, since all their sub-states are failed.
# * Calling the crew with no failed component (f = 0) is allowed only with
#   continuous inspection and a positive delay (Remark rmk:call0); repairing at
#   an inspection needs f >= 1.
#
# Indexing: f = number of failed components (0..n), w = n - f; arrays indexed
# by f use f + 1; states x are integers with bit e-1 = 1 iff component e works.

# Cost parameters of the optimal-policy experiments: setup > 0 (Remark rmk:call0).
const OPT_RATIOS = COST_RATIOS
const OPT_SETUPS = [0.05, 0.1, 0.2, 0.5]
const OPTIMAL_DIR = joinpath(REPO_DIR, "output", "optimal")
# Version of the row files written by `run_optimal_row`; readers refuse other versions.
# 1: tables and r-out-of-n quantities from the survivor-count chain (commit 85c1299).
# 2: from the closed forms in BigFloat, with the chain as a recorded check (engine_check_*).
# 3: as 2, for all systems in one directory (state-based policies for n <= OPT_STATE_MAX_N),
#    with the time of every step in full precision (seconds_rout, seconds_tables,
#    seconds_count, seconds_state_pre, seconds_state, ...), the closed forms computed
#    afresh in every row, and the peak memory of the process (maxrss_mb).
const OPTIMAL_VERSION = 3
# State-based policies are computed for the systems with n <= OPT_STATE_MAX_N (all but
# the two systems with n = 30) by scripts/run_optimal.jl.
const OPT_STATE_MAX_N = 27
# 2^n states. Peak memory of a state-based row is about 33 2^n bytes plus four
# (w+1) x binomial(n, w) Float64 matrices of the largest level: ~0.4 GB for n = 22,
# ~7 GB for n = 26, ~14 GB for n = 27, ~115 GB for n = 30 (see README).
const STATE_MAX_N_HARD = 30

"""
    Costs(n; ratio, setup, c_sys = 1.0)

Repair cost c_cmp(j) = c_set + j c_cmp for j failed components, with
c_cmp = ratio c_sys / n and c_set = setup c_sys, plus c_sys if the system failed.
"""
struct Costs
    c_sys::Float64
    c_set::Float64
    c_cmp::Float64
end

Costs(n::Integer; ratio::Real, setup::Real, c_sys::Real = 1.0) =
    Costs(Float64(c_sys), Float64(setup * c_sys), Float64(ratio * c_sys / n))

repair_cost(c::Costs, f::Integer) = c.c_set + c.c_cmp * f
failure_cost(c::Costs, f::Integer) = c.c_set + c.c_cmp * f + c.c_sys

"Whether the crew may be called (or a repair made) with f failed components."
call_allowed(kind::Symbol, f::Integer) = !(f == 0 && (kind == :none || kind == :insp))

# ---------------------------------------------------------------------------
# After-call and one-period quantities of a residual system with w components
# ---------------------------------------------------------------------------

"""
    CallTables

For w = 1..n and i = 1..w, for an LFMO system with w components whose system
fails at the i-th failure (an i-out-of-w:F system, i.e. T_{i:w}) and a delay
D started at time 0:

- `pf[w, i]`  = P(T_{i:w} <= D)
- `enA[w, i]` = E[N^{(w)}(min(T_{i:w}, D))]
- `enF[w, i]` = E[N^{(w)}(T_{i:w}); T_{i:w} <= D]
- `ell[w, i]` = E[min(T_{i:w}, D)]
- `pN[w, j+1]` = P(N^{(w)}(D) = j), j = 0..w

with D exponential with mean `delay.value`, or constant `delay.value` (also for
an inspection period). For `:none` everything is zero (D = 0).
"""
struct CallTables
    kind::Symbol
    pf::Matrix{Float64}
    enA::Matrix{Float64}
    enF::Matrix{Float64}
    ell::Matrix{Float64}
    pN::Matrix{Float64}
end

function _poisson_tail(λ::Float64)
    λ < 700 || error("uniformization: rate × time = $(λ) too large")
    M = ceil(Int, λ + 10 * sqrt(λ) + 25)
    p = zeros(Float64, M + 1)                 # p[m+1] = P(Poisson(λ) = m)
    p[1] = exp(-λ)
    for m in 1:M
        p[m + 1] = p[m] * λ / m
    end
    tail = zeros(Float64, M + 1)              # tail[m+1] = P(Poisson(λ) > m)
    acc = 0.0
    for m in M:-1:0
        tail[m + 1] = acc
        acc += p[m + 1]
    end
    return p, tail
end

"""
    call_tables(cf::ClosedForms, delay; df = nothing) -> CallTables

The after-call and one-period tables from the closed forms (prop:after call det,
lem:period, lem:pbbs and exps) in BigFloat, rounded to Float64:
pf = 1 - P(D < T_{i:w}), enA = Σ_j j P(N^{(w)}(min(T_{i:w}, D)) = j),
enF = Σ_{j >= i} j P(N^{(w)}(T_{i:w}) = j, T_{j:w} < D), ell = E[min(T_{i:w}, D)]
and pN = P(N^{(w)}(D) = j). `df` may pass `delay_forms(cf, delay)`.
"""
function call_tables(cf::ClosedForms, delay::Delay; df::Union{DelayForms,Nothing} = nothing)
    n = cf.n
    pf = zeros(n, n); enA = zeros(n, n); enF = zeros(n, n); ell = zeros(n, n)
    pN = zeros(n, n + 1)
    delay.kind == :none && return CallTables(:none, pf, enA, enF, ell, pN)
    df === nothing && (df = delay_forms(cf, delay))
    df.kind == delay.kind || error("delay forms are for $(df.kind), not $(delay.kind)")
    setprecision(BigFloat, cf.prec) do
        for w in 1:n
            for i in 1:w
                pf[w, i] = Float64(1 - df.pnot[w][i])
                enA[w, i] = Float64(df.EA[w][i])
                enF[w, i] = Float64(df.T1[w][i])
                ell[w, i] = Float64(df.ell[w][i])
            end
            for j in 0:w
                pN[w, j + 1] = Float64(df.PN[w][j + 1])
            end
        end
    end
    _check_tables(pf, enA, enF, ell, pN, delay)
    return CallTables(delay.kind, pf, enA, enF, ell, pN)
end

"""
    call_tables_chain(R, exits, delay) -> CallTables

Check engine. From the survivor-count chain `(R, exits)` of `transition_rates` for n
components: the residual chain with w working components is the block of counts
f..n, f = n - w. For each (w, i) the chain is frozen at i further failures and
run for the delay from count 0: resolvent for an exponential delay, and
uniformization with rate Ψ(w) (the largest exit rate) for a constant one.
"""
function call_tables_chain(R::Matrix{Float64}, exits::Vector{Float64}, delay::Delay)
    n = length(exits) - 1
    pf = zeros(n, n); enA = zeros(n, n); enF = zeros(n, n); ell = zeros(n, n)
    pN = zeros(n, n + 1)
    delay.kind == :none && return CallTables(:none, pf, enA, enF, ell, pN)
    E0 = zeros(n + 1); F0 = zeros(n + 1); X = zeros(n + 1)
    v = zeros(n + 1); nv = zeros(n + 1)
    for w in 1:n
        f = n - w
        q = exits[f + 1]                       # Ψ(w)
        p, tail = delay.kind == :exp ? (Float64[], Float64[]) : _poisson_tail(q * delay.value)
        for i in 1:w
            fill!(E0, 0.0); fill!(F0, 0.0)
            if delay.kind == :exp
                γ = 1 / delay.value
                X[1] = 1 / (γ + exits[f + 1])
                for c in 1:w
                    acc = 0.0
                    for c1 in 0:(min(c, i) - 1)
                        acc += X[c1 + 1] * R[f + c1 + 1, f + c + 1]
                    end
                    X[c + 1] = acc / (γ + (c < i ? exits[f + c + 1] : 0.0))
                end
                for c in 0:w
                    E0[c + 1] = γ * X[c + 1]
                    F0[c + 1] = c < i ? X[c + 1] : 0.0
                end
            else
                t = delay.value
                fill!(v, 0.0); v[1] = 1.0
                M = length(p) - 1
                for m in 0:M
                    pm = p[m + 1]; tm = tail[m + 1] / q
                    for c in 0:w
                        E0[c + 1] += pm * v[c + 1]
                        F0[c + 1] += tm * v[c + 1]
                    end
                    m == M && break
                    fill!(nv, 0.0)
                    for c in 0:w
                        vc = v[c + 1]
                        vc == 0.0 && continue
                        if c < i
                            nv[c + 1] += vc * (1 - exits[f + c + 1] / q)
                            for c2 in (c + 1):w
                                nv[c2 + 1] += vc * R[f + c + 1, f + c2 + 1] / q
                            end
                        else
                            nv[c + 1] += vc
                        end
                    end
                    v, nv = nv, v
                end
                E0[1] = exp(-exits[f + 1] * t)                 # exact for the starting count
                F0[1] = -expm1(-exits[f + 1] * t) / exits[f + 1]
                for c in i:w
                    F0[c + 1] = 0.0
                end
            end
            a = 0.0; b = 0.0; bf = 0.0; e = 0.0
            for c in 0:w
                b += c * E0[c + 1]
                if c >= i
                    a += E0[c + 1]
                    bf += c * E0[c + 1]
                else
                    e += F0[c + 1]
                end
            end
            pf[w, i] = a; enA[w, i] = b; enF[w, i] = bf; ell[w, i] = e
            if i == w
                for j in 0:w
                    pN[w, j + 1] = E0[j + 1]
                end
            end
        end
    end
    _check_tables(pf, enA, enF, ell, pN, delay)
    return CallTables(delay.kind, pf, enA, enF, ell, pN)
end

function _check_tables(pf, enA, enF, ell, pN, delay::Delay)
    n = size(pf, 1); tol = 1e-9
    for w in 1:n
        abs(sum(pN[w, 1:(w + 1)]) - 1) <= tol || error("call tables: P(N(D) = ·) does not sum to 1 (w = $(w))")
        for i in 1:w
            -tol <= pf[w, i] <= 1 + tol || error("call tables: pf[$(w),$(i)] = $(pf[w, i])")
            -tol <= enF[w, i] <= enA[w, i] + tol || error("call tables: enF > enA at ($(w),$(i))")
            enA[w, i] <= w + tol || error("call tables: enA[$(w),$(i)] = $(enA[w, i]) > w")
            0 < ell[w, i] <= delay.value * (1 + tol) || error("call tables: ell[$(w),$(i)] = $(ell[w, i])")
            # P(T_{i:w} <= D) is nonincreasing in i
            i == 1 || pf[w, i] <= pf[w, i - 1] + tol || error("call tables: pf not nonincreasing in i")
        end
    end
    return true
end

# ---------------------------------------------------------------------------
# Count-based policies
# ---------------------------------------------------------------------------

"""
    CountModel(W, cf::ClosedForms, delay)
    CountModel(W, R, exits, delay)            (check engine: chain tables)
    CountModel(W, R, exits, tables::CallTables)

Everything a count-based pass needs, for a system with working-state counts
`W[f+1]` (signature file), survivor-count chain rates `(R, exits)` and the
after-call tables of the delay (closed forms by default). Per
f with Sbar_f > 0 (index f + 1):

- continuous inspection (`:none`, `:exp`, `:det`): after a call, `pf` =
  P(system fails before the crew arrives), `en` = E[further failures],
  `ell` = E[time until the end of the cycle] (prop:after call det, eq:call det);
- inspections (`:insp`): during one period without repair, `pf` = P(system
  fails), `en` = E[further failures; system fails], `ell` = E[length of the
  period], and `q[f+1][j+1]` = P(next inspection finds f + j failed and the
  system working), j = 0..w-1 (lem:period).
"""
struct CountModel
    n::Int
    kind::Symbol
    Sbar::Vector{Float64}
    R::Matrix{Float64}
    exits::Vector{Float64}
    tables::CallTables
    pf::Vector{Float64}
    en::Vector{Float64}
    ell::Vector{Float64}
    q::Vector{Vector{Float64}}
end

CountModel(W::AbstractVector{<:Integer}, cf::ClosedForms, delay::Delay) =
    CountModel(W, cf.R, cf.exits, call_tables(cf, delay))

CountModel(W::AbstractVector{<:Integer}, R::Matrix{Float64}, exits::Vector{Float64}, delay::Delay) =
    CountModel(W, R, exits, call_tables_chain(R, exits, delay))

function CountModel(W::AbstractVector{<:Integer}, R::Matrix{Float64}, exits::Vector{Float64}, T::CallTables)
    n = length(W) - 1
    size(R) == (n + 1, n + 1) || error("rates are not for n = $(n) components")
    size(T.pf) == (n, n) || error("call tables are not for n = $(n) components")
    delay = (kind = T.kind,)
    Sbar = [Float64(big(W[f + 1]) // binomial(big(n), f)) for f in 0:n]
    pf = zeros(n + 1); en = zeros(n + 1); ell = zeros(n + 1)
    q = [Float64[] for _ in 0:n]
    for f in 0:(n - 1)
        Sbar[f + 1] > 0 || continue
        w = n - f
        for i in 1:w                            # conditional signature s^(f)_i
            si = (Sbar[f + i] - Sbar[f + i + 1]) / Sbar[f + 1]
            pf[f + 1] += si * T.pf[w, i]
            en[f + 1] += si * (delay.kind == :insp ? T.enF[w, i] : T.enA[w, i])
            ell[f + 1] += si * T.ell[w, i]
        end
        if delay.kind == :insp
            q[f + 1] = [T.pN[w, j + 1] * Sbar[f + j + 1] / Sbar[f + 1] for j in 0:(w - 1)]
            # Samaniego: survive the period (and land somewhere) or fail in it
            abs(sum(q[f + 1]) + pf[f + 1] - 1) <= 1e-9 ||
                error("count model: one-period probabilities do not sum to 1 at f = $(f)")
        end
    end
    return CountModel(n, delay.kind, Sbar, R, exits, T, pf, en, ell, q)
end

"""
    count_pass(cm, costs, g; forced = nothing) -> (actions, E C, E T)

One pass of alg:count exp (or its deterministic / no-delay / inspection
versions): backward recursion over f = n-1..0 of the Bellman equation
eq:bellman count (eq:bellman fixed for inspections). `actions[f+1]` for
f = 0..n-1 is 1 (call / repair), 0 (continue / wait) or -1 (Sbar_f = 0).
With `forced = f -> Bool` the actions are fixed where calling is allowed,
which evaluates that policy; e.g. `f -> f >= r` is the r-out-of-n policy.
Ties are resolved by calling, as in the paper.
"""
function count_pass(cm::CountModel, costs::Costs, g::Real; forced = nothing)
    n = cm.n; Sb = cm.Sbar
    c = zeros(n + 1); l = zeros(n + 1)
    act = fill(Int8(-1), n)
    for f in (n - 1):-1:0
        Sb[f + 1] > 0 || continue
        w = n - f
        if cm.kind == :insp
            qf = cm.q[f + 1]
            cw = (costs.c_set + costs.c_sys + costs.c_cmp * f) * cm.pf[f + 1] + costs.c_cmp * cm.en[f + 1]
            lw = cm.ell[f + 1]
            for j in 1:(w - 1)
                qj = qf[j + 1]
                qj > 0 || continue
                cw += qj * c[f + j + 1]
                lw += qj * l[f + j + 1]
            end
            c_cont = cw / (1 - qf[1]); l_cont = lw / (1 - qf[1])
            c_call = repair_cost(costs, f); l_call = 0.0
        else
            c_call = costs.c_set + costs.c_cmp * (f + cm.en[f + 1]) + costs.c_sys * cm.pf[f + 1]
            l_call = cm.ell[f + 1]
            ex = cm.exits[f + 1]
            c_cont = 0.0; l_cont = 1 / ex
            for j in 1:w
                Pj = cm.R[f + 1, f + j + 1] / ex
                ps = Sb[f + j + 1] / Sb[f + 1]
                c_cont += Pj * (ps * c[f + j + 1] + (1 - ps) * failure_cost(costs, f + j))
                l_cont += Pj * ps * l[f + j + 1]
            end
        end
        call = if !call_allowed(cm.kind, f)
            false
        elseif forced !== nothing
            forced(f)::Bool
        else
            c_call - g * l_call <= c_cont - g * l_cont
        end
        act[f + 1] = call ? Int8(1) : Int8(0)
        c[f + 1] = call ? c_call : c_cont
        l[f + 1] = call ? l_call : l_cont
    end
    return act, c[1], l[1]
end

"""
    dinkelbach(pass, g0; maxit = 100, rtol = 1e-12) -> (g, policy, iterations)

alg:dinkelbach. `pass(g)` returns (policy, E C, E T) for a policy minimizing
E[C - g T]; `g0` is the long-term mean cost of some policy (the best
r-out-of-n one, or the count-based optimum for state-based policies). Stops
when F(g) = E C - g E T >= -rtol max(1, |E C|).
"""
function dinkelbach(pass, g0::Real; maxit::Integer = 100, rtol::Real = 1e-12)
    g = Float64(g0)
    for it in 1:maxit
        pol, c, l = pass(g)
        l > 0 || error("Dinkelbach: expected cycle length $(l) <= 0")
        if c - g * l >= -rtol * max(1.0, abs(c))
            return (c / l, pol, it)
        end
        g = c / l
    end
    error("Dinkelbach did not converge in $(maxit) iterations")
end

"""
    is_threshold(actions, kind) -> Bool

Whether a count-based policy is an r-out-of-n policy on the reachable counts
(Sbar_f > 0) where calling is allowed: once it calls, it calls at every larger f.
"""
function is_threshold(act::AbstractVector{Int8}, kind::Symbol)
    seen = false
    for f in 0:(length(act) - 1)
        a = act[f + 1]
        (a < 0 || !call_allowed(kind, f)) && continue
        if a == 1
            seen = true
        elseif seen
            return false
        end
    end
    return true
end

"String of a count policy: 'C' call, '.' continue, '-' unreachable, for f = 0..n-1."
policy_string(act::AbstractVector{Int8}) = String([a < 0 ? '-' : (a == 1 ? 'C' : '.') for a in act])

# ---------------------------------------------------------------------------
# State-based policies
# ---------------------------------------------------------------------------

"""
    StateSpace(sys)

The 2^n states of `sys`: structure function `phi[x+1]`, the states with w
working components in increasing order `levels[w+1]`, and the position
`rank[x+1]` of each state in its level.
"""
struct StateSpace
    n::Int
    phi::BitVector
    levels::Vector{Vector{UInt32}}
    rank::Vector{Int32}
end

function StateSpace(sys::TwoTerminalSystem)
    n = ncomponents(sys)
    n <= STATE_MAX_N_HARD || error("$(sys.name): n = $(n) is too large for state-based policies")
    cs = compile_system(sys)
    N = 1 << n
    phi = falses(N)
    levels = [Vector{UInt32}(undef, binomial(n, w)) for w in 0:n]
    filled = zeros(Int, n + 1)
    rank = Vector{Int32}(undef, N)
    for x in 0:(N - 1)
        phi[x + 1] = works(cs, UInt64(x))
        w = count_ones(x)
        filled[w + 1] += 1
        levels[w + 1][filled[w + 1]] = UInt32(x)
        rank[x + 1] = Int32(filled[w + 1])
    end
    return StateSpace(n, phi, levels, rank)
end

"Number of working states with f failed components, f = 0..n (must equal W)."
working_counts(S::StateSpace) = [count(x -> S.phi[x + 1], S.levels[S.n - f + 1]) for f in 0:S.n]

@inline _lowbit(y::UInt32) = y & (~y + one(UInt32))

"""
    StatePre

Per-state version of the `CountModel` quantities `pf`, `en`, `ell`, with the
per-state signature s^x instead of s^(f) (prop:after call det, lem:period),
indexed by `x + 1`; zero at failed states.
"""
struct StatePre
    pf::Vector{Float64}
    en::Vector{Float64}
    ell::Vector{Float64}
end

"""
    state_precompute(S, cm) -> StatePre

Sweeps the states by number of working components; the counts W^x_i of
def:cond sig are the level sums of lem:level sums with u = Phi (integers).
"""
function state_precompute(S::StateSpace, cm::CountModel)
    n = S.n; T = cm.tables; N = 1 << n
    n == cm.n || error("state space and count model have different n")
    pf = zeros(N); en = zeros(N); ell = zeros(N)
    binf = [Float64(binomial(w, i)) for w in 0:n, i in 0:n]       # binf[w+1, i+1]
    prevW = zeros(Int64, 1, 1)
    for w in 0:n
        L = S.levels[w + 1]
        curW = zeros(Int64, w + 1, length(L))
        for idx in eachindex(L)
            x = L[idx]
            S.phi[x + 1] || continue
            curW[1, idx] = 1
            y = x
            while y != 0
                b = _lowbit(y)
                y ⊻= b
                r = S.rank[(x ⊻ b) + 1]
                @inbounds for j in 1:w
                    curW[j + 1, idx] += prevW[j, r]
                end
            end
            for j in 1:w
                qt, rm = divrem(curW[j + 1, idx], j)
                rm == 0 || error("level sums: W^x not divisible by $(j)")
                curW[j + 1, idx] = qt
            end
            a = 0.0; bb = 0.0; e = 0.0
            for i in 1:w                        # per-state signature s^x_i
                si = curW[i, idx] / binf[w + 1, i] - curW[i + 1, idx] / binf[w + 1, i + 1]
                a += si * T.pf[w, i]
                bb += si * (cm.kind == :insp ? T.enF[w, i] : T.enA[w, i])
                e += si * T.ell[w, i]
            end
            pf[x + 1] = a; en[x + 1] = bb; ell[x + 1] = e
        end
        prevW = curW
    end
    return StatePre(pf, en, ell)
end

"""
    state_pass(S, cm, pre, costs, g; forced = nothing) -> ((calls, working), E C, E T)

One pass of alg:state exp / alg:state det / alg:state fixed. `calls[f+1]` and
`working[f+1]` count, among the working states with f failed components, those
at which the policy calls (repairs). `forced = (x, f) -> Bool` fixes the
actions where calling is allowed (policy evaluation).
"""
function state_pass(S::StateSpace, cm::CountModel, pre::StatePre, costs::Costs, g::Real; forced = nothing)
    n = S.n; kind = cm.kind; T = cm.tables
    binf = [Float64(binomial(w, i)) for w in 0:n, i in 0:n]
    calls = zeros(Int, n + 1); working = zeros(Int, n + 1)
    pj = zeros(n + 1)
    prevC = zeros(1, 1); prevL = zeros(1, 1)
    for w in 0:n
        f = n - w
        L = S.levels[w + 1]; m = length(L)
        curC = Matrix{Float64}(undef, w + 1, m)
        curL = Matrix{Float64}(undef, w + 1, m)
        # probability of moving to one given sub-state with j fewer working components
        if w >= 1
            if kind == :insp
                for j in 0:w
                    pj[j + 1] = T.pN[w, j + 1] / binf[w + 1, j + 1]
                end
            else
                ex = cm.exits[f + 1]
                for j in 1:w
                    pj[j + 1] = cm.R[f + 1, f + j + 1] / ex / binf[w + 1, j + 1]
                end
            end
        end
        for idx in 1:m
            x = L[idx]
            if !S.phi[x + 1]                    # closed form: all sub-states are failed
                for j in 0:w
                    curC[j + 1, idx] = kind == :insp ? 0.0 : binf[w + 1, j + 1] * failure_cost(costs, f + j)
                    curL[j + 1, idx] = 0.0
                end
                continue
            end
            for j in 1:w
                curC[j + 1, idx] = 0.0; curL[j + 1, idx] = 0.0
            end
            y = x
            while y != 0
                b = _lowbit(y)
                y ⊻= b
                r = S.rank[(x ⊻ b) + 1]
                @inbounds for j in 1:w
                    curC[j + 1, idx] += prevC[j, r]
                    curL[j + 1, idx] += prevL[j, r]
                end
            end
            sc = 0.0; sl = 0.0
            for j in 1:w
                curC[j + 1, idx] /= j
                curL[j + 1, idx] /= j
                sc += pj[j + 1] * curC[j + 1, idx]
                sl += pj[j + 1] * curL[j + 1, idx]
            end
            if kind == :insp
                p0 = T.pN[w, 1]
                c_cont = ((costs.c_set + costs.c_sys + costs.c_cmp * f) * pre.pf[x + 1] +
                         costs.c_cmp * pre.en[x + 1] + sc) / (1 - p0)
                l_cont = (pre.ell[x + 1] + sl) / (1 - p0)
                c_call = repair_cost(costs, f); l_call = 0.0
            else
                c_call = costs.c_set + costs.c_cmp * (f + pre.en[x + 1]) + costs.c_sys * pre.pf[x + 1]
                l_call = pre.ell[x + 1]
                c_cont = sc; l_cont = 1 / cm.exits[f + 1] + sl
            end
            call = if !call_allowed(kind, f)
                false
            elseif forced !== nothing
                forced(x, f)::Bool
            else
                c_call - g * l_call <= c_cont - g * l_cont
            end
            working[f + 1] += 1
            calls[f + 1] += call
            curC[1, idx] = call ? c_call : c_cont
            curL[1, idx] = call ? l_call : l_cont
        end
        prevC = curC; prevL = curL
    end
    return (calls, working), prevC[1, 1], prevL[1, 1]
end

"String of a state policy: calls/working for f = 0..n-1, separated by ';'."
state_calls_string(calls, working) = join(("$(calls[f + 1])/$(working[f + 1])" for f in 0:(length(calls) - 2)), ";")

"Number of levels f at which the state-based policy calls at some, but not all, working states."
mixed_levels(calls, working) = count(f -> 0 < calls[f] < working[f], eachindex(calls))

# ---------------------------------------------------------------------------
# One grid row
# ---------------------------------------------------------------------------

"""
    rout_quantities(s, cf::ClosedForms, delay) -> (rs, p_fail, e_nfail, e_trep)

The r-out-of-n quantities for every admissible trigger: r = 0..n with an
exponential or deterministic delay, r = 1..n otherwise (closed forms;
`policy_quantities_all`).
"""
rout_quantities(s, cf::ClosedForms, delay::Delay; df = nothing) = policy_quantities_all(s, cf, delay; df = df)

"""
    rout_quantities_chain(s, R, exits, delay) -> (rs, p_fail, e_nfail, e_trep)

The same with the survivor-count chain engine (check).
"""
function rout_quantities_chain(s, R, exits, delay::Delay)
    n = length(s)
    pf, en, et = policy_quantities_chain(s, R, exits, delay)
    rs = collect(1:n)
    if delay.kind in (:exp, :det)
        pf0, en0, et0 = policy_quantities_r0_chain(s, R, exits, delay)
        pushfirst!(rs, 0); pushfirst!(pf, pf0); pushfirst!(en, en0); pushfirst!(et, et0)
    end
    return rs, pf, en, et
end

const OPTIMAL_HEADER = ["ratio", "setup", "best_r", "ltmc_best_r", "ltmc_corrective",
                        "ltmc_count", "iters_count", "seconds_count", "count_policy", "count_threshold",
                        "ltmc_state", "iters_state", "state_mixed_levels", "state_calls", "seconds_state"]

_seconds(t0::UInt64) = (time_ns() - t0) / 1e9
_maxrss_mb() = Sys.maxrss() / 2^20

"""
    run_optimal_row(row; systems, dir = OPTIMAL_DIR, state_max_n = 14, cache = Dict()) -> Symbol

For one grid row and every cost pair (OPT_RATIOS × OPT_SETUPS): the best
r-out-of-n policy, the optimal count-based policy, and, if n <= `state_max_n`,
the optimal state-based policy. Writes `dir/row_<row_id>.csv` atomically;
returns `:skipped` if it exists. Errors if state <= count <= best r fails.

All production values come from the closed forms in BigFloat: the r-out-of-n
quantities (`rout_quantities`) and the after-call and one-period tables
(`call_tables`), which feed both the count-based and the state-based passes.
The survivor-count chain is run afterwards only as a check, timed separately.

Timings (seconds, full precision; `time_ns`), in the metadata:
`seconds_rout` = closed forms (with the 512-bit rates), delay forms and r-out-of-n
quantities, all computed afresh for this row; `seconds_tables` = after-call tables
and count model, given the delay forms; `seconds_states` = state space (0 if it
was taken from `cache`, see `states_cached`); `seconds_state_pre` = per-state
precomputation; `seconds_chain` = the check. Per cost pair: `seconds_count` and
`seconds_state` = Dinkelbach's method (`iters_*` passes). `maxrss_mb` is the peak
resident memory of the process when the row is finished. Only the state space is
cached between rows (by `cache`), since it depends only on the system.
"""
function run_optimal_row(row; systems::Vector{TwoTerminalSystem}, dir::AbstractString = OPTIMAL_DIR,
                         state_max_n::Integer = 14, cache::Dict = Dict(), verbose::Bool = false)
    path = row_output_path(row.row_id; dir = dir)
    isfile(path) && return :skipped
    t0 = time_ns()
    sys = get_system(systems, row.system)
    n = ncomponents(sys)
    sigfile = read_signature_file(signature_path(sys))
    s = Float64.(sigfile.s)
    model = lfmo_model(row.model, row.beta)
    delay = Delay(row.delay_type, parse(Float64, row.delay))
    # r-out-of-n quantities: closed forms, computed afresh (not cached) so that their
    # time is that of one configuration
    tr = time_ns()
    cf = closed_forms(model, n)
    df = delay.kind == :none ? nothing : delay_forms(cf, delay)
    rs, pf, en, et = rout_quantities(sigfile.s, cf, delay; df = df)
    seconds_rout = _seconds(tr)
    # after-call and one-period tables (closed forms, sharing the delay forms)
    tt = time_ns()
    T = call_tables(cf, delay; df = df)
    cm = CountModel(sigfile.W, cf.R, cf.exits, T)
    seconds_tables = _seconds(tt)
    # independent check: the survivor-count chain engine (Float64), not used below
    tk = time_ns()
    _, cpf, cen, cet = rout_quantities_chain(s, cf.R, cf.exits, delay)
    Tc = call_tables_chain(cf.R, cf.exits, delay)
    check_prob = max(max_diff(pf, cpf; relative = false), max_diff(T.pf, Tc.pf; relative = false),
                     max_diff(T.pN, Tc.pN; relative = false), max_diff(T.enF, Tc.enF; relative = false))
    check_rel = max(max_diff(en, cen; relative = true), max_diff(et, cet; relative = true),
                    max_diff(T.enA, Tc.enA; relative = true), max_diff(T.ell, Tc.ell; relative = true))
    seconds_chain = _seconds(tk)
    check_prob <= 1e-9 && check_rel <= 1e-9 ||
        error("row $(row.row_id): closed forms and chain disagree (probabilities $(check_prob), relative $(check_rel))")
    do_state = n <= state_max_n
    S = pre = nothing
    seconds_states = 0.0; states_cached = false; seconds_state_pre = 0.0
    if do_state
        key = (:states, sys.id)
        if haskey(cache, key)
            S = cache[key]; states_cached = true
        else
            empty!(cache)                       # at most one state space in memory
            GC.gc()
            ts = time_ns()
            S = StateSpace(sys)
            seconds_states = _seconds(ts)
            working_counts(S) == sigfile.W || error("$(sys.name): state space disagrees with signature file")
            cache[key] = S
        end
        tp = time_ns()
        pre = state_precompute(S, cm)
        seconds_state_pre = _seconds(tp)
    end
    tol = 1e-9
    lines = String[]
    for ρ in OPT_RATIOS, σ in OPT_SETUPS
        costs = Costs(n; ratio = ρ, setup = σ)
        Lr = [(costs.c_set + costs.c_cmp * en[i] + costs.c_sys * pf[i]) / et[i] for i in eachindex(rs)]
        ib = argmin(Lr)
        tc = time_ns()
        gc, actc, itc = dinkelbach(g -> count_pass(cm, costs, g), Lr[ib])
        secc = _seconds(tc)
        gc <= Lr[ib] * (1 + tol) || error("row $(row.row_id): count-based optimum $(gc) > best r-out-of-n $(Lr[ib])")
        gs = NaN; its = 0; mixed = ""; callstr = ""; secs = 0.0
        if do_state
            ts = time_ns()
            gs, (calls, working), its = dinkelbach(g -> state_pass(S, cm, pre, costs, g), gc)
            secs = _seconds(ts)
            gs <= gc * (1 + tol) || error("row $(row.row_id): state-based optimum $(gs) > count-based $(gc)")
            mixed = string(mixed_levels(calls, working))
            callstr = state_calls_string(calls, working)
            verbose && println("  ratio $(ρ) setup $(σ): state-based $(its) passes, $(round(secs; digits = 2)) s")
        end
        push!(lines, join((ρ, σ, rs[ib], repr(Lr[ib]), repr(Lr[end]), repr(gc), itc, repr(secc),
                           policy_string(actc), is_threshold(actc, delay.kind), repr(gs), its, mixed, callstr,
                           repr(secs)), ","))
    end
    seconds = _seconds(t0)

    io = IOBuffer()
    println(io, "# SimpleDelayed optimal count-based and state-based policies (sec:optimal)")
    for (k, val) in pairs(row)
        println(io, "# ", k, " = ", val)
    end
    println(io, "# name = ", sys.name)
    println(io, "# n_components = ", n)
    println(io, "# param = ", string(Float64(model.param)))
    println(io, "# mean_batch_n = ", mean_batch(model, n))
    println(io, "# state_based = ", do_state)
    println(io, "# n_working_states = ", sum(sigfile.W))
    println(io, "# seconds = ", repr(seconds))
    println(io, "# seconds_rout = ", repr(seconds_rout))
    println(io, "# seconds_tables = ", repr(seconds_tables))
    println(io, "# seconds_states = ", repr(seconds_states))
    println(io, "# states_cached = ", states_cached)
    println(io, "# seconds_state_pre = ", repr(seconds_state_pre))
    println(io, "# maxrss_mb = ", repr(round(_maxrss_mb(); digits = 1)))
    println(io, "# engine = closed forms, BigFloat ", cf.prec, " bits")
    println(io, "# engine_check_prob = ", repr(check_prob))
    println(io, "# engine_check_rel = ", repr(check_rel))
    println(io, "# seconds_chain = ", repr(seconds_chain))
    println(io, "# optimal_version = ", OPTIMAL_VERSION)
    println(io, "# julia_version = ", VERSION)
    println(io, "# host = ", _hostname())
    println(io, join(OPTIMAL_HEADER, ","))
    for line in lines
        println(io, line)
    end
    atomic_write(path, take!(io))
    return :done
end

"""
    optimal_warmup(systems; state_max_n = OPT_STATE_MAX_N) -> seconds

Runs `run_optimal_row` on four rows of the smallest system (LFMO-C, β = 0.2, no
delay and the three delays with d = 0.1) in a temporary directory, which is then
deleted, so that the timings of the production rows exclude compilation.
"""
function optimal_warmup(systems::Vector{TwoTerminalSystem}; state_max_n::Integer = OPT_STATE_MAX_N)
    t0 = time_ns()
    small = argmin(s -> ncomponents(s), systems).id
    rows = [r for r in read_grid() if r.system == small && r.model == "LFMO-C" && r.beta == "0.2" &&
            (r.delay_type == "none" || r.delay == "0.1")]
    length(rows) == 4 || error("warm-up: expected 4 rows, found $(length(rows))")
    dir = mktempdir()
    try
        for r in rows
            run_optimal_row(r; systems = systems, dir = dir, state_max_n = state_max_n, cache = Dict())
        end
    finally
        rm(dir; recursive = true, force = true)
    end
    return _seconds(t0)
end

"""
    read_optimal_output(path; versions = (OPTIMAL_VERSION,)) -> (meta, header, rows)

Errors unless the file's `optimal_version` is in `versions`. Other versions are
accepted only by scripts/compare_optimal.jl, which compares old and new outputs.
"""
function read_optimal_output(path::AbstractString; versions = (OPTIMAL_VERSION,))
    meta = Dict{String,String}()
    for line in eachline(path)
        startswith(line, "#") || continue
        m = match(r"^#\s*([A-Za-z_]+)\s*=\s*(.*)$", line)
        m === nothing || (meta[m.captures[1]] = strip(m.captures[2]))
    end
    get(meta, "optimal_version", "") in string.(versions) ||
        error("$(path): optimal_version $(get(meta, "optimal_version", "?")), expected one of $(versions)")
    header, rows = read_simple_csv(path)
    expected = meta["optimal_version"] == string(OPTIMAL_VERSION) ? OPTIMAL_HEADER :
               filter(!=("seconds_count"), OPTIMAL_HEADER)
    header == expected || error("$(path): unexpected header")
    return (meta = meta, header = header, rows = rows)
end
