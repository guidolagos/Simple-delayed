# Long-term mean cost (LTMC) of the r-out-of-n repair policy with delay, for
# LFMO lifetimes (Theorems thm:repairdelay and thm:LFMOrepairdelay of the paper).
#
# By the signature decomposition, for each k-out-of-n system (system failure
# time T_{k:n}) we need, with T_rep = min(T_{r:n} + D, T_{k:n}):
#
#   p_fail  = P(T_{k:n} <= T_{r:n} + D)      system failed before the repair started
#   e_nfail = E[N(T_rep)]                    failed components when the repair starts
#   e_trep  = E[T_rep]                       cycle length
#
# and the system values are the s_k-weighted sums. With c_sys the system
# failure cost and c_cmp(j) = j c_cmp,
#
#   LTMC = (c_cmp e_nfail + c_sys p_fail) / e_trep .
#
# Method. The number of failed components is a pure-death (pure-birth in the
# failed count) Markov chain with rates from `transition_rates`. For the
# k-out-of-n system we freeze the chain at every count c >= k (the system has
# failed; the count at failure may exceed k because of simultaneous failures).
# Then:
#   * phase 1 (until the trigger T_{r:n}): forward recursion of visit
#     probabilities and of the hitting distribution of {c >= min(r,k)};
#   * delay phase, from each trigger count c0 in r..k-1:
#       exponential D (rate γ): resolvent (γI - Q)^{-1};
#       deterministic D: exp(Q D) and ∫_0^D exp(Q s) ds by uniformization;
#   * inspection every Δ: the chain observed at inspection times, with
#     G = e_0 (I - exp(QΔ) restricted to counts < min(r,k))^{-1}.
# The generator is upper triangular, so every solve is a substitution, and all
# sums have nonnegative terms (no cancellation). The results agree with the
# closed forms of the paper (thm:LFMO, thm:determinstic, thm:LFMOsch) to about
# 1e-15; see test/reference/.

"""
    Delay(kind, value)

Delay between the trigger T_{r:n} and the repair start. `kind` is one of
`:none` (value ignored), `:exp` (exponential with mean `value`), `:det`
(constant `value`) or `:insp` (inspections every `value` time units; the repair
starts at the first inspection after the trigger).
"""
struct Delay
    kind::Symbol
    value::Float64
    function Delay(kind::Symbol, value::Real)
        kind in (:none, :exp, :det, :insp) || error("unknown delay kind $(kind)")
        kind == :none || value > 0 || error("delay $(kind) needs a positive value, got $(value)")
        return new(kind, kind == :none ? 0.0 : Float64(value))
    end
end

Delay(kind::AbstractString, value::Real) = Delay(Symbol(kind), value)

# Generator of the chain frozen at counts >= k. Q[c+1, c2+1].
function _frozen_generator(R::Matrix{Float64}, exits::Vector{Float64}, k::Int)
    N = size(R, 1)
    Q = zeros(Float64, N, N)
    for c in 0:(k - 1)
        Q[c + 1, c + 1] = -exits[c + 1]
        for c2 in (c + 1):(N - 1)
            Q[c + 1, c2 + 1] = R[c + 1, c2 + 1]
        end
    end
    return Q
end

# Exponential delay with rate γ: X = (γI - Q)^{-1} (rows c0 < k).
# Law of the count at the stopping time: γ X[c0, :]; occupation of transient
# counts: X[c0, c] for c < k.
function _resolvent(Q::Matrix{Float64}, k::Int, γ::Float64)
    N = size(Q, 1)
    X = zeros(Float64, N, N)
    for c0 in 0:(k - 1)
        X[c0 + 1, c0 + 1] = 1 / (γ - Q[c0 + 1, c0 + 1])
        for c in (c0 + 1):(N - 1)
            acc = 0.0
            for c1 in c0:(min(c, k) - 1)
                acc += X[c0 + 1, c1 + 1] * Q[c1 + 1, c + 1]
            end
            X[c0 + 1, c + 1] = acc / (γ - Q[c + 1, c + 1])
        end
    end
    return γ .* X, X
end

# Product of upper-triangular matrices, rows 1..nrows only (other rows: copy of A).
function _triu_mul(A::Matrix{Float64}, B::Matrix{Float64}, nrows::Int)
    N = size(A, 1)
    C = copy(A)
    for i in 1:nrows, j in i:N
        acc = 0.0
        for l in i:j
            acc += A[i, l] * B[l, j]
        end
        C[i, j] = acc
    end
    return C
end

# E = exp(Q t) and F = ∫_0^t exp(Q s) ds by uniformization with rate q ≥ max exit rate.
function _uniformize(Q::Matrix{Float64}, k::Int, q::Float64, t::Float64)
    N = size(Q, 1)
    λ = q * t
    λ < 700 || error("uniformization: q t = $(λ) too large")
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
    P = Matrix{Float64}(undef, N, N)
    for i in 1:N, j in 1:N
        P[i, j] = (i == j ? 1.0 : 0.0) + Q[i, j] / q
    end
    Pm = zeros(Float64, N, N)
    for i in 1:N
        Pm[i, i] = 1.0
    end
    E = zeros(Float64, N, N)
    F = zeros(Float64, N, N)
    for m in 0:M
        for i in 1:k, j in i:N
            E[i, j] += p[m + 1] * Pm[i, j]
            F[i, j] += tail[m + 1] / q * Pm[i, j]
        end
        m < M && (Pm = _triu_mul(Pm, P, k))
    end
    for i in 1:k                               # exact diagonal (triangular generator)
        E[i, i] = exp(Q[i, i] * t)
        F[i, i] = expm1(Q[i, i] * t) / Q[i, i]
    end
    for i in (k + 1):N                         # frozen counts stay put
        E[i, i] = 1.0
    end
    return E, F
end

# Phase 1: visit probabilities v of counts < m = min(r,k), hitting law h of
# counts >= m, and expected hitting time τ.
function _phase1(R::Matrix{Float64}, exits::Vector{Float64}, m::Int)
    N = size(R, 1)
    v = zeros(Float64, m)
    v[1] = 1.0
    for c in 1:(m - 1)
        acc = 0.0
        for c1 in 0:(c - 1)
            acc += v[c1 + 1] * R[c1 + 1, c + 1] / exits[c1 + 1]
        end
        v[c + 1] = acc
    end
    h = zeros(Float64, N)
    for c in m:(N - 1)
        acc = 0.0
        for c1 in 0:(m - 1)
            acc += v[c1 + 1] * R[c1 + 1, c + 1] / exits[c1 + 1]
        end
        h[c + 1] = acc
    end
    τ = 0.0
    for c in 0:(m - 1)
        τ += v[c + 1] / exits[c + 1]
    end
    return h, τ
end

"""
    policy_quantities(s, R, exits, delay) -> (p_fail, e_nfail, e_trep)

For the system with signature `s` (`s[k]`, k = 1..n, as Float64) and the
survivor-count chain `(R, exits)` of `transition_rates`, returns three vectors
indexed by the trigger r = 1..n with P(T_rep = T_fail), E[N(T_rep)] and
E[T_rep] for the r-out-of-n policy with `delay`.
"""
function policy_quantities(s::AbstractVector{<:Real}, R::Matrix{Float64}, exits::Vector{Float64}, delay::Delay)
    n = length(s)
    size(R) == (n + 1, n + 1) && length(exits) == n + 1 || error("rates are not for n = $(n) components")
    pf = zeros(Float64, n); en = zeros(Float64, n); et = zeros(Float64, n)
    N = n + 1
    for k in 1:n
        sk = Float64(s[k])
        sk == 0 && continue
        Q = _frozen_generator(R, exits, k)
        Mst, F = if delay.kind == :exp
            _resolvent(Q, k, 1 / delay.value)
        elseif delay.kind == :det || delay.kind == :insp
            _uniformize(Q, k, exits[1], delay.value)
        else
            Q, Q      # unused
        end
        for r in 1:n
            m = min(r, k)
            a = 0.0; b = 0.0; c_ = 0.0
            if delay.kind == :insp
                G = zeros(Float64, m)          # expected visits at inspection times
                for c in 0:(m - 1)
                    acc = c == 0 ? 1.0 : 0.0
                    for c1 in 0:(c - 1)
                        acc += G[c1 + 1] * Mst[c1 + 1, c + 1]
                    end
                    G[c + 1] = acc / (-expm1(Q[c + 1, c + 1] * delay.value))
                end
                for c1 in 0:(m - 1)
                    for c in k:(N - 1)
                        a += G[c1 + 1] * Mst[c1 + 1, c + 1]
                    end
                    for c in m:(N - 1)
                        b += G[c1 + 1] * c * Mst[c1 + 1, c + 1]
                    end
                    for c in 0:(k - 1)
                        c_ += G[c1 + 1] * F[c1 + 1, c + 1]
                    end
                end
            else
                h, τ = _phase1(R, exits, m)
                c_ = τ
                for c in k:(N - 1)
                    a += h[c + 1]
                end
                if delay.kind == :none || r >= k
                    for c in m:(N - 1)
                        b += c * h[c + 1]
                    end
                else
                    for c in k:(N - 1)
                        b += c * h[c + 1]
                    end
                    for c0 in r:(k - 1)
                        h0 = h[c0 + 1]
                        h0 == 0 && continue
                        for c in k:(N - 1)
                            a += h0 * Mst[c0 + 1, c + 1]
                        end
                        for c in c0:(N - 1)
                            b += h0 * c * Mst[c0 + 1, c + 1]
                        end
                        for c in c0:(k - 1)
                            c_ += h0 * F[c0 + 1, c + 1]
                        end
                    end
                end
            end
            pf[r] += sk * a
            en[r] += sk * b
            et[r] += sk * c_
        end
    end
    return pf, en, et
end

"""
    ltmc(p_fail, e_nfail, e_trep, n; ratio, setup = 0.0, c_sys = 1.0) -> Vector{Float64}

Long-term mean cost for each trigger r, with repair cost
c_cmp(j) = c_setup + j c_cmp for j failed components, where
`ratio` = n c_cmp / c_sys and `setup` = c_setup / c_sys (both relative to the
system failure cost):

    LTMC(r) = (c_setup + c_cmp e_nfail[r] + c_sys p_fail[r]) / e_trep[r].

`setup = 0` is the purely proportional cost c_cmp(j) = j c_cmp.
"""
function ltmc(pf::AbstractVector, en::AbstractVector, et::AbstractVector, n::Integer;
              ratio::Real, setup::Real = 0.0, c_sys::Real = 1.0)
    ratio > 0 || error("ratio must be positive")
    setup >= 0 || error("setup must be nonnegative")
    c_cmp = ratio * c_sys / n
    return (setup * c_sys .+ c_cmp .* en .+ c_sys .* pf) ./ et
end

# Cost parameters of the experiments (post-processing only):
# ratio = n c_cmp/c_sys and setup = c_setup/c_sys. ratio = 1/d for the choices
# c_sys = d n c_cmp, d = 1, 2, 5, 10, of the subsystems paper.
const COST_RATIOS = [1 / 10, 1 / 5, 1 / 2, 2 / 3, 1.0, 4 / 3, 2.0, 4.0, 10.0]
const COST_SETUPS = [0.0, 0.05, 0.1, 0.2, 0.5]
