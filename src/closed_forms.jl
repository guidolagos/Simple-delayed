# Closed forms of Section sec:LFMO of the paper, evaluated in multiple precision.
#
# Lemma lem:pbbs and exps gives the probabilities and expectations of an LFMO
# vector with w components as alternating sums over l = 0..w of binomial
# coefficients times Ψ(l), 1/Ψ(l), exp(-tΨ(l)), 1/(γ+Ψ(l)), ... Theorems thm:LFMO
# (exponential delay), thm:determinstic (constant delay) and thm:LFMOsch
# (inspections) combine them into the three cost-free quantities of the
# r-out-of-n policy for every k-out-of-n system (`policy_quantities`), and
# Proposition prop:after call det and Lemma lem:period into the after-call and
# one-period quantities of the optimal policies (`call_tables`).
#
# Precision. The terms of the sums add up, in absolute value, to about
# C(w,f) 2^f <= 3^w times the result, so in Float64 they lose about w log10(3)
# digits (15 for w = 30; the long-term mean cost of a 30-component system would
# be wrong by up to ~10%). We therefore evaluate every closed form in BigFloat
# with CLOSED_PRECISION = 256 bits, which leaves more than 200 bits for n <= 30,
# and round the final results to Float64.
#
# Check. The survivor-count chain (`policy_quantities_chain`, `call_tables_chain`)
# evaluates the same quantities in Float64 with nonnegative terms only. Every
# row written by `run_row` and `run_optimal_row` records the largest difference
# between the two engines (`engine_check` in the row metadata).
#
# Notation: w = number of components of the (residual) LFMO vector, f, i, j =
# numbers of failed components, r = trigger, k = index of the k-out-of-n system.
# Vectors indexed by a count x start at x = 0, so they use x + 1.

const CLOSED_PRECISION = 256   # bits

"""
    ClosedForms

Delay-independent closed forms of an LFMO model with `n` components, in
BigFloat (`CLOSED_PRECISION` bits):

- `Psi[x+1]` = Ψ(x), x = 0..n;
- `psiP[w, f]` = Ψ(w) P^{(w)}_f (def:P), the rate at which f of w working
  components fail together;
- `pntn[f1, f2]` = P(N^{(n)}(T_{f1:n}) = f2), f1 <= f2 (eq:PNTrn);
- `ETn[f]` = E T_{f:n} (eq:ETrn);

plus the Float64 rates `R`, `exits` of `transition_rates`, used by the
continuation steps of the optimal-policy passes, and the precision `prec` (bits)
used for these and all the derived closed forms.
"""
struct ClosedForms
    model::LFMOModel
    n::Int
    Psi::Vector{BigFloat}
    psiP::Matrix{BigFloat}
    pntn::Matrix{BigFloat}
    ETn::Vector{BigFloat}
    R::Matrix{Float64}
    exits::Vector{Float64}
    prec::Int
end

# Σ_{l=w-g}^{w} C(w,l) C(l,w-g) (-1)^{l-(w-g)} φ(l), with phi[l+1] = φ(l)
# (eq:PNwt, def:H and the inner sums of eq:PNTrn).
function _csum_A(phi::Vector{BigFloat}, w::Int, g::Int)
    s = BigFloat(0)
    for l in (w - g):w
        t = (binomial(w, l) * binomial(l, w - g)) * phi[l + 1]
        s = isodd(l - (w - g)) ? s - t : s + t
    end
    return s
end

# Σ_{l=w-f+1}^{w} C(w,l) C(l-1,w-f) (-1)^{l-(w-f+1)} φ(l) (eq:ETrn, eq:Kclosed,
# eq:EminDT, eq:Ddet_Tkn).
function _csum_E(phi::Vector{BigFloat}, w::Int, f::Int)
    s = BigFloat(0)
    for l in (w - f + 1):w
        t = (binomial(w, l) * binomial(l - 1, w - f)) * phi[l + 1]
        s = isodd(l - (w - f + 1)) ? s - t : s + t
    end
    return s
end

# Ψ(w) P^{(w)}_f = Σ_{l=w-f}^{w} (-1)^{l-(w-f)+1} C(w,f) C(f,w-l) Ψ(l) (def:P).
function _csum_psiP(Psi::Vector{BigFloat}, w::Int, f::Int)
    s = BigFloat(0)
    for l in (w - f):w
        t = binomial(f, w - l) * Psi[l + 1]
        s = isodd(l - (w - f) + 1) ? s - t : s + t
    end
    return binomial(w, f) * s
end

# T[f1, f2] = Σ_{g=0}^{f1-1} Ψ(w-g) P^{(w-g)}_{f2-g} A_g(w; φ) for 1 <= f1 <= f2 <= w,
# and 0 for f1 > f2: eq:PNTrn (φ = 1/Ψ), eq:PNTrn t (φ = (1-e^{-tΨ})/Ψ) and
# eq:PNTrn Dind (φ = 1/(γ+Ψ)).
function _csum_pnt(psiP::Matrix{BigFloat}, phi::Vector{BigFloat}, w::Int)
    T = fill(BigFloat(0), w, w)
    Ag = [_csum_A(phi, w, g) for g in 0:(w - 1)]
    for f2 in 1:w
        acc = BigFloat(0)
        for f1 in 1:f2
            g = f1 - 1
            acc = acc + psiP[w - g, f2 - g] * Ag[g + 1]
            T[f1, f2] = acc
        end
    end
    return T
end

# phi[l+1] = fun(l) for l = 1..n, and `at0` for l = 0 (NaN where φ(0) is never used).
_csum_phi(fun, n::Int; at0 = BigFloat(NaN)) = [l == 0 ? BigFloat(at0) : BigFloat(fun(l)) for l in 0:n]

"""
    closed_forms(model, n; precision = CLOSED_PRECISION) -> ClosedForms
"""
function closed_forms(m::LFMOModel, n::Integer; precision::Integer = CLOSED_PRECISION)
    n = Int(n)
    R, exits = transition_rates(m, n)
    setprecision(BigFloat, precision) do
        Psi = [_psi_big(m.name, m.param, x) for x in 0:n]
        psiP = fill(BigFloat(0), n, n)
        for w in 1:n, f in 1:w
            psiP[w, f] = _csum_psiP(Psi, w, f)
        end
        phi_inv = _csum_phi(l -> 1 / Psi[l + 1], n)
        pntn = _csum_pnt(psiP, phi_inv, n)
        ETn = [_csum_E(phi_inv, n, f) for f in 1:n]
        return ClosedForms(m, n, Psi, psiP, pntn, ETn, R, exits, Int(precision))
    end
end

"""
    DelayForms

Closed forms that depend on the delay, for every w = 1..n and i = 1..w, with D
the delay (exponential with rate `gamma`, or constant `t`; for inspections,
t = Δ):

- `PN[w][j+1]` = P(N^{(w)}(D) = j), j = 0..w (eq:PNwt, mixed over D if exponential);
- `pnt3[w][f1, f2]` = P(N^{(w)}(T_{f1:w}) = f2, T_{f2:w} < D), f1 <= f2
  (eq:PNTrn t, eq:PNTrn Dind);
- `ell[w][i]` = E[min(T_{i:w}, D)] (eq:Kclosed, eq:EminDT);
- `pnot[w][i]` = P(D < T_{i:w}) (eq:Ddet_Tkn; 1 - P(T_{i:w} <= t) for constant D);
- `MA[w][i]`, `EA[w][i]` = Σ_j P(N^{(w)}(min(T_{i:w}, D)) = j) (≈ 1) and
  Σ_j j P(N^{(w)}(min(T_{i:w}, D)) = j), with the probabilities of eq:PNcDdetT
  (exponential) or of thm:determinstic (constant);
- `T1[w][i]` = E[N^{(w)}(T_{i:w}); T_{i:w} < D];
- `H[f+1]` = E Σ_{m >= 0} 1{N^{(n)}(mΔ) = f}, f = 0..n-1 (def:H; inspections only).
"""
struct DelayForms
    kind::Symbol
    gamma::BigFloat
    PN::Vector{Vector{BigFloat}}
    pnt3::Vector{Matrix{BigFloat}}
    ell::Vector{Vector{BigFloat}}
    pnot::Vector{Vector{BigFloat}}
    MA::Vector{Vector{BigFloat}}
    EA::Vector{Vector{BigFloat}}
    T1::Vector{Vector{BigFloat}}
    H::Vector{BigFloat}
end

# P(N^{(w)}(min(T_{i:w}, D)) = j): eq:PNcDdetT for an exponential delay (with
# P(... = 0) = γ/(γ+Ψ(w)), thm:LFMO), and for a constant delay P(N(D) = j) if
# j < i, P(N(T_{i:w}) = j, T_{j:w} < D) otherwise (thm:determinstic).
function _csum_pmin(kind::Symbol, gamma::BigFloat, Psi::Vector{BigFloat}, PN::Vector{BigFloat},
                    T::Matrix{BigFloat}, w::Int, i::Int, j::Int)
    if kind == :exp
        j == 0 && return gamma / (gamma + Psi[w + 1])
        base = T[min(i, j), j]
        return j < i ? base * (gamma / (gamma + Psi[w - j + 1])) : base
    end
    return j < i ? PN[j + 1] : T[i, j]
end

"""
    delay_forms(cf, delay) -> DelayForms

For `delay.kind` in (:exp, :det, :insp); the value of the delay is the mean of
the exponential delay, the constant delay, or the inspection period.
"""
function delay_forms(cf::ClosedForms, delay::Delay)
    delay.kind in (:exp, :det, :insp) || error("delay_forms: no delay forms for $(delay.kind)")
    n = cf.n
    setprecision(BigFloat, cf.prec) do
        Psi = cf.Psi
        v = BigFloat(delay.value)
        kind = delay.kind
        if kind == :exp
            gamma = 1 / v
            phiPN = _csum_phi(l -> gamma / (gamma + Psi[l + 1]), n; at0 = 1)
            phiT = _csum_phi(l -> 1 / (gamma + Psi[l + 1]), n)
        else
            gamma = BigFloat(0)
            phiPN = _csum_phi(l -> exp(-v * Psi[l + 1]), n; at0 = 1)
            phiT = _csum_phi(l -> -expm1(-v * Psi[l + 1]) / Psi[l + 1], n)
        end
        PN = Vector{Vector{BigFloat}}(undef, n); pnt3 = Vector{Matrix{BigFloat}}(undef, n)
        ell = Vector{Vector{BigFloat}}(undef, n); pnot = Vector{Vector{BigFloat}}(undef, n)
        MA = Vector{Vector{BigFloat}}(undef, n); EA = Vector{Vector{BigFloat}}(undef, n)
        T1 = Vector{Vector{BigFloat}}(undef, n)
        for w in 1:n
            pn = [_csum_A(phiPN, w, j) for j in 0:w]
            T = _csum_pnt(cf.psiP, phiT, w)
            el = [_csum_E(phiT, w, i) for i in 1:w]
            po = Vector{BigFloat}(undef, w)
            if kind == :exp
                for i in 1:w
                    po[i] = gamma * el[i]                  # eq:Ddet_Tkn = γ × eq:EminDT
                end
            else
                acc = BigFloat(0)
                for i in 1:w
                    acc = acc + pn[i]                      # Σ_{j < i} P(N(t) = j)
                    po[i] = acc
                end
            end
            ma = Vector{BigFloat}(undef, w); ea = Vector{BigFloat}(undef, w); t1 = Vector{BigFloat}(undef, w)
            for i in 1:w
                sm = BigFloat(0); se = BigFloat(0); st = BigFloat(0)
                for j in 0:w
                    p = _csum_pmin(kind, gamma, Psi, pn, T, w, i, j)
                    sm = sm + p
                    se = se + j * p
                end
                for j in i:w
                    st = st + j * T[i, j]
                end
                ma[i] = sm; ea[i] = se; t1[i] = st
            end
            PN[w] = pn; pnt3[w] = T; ell[w] = el; pnot[w] = po; MA[w] = ma; EA[w] = ea; T1[w] = t1
        end
        H = BigFloat[]
        if kind == :insp
            phiH = _csum_phi(l -> 1 / (-expm1(-v * Psi[l + 1])), n)
            H = [_csum_A(phiH, n, f) for f in 0:(n - 1)]
        end
        return DelayForms(kind, gamma, PN, pnt3, ell, pnot, MA, EA, T1, H)
    end
end

# P(N^{(n)}(T_{r:n}) = l), with T_{0:n} := 0 (thm:LFMO, thm:determinstic for r = 0).
function _pnt_n(cf::ClosedForms, r::Int, l::Int)
    r == 0 && return l == 0 ? BigFloat(1) : BigFloat(0)
    return l >= r ? cf.pntn[r, l] : BigFloat(0)
end

"""
    kout_quantities(cf, df_or_nothing, kind) -> (rs, pf, en, et)

For every trigger r in `rs` and every k = 1..n, the three cost-free quantities
of the r-out-of-n policy for the k-out-of-n system, as BigFloat matrices
indexed [index of r in rs, k]:

    pf = P(T_{k:n} <= T_{r:n} + D),  en = E[N^{(n)}(T_rep)],  et = E[T_rep],

with T_rep = min(T_{r:n} + D, T_{k:n}), from thm:LFMO (`:exp`),
thm:determinstic (`:det`), thm:LFMOsch (`:insp`), and the case D = 0 of
thm:determinstic (`:none`). `rs` is 0..n for `:exp` and `:det`, 1..n otherwise.
The sums over j of eqtheo:PNTrn and eqtheo:PNTrn det are grouped by l:
Σ_j j P(N(T_rep) = j) = Σ_{j >= k} j P(N(T_{r:n}) = j)
                       + Σ_{l=r}^{k-1} P(N(T_{r:n}) = l) (l MA + EA)[n-l][k-l].
"""
function kout_quantities(cf::ClosedForms, df::Union{DelayForms,Nothing}, kind::Symbol)
    n = cf.n
    rs = kind in (:exp, :det) ? collect(0:n) : collect(1:n)
    (kind == :none) == (df === nothing) || error("kout_quantities: delay forms needed iff there is a delay")
    setprecision(BigFloat, cf.prec) do
        pf = fill(BigFloat(0), length(rs), n); en = fill(BigFloat(0), length(rs), n)
        et = fill(BigFloat(0), length(rs), n)
        for (ir, r) in enumerate(rs), k in 1:n
            if r >= k                                   # the system fails before the trigger
                s = BigFloat(0)
                for j in k:n
                    s = s + j * cf.pntn[k, j]
                end
                pf[ir, k] = BigFloat(1); en[ir, k] = s; et[ir, k] = cf.ETn[k]
                continue
            end
            if kind == :none                            # repair at T_{r:n} unless N(T_{r:n}) >= k
                p = BigFloat(0); s = BigFloat(0)
                for j in r:n
                    q = cf.pntn[r, j]
                    s = s + j * q
                    j >= k && (p = p + q)
                end
                pf[ir, k] = p; en[ir, k] = s; et[ir, k] = cf.ETn[r]
            elseif kind == :insp                        # thm:LFMOsch, r < k
                pT = BigFloat(0); e = BigFloat(0); s = BigFloat(0)
                for rp in 0:(r - 1)
                    w = n - rp; h = df.H[rp + 1]
                    s1 = BigFloat(0)
                    for l in r:(k - 1)
                        s1 = s1 + df.PN[w][l - rp + 1]
                    end
                    pT = pT + h * s1
                    e = e + h * df.ell[w][k - rp]
                    s2 = BigFloat(0)
                    for j in r:(k - 1)
                        s2 = s2 + j * df.PN[w][j - rp + 1]
                    end
                    for j in k:n
                        s2 = s2 + j * df.pnt3[w][k - rp, j - rp]
                    end
                    s = s + h * s2
                end
                pf[ir, k] = 1 - pT; en[ir, k] = s; et[ir, k] = e
            else                                        # thm:LFMO / thm:determinstic, r < k
                pT = BigFloat(0); extra = BigFloat(0); s = BigFloat(0)
                for j in k:n
                    s = s + j * _pnt_n(cf, r, j)
                end
                for l in r:(k - 1)
                    q = _pnt_n(cf, r, l)
                    iszero(q) && continue
                    w = n - l; i = k - l
                    pT = pT + q * df.pnot[w][i]
                    extra = extra + q * df.ell[w][i]
                    s = s + q * (l * df.MA[w][i] + df.EA[w][i])
                end
                et0 = r == 0 ? BigFloat(0) : cf.ETn[r]
                pf[ir, k] = 1 - pT; en[ir, k] = s
                et[ir, k] = et0 + (kind == :exp ? pT / df.gamma : extra)
            end
        end
        return rs, pf, en, et
    end
end

"""
    policy_quantities_all(s, cf, delay) -> (rs, p_fail, e_nfail, e_trep)

Signature decomposition (thm:repairdelay, thm:LFMOrepairdelay) of the closed
forms: for every admissible trigger r in `rs` (0..n for exponential and
constant delays, 1..n otherwise), the Float64 values of P(T_rep = T_fail),
E[N(T_rep)] and E[T_rep] for the system with signature `s`. `df` may pass the
`delay_forms(cf, delay)` already computed (they are shared with `call_tables`).
"""
function policy_quantities_all(s::AbstractVector{<:Real}, cf::ClosedForms, delay::Delay;
                               df::Union{DelayForms,Nothing} = nothing)
    n = cf.n
    length(s) == n || error("signature has $(length(s)) entries, closed forms are for n = $(n)")
    if delay.kind == :none
        df = nothing
    elseif df === nothing
        df = delay_forms(cf, delay)
    else
        df.kind == delay.kind || error("delay forms are for $(df.kind), not $(delay.kind)")
    end
    rs, P, E, T = kout_quantities(cf, df, delay.kind)
    setprecision(BigFloat, cf.prec) do
        sb = [BigFloat(x) for x in s]
        pf = zeros(length(rs)); en = zeros(length(rs)); et = zeros(length(rs))
        for ir in eachindex(rs)
            a = BigFloat(0); b = BigFloat(0); c = BigFloat(0)
            for k in 1:n
                iszero(sb[k]) && continue
                a = a + sb[k] * P[ir, k]; b = b + sb[k] * E[ir, k]; c = c + sb[k] * T[ir, k]
            end
            pf[ir] = Float64(a); en[ir] = Float64(b); et[ir] = Float64(c)
        end
        return rs, pf, en, et
    end
end

"""
    policy_quantities(s, cf::ClosedForms, delay) -> (p_fail, e_nfail, e_trep)

As `policy_quantities_all`, for r = 1..n only.
"""
function policy_quantities(s::AbstractVector{<:Real}, cf::ClosedForms, delay::Delay)
    rs, pf, en, et = policy_quantities_all(s, cf, delay)
    i1 = findfirst(==(1), rs)
    return pf[i1:end], en[i1:end], et[i1:end]
end

"""
    policy_quantities_r0(s, cf::ClosedForms, delay) -> (p_fail, e_nfail, e_trep)

As `policy_quantities_all`, for r = 0 only (exponential or constant delay).
"""
function policy_quantities_r0(s::AbstractVector{<:Real}, cf::ClosedForms, delay::Delay)
    delay.kind in (:exp, :det) || error("r = 0 needs continuous inspection with a delay (exp or det), got $(delay.kind)")
    rs, pf, en, et = policy_quantities_all(s, cf, delay)
    return pf[1], en[1], et[1]
end

"""
    max_diff(a, b; relative) -> Float64

Largest difference between two arrays of the same shape: absolute, or relative
to |b| (entries with b = 0 are compared in absolute value).
"""
function max_diff(a::AbstractArray, b::AbstractArray; relative::Bool)
    size(a) == size(b) || error("max_diff: shapes differ")
    m = 0.0
    for (x, y) in zip(a, b)
        d = abs(Float64(x) - Float64(y))
        relative && y != 0 && (d /= abs(Float64(y)))
        m = max(m, d)
    end
    return m
end
