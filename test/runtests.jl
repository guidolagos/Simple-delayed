using Test
using SimpleDelayed

# Build a small system from an edge list.
mksys(edges; directed = false, name = "test") =
    TwoTerminalSystem(0, name, name, directed, "Start", "End", [(String(u), String(v)) for (u, v) in edges])

sig(edges; directed = false) = signature(count_working_states(mksys(edges; directed = directed)))

# Independent structure function: plain depth-first search on adjacency lists.
function naive_works(sys::TwoTerminalSystem, x::UInt64)
    adj = Dict{String,Vector{String}}()
    for (e, (u, v)) in enumerate(sys.edges)
        (x >> (e - 1)) & 1 == 1 || continue
        push!(get!(adj, u, String[]), v)
        sys.directed || push!(get!(adj, v, String[]), u)
    end
    seen = Set([sys.source])
    stack = [sys.source]
    while !isempty(stack)
        u = pop!(stack)
        u == sys.target && return true
        for v in get(adj, u, String[])
            v in seen || (push!(seen, v); push!(stack, v))
        end
    end
    return false
end

# Deterministic pseudo-random 64-bit states (xorshift64), no RNG dependency.
function xorshift_states(seed::UInt64, count::Integer, n::Integer)
    x = seed
    out = UInt64[]
    mask = (UInt64(1) << n) - one(UInt64)
    for _ in 1:count
        x ⊻= x << 13; x ⊻= x >> 7; x ⊻= x << 17
        push!(out, x & mask)
    end
    return out
end

# Number of minimal cut sets: failed sets F such that the system fails, and
# repairing any single component of F makes it work. Exhaustive; small n only.
function count_minimal_cuts(sys::TwoTerminalSystem)
    n = ncomponents(sys)
    cs = compile_system(sys)
    N = UInt64(1) << n
    phi = [works(cs, x) for x in UInt64(0):(N - one(UInt64))]
    full = N - one(UInt64)
    count = 0
    for x in UInt64(0):full
        phi[x + 1] && continue
        failed = full & ~x
        minimal = true
        for e in 0:(n - 1)
            if (failed >> e) & 1 == 1 && !phi[(x | (UInt64(1) << e)) + 1]
                minimal = false
                break
            end
        end
        count += minimal
    end
    return count
end

function read_reference(path)
    ref = Dict{Int,Vector{Int}}()
    for line in eachline(path)
        (startswith(line, "#") || startswith(line, "system") || isempty(strip(line))) && continue
        f = split(line, ',')
        push!(get!(ref, parse(Int, f[1]), Int[]), parse(Int, f[3]))
    end
    return ref
end


function read_ltmc_reference(path)
    header, rows = SimpleDelayed.read_simple_csv(path)
    cases = Dict{Tuple{Int,String,String,String,String},Vector{NTuple{4,Float64}}}()
    for r in rows
        key = (parse(Int, r[1]), r[2], r[3], r[4], r[5])
        push!(get!(cases, key, NTuple{4,Float64}[]),
              (parse(Float64, r[6]), parse(Float64, r[7]), parse(Float64, r[8]), parse(Float64, r[9])))
    end
    return cases
end

# Production engine: closed forms of the paper in BigFloat (closed_forms.jl).
quantities0(sysid, model, beta, kind, value; systems) = begin
    sys = get_system(systems, sysid)
    s = read_signature_file(signature_path(sys)).s
    policy_quantities_r0(s, closed_forms(lfmo_model(model, beta), ncomponents(sys)), Delay(kind, value))
end

# Parametric pass helpers for the optimal-policy tests; cf.R, cf.exits are the chain rates
function setup_row(systems, sysid, model, beta)
    sys = get_system(systems, sysid)
    sf = read_signature_file(signature_path(sys))
    cf = closed_forms(lfmo_model(model, beta), ncomponents(sys))
    return sys, sf, cf
end

quantities(sysid, model, beta, kind, value; systems) = begin
    sys = get_system(systems, sysid)
    s = read_signature_file(signature_path(sys)).s
    policy_quantities(s, closed_forms(lfmo_model(model, beta), ncomponents(sys)), Delay(kind, value))
end

# Check engine: survivor-count chain in Float64 (ltmc.jl).
quantities_chain(sysid, model, beta, kind, value; systems) = begin
    sys = get_system(systems, sysid)
    s = Float64.(read_signature_file(signature_path(sys)).s)
    R, ex = transition_rates(lfmo_model(model, beta), ncomponents(sys))
    rout_quantities_chain(s, R, ex, Delay(kind, value))
end

@testset "SimpleDelayed" begin

    @testset "known signatures" begin
        # series system of 4 components: fails at the first failure
        @test sig([("Start", "a"), ("a", "b"), ("b", "c"), ("c", "End")]) == [1, 0, 0, 0]
        # parallel system of 3 components: fails at the last failure
        @test sig([("Start", "End"), ("Start", "End"), ("Start", "End")]) == [0, 0, 1]
        # component 1 in series with the parallel pair (2, 3)
        @test sig([("Start", "a"), ("a", "End"), ("a", "End")]) == [1 // 3, 2 // 3, 0]
        # bridge system (Samaniego 2007): s = (0, 1/5, 3/5, 1/5, 0)
        bridge = [("Start", "1"), ("Start", "2"), ("1", "2"), ("1", "End"), ("2", "End")]
        @test sig(bridge) == [0, 1 // 5, 3 // 5, 1 // 5, 0]
        # directedness is honoured: an arc End -> Start is irrelevant when directed
        @test sig([("Start", "End"), ("End", "Start")]; directed = true) == [1 // 2, 1 // 2]
        @test sig([("Start", "End"), ("End", "Start")]; directed = false) == [0, 1]
        # dangling component (semi-coherent, irrelevant component)
        @test sig([("Start", "a"), ("a", "End"), ("b", "a")]) == [2 // 3, 1 // 3, 0]
    end

    @testset "state convention and counts" begin
        bridge = mksys([("Start", "1"), ("Start", "2"), ("1", "2"), ("1", "End"), ("2", "End")])
        @test works(bridge, state_from_working(5, [1, 4]))          # path Start-1-End
        @test works(bridge, state_from_working(5, [2, 3, 4]))       # path Start-2-1-End
        @test !works(bridge, state_from_working(5, [1, 2, 3]))      # nothing reaches End
        @test count_working_states(bridge) == [1, 5, 8, 2, 0, 0]
        # result independent of the chunking of the state space
        @test count_working_states(bridge; nchunks = 1) == count_working_states(bridge; nchunks = 7)
        @test_throws ErrorException signature([0, 0, 0])            # Φ(1) = 0 is rejected
    end

    systems = load_systems("yeh2024")

    @testset "Yeh et al. (2024) data files" begin
        @test [s.id for s in systems] == collect(1:20)
        @test [ncomponents(s) for s in systems] ==
              [5, 8, 8, 9, 12, 14, 21, 13, 12, 12, 12, 13, 30, 26, 14, 21, 27, 22, 30, 24]
        @test [nnodes(s) for s in systems] ==
              [4, 6, 5, 6, 9, 7, 11, 9, 8, 8, 7, 8, 16, 21, 9, 10, 18, 13, 20, 16]
        @test !any(s.directed for s in systems)     # all undirected, as in Yeh et al.
        # CSV edge lists coincide with the raw file (networks 21 and 22 of the raw
        # file are dropped; see data/systems/README.md)
        raw = parse_yeh_raw(joinpath(SimpleDelayed.SYSTEMS_DIR, "yeh2024_raw.txt"))
        @test [r.id for r in raw] == collect(1:22)
        for s in systems
            r = raw[s.id]
            @test (r.id, r.label, r.edges) == (s.id, s.label, s.edges)
        end
    end

    @testset "minimal cuts vs. Table 9 of Yeh et al. (2024), systems with n <= 16" begin
        # Table 9 values, except systems 6 (Table 9: 25) and 10 (Table 9 repeats
        # the row of system 15); see data/systems/README.md.
        expected = Dict(1 => 4, 2 => 9, 3 => 8, 4 => 9, 5 => 28, 6 => 26, 8 => 24,
                        9 => 19, 10 => 20, 11 => 20, 12 => 29, 15 => 25)
        for sys in systems
            ncomponents(sys) <= 16 || continue
            @test count_minimal_cuts(sys) == expected[sys.id]
        end
    end

    @testset "structure function vs. naive DFS (random states, all Yeh systems)" begin
        for sys in systems
            cs = compile_system(sys)
            states = xorshift_states(UInt64(0x9E3779B97F4A7C15) + UInt64(sys.id), 3000, ncomponents(sys))
            @test all(works(cs, x) == naive_works(sys, x) for x in states)
        end
    end

    @testset "counts vs. independent C++ reference (Yeh systems with n <= 22)" begin
        ref = read_reference(joinpath(@__DIR__, "reference", "yeh2024_W_reference.csv"))
        @test sort(collect(keys(ref))) == collect(1:20)
        for sys in systems
            ncomponents(sys) <= 22 || continue
            W = count_working_states(sys)
            @test W == ref[sys.id]
            @test sum(signature(W)) == 1
        end
    end

    @testset "signature file round trip" begin
        sys = get_system(systems, 1)
        W = count_working_states(sys)
        dir = mktempdir()
        path = joinpath(dir, "bridge.csv")
        write_signature_file(path, sys, W; seconds = 0.0, threads = 1)
        r = read_signature_file(path)
        @test r.W == W
        @test r.s == [0, 1 // 5, 3 // 5, 1 // 5, 0]
        @test r.meta["name"] == "yeh2024_01"
        @test readdir(dir) == ["bridge.csv"]          # no leftover temporary file
    end

    @testset "LFMO models: calibration by the pairwise beta factor" begin
        for name in MODELS[2:end], b in ("0.01", "0.05", "0.1", "0.2", "0.5")
            m = lfmo_model(name, b)
            @test abs(psi(m, 1) - 1) < 1e-60
            @test abs(2 - psi(m, 2) - parse(BigFloat, b)) < 1e-60
            @test psi(m, 0) == 0
        end
        # closed forms of the parameter
        @test Float64(lfmo_model("LFMO-C", "0.2").param) ≈ 0.2 / 1.8
        @test Float64(lfmo_model("BFR", "0.2").param) == 0.2
        @test Float64(lfmo_model("LSh", "0.2").param) == 0.2
        @test Float64(lfmo_model("LFMO-S", "0.5").param) ≈ log2(1.5)
        # beta = 1/2 gives eta = (sqrt(5)-1)/2 for Gamma and 1/sqrt(12) for inverse Gaussian
        @test Float64(lfmo_model("LFMO-G", "0.5").param) ≈ (sqrt(5) - 1) / 2
        @test Float64(lfmo_model("LFMO-IG", "0.5").param) ≈ 1 / sqrt(12)
        # the subsystems-paper parameters: 10x/(x+9), log(1+x(e-1)), sqrt(2x+1/4)-1/2, sqrt(x)
        @test Float64(lfmo_model("LFMO-C", string(2 - 20 / 11)).param) ≈ 0.1
        # inverse Gaussian cannot exceed beta = 2 - sqrt(2)
        @test_throws ErrorException lfmo_model("LFMO-IG", "0.6")
        @test_throws ErrorException lfmo_model("Exp", "0.1")
        @test_throws ErrorException lfmo_model("LFMO-C", "0")
    end

    @testset "transition rates of the survivor-count chain" begin
        n = 30
        for name in MODELS, b in (name == "Exp" ? ("0",) : ("0.01", "0.5"))
            m = lfmo_model(name, b)
            R, ex = transition_rates(m, n)
            @test all(>=(0), R)
            for c in 0:(n - 1)
                @test isapprox(sum(R[c + 1, :]), ex[c + 1]; rtol = 1e-13)
                @test ex[c + 1] ≈ Float64(psi(m, n - c))
                # mean batch = w Ψ(1) / Ψ(w)
                w = n - c
                @test isapprox(sum((f) * R[c + 1, c + f + 1] for f in 1:w) / ex[c + 1], mean_batch(m, w); rtol = 1e-12)
            end
        end
        # iid exponential: single failures only
        R, ex = transition_rates(lfmo_model("Exp", "0"), 6)
        @test all(R[c + 1, c2 + 1] == (c2 == c + 1 ? 6 - c : 0) for c in 0:5, c2 in 0:6)
        # lethal shock: drift (1-θ) per component plus killing θ of all survivors
        R, ex = transition_rates(lfmo_model("LSh", "0.2"), 5)
        @test R[1, 2] ≈ 5 * 0.8 && R[1, 6] ≈ 0.2 && all(abs(R[1, j]) < 1e-15 for j in 3:5)
        # BFR: binomial batches, shock rate 1/p
        R, ex = transition_rates(lfmo_model("BFR", "0.2"), 4)
        @test all(R[1, f + 1] ≈ binomial(4, f) * 0.2^f * 0.8^(4 - f) / 0.2 for f in 1:4)
    end

    @testset "LTMC quantities: known cases (iid exponential, no delay)" begin
        # bridge system (system 1): s = (0, 1/5, 3/5, 1/5, 0); E T_{k:5} = sum_{i=5-k+1}^{5} 1/i
        pf, en, et = quantities(1, "Exp", "0", "none", 0; systems = systems)
        ET = cumsum([1 / i for i in 5:-1:1])
        @test pf ≈ [0, 0.2, 0.8, 1, 1]
        # single failures: N(T_rep) = min(r, K) with P(K = k) = s_k
        @test en ≈ [1, 2, 2.8, 3, 3]
        @test et[1] ≈ ET[1]
        @test et[5] ≈ 0.2 * ET[2] + 0.6 * ET[3] + 0.2 * ET[4]
        # series system of n components: fails at the first failure, whatever r and delay
        sys = mksys([("Start", "a"), ("a", "b"), ("b", "End")])
        s = signature(count_working_states(sys))
        cf = closed_forms(lfmo_model("LFMO-G", "0.2"), 3)
        for d in (Delay(:none, 0), Delay(:exp, 0.3), Delay(:det, 0.3), Delay(:insp, 0.3))
            pf, en, et = policy_quantities(s, cf, d)
            @test all(≈(1), pf)
            @test all(≈(1 / cf.exits[1]), et)
            pf, en, et = policy_quantities_chain(Float64.(s), cf.R, cf.exits, d)
            @test all(≈(1), pf)
            @test all(≈(1 / cf.exits[1]), et)
        end
    end

    @testset "LTMC quantities vs. closed forms of the paper (test/reference)" begin
        ref = read_ltmc_reference(joinpath(@__DIR__, "reference", "ltmc_reference.csv"))
        @test length(ref) == 43
        @test count(((k, v),) -> v[1][1] == 0, ref) == 20      # r = 0 for the exp and det cases
        for ((sysid, model, beta, kind, value), vals) in ref
            # the chain engine against the same independent reference
            rsc, pfc, enc, etc = quantities_chain(sysid, model, beta, kind, parse(Float64, value); systems = systems)
            for (r, p, e, t) in vals
                i = findfirst(==(Int(r)), rsc)
                @test isapprox(pfc[i], p; rtol = 1e-10, atol = 1e-13)
                @test isapprox(enc[i], e; rtol = 1e-10)
                @test isapprox(etc[i], t; rtol = 1e-10)
            end
            pf, en, et = quantities(sysid, model, beta, kind, parse(Float64, value); systems = systems)
            for (r, p, e, t) in vals
                r = Int(r)
                if r == 0
                    pf0, en0, et0 = quantities0(sysid, model, beta, kind, parse(Float64, value); systems = systems)
                    @test isapprox(pf0, p; rtol = 1e-10, atol = 1e-13)
                    @test isapprox(en0, e; rtol = 1e-10)
                    @test isapprox(et0, t; rtol = 1e-10)
                else
                    @test isapprox(pf[r], p; rtol = 1e-10, atol = 1e-13)
                    @test isapprox(en[r], e; rtol = 1e-10)
                    @test isapprox(et[r], t; rtol = 1e-10)
                end
            end
        end
    end

    @testset "closed forms vs. survivor-count chain (all delays, n up to 30)" begin
        # (system, n): 1 (5), 5 (12), 7 (21), 13 (30), 19 (30)
        for (sysid, model, beta) in ((1, "LFMO-C", "0.2"), (5, "BFR", "0.01"), (7, "LFMO-IG", "0.5"),
                                     (13, "LFMO-S", "0.2"), (19, "LSh", "0.5"), (13, "Exp", "0"),
                                     (19, "LFMO-G", "0.05"), (13, "BFR", "0.2"))
            sys, sf, cf = setup_row(systems, sysid, model, beta)
            for d in (Delay(:none, 0), Delay(:exp, 0.01), Delay(:exp, 1.0), Delay(:det, 0.01), Delay(:det, 1.0),
                      Delay(:insp, 0.01), Delay(:insp, 1.0))
                rs, pf, en, et = rout_quantities(sf.s, cf, d)
                rsc, pfc, enc, etc = rout_quantities_chain(Float64.(sf.s), cf.R, cf.exits, d)
                @test rs == rsc
                @test max_diff(pf, pfc; relative = false) < 1e-12
                @test max_diff(en, enc; relative = true) < 1e-12
                @test max_diff(et, etc; relative = true) < 1e-12
                T = call_tables(cf, d); Tc = call_tables_chain(cf.R, cf.exits, d)
                for (a, b) in ((T.pf, Tc.pf), (T.pN, Tc.pN), (T.enF, Tc.enF))
                    @test max_diff(a, b; relative = false) < 1e-12
                end
                for (a, b) in ((T.enA, Tc.enA), (T.ell, Tc.ell))
                    @test max_diff(a, b; relative = true) < 1e-12
                end
            end
        end
        # small n, every model family and beta: closed forms and chain agree to ~1e-14
        for sysid in (1, 4), (model, beta) in (("Exp", "0"), ("LFMO-C", "0.01"), ("LFMO-G", "0.5"),
                                               ("LFMO-IG", "0.2"), ("LFMO-S", "0.05"), ("BFR", "0.1"),
                                               ("LSh", "0.2"))
            sys, sf, cf = setup_row(systems, sysid, model, beta)
            for d in (Delay(:none, 0), Delay(:exp, 0.03), Delay(:exp, 1.0), Delay(:det, 0.03), Delay(:det, 1.0),
                      Delay(:insp, 0.03), Delay(:insp, 1.0))
                rs, pf, en, et = rout_quantities(sf.s, cf, d)
                _, pfc, enc, etc = rout_quantities_chain(Float64.(sf.s), cf.R, cf.exits, d)
                @test max_diff(pf, pfc; relative = false) < 1e-13
                @test max_diff(en, enc; relative = true) < 1e-13
                @test max_diff(et, etc; relative = true) < 1e-13
                T = call_tables(cf, d); Tc = call_tables_chain(cf.R, cf.exits, d)
                for (a, b) in ((T.pf, Tc.pf), (T.pN, Tc.pN), (T.enF, Tc.enF))
                    @test max_diff(a, b; relative = false) < 1e-13
                end
                for (a, b) in ((T.enA, Tc.enA), (T.ell, Tc.ell))
                    @test max_diff(a, b; relative = true) < 1e-13
                end
            end
        end
        # the precision is enough: 256 and 512 bits agree far beyond Float64 (n = 30)
        m = lfmo_model("BFR", "0.2")
        cf256 = closed_forms(m, 30); cf512 = closed_forms(m, 30; precision = 512)
        @test cf256.prec == CLOSED_PRECISION == 256 && cf512.prec == 512
        for d in (Delay(:exp, 0.01), Delay(:det, 0.01), Delay(:insp, 0.01))
            rs, P, E, T = kout_quantities(cf256, delay_forms(cf256, d), d.kind)
            _, P2, E2, T2 = kout_quantities(cf512, delay_forms(cf512, d), d.kind)
            @test maximum(abs.(P .- P2)) < 1e-50
            @test maximum(abs.(E .- E2) ./ E2) < 1e-50
            @test maximum(abs.(T .- T2) ./ T2) < 1e-50
        end
        # ... whereas in Float64 the same closed forms are useless for n = 30: the
        # alternating sum for P(N(t) = f) (eq:PNwt) is off by orders of magnitude
        Psi = [Float64(psi(m, x)) for x in 0:30]
        f, w, t = 25, 30, 0.01
        naive = sum(binomial(w, l) * binomial(l, w - f) * (-1)^(l - (w - f)) * exp(-t * Psi[l + 1]) for l in (w - f):w)
        exact = setprecision(BigFloat, 256) do
            Pb = [psi(m, x) for x in 0:30]
            Float64(sum(binomial(w, l) * binomial(l, w - f) * (-1)^(l - (w - f)) * exp(-t * Pb[l + 1]) for l in (w - f):w))
        end
        @test exact > 0
        @test abs(naive - exact) > 100 * exact
    end

    @testset "LTMC quantities: vanishing delays recover the no-delay case" begin
        for (sysid, model, beta) in ((5, "BFR", "0.5"), (13, "LFMO-S", "0.2"), (19, "LSh", "0.1"))
            base = quantities(sysid, model, beta, "none", 0; systems = systems)
            for (kind, d, tol) in (("exp", 1e-9, 1e-6), ("det", 1e-9, 1e-6), ("insp", 1e-8, 1e-5))
                q = quantities(sysid, model, beta, kind, d; systems = systems)
                @test all(isapprox.(q[1], base[1]; atol = tol))
                @test all(isapprox.(q[2], base[2]; rtol = tol))
                @test all(isapprox.(q[3], base[3]; rtol = tol))
            end
            # longer delays: more system failures and more failed components at the repair
            short = quantities(sysid, model, beta, "det", 0.01; systems = systems)
            long = quantities(sysid, model, beta, "det", 0.3; systems = systems)
            @test all(long[1] .>= short[1] .- 1e-12)
            @test all(long[2] .>= short[2] .- 1e-12)
        end
    end

    @testset "LTMC and cost ratio" begin
        pf, en, et = [0.5, 1.0], [1.0, 2.0], [1.0, 2.0]
        # ratio = n c_cmp/c_sys = 1, n = 2, c_sys = 1: c_cmp = 1/2
        @test ltmc(pf, en, et, 2; ratio = 1.0) ≈ [(0.5 * 1 + 0.5) / 1, (0.5 * 2 + 1) / 2]
        # plus a setup cost of 0.1 c_sys per repair visit
        @test ltmc(pf, en, et, 2; ratio = 1.0, setup = 0.1) ≈ [(0.1 + 0.5 * 1 + 0.5) / 1, (0.1 + 0.5 * 2 + 1) / 2]
        # iid exponential, proportional cost, no delay: repairing at every failure costs
        # c_cmp per failure, i.e. LTMC(1) = n c_cmp / (component mean life) = ratio
        q = quantities(5, "Exp", "0", "none", 0; systems = systems)
        @test ltmc(q..., 12; ratio = 4.0)[1] ≈ 4.0
        @test ltmc(q..., 12; ratio = 0.25)[1] ≈ 0.25
        # ratio scales the component term: c_cmp = ratio c_sys / n
        @test ltmc([0.0], [3.0], [1.0], 6; ratio = 2.0) ≈ [3.0 * 2.0 / 6]
    end

    @testset "experiment grid and row output" begin
        grid = read_grid()
        @test length(grid) == 20 * 31 * 16
        @test all(Tuple(a) == Tuple(b) for (a, b) in zip(grid, make_grid()))
        @test [g.row_id for g in grid] == collect(1:length(grid))
        dir = mktempdir()
        cache = Dict()
        row = grid[findfirst(g -> g.system == 1 && g.model == "LFMO-C" && g.beta == "0.2" &&
                                  g.delay_type == "exp" && g.delay == "0.1", grid)]
        @test run_row(row; systems = systems, dir = dir, cache = cache) == :done
        @test run_row(row; systems = systems, dir = dir, cache = cache) == :skipped
        out = read_row_output(row_output_path(row.row_id; dir = dir))
        ref = read_ltmc_reference(joinpath(@__DIR__, "reference", "ltmc_reference.csv"))[(1, "LFMO-C", "0.2", "exp", "0.1")]
        @test out.r == 0:5                                 # r = 0 with an exponential delay
        @test isapprox(out.p_fail, [v[2] for v in ref]; rtol = 1e-10)
        @test out.meta["model"] == "LFMO-C"
        @test out.meta["engine_version"] == string(ENGINE_VERSION)
        @test parse(Float64, out.meta["engine_check_prob"]) < 1e-12
        @test parse(Float64, out.meta["engine_check_rel"]) < 1e-12
        @test length(readdir(dir)) == 1                    # no leftover temporary file
        # inspections: no r = 0
        row = grid[findfirst(g -> g.system == 1 && g.model == "Exp" && g.delay_type == "insp", grid)]
        @test run_row(row; systems = systems, dir = dir, cache = cache) == :done
        @test read_row_output(row_output_path(row.row_id; dir = dir)).r == 1:5
    end

    @testset "optimal policies: after-call quantities" begin
        for (sysid, model, beta) in ((1, "LFMO-C", "0.2"), (5, "BFR", "0.5"), (13, "LFMO-S", "0.05"))
            sys, sf, cf = setup_row(systems, sysid, model, beta)
            R, ex = cf.R, cf.exits
            n = ncomponents(sys); s = Float64.(sf.s)
            for d in (Delay(:exp, 0.3), Delay(:det, 0.1), Delay(:exp, 3.0), Delay(:det, 1.0))
                cm = CountModel(sf.W, cf, d)
                # a call at f = 0 is the r = 0 policy of the LTMC engine (different code path)
                pf0, en0, et0 = policy_quantities_r0(sf.s, cf, d)
                @test isapprox(cm.pf[1], pf0; rtol = 1e-11, atol = 1e-14)
                @test isapprox(cm.en[1], en0; rtol = 1e-11)
                @test isapprox(cm.ell[1], et0; rtol = 1e-11)
                d.kind == :exp || continue
                # exponential delay: the signature form agrees with the recursion eq:ccall exp
                γ = 1 / d.value
                Pf = zeros(n + 1); EN = zeros(n + 1); Ll = zeros(n + 1)
                for f in (n - 1):-1:0
                    cm.Sbar[f + 1] > 0 || continue
                    w = n - f; Ψw = ex[f + 1]
                    a = 0.0; b = 0.0; c = 0.0
                    for j in 1:w
                        Pj = R[f + 1, f + j + 1] / Ψw
                        ps = cm.Sbar[f + j + 1] / cm.Sbar[f + 1]
                        a += Pj * (ps * Pf[f + j + 1] + (1 - ps))
                        b += Pj * (ps * EN[f + j + 1] + (1 - ps) * (f + j))
                        c += Pj * ps * Ll[f + j + 1]
                    end
                    Pf[f + 1] = Ψw * a / (γ + Ψw)
                    EN[f + 1] = (γ * f + Ψw * b) / (γ + Ψw)
                    Ll[f + 1] = (1 + Ψw * c) / (γ + Ψw)
                    @test isapprox(cm.pf[f + 1], Pf[f + 1]; rtol = 1e-10, atol = 1e-14)
                    @test isapprox(f + cm.en[f + 1], EN[f + 1]; rtol = 1e-10)
                    @test isapprox(cm.ell[f + 1], Ll[f + 1]; rtol = 1e-10)
                end
            end
            # inspections: one period either ends at an inspection or with a system failure
            cm = CountModel(sf.W, cf, Delay(:insp, 0.3))
            for f in 0:(n - 1)
                cm.Sbar[f + 1] > 0 && @test isapprox(sum(cm.q[f + 1]) + cm.pf[f + 1], 1; atol = 1e-12)
            end
        end
    end

    @testset "optimal policies: forced thresholds reproduce the r-out-of-n engine" begin
        for (sysid, model, beta) in ((3, "LFMO-G", "0.2"), (5, "LSh", "0.5"), (11, "BFR", "0.2"))
            sys, sf, cf = setup_row(systems, sysid, model, beta)
            n = ncomponents(sys)
            S = StateSpace(sys)
            @test working_counts(S) == sf.W
            for d in (Delay(:none, 0), Delay(:exp, 0.3), Delay(:det, 0.1), Delay(:insp, 0.3))
                rs, pf, en, et = rout_quantities(sf.s, cf, d)
                @test rs == (d.kind in (:exp, :det) ? collect(0:n) : collect(1:n))
                cm = CountModel(sf.W, cf, d)
                pre = state_precompute(S, cm)
                costs = Costs(n; ratio = 1.0, setup = 0.1)
                for (i, r) in enumerate(rs)
                    L = (costs.c_set + costs.c_cmp * en[i] + costs.c_sys * pf[i]) / et[i]
                    _, c, l = count_pass(cm, costs, 0.0; forced = f -> f >= r)
                    @test isapprox(c / l, L; rtol = 1e-11)
                    _, c, l = state_pass(S, cm, pre, costs, 0.0; forced = (x, f) -> f >= r)
                    @test isapprox(c / l, L; rtol = 1e-11)
                end
            end
        end
    end

    @testset "optimal policies vs. independent reference (test/reference/optimal_reference.csv)" begin
        header, rows = SimpleDelayed.read_simple_csv(joinpath(@__DIR__, "reference", "optimal_reference.csv"))
        col(name) = findfirst(==(name), header)
        @test length(rows) == 64
        cache = Dict()
        for r in rows
            sysid = parse(Int, r[col("system")])
            sys, sf, cf = setup_row(systems, sysid, r[col("model")], r[col("beta")])
            n = ncomponents(sys)
            d = Delay(r[col("delay_type")], parse(Float64, r[col("delay")]))
            costs = Costs(n; ratio = parse(Float64, r[col("ratio")]), setup = parse(Float64, r[col("setup")]))
            rs, pf, en, et = rout_quantities(sf.s, cf, d)
            L = (costs.c_set .+ costs.c_cmp .* en .+ costs.c_sys .* pf) ./ et
            ib = argmin(L)
            gr = parse(Float64, r[col("ltmc_best_r")])
            @test isapprox(L[ib], gr; rtol = 1e-11)
            @test isapprox(L[findfirst(==(parse(Int, r[col("best_r")])), rs)], gr; rtol = 1e-11)
            cm = CountModel(sf.W, cf, d)
            gc, act, _ = dinkelbach(g -> count_pass(cm, costs, g), L[ib])
            @test isapprox(gc, parse(Float64, r[col("ltmc_count")]); rtol = 1e-11)
            refpol = r[col("count_policy")]
            if policy_string(act) != refpol           # only a numerical tie could explain this
                _, c, l = count_pass(cm, costs, 0.0; forced = f -> refpol[f + 1] == 'C')
                @test isapprox(c / l, gc; rtol = 1e-11)
            end
            S = get!(cache, sysid) do
                StateSpace(sys)
            end
            pre = state_precompute(S, cm)
            gs, _, _ = dinkelbach(g -> state_pass(S, cm, pre, costs, g), gc)
            @test isapprox(gs, parse(Float64, r[col("ltmc_state")]); rtol = 1e-11)
        end
    end

    @testset "optimal policies: [Exp] with continuous inspection is a threshold policy" begin
        for sysid in (2, 6, 9)
            sys, sf, cf = setup_row(systems, sysid, "Exp", "0")
            n = ncomponents(sys)
            for d in (Delay(:none, 0), Delay(:exp, 0.01), Delay(:exp, 1.0), Delay(:det, 0.01), Delay(:det, 1.0))
                cm = CountModel(sf.W, cf, d)
                rs, pf, en, et = rout_quantities(sf.s, cf, d)
                for ρ in OPT_RATIOS, σ in OPT_SETUPS
                costs = Costs(n; ratio = ρ, setup = σ)
                L = (costs.c_set .+ costs.c_cmp .* en .+ costs.c_sys .* pf) ./ et
                gc, act, _ = dinkelbach(g -> count_pass(cm, costs, g), minimum(L))
                @test is_threshold(act, d.kind)
                @test isapprox(gc, minimum(L); rtol = 1e-11)    # hence no better than the best r
                end
            end
        end
        @test !is_threshold(Int8[0, 0, 1, 0, -1], :exp)
        @test is_threshold(Int8[0, 1, -1, 1, -1], :exp)
        @test is_threshold(Int8[1, 0, 1, 1], :insp)               # f = 0 is not a decision at inspections
    end

    @testset "optimal policies: grid rows" begin
        grid = read_grid()
        dir = mktempdir()
        cache = Dict()
        for (sysid, model, kind) in ((1, "BFR", "det"), (3, "LFMO-S", "insp"), (13, "LFMO-C", "exp"))
            row = grid[findfirst(g -> g.system == sysid && g.model == model && g.beta == "0.2" &&
                                      g.delay_type == kind && g.delay == "0.1", grid)]
            @test run_optimal_row(row; systems = systems, dir = dir, cache = cache) == :done
            @test run_optimal_row(row; systems = systems, dir = dir, cache = cache) == :skipped
            res = read_optimal_output(row_output_path(row.row_id; dir = dir))
            @test res.header == OPTIMAL_HEADER
            @test res.meta["optimal_version"] == string(OPTIMAL_VERSION) == "3"
            @test length(res.rows) == length(OPT_RATIOS) * length(OPT_SETUPS)
            @test res.meta["state_based"] == string(sysid != 13)
            @test parse(Float64, res.meta["engine_check_prob"]) < 1e-12
            @test parse(Float64, res.meta["engine_check_rel"]) < 1e-12
            for k in ("seconds", "seconds_rout", "seconds_tables", "seconds_states", "seconds_state_pre",
                      "seconds_chain", "maxrss_mb")
                @test parse(Float64, res.meta[k]) >= 0
            end
            @test parse(Float64, res.meta["seconds_rout"]) > 0
            @test res.meta["states_cached"] in ("true", "false")
            @test parse(Int, res.meta["n_working_states"]) == sum(read_signature_file(signature_path(get_system(systems, sysid))).W)
            c(name) = findfirst(==(name), res.header)
            # the row holds exactly the closed-form values (the chain is only a check)
            sys, sf, cf = setup_row(systems, sysid, model, "0.2")
            n = ncomponents(sys)
            d = Delay(kind, 0.1)
            rs, pf, en, et = rout_quantities(sf.s, cf, d)
            cm = CountModel(sf.W, cf, d)
            S = sysid == 13 ? nothing : StateSpace(sys)
            pre = S === nothing ? nothing : state_precompute(S, cm)
            for r in res.rows
                costs = Costs(n; ratio = parse(Float64, r[c("ratio")]), setup = parse(Float64, r[c("setup")]))
                L = (costs.c_set .+ costs.c_cmp .* en .+ costs.c_sys .* pf) ./ et
                @test r[c("ltmc_best_r")] == repr(minimum(L))
                gc, _, itc = dinkelbach(g -> count_pass(cm, costs, g), minimum(L))
                @test r[c("ltmc_count")] == repr(gc) && r[c("iters_count")] == string(itc)
                @test parse(Float64, r[c("seconds_count")]) >= 0
                if S !== nothing
                    gs, _, _ = dinkelbach(g -> state_pass(S, cm, pre, costs, g), gc)
                    @test r[c("ltmc_state")] == repr(gs)
                end
            end
            for r in res.rows
                gr, gc, gs = parse.(Float64, (r[c("ltmc_best_r")], r[c("ltmc_count")], r[c("ltmc_state")]))
                @test gc <= gr * (1 + 1e-12)
                sysid == 13 ? (@test isnan(gs)) : (@test gs <= gc * (1 + 1e-12))
                @test length(r[c("count_policy")]) == ncomponents(get_system(systems, sysid))
            end
        end
        @test length(readdir(dir)) == 3
        # readers refuse other versions; compare_optimal.jl may read them explicitly
        lines = String.(split(read(joinpath(dir, first(readdir(dir))), String), '\n'))
        for (i, l) in enumerate(lines)
            l == "# optimal_version = 3" && (lines[i] = "# optimal_version = 2")
        end
        hdr = findfirst(l -> startswith(l, "ratio,"), lines)
        for i in hdr:length(lines)                    # drop the 8th field, seconds_count (version 3)
            isempty(lines[i]) && continue
            f = split(lines[i], ','); deleteat!(f, 8); lines[i] = join(f, ',')
        end
        p = joinpath(mktempdir(), "row_00001.csv")
        write(p, join(lines, '\n'))
        @test_throws ErrorException read_optimal_output(p)
        @test read_optimal_output(p; versions = (2, 3)).header == filter(!=("seconds_count"), OPTIMAL_HEADER)
        # warm-up leaves nothing behind and returns a time
        @test optimal_warmup(systems; state_max_n = 14) > 0
    end
end
