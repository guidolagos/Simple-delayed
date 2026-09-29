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

function read_reference(path)
    ref = Dict{Int,Vector{Int}}()
    for line in eachline(path)
        (startswith(line, "#") || startswith(line, "system") || isempty(strip(line))) && continue
        f = split(line, ',')
        push!(get!(ref, parse(Int, f[1]), Int[]), parse(Int, f[3]))
    end
    return ref
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
        @test length(systems) == 22
        @test [ncomponents(s) for s in systems] ==
              [5, 8, 8, 9, 12, 14, 21, 13, 12, 12, 12, 13, 30, 26, 14, 21, 27, 22, 30, 24, 30, 23]
        @test [s.id for s in systems if s.directed] == [6, 7]
        raw = parse_yeh_raw(joinpath(SimpleDelayed.SYSTEMS_DIR, "yeh2024_raw.txt"))
        @test length(raw) == length(systems)
        for (r, s) in zip(raw, systems)
            @test (r.id, r.label, r.directed, r.edges) == (s.id, s.label, s.directed, s.edges)
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
        @test sort(collect(keys(ref))) == collect(1:22)
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
end
