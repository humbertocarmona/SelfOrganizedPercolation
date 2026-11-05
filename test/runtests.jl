using SelfOrganizedPercolation
using SelfOrganizedPercolation: State, resolve!, EMPTY, heap_push!, heap_pop!, Threshold
using Random, Test

@testset "lattice geometry" begin
    for (lat, degree) in ((SquareLattice(7), 4), (TriangularLattice(7), 6), (TriangularLattice(8), 6))
        W, nx = lat.width, 9
        for i in 1:W*nx
            nbrs = filter(!=(0), collect(neighbors(lat, i, nx, true)))
            @test allunique(nbrs)
            @test i ∉ nbrs
            for j in nbrs
                @test i in neighbors(lat, j, nx, true)   # symmetric
            end
            x = (i - 1) ÷ W
            0 < x < nx - 1 && @test length(nbrs) == degree
        end
        for i in 1:W*nx, j in neighbors(lat, i, nx, false)
            j == 0 && continue
            @test i in neighbors(lat, j, nx, false)
        end
    end
    @test neighbors(SquareLattice(5), 1, 3, false) == (6, 0, 2, 0)
    @test neighbors(SquareLattice(5), 1, 3, true) == (6, 5, 2, 0)
end

@testset "threshold heap" begin
    h = Threshold[]
    xs = rand(Xoshiro(1), 1000)
    foreach(x -> heap_push!(h, (x, 1)), xs)
    @test [heap_pop!(h)[1] for _ in 1:1000] == sort(xs)
end

@testset "feedback and bookkeeping" begin
    m = SOPModel(SquareLattice(64); percolation=SitePercolation(), seeds=20, c=0.064, f_T=10 / 64, p0=0.6, tmax=200)
    r = simulate(m; rng=Xoshiro(3))
    for t in 1:size(r.p, 1)-1
        @test r.p[t+1] ≈ clamp(r.p[t] + 0.064 * (10 / 64 - r.f[t]), 0, 1)
    end
    @test all(0 .<= r.p .<= 1)
    @test r.mass[1] == 20 + round(Int, sum(r.f) * 64)
end

@testset "time-dependent target" begin
    m = SOPModel(SquareLattice(64); seeds=20, c=0.064, f_T=PowerLawTarget(10 / 64, 0.5), p0=0.7, tmax=50)
    r = simulate(m; rng=Xoshiro(4))
    @test r.p[6] ≈ clamp(r.p[5] + 0.064 * (10 / 64 * sqrt(4) - r.f[5]), 0, 1)
end

@testset "reproducibility" begin
    m = SOPModel(SquareLattice(128); seeds=20, c=0.128, f_T=30 / 128, tmax=300)
    a = ensemble(m, 8; seed=9, ntasks=1)
    b = ensemble(m, 8; seed=9, ntasks=4)
    @test all(x.p == y.p && x.f == y.f for (x, y) in zip(a, b))
end

@testset "limits of static percolation" begin
    for perc in (SitePercolation(), BondPercolation())
        full = SOPModel(SquareLattice(32); percolation=perc, seeds=32, c=0.0, p0=1.0)
        r = simulate(full; rng=Xoshiro(1), keep_sites=true)
        @test r.status === :spanned
        @test all(!=(EMPTY), r.sites)
        none = SOPModel(SquareLattice(32); percolation=perc, seeds=32, c=0.0, p0=0.0)
        @test simulate(none; rng=Xoshiro(1)).status === :extinct
    end
end

@testset "site is tested once" begin
    # A site next to two active seeds gets one chance per species, not two.
    m = SOPModel(SquareLattice(3); percolation=SitePercolation(), seeds=2, c=0.0, p0=0.5, tmax=1)
    hits = trials = 0
    for seed in 1:20000
        rng = Xoshiro(seed)
        st = SelfOrganizedPercolation.initial_state(m, rng)
        wall = findall(!=(EMPTY), st.sites[1:3])
        wall == [1, 3] || continue
        SelfOrganizedPercolation.propose!(st, m, 1, rng)
        trials += 1
        hits += any(c -> c[1] == 2, st.candidates)
    end
    @test abs(hits / trials - 0.5) < 0.03
end

@testset "uniform choice among candidate bonds" begin
    wins = 0
    rng = Xoshiro(5)
    n = 30000
    st = State(UInt8[0], UInt8[], 1, [Int[], Int[]], [Int[], Int[]], Tuple{Int,Int}[],
        [Threshold[], Threshold[]], [1.0, 1.0], [0, 0], [false, false], 0)
    for _ in 1:n
        st.sites[1] = EMPTY
        append!(st.candidates, [(1, 1), (1, 1), (1, 2)])
        resolve!(st, 1, rng)
        wins += st.sites[1] == 0x01
        foreach(empty!, st.next)
    end
    @test abs(wins / n - 2 / 3) < 0.015
end

# Every occupied site must connect to a seed through sites of its own species.
function connected_to_wall(sites, lat, periodic, s)
    W, nx = size(sites)
    seen = falses(length(sites))
    stack = [i for i in 1:W if sites[i] == s]
    seen[stack] .= true
    while !isempty(stack)
        for j in neighbors(lat, pop!(stack), nx, periodic)
            j != 0 && sites[j] == s && !seen[j] && (seen[j] = true; push!(stack, j))
        end
    end
    return all(i -> sites[i] != s || seen[i], eachindex(sites))
end

@testset "immiscible fluids" begin
    for lat in (SquareLattice(96), TriangularLattice(96)),
        perc in (SitePercolation(), BondPercolation()),
        rej in (Permanent(), Persistent()), ns in (2, 3)

        m = SOPModel(lat; percolation=perc, rejection=rej, species=ns, seeds=8, c=0.192,
            f_T=15 / 96, p0=0.8, periodic=true, right=GrowingRight(64), tmax=400)
        r = simulate(m; rng=Xoshiro(11), keep_sites=true)
        @test size(r.p, 2) == ns
        for s in 1:ns
            @test count(==(UInt8(s)), r.sites) == r.mass[s]
            @test connected_to_wall(r.sites, lat, true, UInt8(s))
        end
    end
    # Identical fluids share the medium evenly on average.
    m = SOPModel(SquareLattice(128); species=2, seeds=30, c=0.128, f_T=30 / 128, p0=0.7, tmax=400)
    masses = ensemble(m, 200; f=r -> r.mass)
    @test abs(sum(x[1] for x in masses) / sum(x[2] for x in masses) - 1) < 0.03
end

@testset "competing seeds with a smaller active fraction" begin
    m = SOPModel(SquareLattice(100); species=2, seeds=40, active_fraction=0.5, c=0.0, p0=0.0, tmax=1)
    st = SelfOrganizedPercolation.initial_state(m, Xoshiro(1))
    @test length.(st.fronts) == [20, 20]
    @test count(==(0x01), st.sites) == 40
end

@testset "right boundary" begin
    kw = (; percolation=SitePercolation(), seeds=20, c=0.064, f_T=60 / 64, p0=0.8)
    fixed = SOPModel(SquareLattice(64); right=FixedRight(40), kw...)
    r = simulate(fixed; rng=Xoshiro(2))
    @test r.status === :spanned
    @test r.length == 40

    # Moving the wall must not change the dynamics, and the front never touches it.
    tmax = size(r.p, 1) - 1
    grown = SOPModel(SquareLattice(64); right=GrowingRight(4), tmax=tmax, kw...)
    g = simulate(grown; rng=Xoshiro(2), keep_sites=true)
    @test g.p ≈ r.p[1:tmax, :]
    @test g.status === :tmax
    @test size(g.sites, 2) >= 40

    capped = SOPModel(SquareLattice(64); right=GrowingRight(4; maxlength=16), kw...)
    @test simulate(capped; rng=Xoshiro(2)).status === :maxlength
end

@testset "persistent front survives where the permanent one dies" begin
    kw = (; seeds=5, c=0.256, f_T=2 / 128, p0=0.4, tmax=4000, right=GrowingRight(64))
    alive(rej) = count(r -> r.status !== :extinct,
        ensemble(SOPModel(SquareLattice(128); rejection=rej, kw...), 40))
    @test alive(Persistent()) == 40
    @test alive(Permanent()) < 40
end

@testset "no type instability in the kernel" begin
    for lat in (SquareLattice(32), TriangularLattice(32)), rej in (Permanent(), Persistent()),
        perc in (SitePercolation(), BondPercolation())

        m = SOPModel(lat; percolation=perc, rejection=rej, seeds=4, f_T=0.25, p0=0.7, tmax=20)
        @test (@inferred simulate(m; rng=Xoshiro(1))) isa SOPResult
    end
end
