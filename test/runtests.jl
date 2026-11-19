using SelfOrganizedPercolation
using SelfOrganizedPercolation: State, resolve!, EMPTY, heap_push!, heap_pop!, Threshold, threshold
using Random, Test

@testset "lattice geometry" begin
    for (lat, degree) in ((SquareLattice(7), 4), (TriangularLattice(7), 6), (TriangularLattice(8), 6),
        (CubicLattice(4, 5), 6), (CubicLattice(3), 6))
        W, nx = SelfOrganizedPercolation.transverse_sites(lat), 9
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
    # Cubic: y fastest, then z, then x. Site (x, y, z) = (1, 2, 3) on a 4 × 5 wall.
    c = CubicLattice(4, 5)
    i = 20 * 1 + 4 * 3 + 2 + 1
    @test neighbors(c, i, 3, false) == (i + 20, i - 1, i + 1, i - 4, i + 4, i - 20)
    @test neighbors(c, 1, 3, true) == (21, 4, 2, 17, 5, 0)
    @test neighbors(c, 1, 3, false) == (21, 0, 2, 0, 5, 0)
    @test_throws ArgumentError SOPModel(CubicLattice(2, 8); seeds=4)
    @test_throws ArgumentError IPModel(CubicLattice(8); seeds=8, layout=BlockSeeds(2))
end

@testset "threshold heap" begin
    h = Threshold[]
    xs = rand(Xoshiro(1), 1000)
    foreach(x -> heap_push!(h, (x, 1)), xs)
    @test [heap_pop!(h)[1] for _ in 1:1000] == sort(xs)
end

@testset "medium thresholds" begin
    m = Medium(7)
    xs = [threshold(m, j) for j in 1:200_000]
    @test all(0 .<= xs .< 1)
    @test xs == [threshold(Medium(7), j) for j in 1:200_000]        # deterministic
    @test abs(sum(xs) / length(xs) - 1 / 2) < 0.003
    @test abs(sum(abs2, xs .- 1 / 2) / length(xs) - 1 / 12) < 0.001
    counts = zeros(Int, 20)
    foreach(x -> counts[floor(Int, 20x)+1] += 1, xs)
    @test all(c -> abs(c / 10_000 - 1) < 0.05, counts)              # flat histogram
    @test abs(sum((xs[1:end-1] .- 0.5) .* (xs[2:end] .- 0.5)) / length(xs)) < 0.002   # neighbors uncorrelated
    @test count(j -> threshold(Medium(8), j) == threshold(m, j), 1:1000) == 0
    @test all(threshold(m, i, i + 3) == threshold(m, i + 3, i) for i in 1:1000)   # bonds are symmetric
    bonds = [threshold(m, i, i + 1) for i in 1:200_000]
    @test abs(sum(bonds) / length(bonds) - 1 / 2) < 0.003
    @test abs(sum((bonds .- 0.5) .* (xs .- 0.5)) / length(xs)) < 0.002         # unrelated to sites
end

# Sites joined to the seeds by elements with threshold below p: what Persistent must
# fill at constant p, whatever the order of invasion.
function static_cluster(lat, W, nx, seeds, medium, perc, p)
    filled = falses(W * nx)
    filled[seeds] .= true
    stack = copy(seeds)
    while !isempty(stack)
        i = pop!(stack)
        for j in neighbors(lat, i, nx, false)
            (j == 0 || filled[j]) && continue
            x = perc isa SitePercolation ? threshold(medium, j) : threshold(medium, i, j)
            x < p && (filled[j] = true; push!(stack, j))
        end
    end
    return reshape(filled, W, nx)
end

@testset "persistent rule reads the medium" begin
    for lat in (SquareLattice(64), TriangularLattice(64)),
        (perc, p) in ((SitePercolation(), 0.47), (BondPercolation(), 0.32))   # below p_c

        m = SOPModel(lat; percolation=perc, rejection=Persistent(), right=FixedRight(48),
            seeds=16, layout=BlockSeeds(4), c=0.0, p0=p, tmax=3000)
        medium = Medium(42)
        a = simulate(m; rng=Xoshiro(1), medium, keep_sites=true)
        b = simulate(m; rng=Xoshiro(2), medium, keep_sites=true)
        @test a.status === :tmax                       # finite cluster, never spanned
        @test a.sites == b.sites                       # the medium alone sets the cluster
        seeds = findall(!=(EMPTY), a.sites[:, 1])
        @test (a.sites .!= EMPTY) == static_cluster(lat, 64, 48, seeds, medium, perc, p)
        c = simulate(m; rng=Xoshiro(1), medium=Medium(43), keep_sites=true)
        @test c.sites != a.sites
    end
    # By default the seed fixes the medium, so models with the same seed share it.
    m = SOPModel(SquareLattice(32); rejection=Persistent(), seeds=8, p0=0.7, tmax=50)
    @test simulate(m; rng=Xoshiro(5)).p == simulate(m; rng=Xoshiro(5)).p
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

@testset "ratchet feedback" begin
    # Eq. (1) of the SOP-IP note: only the positive part of c[f_T - f] is added.
    for perc in (SitePercolation(), BondPercolation()), rej in (Permanent(), Persistent())
        m = SOPModel(SquareLattice(64); percolation=perc, rejection=rej, seeds=20, c=0.064,
            f_T=10 / 64, p0=0.6, tmax=200, feedback=Ratchet())
        r = simulate(m; rng=Xoshiro(3))
        for t in 1:size(r.p, 1)-1
            @test r.p[t+1] ≈ min(r.p[t] + 0.064 * max(10 / 64 - r.f[t], 0), 1)
        end
        @test all(diff(r.p[:, 1]) .>= 0)                    # p never decreases
        @test any(r.f .> 10 / 64)                           # ...although the front overshoots
        @test all(0 .<= r.p .<= 1)
    end
    # The default is the two-sided rule, and it differs from the ratchet.
    kw = (; seeds=20, c=0.064, f_T=10 / 64, p0=0.6, tmax=200)
    @test SOPModel(SquareLattice(64); kw...).feedback === TwoSided()
    two = simulate(SOPModel(SquareLattice(64); kw...); rng=Xoshiro(3))
    @test two.p == simulate(SOPModel(SquareLattice(64); kw..., feedback=TwoSided()); rng=Xoshiro(3)).p
    @test any(diff(two.p[:, 1]) .< 0)
    # Several species: each has its own ratchet.
    m = SOPModel(SquareLattice(64); species=2, seeds=10, c=0.1, f_T=5 / 64, p0=[0.5, 0.7],
        tmax=100, feedback=Ratchet())
    r = simulate(m; rng=Xoshiro(6))
    @test all(diff(r.p; dims=1) .>= 0)
    # Saturates at 1 instead of overshooting.
    r = simulate(SOPModel(SquareLattice(32); seeds=4, c=10.0, f_T=0.5, p0=0.95, tmax=5,
        feedback=Ratchet()); rng=Xoshiro(1))
    @test maximum(r.p) == 1.0
end

@testset "ratchet extinction clusters are invasion percolation ponds" begin
    # With Persistent rejection, p0 = 0 and f_T < 1/A, the SOP cluster at every
    # extinction is the IP cluster of the same size in the same medium (Proposition 2
    # of the SOP-IP note).
    for lat in (SquareLattice(32), TriangularLattice(32)), perc in (SitePercolation(), BondPercolation())
        W = SelfOrganizedPercolation.transverse_sites(lat)
        medium = Medium(11)
        sop = SOPModel(lat; percolation=perc, rejection=Persistent(), right=FixedRight(40),
            seeds=8, p0=0.0, c=1.0, f_T=0.5 / W, tmax=100_000, feedback=Ratchet())
        ip = IPModel(lat; percolation=perc, right=FixedRight(40), seeds=8)
        a = simulate(sop; rng=Xoshiro(2), medium, keep_times=true)
        b = simulate(ip; rng=Xoshiro(2), medium, keep_times=true)
        @test a.status === :spanned
        nseeds = count(==(0), b.times)
        extinctions = 0
        for k in 2:size(a.f, 1)
            a.f[k, 1] == 0 && a.f[k-1, 1] > 0 || continue   # the front has just died
            cluster = 0 .<= a.times .<= k
            n = count(cluster) - nseeds
            @test cluster == (0 .<= b.times .<= n)
            @test b.x[n+1] >= a.p[k, 1]                     # the next IP threshold is above p
            extinctions += 1
        end
        @test extinctions > 3
        # The two-sided rule does not give the same guarantee: p falls inside a pond.
        two = SOPModel(lat; percolation=perc, rejection=Persistent(), right=FixedRight(40),
            seeds=8, p0=0.0, c=1.0, f_T=0.5 / W, tmax=100_000)
        @test any(diff(simulate(two; rng=Xoshiro(2), medium).p[:, 1]) .< 0)
    end
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

@testset "ensemble with mixed return types" begin
    m = SOPModel(SquareLattice(32); seeds=4, f_T=0.25, p0=0.7, tmax=20)
    r = ensemble(m, 6; f=x -> isodd(length(x.p)) ? nothing : 1, ntasks=1)
    @test length(r) == 6
    mixed = ensemble(m, 4; f=x -> x.mass[1] > 0 ? (x.mass[1] % 2 == 0 ? nothing : 1) : nothing)
    @test length(mixed) == 4
    @test ensemble(m, 3; f=x -> 1.0) isa Vector{Float64}
end

@testset "seed layouts" begin
    m = SOPModel(SquareLattice(64); seeds=32, layout=BlockSeeds(4), c=0.0, p0=0.0, tmax=1)
    st = SelfOrganizedPercolation.initial_state(m, Xoshiro(1))
    wall = findall(!=(EMPTY), st.sites[1:64])
    @test wall == vcat(5:12, 21:28, 37:44, 53:60)          # equal blocks and gaps
    @test wall[1] > 1 && wall[end] < 64                    # away from the periodic edge
    @test_throws ArgumentError SOPModel(SquareLattice(64); seeds=30, layout=BlockSeeds(4))
    @test_throws ArgumentError SOPModel(SquareLattice(64); species=2, seeds=8, layout=BlockSeeds(4))
    # The default layout gives the same run as before the option existed.
    a = simulate(SOPModel(SquareLattice(64); seeds=20, p0=0.7, tmax=50); rng=Xoshiro(1))
    b = simulate(SOPModel(SquareLattice(64); seeds=20, p0=0.7, tmax=50, layout=RandomSeeds()); rng=Xoshiro(1))
    @test a.p == b.p
end

@testset "arrival times" begin
    for rej in (Permanent(), Persistent()), right in (FixedRight(40), GrowingRight(8))
        m = SOPModel(SquareLattice(64); rejection=rej, right, seeds=20, c=0.1, f_T=0.3, p0=0.8, tmax=300)
        plain = simulate(m; rng=Xoshiro(3), keep_sites=true)
        r = simulate(m; rng=Xoshiro(3), keep_sites=true, keep_times=true)
        @test r.p == plain.p && r.sites == plain.sites          # recording changes nothing
        @test size(r.times) == size(r.sites)
        @test all((r.times .>= 0) .== (r.sites .!= 0))
        @test count(==(0), r.times) == 20
        for t in 1:size(r.f, 1)
            @test count(==(t), r.times) == round(Int, r.f[t, 1] * 64)
        end
        # Every filled site has a neighbor filled exactly one step earlier, unless it
        # was reached through a throat that waited, which only memory allows.
        if rej isa Permanent
            W, nx = size(r.times)
            for i in eachindex(r.times)
                r.times[i] > 0 || continue
                @test any(j -> j != 0 && r.times[j] == r.times[i] - 1, neighbors(m.lattice, i, nx, false))
            end
        end
    end
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
        [Threshold[], Threshold[]], [1.0, 1.0], [0, 0], [false, false], 0, Int32[], Int32(0), Medium(1),
        SelfOrganizedPercolation.Trapper(0))
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
    nx = size(sites)[end]
    W = length(sites) ÷ nx
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
    for lat in (SquareLattice(96), TriangularLattice(96), CubicLattice(12, 8)),
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
    for lat in (SquareLattice(32), TriangularLattice(32), CubicLattice(6)), rej in (Permanent(), Persistent()),
        perc in (SitePercolation(), BondPercolation()), fb in (TwoSided(), Ratchet())

        m = SOPModel(lat; percolation=perc, rejection=rej, seeds=4, f_T=0.25, p0=0.7, tmax=20, feedback=fb)
        @test (@inferred simulate(m; rng=Xoshiro(1))) isa SOPResult
    end
end

# Invasion percolation by the simple O(MN) algorithm: scan the whole perimeter for
# the smallest threshold at every step. Returns the sites in invasion order.
function brute_invasion(lat, W, nx, periodic, seeds, medium, perc, nmax)
    filled = falses(W * nx)
    filled[seeds] .= true
    order = Int[]
    while length(order) < nmax
        best, site = Inf, 0
        for i in findall(filled), j in neighbors(lat, i, nx, periodic)
            (j == 0 || filled[j]) && continue
            x = perc isa SitePercolation ? threshold(medium, j) : threshold(medium, i, j)
            x < best && ((best, site) = (x, j))
        end
        site == 0 && break
        filled[site] = true
        push!(order, site)
        (site - 1) ÷ W == nx - 1 && break
    end
    return order
end

@testset "invasion percolation matches the simple algorithm" begin
    for lat in (SquareLattice(16), TriangularLattice(16), CubicLattice(5, 4)),
        perc in (SitePercolation(), BondPercolation()), periodic in (false, true)

        m = IPModel(lat; percolation=perc, periodic, seeds=5, right=FixedRight(24))
        medium = Medium(3)
        r = simulate(m; rng=Xoshiro(1), medium, keep_sites=true, keep_times=true)
        @test r.status === :spanned
        seeds = findall(==(0), vec(r.times))
        order = sortperm(vec(r.times))[length(seeds)+count(<(0), r.times)+1:end]
        W = SelfOrganizedPercolation.transverse_sites(lat)
        @test order == brute_invasion(lat, W, 24, periodic, seeds, medium, perc, typemax(Int))
        @test r.mass[1] == 5 + length(r.x) == count(!=(EMPTY), r.sites)
        @test r.length == 24
    end
end

@testset "invasion percolation shares the medium with SOP" begin
    # IP fills every element below p reachable from the seeds before any element above
    # it, so up to its first threshold above p it has filled exactly the cluster that
    # Persistent SOP fills at constant p in the same medium.
    for lat in (SquareLattice(64), TriangularLattice(64)),
        (perc, p) in ((SitePercolation(), 0.47), (BondPercolation(), 0.32))

        medium = Medium(42)
        sop = SOPModel(lat; percolation=perc, rejection=Persistent(), right=FixedRight(48),
            seeds=16, layout=BlockSeeds(4), c=0.0, p0=p, tmax=3000)
        ip = IPModel(lat; percolation=perc, right=FixedRight(48), seeds=16, layout=BlockSeeds(4))
        a = simulate(sop; rng=Xoshiro(1), medium, keep_sites=true)
        b = simulate(ip; rng=Xoshiro(1), medium, keep_times=true)
        nstar = findfirst(>=(p), b.x) - 1
        @test (a.sites .!= EMPTY) == (0 .<= b.times .<= nstar)
    end
    # With the same seed, both models put the seeds in the same places.
    sop = SOPModel(SquareLattice(64); rejection=Persistent(), seeds=20, p0=0.0, tmax=1)
    ip = IPModel(SquareLattice(64); seeds=20, nmax=5)
    a = simulate(sop; rng=Xoshiro(8), keep_sites=true)
    b = simulate(ip; rng=Xoshiro(8), keep_times=true)
    @test findall(==(0), b.times[:, 1]) == findall(!=(EMPTY), a.sites[:, 1])
end

@testset "invasion percolation with several fluids" begin
    for lat in (SquareLattice(64), TriangularLattice(64), CubicLattice(10, 8)), perc in (SitePercolation(), BondPercolation())
        m = IPModel(lat; percolation=perc, species=3, seeds=8, periodic=true, right=GrowingRight(16), nmax=6000)
        r = simulate(m; rng=Xoshiro(11), keep_sites=true, keep_times=true)
        @test r.status === :nmax
        @test length(r.x) == 6000
        for s in 1:3
            @test count(==(UInt8(s)), r.sites) == r.mass[s] == 8 + count(==(UInt8(s)), r.species)
            @test connected_to_wall(r.sites, lat, true, UInt8(s))
        end
        @test sort(filter(>(0), vec(r.times))) == 1:6000
    end
    # Fluids are picked uniformly, so identical fluids share the medium evenly.
    m = IPModel(SquareLattice(128); species=2, seeds=30, periodic=true, right=GrowingRight(32), nmax=4000)
    masses = ensemble(m, 100; f=r -> r.mass)
    @test abs(sum(x[1] for x in masses) / sum(x[2] for x in masses) - 1) < 0.05
    # The medium is fixed while rng changes which fluid moves.
    medium = Medium(1)
    a = simulate(m; rng=Xoshiro(1), medium, keep_sites=true)
    b = simulate(m; rng=Xoshiro(2), medium, keep_sites=true)
    @test a.sites != b.sites
end

@testset "invasion percolation boundaries" begin
    m = IPModel(SquareLattice(32); percolation=SitePercolation(), seeds=32, right=GrowingRight(4), nmax=20000)
    r = simulate(m; rng=Xoshiro(1), keep_sites=true)
    @test r.status === :nmax
    @test size(r.sites, 2) > r.length                      # the wall stays ahead
    capped = IPModel(SquareLattice(32); seeds=32, right=GrowingRight(4; maxlength=16))
    c = simulate(capped; rng=Xoshiro(1), keep_sites=true)
    @test c.status === :maxlength
    @test c.length <= 15 && size(c.sites, 2) == 16
    none = IPModel(SquareLattice(32); seeds=4, active_fraction=0.0)
    @test simulate(none; rng=Xoshiro(1)).status === :extinct
    quiet = simulate(IPModel(SquareLattice(32); seeds=4); rng=Xoshiro(1), keep_history=false)
    @test isempty(quiet.x) && quiet.status === :spanned
    for lat in (SquareLattice(32), TriangularLattice(32), CubicLattice(6)), perc in (SitePercolation(), BondPercolation())
        m = IPModel(lat; percolation=perc, species=2, seeds=4, nmax=200)
        @test (@inferred simulate(m; rng=Xoshiro(1))) isa IPResult
    end
end

@testset "cubic lattice" begin
    # Results keep the wall's two dimensions, and the default length is ny.
    m = SOPModel(CubicLattice(8, 6); seeds=10, p0=0.6, tmax=20)
    r = simulate(m; rng=Xoshiro(1), keep_sites=true, keep_times=true)
    @test ndims(r.sites) == 3 && size(r.sites)[1:2] == (8, 6)
    @test size(r.times) == size(r.sites)
    @test m.right.length == 8
    @test count(!=(EMPTY), r.sites) == r.mass[1]

    # IP in the same medium as Persistent SOP, up to its first threshold above p.
    for (perc, p) in ((SitePercolation(), 0.28), (BondPercolation(), 0.22))   # below p_c
        medium = Medium(9)
        sop = SOPModel(CubicLattice(16); percolation=perc, rejection=Persistent(), right=FixedRight(24),
            seeds=40, c=0.0, p0=p, tmax=2000)
        ip = IPModel(CubicLattice(16); percolation=perc, right=FixedRight(24), seeds=40)
        a = simulate(sop; rng=Xoshiro(2), medium, keep_sites=true)
        b = simulate(ip; rng=Xoshiro(2), medium, keep_times=true)
        @test a.status === :tmax
        nstar = findfirst(>=(p), b.x) - 1
        @test (a.sites .!= EMPTY) == (0 .<= b.times .<= nstar)
    end

    # Static percolation thresholds: the largest invaded threshold approaches p_c.
    for (perc, pc) in ((SitePercolation(), 0.3116), (BondPercolation(), 0.2488))
        m = IPModel(CubicLattice(48); percolation=perc, periodic=true, seeds=48^2)
        r = simulate(m; rng=Xoshiro(3))
        @test r.status === :spanned
        @test abs(maximum(r.x) - pc) < 0.02
    end
end

# Menisci on the front after a run, counted from scratch: perimeter sites for site
# percolation, (contributor, empty neighbor) bonds for bond percolation.
function brute_menisci(sites, lat, nx, periodic, contributors, perc)
    empties = Int[]
    for i in contributors, j in neighbors(lat, i, nx, periodic)
        j != 0 && sites[j] == EMPTY && push!(empties, j)
    end
    return perc isa SitePercolation ? length(unique(empties)) : length(empties)
end

@testset "capacitive volume: invasion order is that of IP" begin
    for lat in (SquareLattice(32), TriangularLattice(32), CubicLattice(6)),
        perc in (SitePercolation(), BondPercolation()), periodic in (false, true)

        kw = (percolation=perc, periodic, seeds=6, right=FixedRight(32))
        ip = simulate(IPModel(lat; kw...); rng=Xoshiro(4))
        cv = simulate(CVModel(lat; compliance=50.0, meniscus=0.3, kw...); rng=Xoshiro(4))
        @test cv.x == ip.x
        @test cv.status === ip.status === :spanned
        @test cv.mass == ip.mass[1]
    end
end

@testset "capacitive volume: menisci on the front" begin
    for lat in (SquareLattice(16), TriangularLattice(16)),
        perc in (SitePercolation(), BondPercolation()), frac in (1.0, 0.5)

        medium = Medium(5)
        kw = (percolation=perc, periodic=true, seeds=8, active_fraction=frac,
            right=GrowingRight(8))
        full = simulate(CVModel(lat; compliance=10.0, meniscus=1.0, nmax=300, kw...);
            rng=Xoshiro(2), medium)
        st = SelfOrganizedPercolation.initial_state(IPModel(lat; kw...), Xoshiro(2),
            medium, false)
        for k in (1, 2, 17, 150, 300)
            r = simulate(CVModel(lat; compliance=10.0, nmax=k - 1, kw...);
                rng=Xoshiro(2), medium, keep_sites=true, keep_times=true)
            contributors = vcat(st.active, findall(>(0), vec(r.times)))
            nx = size(r.sites, 2)
            @test full.nf[k] == brute_menisci(vec(r.sites), lat, nx, true,
                contributors, perc)
        end
    end
end

@testset "capacitive volume: volume balance and bursts" begin
    lat = SquareLattice(64)
    # Constant compliance: C (p - p0) = rate t - (pores invaded), at every invasion.
    cv = CVModel(lat; compliance=200.0, rate=3.0, p0=0.1, seeds=64, periodic=true)
    r = simulate(cv; rng=Xoshiro(7))
    k = 1:length(r.x)
    @test all(@. isapprox(200.0 * (r.p - 0.1), 3.0 * r.t - k; atol=1e-8))
    # With menisci, step by step.
    cm = CVModel(lat; compliance=20.0, meniscus=0.5, rate=2.0, seeds=64)
    s = simulate(cm; rng=Xoshiro(7))
    prev, tprev = 0.0, 0.0
    for n in eachindex(s.x)
        C = 20.0 + 0.5 * s.nf[n]
        @test s.t[n] - tprev ≈ C * max(s.x[n] - prev, 0.0) / 2.0
        @test s.p[n] ≈ max(prev, s.x[n]) - 1 / C
        prev, tprev = s.p[n], s.t[n]
    end
    b = bursts(s)
    @test sum(b.size) == length(s.x)
    @test issorted(b.start) && all(>(0), b.drop)
    @test b.peak == s.x[b.start]
    @test all(diff(b.time) .> 0)
    # A stiff system cannot sustain a burst: one pore per burst.
    stiff = simulate(CVModel(lat; compliance=0.01, seeds=64); rng=Xoshiro(7))
    @test all(==(1), bursts(stiff).size)
    # A soft one has bursts of many pores, and its peak pressure stays near p_c.
    soft = simulate(CVModel(lat; compliance=5000.0, seeds=64, periodic=true);
        rng=Xoshiro(7))
    @test maximum(bursts(soft).size) > 20
    @test abs(maximum(soft.x) - 0.5) < 0.05     # bond p_c on the square lattice
end

@testset "capacitive volume: arguments, stability, ensemble" begin
    lat = SquareLattice(16)
    @test_throws ArgumentError CVModel(lat; compliance=1.0, species=2)
    @test_throws ArgumentError CVModel(lat; compliance=0.0)
    @test_throws ArgumentError CVModel(lat; compliance=1.0, rate=0.0)
    for l in (SquareLattice(16), CubicLattice(5)), perc in (SitePercolation(), BondPercolation())
        m = CVModel(l; compliance=5.0, meniscus=0.2, percolation=perc, seeds=4)
        @test (@inferred simulate(m; rng=Xoshiro(1))) isa CVResult
    end
    quiet = simulate(CVModel(lat; compliance=5.0, seeds=4); keep_history=false)
    @test isempty(quiet.x) && quiet.status === :spanned
    sizes = ensemble(CVModel(lat; compliance=5.0, seeds=4), 4; f=r -> length(r.x))
    @test length(sizes) == 4
end

# Empty sites that reach the last column through empty sites.
function open_to_outlet(sites, lat, nx, periodic)
    W = length(sites) ÷ nx
    seen = falses(length(sites))
    stack = [i for i in W*(nx-1)+1:W*nx if sites[i] == EMPTY]
    seen[stack] .= true
    while !isempty(stack)
        for j in neighbors(lat, pop!(stack), nx, periodic)
            j != 0 && sites[j] == EMPTY && !seen[j] && (seen[j] = true; push!(stack, j))
        end
    end
    return seen
end

# Invasion percolation with trapping by brute force: after every invasion, every
# empty site that no longer reaches the outlet is trapped.
function brute_trapping(lat, W, nx, periodic, seeds, medium, perc)
    sites = zeros(UInt8, W * nx)
    sites[seeds] .= 0x01
    order = Int[]
    while true
        best, site = Inf, 0
        for i in findall(==(0x01), sites), j in neighbors(lat, i, nx, periodic)
            (j == 0 || sites[j] != EMPTY) && continue
            x = perc isa SitePercolation ? threshold(medium, j) : threshold(medium, i, j)
            x < best && ((best, site) = (x, j))
        end
        site == 0 && break
        sites[site] = 0x01
        push!(order, site)
        (site - 1) ÷ W == nx - 1 && break
        reach = open_to_outlet(sites, lat, nx, periodic)
        for u in eachindex(sites)
            sites[u] == EMPTY && !reach[u] && (sites[u] = TRAPPED)
        end
    end
    return order, sites
end

# The trapping is exact: every empty site reaches the outlet, and no trapped
# site touches an empty one.
function trapping_exact(sites, lat, periodic)
    nx = size(sites)[end]
    v = vec(sites)
    reach = open_to_outlet(v, lat, nx, periodic)
    all(u -> v[u] != EMPTY || reach[u], eachindex(v)) || return false
    for u in eachindex(v)
        v[u] == TRAPPED || continue
        any(j -> j != 0 && v[j] == EMPTY, neighbors(lat, u, nx, periodic)) && return false
    end
    return true
end

@testset "invasion percolation with trapping matches brute force" begin
    ntrapped = 0
    for lat in (SquareLattice(14), TriangularLattice(14), CubicLattice(5, 4)),
        perc in (SitePercolation(), BondPercolation()), periodic in (false, true)

        m = IPModel(lat; percolation=perc, periodic, seeds=4, right=FixedRight(20),
            trapping=true)
        medium = Medium(8)
        r = simulate(m; rng=Xoshiro(3), medium, keep_sites=true, keep_times=true)
        seeds = findall(==(0), vec(r.times))
        W = SelfOrganizedPercolation.transverse_sites(lat)
        order, sites = brute_trapping(lat, W, 20, periodic, seeds, medium, perc)
        invaded = sortperm(vec(r.times))[length(seeds)+count(<(0), r.times)+1:end]
        @test invaded == order
        @test vec(r.sites) == sites
        ntrapped += count(==(TRAPPED), r.sites)
    end
    @test ntrapped > 100
end

@testset "trapping is exact in every model" begin
    for lat in (SquareLattice(48), TriangularLattice(48), CubicLattice(8)),
        perc in (SitePercolation(), BondPercolation()), periodic in (false, true)

        kw = (percolation=perc, periodic, seeds=6, right=GrowingRight(8), trapping=true)
        ip = simulate(IPModel(lat; nmax=4000, species=2, kw...); rng=Xoshiro(2),
            keep_sites=true)
        @test trapping_exact(ip.sites, lat, periodic)
        @test count(==(TRAPPED), ip.sites) > 0
        for rej in (Permanent(), Persistent())
            sop = simulate(SOPModel(lat; rejection=rej, species=2, f_T=0.1, c=0.1,
                p0=0.6, tmax=300, kw...); rng=Xoshiro(2), keep_sites=true)
            @test trapping_exact(sop.sites, lat, periodic)
            for s in 1:2
                @test count(==(UInt8(s)), sop.sites) == sop.mass[s]
            end
        end
        cv = simulate(CVModel(lat; compliance=50.0, meniscus=0.5, nmax=4000, kw...);
            rng=Xoshiro(2), keep_sites=true)
        @test trapping_exact(cv.sites, lat, periodic)
        @test cv.x == simulate(IPModel(lat; nmax=4000, kw...); rng=Xoshiro(2)).x
    end
    # Off by default, and nothing is trapped then.
    r = simulate(IPModel(SquareLattice(32); seeds=8); rng=Xoshiro(1), keep_sites=true)
    @test !any(==(TRAPPED), r.sites)
    for lat in (SquareLattice(32), CubicLattice(6)), perc in (SitePercolation(), BondPercolation())
        @test (@inferred simulate(SOPModel(lat; percolation=perc, seeds=4, tmax=30,
            rejection=Persistent(), trapping=true); rng=Xoshiro(1))) isa SOPResult
        @test (@inferred simulate(IPModel(lat; percolation=perc, seeds=4,
            trapping=true); rng=Xoshiro(1))) isa IPResult
    end
end

@testset "capacitive volume with trapping: menisci and a steady front" begin
    for lat in (SquareLattice(16), TriangularLattice(16)),
        perc in (SitePercolation(), BondPercolation())

        medium = Medium(5)
        kw = (percolation=perc, periodic=true, seeds=8, right=GrowingRight(8),
            trapping=true)
        full = simulate(CVModel(lat; compliance=10.0, meniscus=1.0, nmax=600, kw...);
            rng=Xoshiro(2), medium)
        st = SelfOrganizedPercolation.initial_state(IPModel(lat; kw...), Xoshiro(2),
            medium, false)
        for k in (2, 40, 300, 600)
            r = simulate(CVModel(lat; compliance=10.0, nmax=k - 1, kw...);
                rng=Xoshiro(2), medium, keep_sites=true, keep_times=true)
            contributors = vcat(st.active, findall(>(0), vec(r.times)))
            @test full.nf[k] == brute_menisci(vec(r.sites), lat, size(r.sites, 2),
                true, contributors, perc)
        end
    end
    # Furuberg et al.: without trapping the front keeps growing; with it, it settles.
    lat = SquareLattice(32)
    kw = (seeds=32, periodic=true, right=GrowingRight(16), nmax=30000)
    free = simulate(CVModel(lat; compliance=0.0, meniscus=1.0, kw...); rng=Xoshiro(1))
    trap = simulate(CVModel(lat; compliance=0.0, meniscus=1.0, trapping=true, kw...);
        rng=Xoshiro(1))
    @test free.nf[end] > 1.6 * free.nf[end÷2]
    @test trap.nf[end] < 1.3 * trap.nf[end÷2]
end
