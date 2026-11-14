# Invasion percolation: cost per invasion and the mass at breakthrough, as in Masson
# & Pride (2014): fluid enters through the whole left wall, the sides are periodic,
# and the run stops when the fluid first reaches the right wall.
using SelfOrganizedPercolation, Random, Statistics, Printf

wall(lat, L) = lat isa Type{CubicLattice} ? L^2 : L

function breakthrough(lat, perc, L, n)
    m = IPModel(lat(L); percolation=perc, periodic=true, seeds=wall(lat, L), right=FixedRight(L))
    simulate(IPModel(lat(16); percolation=perc, periodic=true, seeds=wall(lat, 16)); rng=Xoshiro(0))   # compile
    M, T = Float64[], Float64[]
    for i in 1:n
        t = @elapsed r = simulate(m; rng=Xoshiro(i), keep_history=false)
        r.status === :spanned || error("run $i did not span")
        push!(M, sum(r.mass) - wall(lat, L))
        push!(T, t)
    end
    return mean(M), mean(T)
end

function report(name, lat, perc, sizes, D0)
    println("\n", name, ": mass M at breakthrough and CPU time")
    Ls, Ms = Int[], Float64[]
    for L in sizes
        M, T = breakthrough(lat, perc, L, L <= sizes[3] ? 16 : 4)
        push!(Ls, L); push!(Ms, M)
        @printf("  L = %5d  M = %.3e  T = %7.3f s  T/(M log2 M) = %.1f ns\n", L, M, T, 1e9T / (M * log2(M)))
    end
    D = cov(log.(Ls), log.(Ms)) / var(log.(Ls))
    @printf("  M ~ L^D with D = %.3f (no trapping: %s)\n", D, D0)
end

L2 = (256, 512, 1024, 2048, 4096, 8192)
report("square site", SquareLattice, SitePercolation(), L2, "91/48 = 1.896")
report("square bond", SquareLattice, BondPercolation(), L2, "91/48 = 1.896")
report("triangular site", TriangularLattice, SitePercolation(), L2, "91/48 = 1.896")
# Masson & Pride (2014) find M = 0.407 N^0.852 on cubes of N = L^3 sites, that is D = 2.56.
L3 = (32, 64, 128, 256, 384)
report("cubic site", CubicLattice, SitePercolation(), L3, "2.53; 2.56 in Masson & Pride")
report("cubic bond", CubicLattice, BondPercolation(), L3, "2.53")
