# Invasion percolation: cost per invasion and the mass at breakthrough, as in Masson
# & Pride (2014): fluid enters through the whole left wall, the sides are periodic,
# and the run stops when the fluid first reaches the right wall.
using SelfOrganizedPercolation, Random, Statistics, Printf

function breakthrough(lat, perc, L, n)
    m = IPModel(lat(L); percolation=perc, periodic=true, seeds=L, right=FixedRight(L))
    simulate(IPModel(lat(64); percolation=perc, periodic=true, seeds=64); rng=Xoshiro(0))   # compile
    M, T = Float64[], Float64[]
    for i in 1:n
        t = @elapsed r = simulate(m; rng=Xoshiro(i), keep_history=false)
        r.status === :spanned || error("run $i did not span")
        push!(M, sum(r.mass) - L)
        push!(T, t)
    end
    return mean(M), mean(T)
end

for (name, lat, perc) in (("square site", SquareLattice, SitePercolation()),
        ("square bond", SquareLattice, BondPercolation()),
        ("triangular site", TriangularLattice, SitePercolation()))
    println("\n", name, ": mass M at breakthrough and CPU time")
    Ls, Ms = Int[], Float64[]
    for L in (256, 512, 1024, 2048, 4096, 8192)
        M, T = breakthrough(lat, perc, L, L <= 1024 ? 16 : 4)
        push!(Ls, L); push!(Ms, M)
        @printf("  L = %5d  M = %.3e  T = %7.3f s  T/(M log2 M) = %.1f ns\n", L, M, T, 1e9T / (M * log2(M)))
    end
    D = cov(log.(Ls), log.(Ms)) / var(log.(Ls))
    @printf("  M ~ L^D with D = %.3f (no trapping: 91/48 = 1.896)\n", D)
end
