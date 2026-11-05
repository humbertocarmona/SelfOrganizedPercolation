# Throughput of one realization and ensemble scaling with the number of tasks.
using SelfOrganizedPercolation, Random, Printf

function rate(m; n=1)
    simulate(m; rng=Xoshiro(1))   # compile first
    t = @elapsed r = simulate(m; rng=Xoshiro(2))
    return sum(r.mass), t
end
for (name, lat) in (("square", SquareLattice), ("triangular", TriangularLattice)),
    (pn, perc) in (("site", SitePercolation()), ("bond", BondPercolation())),
    (rn, rej) in (("permanent", Permanent()), ("persistent", Persistent()))

    W = 16384
    m = SOPModel(lat(W); percolation=perc, rejection=rej, periodic=true, seeds=W ÷ 5, p0=0.7,
        k=0.05 / W, target=0.1W, right=GrowingRight(1024), tmax=3000)
    mass, t = rate(m)
    @printf("%-10s %-4s %-10s  W=%d  %.1f M occupied sites/s  (%.2f s)\n", name, pn, rn, W, mass / t / 1e6, t)
end
m = SOPModel(SquareLattice(2048); percolation=SitePercolation(), periodic=true, seeds=400, p0=0.7,
    k=0.05 / 2048, target=205, right=GrowingRight(1024), tmax=4000)
ensemble(m, 2; ntasks=1)
base = 0.0
for n in (1, 2, 4, 8)
    n > Threads.nthreads() && continue
    t = @elapsed ensemble(m, 64; ntasks=n)
    n == 1 && (global base = t)
    @printf("ensemble of 64, %d task(s): %.2f s, speedup %.2f\n", n, t, base / t)
end
