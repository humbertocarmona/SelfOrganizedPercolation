# Checks of what the 1997 paper does not cover: known thresholds on other lattices,
# the Persistent rule against the independent C simulator, and two immiscible fluids.
using SelfOrganizedPercolation, Random, Statistics, Printf

const W = 4096
function plateau(lat, perc, pc, name; n=16)
    # Small target fraction, periodic sides and a right wall that follows the front.
    for fT in (0.2, 0.1, 0.05)
        m = SOPModel(lat; percolation=perc, periodic=true, seeds=W ÷ 5, p0=pc + 0.05,
            c=0.05, f_T=fT, right=GrowingRight(1024), tmax=8000)
        res = ensemble(m, n; f=r -> (r.status, stationary_mean(@view r.p[:, 1]; burn=size(r.p, 1) ÷ 2)))
        alive = [v for (s, v) in res if s !== :extinct]
        @printf("%-22s f_T = %.2f  p* = %.4f ± %.4f   p_c = %.4f   (%d of %d survive)\n", name, fT, mean(alive),
            std(alive) / sqrt(length(alive)), pc, length(alive), n)
    end
end
plateau(SquareLattice(W), SitePercolation(), 0.59274, "square site")
plateau(SquareLattice(W), BondPercolation(), 0.5, "square bond")
plateau(TriangularLattice(W), SitePercolation(), 0.5, "triangular site")
plateau(TriangularLattice(W), BondPercolation(), 2sin(π / 18), "triangular bond")

# Cubic lattice: wall of 128 × 128 sites, p_c from Wang et al., PRE 87, 052107 (2013).
for (perc, pc, name) in ((SitePercolation(), 0.3116, "cubic site"), (BondPercolation(), 0.2488, "cubic bond"))
    for fT in (0.2, 0.1, 0.05)
        m = SOPModel(CubicLattice(128); percolation=perc, periodic=true, seeds=128^2 ÷ 5, p0=pc + 0.05,
            c=0.05, f_T=fT, right=GrowingRight(256), tmax=3000)
        res = ensemble(m, 16; f=r -> (r.status, stationary_mean(@view r.p[:, 1])))
        alive = [v for (s, v) in res if s !== :extinct]
        @printf("%-22s f_T = %.2f  p* = %.4f ± %.4f   p_c = %.4f   (%d of %d survive)\n", name, fT, mean(alive),
            std(alive) / sqrt(length(alive)), pc, length(alive), 16)
    end
end

println("\nPersistent rule, bond, L = 2048, H = 2048, L/20 seeds. C simulator in brackets.")
for (fT, c, cmean, cpeak) in ((0.05, 0.05, 0.4881, 0.5014), (0.2, 0.05, 0.4943, 0.5058), (0.05, 0.2, 0.4793, 0.5015))
    m = SOPModel(SquareLattice(2048); rejection=Persistent(), seeds=102, p0=1.0, c=c, f_T=fT)
    res = ensemble(m, 12; f=r -> begin
        # As in the C code: average once the front has passed a quarter of the lattice.
        P = @view r.p[:, 1]
        tail = @view P[size(r.p, 1) ÷ 3:end]
        (mean(tail), maximum(tail))
    end)
    @printf("f_T = %.2f c = %.2f  mean %.4f [%.4f]  peak %.4f [%.4f]\n", fT, c, mean(first, res), cmean, mean(last, res), cpeak)
end

println("\nTwo identical immiscible fluids, bond, square, L = 2048, periodic, f_T = 0.1 each")
m = SOPModel(SquareLattice(2048); species=2, periodic=true, seeds=205, p0=0.7, c=0.05, f_T=205 / 2048,
    right=GrowingRight(1024), tmax=8000)
res = ensemble(m, 32; f=r -> (r.status, vec(mean(@view(r.p[size(r.p, 1)÷2:end, :]); dims=1)), r.mass))
for (i, (s, ps, mass)) in enumerate(res[1:4])
    @printf("run %d: %-8s p* = (%.4f, %.4f)  mass = (%d, %d)\n", i, s, ps..., mass...)
end
println("statuses: ", Dict(s => count(r -> r[1] === s, res) for s in unique(first.(res))))
