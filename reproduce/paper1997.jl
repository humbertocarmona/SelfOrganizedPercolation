# Reproduces Figs. 1-4 of Alencar, Andrade & Lucena, PRE 56, R2379 (1997).
#
#   julia -t auto --project=reproduce -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#   julia -t auto --project=reproduce reproduce/paper1997.jl
#
# The paper seeds the bottom line and grows upward. Here the fluid enters at the
# left wall and grows to the right, which is the same model rotated.
using SelfOrganizedPercolation, CairoMakie, Random, Statistics, Printf

const OUT = joinpath(@__DIR__, "output")
mkpath(OUT)
const BLUE, ORANGE, AQUA = colorant"#2a78d6", colorant"#eb6834", colorant"#1baf7a"
const COLORS = [BLUE, ORANGE, AQUA]
const INK, MUTED = colorant"#0b0b0b", colorant"#52514e"

set_theme!(Theme(fontsize=15, Axis=(xgridvisible=false, ygridcolor=(:black, 0.07),
    topspinevisible=false, rightspinevisible=false, spinewidth=1, xtickwidth=1, ytickwidth=1,
    xticksize=4, yticksize=4, leftspinecolor=MUTED, bottomspinecolor=MUTED,
    xticklabelcolor=MUTED, yticklabelcolor=MUTED, titlealign=:left, titlefont=:regular,
    titlecolor=INK), Lines=(linewidth=2,)))

const L = 2000
const TBURN = 300
const NENS = 64

# Many runs die early, and the paper can only show ones that survived. Single-run
# figures therefore use the first surviving seed, and the extinct fraction is reported.
function survivor(model, first_seed)
    for seed in first_seed:first_seed+10_000
        r = simulate(model; rng=Xoshiro(seed), keep_sites=false)
        r.status !== :extinct && return r
    end
    error("no surviving run")
end

function ensemble_row(label, model, reference; n=NENS)
    means = Float64[]
    extinct = 0
    results = ensemble(model, n; f=r -> (r.status, r.status === :extinct ? NaN : mean(@view r.p[TBURN+1:end, 1])))
    for (status, value) in results
        status === :extinct ? (extinct += 1) : push!(means, value)
    end
    m, se = mean(means), std(means) / sqrt(length(means))
    @printf("%-36s paper %-6s this work %.4f ± %.4f  (%d of %d extinct)\n", label, reference, m, se, extinct, n)
    return (; label, reference, mean=m, se, extinct, n)
end

# First minimum and following maximum of p, averaged over surviving runs.
function transient_row(label, model, reference; n=200)
    results = ensemble(model, n; f=r -> begin
        P = @view r.p[:, 1]
        w = min(length(P), 250)
        low = argmin(@view P[1:min(w, 60)])
        high = low - 1 + argmax(@view P[low:w])
        (r.status === :extinct, P[low], low - 1, P[high], high - 1)
    end)
    alive = filter(x -> !x[1], results)
    @printf("%-36s paper %-24s this work min %.3f at t = %3.0f, peak %.3f at t = %3.0f\n", label, reference,
        mean(x[2] for x in alive), mean(x[3] for x in alive), mean(x[4] for x in alive), mean(x[5] for x in alive))
end

site(; p0, k=1e-4, target=200, tmax=1000, kw...) = SOPModel(SquareLattice(L);
    percolation=SitePercolation(), p0, k, target, tmax, seeds=100, kw...)
bond(; p0, k=1e-4, target=200, tmax=1000, kw...) = SOPModel(SquareLattice(L);
    percolation=BondPercolation(), p0, k, target, tmax, seeds=100, kw...)

# ---- Fig. 1: cluster, bond percolation, L = 1000, p0 = 1, k = 1e-5, N_T = 200
fig1 = SOPModel(SquareLattice(1000); percolation=BondPercolation(), seeds=100, p0=1.0,
    k=1e-5, target=200)
r1 = simulate(fig1; rng=Xoshiro(1), keep_sites=true)
f = Figure(size=(760, 700), backgroundcolor=:white)
ax = Axis(f[1, 1]; aspect=DataAspect(), xlabel="x, injection direction", ylabel="y",
    title="Fig. 1: bond percolation, L = 1000, p₀ = 1, k = 10⁻⁵, N_T = 200, 100 seeds")
heatmap!(ax, 1:1000, 1:1000, permutedims(Float64.(r1.sites .!= 0));
    colormap=[colorant"#ffffff", colorant"#0b0b0b"], colorrange=(0, 1), interpolate=false)
save(joinpath(OUT, "fig1_cluster.png"), f; px_per_unit=2)
println("Fig. 1: ", r1.status, " after ", size(r1.p, 1), " steps, mass ", r1.mass[1])

# ---- Fig. 2: p(t) for several p0
function fig2panel!(slot, title, pc, ps, build)
    ax = Axis(slot; title, xlabel="time", ylabel="p", limits=(0, 1000, 0.2, 1.0))
    hlines!(ax, [pc]; color=MUTED, linestyle=:dash, linewidth=1, label="p_c = $pc")
    for (i, p0) in enumerate(ps)
        r = simulate(build(p0); rng=Xoshiro(10 + i))
        lines!(ax, 0:size(r.p, 1)-1, r.p[:, 1]; color=COLORS[i], label="p₀ = $p0")
    end
    axislegend(ax; framevisible=false, position=:rt)
    return ax
end
f = Figure(size=(760, 760), backgroundcolor=:white)
fig2panel!(f[1, 1], "Fig. 2(a): site percolation, L = 2000, k = 10⁻⁴, N_T = 200", 0.59274,
    (0.4, 0.8, 1.0), p0 -> site(; p0))
fig2panel!(f[2, 1], "Fig. 2(b): bond percolation, L = 2000, k = 10⁻⁴, N_T = 200", 0.5,
    (0.3, 0.7, 1.0), p0 -> bond(; p0))
save(joinpath(OUT, "fig2_probability.png"), f; px_per_unit=2)

# ---- Fig. 3: N(t). The paper does not state the lattice or p0 for this figure:
# site percolation, L = 2000 and p0 = 1 are assumptions.
cases = ((300, 1e-3, "(a) N_T = 300, k = 10⁻³"), (200, 1e-4, "(b) N_T = 200, k = 10⁻⁴"),
    (100, 1e-5, "(c) N_T = 100, k = 10⁻⁵"))
f = Figure(size=(760, 460), backgroundcolor=:white)
ax = Axis(f[1, 1]; title="Fig. 3: active sites (site percolation, L = 2000, p₀ = 1 assumed)",
    xlabel="time", ylabel="N", limits=(0, 2000, 0, 900))
for (i, (NT, k, label)) in enumerate(cases)
    r = survivor(site(; p0=1.0, k, target=NT, tmax=2000), 20 + i)
    lines!(ax, 0:size(r.N, 1)-1, r.N[:, 1]; color=COLORS[i], linewidth=1, label)
    println("Fig. 3 ", label, ": ", r.status, ", steps ", size(r.N, 1),
        ", std(N) ", round(std(r.N[min(TBURN, end ÷ 2):end, 1]); digits=1))
end
axislegend(ax; framevisible=false, position=:rt)
save(joinpath(OUT, "fig3_active_sites.png"), f; px_per_unit=2)

# ---- Fig. 4: site percolation with N_T(t) = 10 t^(1/2), k = 1e-4, p0 = 0.8
m4 = site(; p0=0.8, target=PowerLawTarget(10, 0.5), tmax=2000)
r4 = survivor(m4, 30)
f = Figure(size=(760, 760), backgroundcolor=:white)
ax = Axis(f[1, 1]; title="Fig. 4(a): active sites against the target N_T(t) = 10 t^(1/2)",
    xlabel="time", ylabel="N", limits=(0, 2000, 0, 600))
t = 0:size(r4.N, 1)-1
lines!(ax, t, r4.N[:, 1]; color=BLUE, linewidth=1, label="N(t)")
lines!(ax, t, 10 .* sqrt.(t); color=INK, linewidth=2, label="N_T(t)")
axislegend(ax; framevisible=false, position=:lt)
ax = Axis(f[2, 1]; title="Fig. 4(b): occupation probability", xlabel="time", ylabel="p",
    limits=(0, 2000, 0.4, 0.9))
hlines!(ax, [0.59274]; color=MUTED, linestyle=:dash, linewidth=1, label="p_c = 0.59274")
lines!(ax, t, r4.p[:, 1]; color=BLUE, label="p(t)")
axislegend(ax; framevisible=false, position=:rt)
save(joinpath(OUT, "fig4_power_law_target.png"), f; px_per_unit=2)

# ---- Quantitative comparison: ensemble means after the transient
println("\nStationary ⟨p⟩, t = $TBURN to tmax, $NENS realizations")
rows = [
    ensemble_row("Fig. 2(a) site, p0 = 0.4", site(; p0=0.4), "≈0.59"),
    ensemble_row("Fig. 2(a) site, p0 = 0.8", site(; p0=0.8), "≈0.59"),
    ensemble_row("Fig. 2(a) site, p0 = 1.0", site(; p0=1.0), "≈0.59"),
    ensemble_row("Fig. 2(b) bond, p0 = 0.3", bond(; p0=0.3), "≈0.5"),
    ensemble_row("Fig. 2(b) bond, p0 = 0.7", bond(; p0=0.7), "≈0.5"),
    ensemble_row("Fig. 2(b) bond, p0 = 1.0", bond(; p0=1.0), "≈0.5"),
    ensemble_row("Fig. 3(a) site, N_T = 300, k = 1e-3", site(; p0=1.0, k=1e-3, target=300, tmax=2000), "n/a"),
    ensemble_row("Fig. 3(b) site, N_T = 200, k = 1e-4", site(; p0=1.0, k=1e-4, target=200, tmax=2000), "n/a"),
    ensemble_row("Fig. 3(c) site, N_T = 100, k = 1e-5", site(; p0=1.0, k=1e-5, target=100, tmax=2000), "n/a"; n=400),
    ensemble_row("Fig. 4 site, N_T = 10√t", site(; p0=0.8, target=PowerLawTarget(10, 0.5), tmax=2000), "≈0.59"; n=400),
]
open(joinpath(OUT, "stationary_p.csv"), "w") do io
    println(io, "case;paper;mean;standard_error;extinct;runs")
    foreach(r -> @printf(io, "%s;%s;%.5f;%.5f;%d;%d\n", r.label, r.reference, r.mean, r.se, r.extinct, r.n), rows)
end

# ---- Transients. Paper values are read off the published figures by eye.
println("\nTransient of p, surviving runs")
power = PowerLawTarget(10, 0.5)
transient_row("Fig. 4(b) site, stated parameters", site(; p0=0.8, target=power, tmax=2000), "min ≈0.42, peak ≈0.82"; n=400)
transient_row("Fig. 2(a) site, p0 = 0.4, k = 1e-4", site(; p0=0.4), "peak ≈0.84 at t ≈ 70")
transient_row("Fig. 2(b) bond, p0 = 0.3, k = 1e-4", bond(; p0=0.3), "peak ≈0.92 at t ≈ 60")
transient_row("  site, p0 = 0.4, k = 5e-5 (diagnostic)", site(; p0=0.4, k=5e-5), "peak ≈0.84 at t ≈ 70")
transient_row("  bond, p0 = 0.3, k = 5e-5 (diagnostic)", bond(; p0=0.3, k=5e-5), "peak ≈0.92 at t ≈ 60")
