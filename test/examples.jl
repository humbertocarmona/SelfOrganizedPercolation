# Examples of what the package can run. Each example is a function that prints a
# short summary and, where a picture helps, saves a figure to test/output/examples/.
#
#     julia -t auto --project=reproduce test/examples.jl           # all examples
#     julia -t auto --project=reproduce test/examples.jl 3 7       # only examples 3 and 7
#
# Including the file only defines the functions, so one example can be run from the
# REPL with, say, `include("test/examples.jl"); ex_fluids()`.
using SelfOrganizedPercolation, CairoMakie, Random, Statistics, Printf

const OUT = joinpath(@__DIR__, "output", "examples")
const FLUID_COLORS = [colorant"#ffffff", colorant"#2a78d6", colorant"#eb6834", colorant"#1baf7a",
    colorant"#8a5cc2", colorant"#c9a227"]

# Sites colored by fluid (white is empty), x along the horizontal axis.
function plot_sites(sites, title; z=nothing)
    A = ndims(sites) == 3 ? sites[:, something(z, size(sites, 2) ÷ 2 + 1), :] : sites
    nfluids = max(Int(maximum(A)), 1)
    f = Figure(size=(700, 620), backgroundcolor=:white)
    ax = Axis(f[1, 1]; aspect=DataAspect(), title, titlealign=:left, titlefont=:regular,
        xlabel="x, injection direction", ylabel="y")
    heatmap!(ax, permutedims(Int.(A)); colormap=FLUID_COLORS[1:nfluids+1], colorrange=(0, nfluids),
        interpolate=false)
    return f
end

# Arrival times, white where nothing arrived.
function plot_times(times, title; z=nothing)
    A = ndims(times) == 3 ? times[:, something(z, size(times, 2) ÷ 2 + 1), :] : times
    f = Figure(size=(760, 620), backgroundcolor=:white)
    ax = Axis(f[1, 1]; aspect=DataAspect(), title, titlealign=:left, titlefont=:regular,
        xlabel="x, injection direction", ylabel="y")
    hm = heatmap!(ax, permutedims(map(t -> t < 0 ? NaN32 : Float32(t), A)); colormap=:viridis,
        nan_color=:white, interpolate=false)
    Colorbar(f[1, 2], hm; label="arrival time")
    return f
end

function savefig(name, f)
    mkpath(OUT)
    file = joinpath(OUT, name)
    save(file, f; px_per_unit=1)
    println("   figure: ", file)
end

# ---------------------------------------------------------------------------- SOP

"1. The original SOP of the 1997 paper: site percolation, one fluid, permanent rejection."
function ex_basic()
    m = SOPModel(SquareLattice(1000); percolation=SitePercolation(), seeds=100, p0=1.0,
        c=0.1, f_T=0.2, tmax=2000)
    r = simulate(m; rng=Xoshiro(1), keep_sites=true)
    @printf("   status %s after %d steps, mass %d, p* = %.4f (p_c = 0.5927)\n", r.status,
        size(r.p, 1), r.mass[1], stationary_mean(@view r.p[:, 1]))
    f = Figure(size=(1300, 560), backgroundcolor=:white)
    ax = Axis(f[1, 1]; aspect=DataAspect(), title="cluster", xlabel="x", ylabel="y")
    heatmap!(ax, permutedims(Float64.(r.sites .!= 0)); colormap=[:white, :black], interpolate=false)
    ax2 = Axis(f[1, 2]; title="p(t) and f(t)", xlabel="t")
    lines!(ax2, r.p[:, 1]; label="p(t)")
    lines!(ax2, r.f[:, 1]; label="f(t)")
    hlines!(ax2, [0.5927]; color=:gray, linestyle=:dash, label="p_c")
    axislegend(ax2; position=:rt)
    savefig("01_basic.png", f)
end

"2. A time-dependent target f_T(t) = a t^α (Fig. 4 of the paper), keeping the first survivor."
function ex_target()
    m = SOPModel(SquareLattice(1000); percolation=SitePercolation(), seeds=100, p0=0.8, c=0.1,
        f_T=PowerLawTarget(0.01, 0.5), tmax=2000)
    # Starting from f_T(0) = 0, many fronts die before the target grows; skip those.
    seed = findfirst(k -> simulate(m; rng=Xoshiro(k)).status !== :extinct, 1:1000)
    r = simulate(m; rng=Xoshiro(seed))
    @printf("   %d runs died first; survivor %s after %d steps, p at t = 500, 1000, end: %.4f %.4f %.4f (p_c = 0.5927)\n",
        seed - 1, r.status, size(r.p, 1), r.p[min(500, end), 1], r.p[min(1000, end), 1], r.p[end, 1])
    f = Figure(size=(800, 420), backgroundcolor=:white)
    ax = Axis(f[1, 1]; xlabel="t", title="f(t) against f_T(t) = 0.01 √t", titlealign=:left, titlefont=:regular)
    lines!(ax, r.f[:, 1]; label="f(t)")
    lines!(ax, 0:size(r.f, 1)-1, t -> 0.01sqrt(t); label="f_T(t)")
    axislegend(ax; position=:lt)
    savefig("02_target.png", f)
end

"3. Permanent and Persistent rejection: with memory, a front with a small target survives."
function ex_rejection()
    for rej in (Permanent(), Persistent())
        m = SOPModel(SquareLattice(512); rejection=rej, seeds=20, p0=0.5, c=0.1, f_T=0.02,
            right=GrowingRight(64), tmax=3000)
        res = ensemble(m, 32; f=r -> (r.status, r.length))
        alive = count(x -> x[1] !== :extinct, res)
        @printf("   %-12s %2d of 32 fronts survive; mean penetration %.0f columns\n",
            nameof(typeof(rej)), alive, mean(last, res))
    end
end

"4. Ensembles: the stationary p* over many runs, and how many of them survive."
function ex_ensemble()
    for fT in (0.2, 0.1, 0.05)
        m = SOPModel(SquareLattice(1024); percolation=SitePercolation(), periodic=true, seeds=200,
            p0=0.65, c=0.05, f_T=fT, right=GrowingRight(256), tmax=3000)
        res = ensemble(m, 16; f=r -> (r.status, stationary_mean(@view r.p[:, 1])))
        ps = [p for (s, p) in res if s !== :extinct]
        @printf("   f_T = %.2f  p* = %.4f ± %.4f  (%d of 16 survive)\n", fT, mean(ps),
            std(ps) / sqrt(length(ps)), length(ps))
    end
end

"5. Several immiscible fluids competing for the medium, with different starting p0."
function ex_fluids()
    m = SOPModel(SquareLattice(600); percolation=BondPercolation(), rejection=Persistent(), species=3,
        seeds=40, p0=[0.55, 0.6, 0.65], c=0.05, f_T=0.05, periodic=true, right=FixedRight(600))
    r = simulate(m; rng=Xoshiro(5), keep_sites=true)
    @printf("   status %s after %d steps; masses %s; final p %s\n", r.status, size(r.p, 1),
        r.mass, round.(r.p[end, :]; digits=4))
    savefig("05_fluids.png", plot_sites(r.sites, "three fluids, bond SOP (blue, orange, green)"))
end

"6. Seed layouts: random seeds against equal blocks along the wall, with arrival times."
function ex_layouts()
    for (name, layout) in (("random", RandomSeeds()), ("blocks", BlockSeeds(4)))
        m = SOPModel(SquareLattice(512); percolation=SitePercolation(), seeds=128, layout,
            p0=0.7, c=0.1, f_T=0.1, right=FixedRight(512))
        r = simulate(m; rng=Xoshiro(6), keep_times=true)
        @printf("   %-7s status %s, mass %d, last arrival %d\n", name, r.status, r.mass[1], maximum(r.times))
        savefig("06_layout_$name.png", plot_times(r.times, "SOP from $name seeds, colored by arrival step"))
    end
end

"7. Inactive seeds: only part of each fluid's wall sites inject; the rest only block."
function ex_active()
    m = SOPModel(SquareLattice(512); rejection=Persistent(), species=2, seeds=100, active_fraction=0.2, p0=0.7,
        c=0.05, f_T=0.05, periodic=true, right=GrowingRight(128), tmax=2000)
    r = simulate(m; rng=Xoshiro(7), keep_sites=true)
    @printf("   status %s; masses %s; front reached column %d\n", r.status, r.mass, r.length)
    savefig("07_active.png", plot_sites(r.sites, "two fluids, 20% of the seeds active"))
end

# ---------------------------------------------------------------------------- IP

"8. Invasion percolation: the largest invaded threshold approaches p_c."
function ex_ip()
    for (name, perc, pc) in (("site", SitePercolation(), 0.5927), ("bond", BondPercolation(), 0.5))
        m = IPModel(SquareLattice(1024); percolation=perc, periodic=true, seeds=1024)
        t = @elapsed r = simulate(m; rng=Xoshiro(8), keep_times=true)
        @printf("   %s: %s after %d invasions in %.2f s; largest threshold %.4f (p_c = %.4f)\n",
            name, r.status, length(r.x), t, maximum(r.x), pc)
        name == "site" && savefig("08_ip.png", plot_times(r.times, "site IP, colored by invasion order"))
    end
end

"9. Invasion percolation with three fluids, each step moving one fluid picked at random."
function ex_ip_fluids()
    m = IPModel(SquareLattice(512); percolation=SitePercolation(), species=3, seeds=30,
        periodic=true, right=GrowingRight(64), nmax=100_000)
    r = simulate(m; rng=Xoshiro(9), keep_sites=true)
    @printf("   status %s; masses %s; invasions per fluid %s\n", r.status, r.mass,
        [count(==(UInt8(s)), r.species) for s in 1:3])
    savefig("09_ip_fluids.png", plot_sites(r.sites, "three fluids, site IP"))
end

"10. The same medium in both models: SOP with a small target follows IP closely."
function ex_shared_medium()
    medium = Medium(10)
    sop = SOPModel(SquareLattice(512); percolation=SitePercolation(), rejection=Persistent(),
        periodic=true, seeds=512, p0=0.65, c=0.05, f_T=0.02, right=FixedRight(512))
    ip = IPModel(SquareLattice(512); percolation=SitePercolation(), periodic=true, seeds=512)
    s = simulate(sop; rng=Xoshiro(1), medium, keep_sites=true)
    a = simulate(ip; rng=Xoshiro(1), medium, keep_sites=true)
    A, B = s.sites .!= 0, a.sites .!= 0
    @printf("   SOP mass %d, IP mass %d, overlap %.3f (see test/compare_SOP_IP.jl)\n",
        s.mass[1], a.mass[1], count(A .& B) / count(A .| B))
    # Fixing the medium while rng changes: only the competition between fluids differs.
    two = IPModel(SquareLattice(256); species=2, seeds=40, periodic=true, right=GrowingRight(64), nmax=20_000)
    x, y = (simulate(two; rng=Xoshiro(k), medium, keep_sites=true) for k in 1:2)
    @printf("   two fluids, same medium, different rng: %.3f of the sites hold the same fluid\n",
        count(x.sites .== y.sites) / length(x.sites))
end

# ---------------------------------------------------------------------------- 3D

"11. Three dimensions: SOP and IP on a simple cubic lattice."
function ex_cubic()
    m = SOPModel(CubicLattice(96); percolation=SitePercolation(), periodic=true, seeds=96^2 ÷ 5,
        p0=0.4, c=0.05, f_T=0.05, right=GrowingRight(128), tmax=1500)
    r = simulate(m; rng=Xoshiro(11))
    @printf("   SOP cubic site: %s, p* = %.4f (p_c = 0.3116)\n", r.status, stationary_mean(@view r.p[:, 1]))
    ip = IPModel(CubicLattice(128); percolation=SitePercolation(), periodic=true, seeds=128^2)
    t = @elapsed a = simulate(ip; rng=Xoshiro(11), keep_sites=true)
    @printf("   IP cubic site: %s after %d invasions in %.2f s, sites array %s\n", a.status,
        length(a.x), t, size(a.sites))
    savefig("11_cubic_ip.png", plot_sites(a.sites, "site IP on a 128³ cube, middle plane"))
end

"12. Large runs: keep only what you need (no history, no lattice)."
function ex_large()
    m = IPModel(SquareLattice(8192); percolation=SitePercolation(), periodic=true, seeds=8192)
    t = @elapsed r = simulate(m; rng=Xoshiro(12), keep_history=false)
    @printf("   IP on 8192² (67 M sites): %d invasions in %.1f s\n", sum(r.mass) - 8192, t)
end

const EXAMPLES = [ex_basic, ex_target, ex_rejection, ex_ensemble, ex_fluids, ex_layouts,
    ex_active, ex_ip, ex_ip_fluids, ex_shared_medium, ex_cubic, ex_large]

# The first line of a docstring, as written (rendering would turn p_c into italics).
function title(f)
    docstr = first(values(Base.Docs.meta(@__MODULE__)[Base.Docs.Binding(@__MODULE__, nameof(f))].docs))
    return strip(first(split(docstr.text[1], '\n')))
end

function main(args=ARGS)
    which = isempty(args) ? eachindex(EXAMPLES) : parse.(Int, args)
    for i in which
        println("\n", title(EXAMPLES[i]))
        EXAMPLES[i]()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
