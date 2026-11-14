# Persistent SOP and invasion percolation on the same quenched medium.
#
#     julia -t auto --project=reproduce test/compare_SOP_IP.jl            # table and figures
#     julia -t auto --project=reproduce test/compare_SOP_IP.jl --movie    # also a movie
#     julia -t auto --project=reproduce test/compare_SOP_IP.jl --movie --frames --every=5
#
# `--movie` writes an mp4 with one frame per SOP step (every `--every` steps), and
# `--frames` also saves each frame as a PNG. Output goes to test/output/. Including
# this file (instead of running it) only defines the functions, so the plotting can
# be used on other runs:
#
#     include("test/compare_SOP_IP.jl")
#     s, a = paired_runs(SitePercolation(), 0.05, 1)
#     save("pair.png", plot_pair(s, a))
#     movie(s, a, "pair.mp4"; frames_dir="frames", every=2)
#
# One fluid enters through the whole left wall of an L × L square lattice with
# periodic sides; both models stop when it first reaches the right wall. For each
# medium and each SOP target f_T we compare:
#   overlap   Jaccard overlap |A ∩ B| / |A ∪ B| of the invaded sites (seeds excluded),
#             after cutting the larger cluster, in its own order of arrival, to the
#             mass of the smaller, so both have about the same size;
#   baseline  the same overlap with IP run on a different medium, which is what
#             independent disorder gives;
#   order     Spearman correlation between SOP arrival time and IP invasion order
#             on the sites both invaded;
#   p*        SOP stationary p (second half of the run), against x_max, the largest
#             threshold IP had to invade to break through.
using SelfOrganizedPercolation, CairoMakie, Random, Statistics, Printf

const L = 512
const NMEDIA = 16
const PC = Dict(SitePercolation => 0.592746, BondPercolation => 0.5)

ranks(v) = invperm(sortperm(v))
spearman(a, b) = cor(ranks(a), ranks(b))
jaccard(a, b) = count(a .& b) / count(a .| b)

# Sites that arrived first, `k` of them (more for SOP, which fills many per step).
function first_arrivals(times, k)
    t = sort!(times[times .> 0])
    return 0 .< times .<= t[min(k, length(t))]
end
function overlap(t1, t2)
    k = min(count(>(0), t1), count(>(0), t2))
    return jaccard(first_arrivals(t1, k), first_arrivals(t2, k))
end

function models(perc, fT)
    sop = SOPModel(SquareLattice(L); percolation=perc, rejection=Persistent(), periodic=true,
        seeds=L, p0=PC[typeof(perc)] + 0.05, c=0.05, f_T=fT, right=FixedRight(L))
    ip = IPModel(SquareLattice(L); percolation=perc, periodic=true, seeds=L, right=FixedRight(L))
    return sop, ip
end

"SOP and IP on medium `i`, run to breakthrough with arrival times. Seeds fill the
whole wall, so they do not depend on `rng`; only the medium matters."
function paired_runs(perc, fT, i; medium=Medium(i))
    sop, ip = models(perc, fT)
    s = simulate(sop; rng=Xoshiro(i), medium, keep_times=true)
    a = simulate(ip; rng=Xoshiro(i), medium, keep_times=true)
    return s, a
end

function compare(perc, fT, i)
    s, a = paired_runs(perc, fT, i)
    _, b = paired_runs(perc, fT, i; medium=Medium(10_000 + i))
    s.status === a.status === b.status === :spanned || error("medium $i did not span")

    same, base = overlap(s.times, a.times), overlap(s.times, b.times)
    both = (s.times .> 0) .& (a.times .> 0)
    order = spearman(s.times[both], a.times[both])
    pstar = mean(@view s.p[size(s.p, 1)÷2:end, 1])
    return (; same, base, order, pstar, xmax=maximum(a.x), msop=sum(s.mass) - L, mip=length(a.x),
        steps=size(s.p, 1))
end

function table()
    for perc in (SitePercolation(), BondPercolation())
        name = perc isa SitePercolation ? "site" : "bond"
        @printf("\nsquare %s, L = %d, %d media, p_c = %.4f\n", name, L, NMEDIA, PC[typeof(perc)])
        @printf("  %-5s  %-15s  %-15s  %-15s  %-7s  %-7s  %-11s  %-11s  %s\n", "f_T", "overlap", "baseline",
            "order", "p*", "x_max", "SOP mass", "IP mass", "SOP steps")
        for fT in (0.4, 0.2, 0.1, 0.05, 0.02)
            rs = Vector{Any}(undef, NMEDIA)
            Threads.@threads for i in 1:NMEDIA
                rs[i] = compare(perc, fT, i)
            end
            se(v) = std(v) / sqrt(length(v))
            f(k) = [getproperty(r, k) for r in rs]
            @printf("  %.2f   %.3f ± %.3f    %.3f ± %.3f    %.3f ± %.3f    %.4f   %.4f   %.3e   %.3e   %.0f\n",
                fT, mean(f(:same)), se(f(:same)), mean(f(:base)), se(f(:base)), mean(f(:order)), se(f(:order)),
                mean(f(:pstar)), mean(f(:xmax)), mean(f(:msop)), mean(f(:mip)), mean(f(:steps)))
        end
    end
end

# ---------------------------------------------------------------------------- plots

const BLUE, ORANGE = colorant"#2a78d6", colorant"#eb6834"
const INK, MUTED, EMPTY_COLOR = colorant"#0b0b0b", colorant"#52514e", colorant"#ffffff"
const BOTH = colorant"#3a3a38"

# Arrival times as an x-by-y matrix, so x (injection) runs along the horizontal axis.
# A cubic lattice is cut at the transverse plane `z`.
xy(times::AbstractMatrix; z=1) = permutedims(times)
xy(times::AbstractArray{<:Any,3}; z=size(times, 2) ÷ 2 + 1) = permutedims(times[:, z, :])

"Arrival time below which `m` sites (seeds excluded) had been invaded."
function time_for_mass(times, m)
    t = sort!(times[times .> 0])
    return isempty(t) ? 0 : t[clamp(m, 1, length(t))]
end

"Sites invaded by time `t`, colored by arrival time (NaN, drawn white, for the rest)."
arrivals(times, t) = map(x -> 0 <= x <= t ? Float32(x) : NaN32, times)

"1 where both models invaded, 2 SOP only, 3 IP only, NaN neither."
function overlap_map(ts, ta, t_sop, t_ip)
    return map(ts, ta) do x, y
        a, b = 0 <= x <= t_sop, 0 <= y <= t_ip
        a && b ? 1.0f0 : a ? 2.0f0 : b ? 3.0f0 : NaN32
    end
end

"""
    plot_configuration(result; upto=nothing, z=nothing, title="", colorrange=nothing)

One run (`SOPResult` or `IPResult` kept with `keep_times=true`): the invaded sites
colored by arrival time, up to arrival time `upto` (all by default). A cubic lattice
is shown as the plane `z`.
"""
function plot_configuration(r; upto=nothing, z=nothing, title="", colorrange=nothing)
    times = z === nothing ? xy(r.times) : xy(r.times; z)
    T = upto === nothing ? maximum(times) : upto
    f = Figure(size=(620, 560), backgroundcolor=:white)
    ax = Axis(f[1, 1]; aspect=DataAspect(), title, xlabel="x, injection direction", ylabel="y")
    hm = heatmap!(ax, arrivals(times, T); colormap=:viridis, nan_color=EMPTY_COLOR,
        colorrange=something(colorrange, (0, max(T, 1))), interpolate=false)
    Colorbar(f[1, 2], hm; label="arrival time")
    return f
end

"""
    plot_pair(sop, ip; step=nothing)

SOP at `step` (its final step by default) beside IP at the same invaded mass, and a
third panel with the sites they share. Both results need `keep_times=true`.
"""
function plot_pair(s, a; step=nothing)
    fig, update! = pair_figure(s, a)
    update!(something(step, size(s.p, 1)))
    return fig
end

# Builds the three panels once and returns the figure with a function that redraws
# them for a given SOP step, so a movie does not rebuild the figure for each frame.
function pair_figure(s, a)
    ts, ta = xy(s.times), xy(a.times)
    nsteps, nip = size(s.p, 1), length(a.x)
    sop_t = Observable(arrivals(ts, 0))
    ip_t = Observable(arrivals(ta, 0))
    both = Observable(overlap_map(ts, ta, 0, 0))
    header = Observable("")

    fig = Figure(size=(1500, 560), backgroundcolor=:white)
    Label(fig[0, 1:3], header; tellwidth=false, halign=:left, color=INK, fontsize=16)
    kw = (; aspect=DataAspect(), xlabel="x, injection direction", titlealign=:left, titlefont=:regular)
    ax1 = Axis(fig[1, 1]; title="Persistent SOP (color: step)", ylabel="y", kw...)
    ax2 = Axis(fig[1, 2]; title="IP at the same mass (color: invasion order)", kw...)
    ax3 = Axis(fig[1, 3]; title="shared sites", kw...)
    heatmap!(ax1, sop_t; colormap=:viridis, nan_color=EMPTY_COLOR, colorrange=(0, nsteps), interpolate=false)
    heatmap!(ax2, ip_t; colormap=:viridis, nan_color=EMPTY_COLOR, colorrange=(0, max(nip, 1)), interpolate=false)
    heatmap!(ax3, both; colormap=[BOTH, BLUE, ORANGE], nan_color=EMPTY_COLOR, colorrange=(1, 3),
        interpolate=false)
    Legend(fig[2, 3], [PolyElement(color=c) for c in (BOTH, BLUE, ORANGE)], ["both", "SOP only", "IP only"];
        orientation=:horizontal, framevisible=false)
    hideydecorations!(ax2; label=false); hideydecorations!(ax3; label=false)
    linkaxes!(ax1, ax2, ax3)

    function update!(step)
        mass = count(x -> 0 < x <= step, s.times)
        tip = time_for_mass(a.times, mass)          # IP invasion index with that mass
        sop_t[] = arrivals(ts, step)
        ip_t[] = arrivals(ta, tip)
        both[] = overlap_map(ts, ta, step, tip)
        shared = count(==(1.0f0), both[])
        total = count(!isnan, both[])
        p = s.p[min(step, nsteps), 1]
        header[] = @sprintf("SOP step %d of %d,  p = %.4f    invaded mass %d (IP: %d of %d)    overlap %.3f",
            step, nsteps, p, mass, min(mass, nip), nip, shared / max(total, 1))
        return fig
    end
    return fig, update!
end

"""
    movie(sop, ip, path; every=1, framerate=30, frames_dir=nothing)

Writes an mp4 of [`plot_pair`](@ref) with one frame every `every` SOP steps, always
ending at the final step. With `frames_dir`, each frame is also saved as a PNG,
`frame_00001.png` and so on, ready for ffmpeg or any other movie tool.
"""
function movie(s, a, path; every=1, framerate=30, frames_dir=nothing)
    fig, update! = pair_figure(s, a)
    steps = unique([1:every:size(s.p, 1); size(s.p, 1)])
    frames_dir === nothing || mkpath(frames_dir)
    record(fig, path, enumerate(steps); framerate) do (k, step)
        update!(step)
        frames_dir === nothing || save(joinpath(frames_dir, @sprintf("frame_%05d.png", k)), fig)
    end
    return path
end

# ---------------------------------------------------------------------------- run

function main(args=ARGS)
    flag(name) = any(==("--" * name), args)
    i = findfirst(startswith("--every="), args)
    every = i === nothing ? 1 : parse(Int, args[i][length("--every=")+1:end])
    out = joinpath(@__DIR__, "output")
    mkpath(out)

    table()
    # Figures for one medium: far from and close to IP.
    for perc in (SitePercolation(), BondPercolation()), fT in (0.4, 0.02)
        name = perc isa SitePercolation ? "site" : "bond"
        s, a = paired_runs(perc, fT, 1)
        file = joinpath(out, @sprintf("pair_%s_fT%.2f.png", name, fT))
        save(file, plot_pair(s, a); px_per_unit=1)
        println("wrote ", file)
    end
    if flag("movie")
        s, a = paired_runs(SitePercolation(), 0.05, 1)
        file = joinpath(out, "pair_site_fT0.05.mp4")
        frames = flag("frames") ? joinpath(out, "frames_site_fT0.05") : nothing
        movie(s, a, file; every, frames_dir=frames)
        println("wrote ", file, frames === nothing ? "" : " and the frames in " * frames)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
