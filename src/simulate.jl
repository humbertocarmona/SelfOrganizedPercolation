const EMPTY = 0x00

"""
Outcome of one run. `p[t+1, s]` is the probability used by species `s` at time
`t`, and `f[t+1, s]` is the number of sites it occupied in that step per unit wall
area, `N/L^(d-1)`.

`status` is `:extinct`, `:spanned`, `:tmax` or `:maxlength`. `sites` is the
lattice (`0` empty, `s` for species `s`) with `keep_sites=true`. `times` is the
arrival time of each site with `keep_times=true`: `t + 1` for a site filled at
time `t`, `0` for a seed, and `-1` for a site never filled.
"""
struct SOPResult
    p::Matrix{Float64}
    f::Matrix{Float64}
    status::Symbol
    mass::Vector{Int}
    length::Int
    sites::Union{Nothing,Matrix{UInt8}}
    times::Union{Nothing,Matrix{Int32}}
end

mutable struct State
    sites::Vector{UInt8}
    tested::Vector{UInt8}          # site percolation only: one bit per species
    nx::Int
    fronts::Vector{Vector{Int}}
    next::Vector{Vector{Int}}
    candidates::Vector{Tuple{Int,Int}}
    heaps::Vector{Vector{Threshold}}
    p::Vector{Float64}
    mass::Vector{Int}
    reached::Vector{Bool}
    xmax::Int                      # rightmost occupied column, 0-based
    times::Vector{Int32}           # arrival times; empty unless they are requested
    step::Int32                    # value recorded for the sites filled in this step
    medium::Medium                 # thresholds of the Persistent rule
end

wall_order(::RandomSeeds, m, W, rng) = randperm(rng, W)

function wall_order(l::BlockSeeds, m, W, rng)
    width = m.seeds ÷ l.blocks
    cell = W / l.blocks
    return [round(Int, (i - 0.5) * cell - width / 2) + j for i in 1:l.blocks for j in 1:width]
end

function initial_state(m::SOPModel, rng; keep_times::Bool=false, medium::Medium=Medium(rand(rng, UInt64)))
    W = transverse_sites(m.lattice)
    nx = initial_length(m.right)
    sites = zeros(UInt8, W * nx)
    times = keep_times ? fill(Int32(-1), W * nx) : Int32[]
    tested = m.percolation isa SitePercolation ? zeros(UInt8, W * nx) : UInt8[]
    fronts = [Int[] for _ in 1:m.species]
    wall = wall_order(m.layout, m, W, rng)
    nactive = round(Int, m.active_fraction * m.seeds)
    for s in 1:m.species, n in 1:m.seeds
        site = wall[(s-1)*m.seeds+n]
        sites[site] = UInt8(s)
        keep_times && (times[site] = 0)
        n <= nactive && push!(fronts[s], site)
    end
    return State(sites, tested, nx, fronts, [Int[] for _ in 1:m.species], Tuple{Int,Int}[],
        [Threshold[] for _ in 1:m.species], copy(m.p0), fill(m.seeds, m.species),
        fill(false, m.species), 0, times, Int32(0), medium)
end

function extend!(st, W::Int, nx::Int)
    old = length(st.sites)
    new = W * nx
    resize!(st.sites, new)
    fill!(view(st.sites, old+1:new), EMPTY)
    if !isempty(st.tested)
        resize!(st.tested, new)
        fill!(view(st.tested, old+1:new), 0x00)
    end
    if !isempty(st.times)
        resize!(st.times, new)
        fill!(view(st.times, old+1:new), Int32(-1))
    end
    st.nx = nx
    return st
end

# Marks site j as seen by species s. False if that species has already seen it.
@inline first_contact!(st, ::BondPercolation, j, s) = true
@inline function first_contact!(st, ::SitePercolation, j, s)
    bit = 0x01 << (s - 1)
    @inbounds st.tested[j] & bit != 0x00 && return false
    @inbounds st.tested[j] |= bit
    return true
end

function propose!(st::State, m::SOPModel{<:Any,<:Any,Permanent}, s::Int, rng)
    p = st.p[s]
    @inbounds for i in st.fronts[s]
        for j in neighbors(m.lattice, i, st.nx, m.periodic)
            (j == 0 || st.sites[j] != EMPTY) && continue
            first_contact!(st, m.percolation, j, s) || continue
            rand(rng) < p && push!(st.candidates, (j, s))
        end
    end
    return nothing
end

function propose!(st::State, m::SOPModel{<:Any,<:Any,Persistent}, s::Int, rng)
    p = st.p[s]
    heap = st.heaps[s]
    @inbounds for i in st.fronts[s]
        for j in neighbors(m.lattice, i, st.nx, m.periodic)
            (j == 0 || st.sites[j] != EMPTY) && continue
            first_contact!(st, m.percolation, j, s) || continue
            heap_push!(heap, (threshold(st.medium, m.percolation, i, j), j))
        end
    end
    @inbounds while !isempty(heap) && heap[1][1] < p
        _, j = heap_pop!(heap)
        st.sites[j] == EMPTY && push!(st.candidates, (j, s))
    end
    return nothing
end

"""
Gives each contested site to one candidate chosen uniformly at random: after a
shuffle, the first candidate for a site claims it. A site with two candidates
from one fluid and one from another therefore goes to the first fluid with
probability 2/3.
"""
function resolve!(st::State, W::Int, rng)
    length(st.next) > 1 && shuffle!(rng, st.candidates)
    @inbounds for (j, s) in st.candidates
        st.sites[j] == EMPTY || continue
        st.sites[j] = UInt8(s)
        isempty(st.times) || (st.times[j] = st.step)
        push!(st.next[s], j)
        st.xmax = max(st.xmax, (j - 1) ÷ W)
    end
    empty!(st.candidates)
    return nothing
end

# Before each step the lattice must have a column to the right of the front.
prepare_right!(st, ::FixedRight, W) = :running
function prepare_right!(st, r::GrowingRight, W)
    st.xmax + 2 <= st.nx && return :running
    st.nx >= r.maxlength && return :maxlength
    extend!(st, W, min(2st.nx, r.maxlength))
    return :running
end

can_advance(st, m::SOPModel{<:Any,<:Any,Permanent}, s) = !isempty(st.fronts[s])
can_advance(st, m::SOPModel{<:Any,<:Any,Persistent}, s) =
    !isempty(st.fronts[s]) || !isempty(st.heaps[s])

function finished(st, m::SOPModel, ::FixedRight)
    all(s -> st.reached[s] || !can_advance(st, m, s), 1:m.species) || return :running
    return any(st.reached) ? :spanned : :extinct
end
finished(st, m::SOPModel, ::GrowingRight) =
    any(s -> can_advance(st, m, s), 1:m.species) ? :running : :extinct

function evolve!(st::State, m::SOPModel, rng, phist, Nhist)
    W = transverse_sites(m.lattice)
    ns = m.species
    t = 0
    status = :tmax
    while t < m.tmax
        status = prepare_right!(st, m.right, W)
        status === :running || break
        st.step = Int32(t + 1)
        for s in 1:ns
            push!(phist, st.p[s])
            propose!(st, m, s, rng)
        end
        resolve!(st, W, rng)
        target = m.f_T(t)
        for s in 1:ns
            front = st.next[s]
            N = length(front)
            push!(Nhist, N)
            st.mass[s] += N
            st.p[s] = clamp(st.p[s] + m.c * (target - N / W), 0.0, 1.0)
            if !st.reached[s]
                @inbounds for j in front
                    (j - 1) ÷ W == st.nx - 1 && (st.reached[s] = true; break)
                end
            end
            st.fronts[s], st.next[s] = front, st.fronts[s]
            empty!(st.next[s])
        end
        t += 1
        status = finished(st, m, m.right)
        status === :running || break
        status = :tmax
    end
    return status
end

"""
    simulate(model; rng=Xoshiro(), medium=Medium(rand(rng, UInt64)), keep_sites=false, keep_times=false)

Runs one realization and returns a [`SOPResult`](@ref).

The `Persistent` rule reads its thresholds from `medium`. By default the medium is
the first draw from `rng`, so a seed fixes the medium, and two models run with the
same seed share it. Pass `medium` to keep the disorder fixed while `rng` varies the
seed positions and the choice among competing fluids.
"""
function simulate(m::SOPModel; rng=Xoshiro(), medium::Medium=Medium(rand(rng, UInt64)),
        keep_sites::Bool=false, keep_times::Bool=false)
    st = initial_state(m, rng; keep_times, medium)
    phist = Float64[]
    Nhist = Int[]
    status = evolve!(st, m, rng, phist, Nhist)
    ns = m.species
    p = permutedims(reshape(phist, ns, :))
    W = transverse_sites(m.lattice)
    f = permutedims(reshape(Nhist ./ W, ns, :))
    sites = keep_sites ? reshape(st.sites, W, st.nx) : nothing
    times = keep_times ? reshape(st.times, W, st.nx) : nothing
    return SOPResult(p, f, status, st.mass, st.xmax + 1, sites, times)
end
