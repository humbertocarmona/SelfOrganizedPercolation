"""
    IPModel(lattice; kwargs...)

Invasion percolation without trapping (Wilkinson & Willemsen 1983), with one or more
immiscible fluids injected from the left wall. Each step picks one fluid uniformly at
random among those that can still advance, and that fluid invades the element with the
smallest threshold on its perimeter: a site for `SitePercolation`, a bond and the site
behind it for `BondPercolation`. Thresholds come from a [`Medium`](@ref), the same one
the `Persistent` rule of [`SOPModel`](@ref) reads.

The perimeter of each fluid is a binary heap, as in Masson & Pride, Transp. Porous Med.
102, 301 (2014), so `M` invasions take `O(M log M)` time. Thresholds are hashed rather
than stored, so the lattice costs one byte per site (two for site percolation) and the
heap 16 bytes per perimeter element.

# Keywords
- `percolation`: `SitePercolation()` or `BondPercolation()` (default).
- `right`: `FixedRight(L)` (default, `L` = lattice width, `ny` for
  `CubicLattice`) stops when a fluid reaches the
  last column; `GrowingRight(L0)` moves the wall ahead of the front.
- `periodic`, `species`, `seeds`, `layout`, `active_fraction`: as in `SOPModel`.
- `nmax`: maximum number of invasions.
"""
struct IPModel{L<:Lattice,P<:Percolation,B,S}
    lattice::L
    percolation::P
    right::B
    periodic::Bool
    species::Int
    seeds::Int
    active_fraction::Float64
    nmax::Int
    layout::S
end

function IPModel(lattice::Lattice; percolation=BondPercolation(),
        right=FixedRight(default_length(lattice)), periodic=false, species=1, seeds=100,
        active_fraction=1.0, nmax=typemax(Int), layout=RandomSeeds())
    check_injection(lattice, percolation, right, species, seeds, layout, active_fraction)
    nmax >= 0 || throw(ArgumentError("nmax must be non-negative"))
    return IPModel(lattice, percolation, right, periodic, species, seeds,
        Float64(active_fraction), Int(nmax), layout)
end

"""
Outcome of one invasion percolation run. `x[n]` is the threshold of the `n`-th
invaded element and `species[n]` the fluid that invaded it; both are empty unless
`keep_history=true`.

`status` is `:spanned` (a fluid reached the last column of a `FixedRight` lattice),
`:extinct` (no fluid has a perimeter left), `:nmax` or `:maxlength`. `sites` is the
lattice with `keep_sites=true`, of size `(W, L)` or `(ny, nz, L)`. `times` with `keep_times=true` gives the invasion
order: `n` for the site invaded at step `n`, `0` for a seed and `-1` for a site
never invaded.
"""
struct IPResult
    x::Vector{Float64}
    species::Vector{UInt8}
    status::Symbol
    mass::Vector{Int}
    length::Int
    sites::Union{Nothing,Array{UInt8}}
    times::Union{Nothing,Array{Int32}}
end

mutable struct IPState
    sites::Vector{UInt8}
    tested::Vector{UInt8}          # site percolation only: one bit per species
    times::Vector{Int32}           # invasion order; empty unless requested
    nx::Int
    xmax::Int                      # rightmost occupied column, 0-based
    heaps::Vector{Vector{Threshold}}
    mass::Vector{Int}
    medium::Medium
end

# Puts the empty neighbors of site `i` on the perimeter of species `s`. In site
# percolation a site enters each perimeter once; in bond percolation each bond does,
# so a site can be in a heap several times.
@inline function add_perimeter!(st::IPState, m::IPModel, s::Int, i::Int)
    heap = st.heaps[s]
    @inbounds for j in neighbors(m.lattice, i, st.nx, m.periodic)
        (j == 0 || st.sites[j] != EMPTY) && continue
        first_contact!(st, m.percolation, j, s) || continue
        heap_push!(heap, (threshold(st.medium, m.percolation, i, j), j))
    end
    return nothing
end

function initial_state(m::IPModel, rng, medium::Medium, keep_times::Bool)
    W = transverse_sites(m.lattice)
    nx = initial_length(m.right)
    st = IPState(zeros(UInt8, W * nx),
        m.percolation isa SitePercolation ? zeros(UInt8, W * nx) : UInt8[],
        keep_times ? fill(Int32(-1), W * nx) : Int32[],
        nx, 0, [Threshold[] for _ in 1:m.species], fill(m.seeds, m.species), medium)
    wall = wall_order(m.layout, m, W, rng)
    for s in 1:m.species, n in 1:m.seeds
        site = wall[(s-1)*m.seeds+n]
        st.sites[site] = UInt8(s)
        keep_times && (st.times[site] = 0)
    end
    # Perimeters only after every seed is placed, so no seed enters a heap.
    nactive = round(Int, m.active_fraction * m.seeds)
    for s in 1:m.species, n in 1:nactive
        add_perimeter!(st, m, s, wall[(s-1)*m.seeds+n])
    end
    return st
end

# Smallest perimeter element of species `s` whose site is still empty, or `nothing`.
# Elements whose site another fluid (or another bond) filled are dropped here.
@inline function next_element!(heap::Vector{Threshold}, sites)
    @inbounds while !isempty(heap)
        x, j = heap[1]
        heap_pop!(heap)
        sites[j] == EMPTY && return (x, j)
    end
    return nothing
end

# Makes room for the right neighbor of column `col` before a site there is invaded.
room_right!(st, ::FixedRight, W, col) = true
function room_right!(st, r::GrowingRight, W, col)
    col + 2 <= st.nx && return true
    st.nx >= r.maxlength && return false
    extend!(st, W, min(2st.nx, r.maxlength))
    return true
end

spans(::FixedRight, st, col) = col == st.nx - 1
spans(::GrowingRight, st, col) = false

function invade!(st::IPState, m::IPModel, rng, xs, who)
    W = transverse_sites(m.lattice)
    active = [s for s in 1:m.species if !isempty(st.heaps[s])]
    keep = !isempty(st.times)
    n = 0
    while n < m.nmax
        isempty(active) && return :extinct
        k = length(active) == 1 ? 1 : rand(rng, 1:length(active))
        s = active[k]
        next = next_element!(st.heaps[s], st.sites)
        if next === nothing
            deleteat!(active, k)
            continue
        end
        x, j = next
        col = (j - 1) ÷ W
        if !room_right!(st, m.right, W, col)
            heap_push!(st.heaps[s], next)   # leave the lattice as it was
            return :maxlength
        end
        n += 1
        @inbounds st.sites[j] = UInt8(s)
        if keep
            n <= typemax(Int32) || throw(OverflowError("invasion order exceeds Int32; use keep_times=false"))
            @inbounds st.times[j] = Int32(n)
        end
        st.mass[s] += 1
        st.xmax = max(st.xmax, col)
        xs === nothing || (push!(xs, x); push!(who, UInt8(s)))
        spans(m.right, st, col) && return :spanned
        add_perimeter!(st, m, s, j)
    end
    return :nmax
end

"""
    simulate(model::IPModel; rng=Xoshiro(), medium=Medium(rand(rng, UInt64)),
             keep_history=true, keep_sites=false, keep_times=false)

Runs invasion percolation and returns an [`IPResult`](@ref). As for `SOPModel`, the
default medium is the first draw from `rng`, and the seeds are placed with the next
draws, so an `IPModel` and an `SOPModel` run with the same seed share the medium and,
for equal `seeds` and `layout`, the seed positions. `rng` also picks the invading
fluid when there are several. Turn off `keep_history` on very large runs: it costs 9
bytes per invasion.
"""
function simulate(m::IPModel; rng=Xoshiro(), medium::Medium=Medium(rand(rng, UInt64)),
        keep_history::Bool=true, keep_sites::Bool=false, keep_times::Bool=false)
    st = initial_state(m, rng, medium, keep_times)
    xs, who = Float64[], UInt8[]
    status = keep_history ? invade!(st, m, rng, xs, who) : invade!(st, m, rng, nothing, nothing)
    shape = lattice_shape(m.lattice, st.nx)
    sites = keep_sites ? reshape(st.sites, shape) : nothing
    times = keep_times ? reshape(st.times, shape) : nothing
    return IPResult(xs, who, status, st.mass, st.xmax + 1, sites, times)
end
