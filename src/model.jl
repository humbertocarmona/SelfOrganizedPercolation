abstract type Percolation end
"Each site is tested at most once per species, as in the site case of the 1997 paper."
struct SitePercolation <: Percolation end
"Each bond from an active site is tested independently."
struct BondPercolation <: Percolation end

abstract type Rejection end
"A failed test is final. This is the original SOP rule."
struct Permanent <: Rejection end
"""
A rejected element keeps a fixed threshold `x` and is occupied whenever
`x < p(t)`. This is the drainage rule: a throat can open later if pressure rises.
"""
struct Persistent <: Rejection end

"The lattice has a fixed `length`; a species that reaches it has spanned."
struct FixedRight
    length::Int
end

"""
    GrowingRight(initial; maxlength=typemax(Int))

The right boundary moves ahead of the front, so the front never reaches a wall.
"""
struct GrowingRight
    initial::Int
    maxlength::Int
end
GrowingRight(initial::Integer; maxlength::Integer=typemax(Int)) = GrowingRight(initial, maxlength)

"Seeds at random positions on the left wall."
struct RandomSeeds end

"""
    BlockSeeds(blocks)

Seeds in `blocks` equal blocks, evenly spaced along the left wall. The blocks are
centered in equal cells, so none wraps across a periodic boundary. With as many seeds
as half the wall, blocks and gaps alternate like a column of a chessboard.
"""
struct BlockSeeds
    blocks::Int
end

initial_length(r::FixedRight) = r.length
initial_length(r::GrowingRight) = r.initial

"Time-dependent target `f_T(t) = a t^α`, the reduced form of Eq. (4) of the 1997 paper."
struct PowerLawTarget
    a::Float64
    α::Float64
end
(f::PowerLawTarget)(t) = f.a * Float64(t)^f.α

"Constant target fraction of active sites."
struct ConstantTarget
    value::Float64
end
(f::ConstantTarget)(t) = f.value

_target(x::Real) = ConstantTarget(Float64(x))
_target(f) = f

# Checks shared by every model that injects fluids from the left wall.
function check_injection(lattice, percolation, right, species, seeds, layout, active_fraction)
    W = transverse_sites(lattice)
    W >= 3 || throw(ArgumentError("the transverse size must be at least 3"))
    1 <= species <= 254 || throw(ArgumentError("species must be between 1 and 254"))
    percolation isa SitePercolation && species > 8 &&
        throw(ArgumentError("site percolation supports at most 8 species"))
    species * seeds <= W || throw(ArgumentError("the left wall has too few sites for the seeds"))
    if layout isa BlockSeeds
        species == 1 || throw(ArgumentError("BlockSeeds supports one species"))
        seeds % layout.blocks == 0 || throw(ArgumentError("seeds must divide evenly into the blocks"))
    end
    0 <= active_fraction <= 1 || throw(ArgumentError("active_fraction must be in [0, 1]"))
    initial_length(right) >= 2 || throw(ArgumentError("the initial length must be at least 2"))
    return nothing
end

"""
    SOPModel(lattice; kwargs...)

Active sites are measured in units of the wall surface area `A = L^(d-1)`, the value of
`transverse_sites(lattice)`: `f(t) = N(t)/A`. The parameters of the 1997 paper
convert as `c = k A` and `f_T = N_T / A`.

# Keywords
- `percolation`: `SitePercolation()` or `BondPercolation()` (default).
- `rejection`: `Permanent()` (default) or `Persistent()`.
- `right`: `FixedRight(L)` (default, `L` = lattice width) or `GrowingRight(L0)`.
- `periodic`: periodic transverse boundaries (default `false`, as in the paper).
- `species`: number of immiscible fluids.
- `seeds`: inhabited sites per species on the left wall.
- `layout`: `RandomSeeds()` (default) or `BlockSeeds(blocks)`.
- `active_fraction`: fraction of those seeds that start active.
- `p0`: initial probability; a number or one value per species.
- `c`: feedback strength in `p(t+1) = p(t) + c[f_T - f(t)]`.
- `f_T`: target fraction of active sites, a number or a callable of `t`.
- `tmax`: maximum number of time steps.
"""
struct SOPModel{L<:Lattice,P<:Percolation,R<:Rejection,B,T,S}
    lattice::L
    percolation::P
    rejection::R
    right::B
    periodic::Bool
    species::Int
    seeds::Int
    active_fraction::Float64
    p0::Vector{Float64}
    c::Float64
    f_T::T
    tmax::Int
    layout::S
end

function SOPModel(lattice::Lattice; percolation=BondPercolation(), rejection=Permanent(),
        right=FixedRight(transverse_sites(lattice)), periodic=false, species=1, seeds=100,
        active_fraction=1.0, p0=1.0, c=0.05, f_T=0.1, tmax=typemax(Int), layout=RandomSeeds())
    check_injection(lattice, percolation, right, species, seeds, layout, active_fraction)
    probabilities = p0 isa Real ? fill(Float64(p0), species) : Float64.(collect(p0))
    length(probabilities) == species || throw(ArgumentError("p0 needs one value per species"))
    all(p -> 0 <= p <= 1, probabilities) || throw(ArgumentError("p0 must be in [0, 1]"))
    return SOPModel(lattice, percolation, rejection, right, periodic, species, seeds,
        Float64(active_fraction), probabilities, Float64(c), _target(f_T), tmax, layout)
end
