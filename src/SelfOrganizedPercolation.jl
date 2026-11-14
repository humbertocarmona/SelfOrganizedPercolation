"""
Self-organized percolation (Alencar, Andrade & Lucena, PRE 56, R2379, 1997) and
invasion percolation, with several immiscible fluids injected from the left wall.

Coordinates: `x` is the injection direction (the last, growable dimension) and
`y` is transverse. Sites are stored column by column, so extending the lattice
to the right only appends memory.
"""
module SelfOrganizedPercolation

using Random

export SquareLattice, TriangularLattice
export SitePercolation, BondPercolation
export Permanent, Persistent
export FixedRight, GrowingRight
export PowerLawTarget
export RandomSeeds, BlockSeeds
export Medium
export SOPModel, SOPResult, IPModel, IPResult
export simulate, ensemble, stationary_mean, neighbors

include("lattices.jl")
include("model.jl")
include("medium.jl")
include("heap.jl")
include("simulate.jl")
include("invasion.jl")
include("ensemble.jl")

end
