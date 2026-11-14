"""
    Medium(seed)

Quenched disorder: a threshold `x` uniform in `[0, 1)` for every site and every bond
of an unbounded lattice. Thresholds are a hash of `seed` and the element, not draws
from a stream, so an element has the same threshold in any model, at any time, and
however the lattice grows. Runs given the same `Medium` see the same porous medium.

A bond takes the same threshold from both ends. Site and bond thresholds are
unrelated, and all fluids share one medium.
"""
struct Medium
    seed::UInt64
end
Medium(seed::Integer) = Medium(seed % UInt64)

# SplitMix64 finalizer (Steele, Lea & Flood 2014): a bijection of UInt64 that turns
# consecutive inputs into independent-looking outputs.
@inline function splitmix(z::UInt64)
    z += 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

# The top 53 bits as a Float64 in [0, 1).
@inline unit(z::UInt64) = Float64(z >> 11) * 0x1.0p-53

"""
    threshold(medium, j)
    threshold(medium, i, j)

Threshold of site `j`, or of the bond between sites `i` and `j`, from their
linear indices. Linear indices do not change when the lattice grows to the right.
"""
@inline threshold(m::Medium, j::Integer) = unit(splitmix(splitmix(m.seed) ⊻ (j % UInt64)))
@inline function threshold(m::Medium, i::Integer, j::Integer)
    a, b = minmax(i, j)
    return unit(splitmix(splitmix(splitmix(m.seed) ⊻ (a % UInt64)) ⊻ (b % UInt64)))
end

# The element a front site `i` tests when it reaches the empty site `j`.
@inline threshold(m::Medium, ::SitePercolation, i, j) = threshold(m, j)
@inline threshold(m::Medium, ::BondPercolation, i, j) = threshold(m, i, j)
