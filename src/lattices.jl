"""Lattice geometry. A new lattice needs a `width` field or `transverse_sites`,
plus a `neighbors` method; the simulation kernel uses nothing else."""
abstract type Lattice end

"Square lattice with `width` sites along the transverse direction."
struct SquareLattice <: Lattice
    width::Int
end

"""
Triangular lattice stored on a rectangle. Columns of constant `x` have odd
columns shifted by half a spacing, so the injection wall stays straight.
"""
struct TriangularLattice <: Lattice
    width::Int
end

transverse_sites(l::Lattice) = l.width

@inline function _wrap(y, W, periodic)
    0 <= y < W && return y
    periodic ? mod(y, W) : -1
end

"""
    neighbors(lattice, i, nx, periodic)

Linear indices of the neighbors of site `i` on a lattice of length `nx`, as a
tuple. Neighbors outside the lattice are returned as 0.
"""
@inline function neighbors(l::SquareLattice, i::Int, nx::Int, periodic::Bool)
    W = l.width
    x, y = divrem(i - 1, W)
    xm = x > 0 ? i - W : 0
    xp = x + 1 < nx ? i + W : 0
    ym = _wrap(y - 1, W, periodic)
    yp = _wrap(y + 1, W, periodic)
    base = W * x + 1
    return (xp, ym < 0 ? 0 : base + ym, yp < 0 ? 0 : base + yp, xm)
end

@inline function neighbors(l::TriangularLattice, i::Int, nx::Int, periodic::Bool)
    W = l.width
    x, y = divrem(i - 1, W)
    base = W * x + 1
    ym = _wrap(y - 1, W, periodic)
    yp = _wrap(y + 1, W, periodic)
    # Second neighbor in the adjacent columns, set by the half-spacing shift.
    y2 = _wrap(isodd(x) ? y + 1 : y - 1, W, periodic)
    hasprev = x > 0
    hasnext = x + 1 < nx
    return (
        hasnext ? i + W : 0,
        hasnext && y2 >= 0 ? base + W + y2 : 0,
        ym < 0 ? 0 : base + ym,
        yp < 0 ? 0 : base + yp,
        hasprev ? i - W : 0,
        hasprev && y2 >= 0 ? base - W + y2 : 0,
    )
end
