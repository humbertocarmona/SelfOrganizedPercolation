"""Lattice geometry. A new lattice needs a `width` field or `transverse_sites`,
plus a `neighbors` method; the simulation kernel uses nothing else. A lattice whose
wall is not a line also defines `default_length` and `lattice_shape`."""
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

"""
    CubicLattice(ny, nz)
    CubicLattice(W)

Simple cubic lattice with an `ny × nz` injection wall (`W × W` for one argument).
Sites are stored with `y` fastest, then `z`, then `x`, so a column of the
two-dimensional lattices becomes a wall layer here. `periodic` wraps both `y` and `z`.
"""
struct CubicLattice <: Lattice
    ny::Int
    nz::Int
end
CubicLattice(W::Integer) = CubicLattice(W, W)

transverse_sites(l::Lattice) = l.width
transverse_sites(l::CubicLattice) = l.ny * l.nz

"Default length of a `FixedRight` lattice: as long as it is wide."
default_length(l::Lattice) = l.width
default_length(l::CubicLattice) = l.ny

"Shape of the arrays returned with `keep_sites` and `keep_times`, `x` last."
lattice_shape(l::Lattice, nx) = (transverse_sites(l), nx)
lattice_shape(l::CubicLattice, nx) = (l.ny, l.nz, nx)

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

@inline function neighbors(l::CubicLattice, i::Int, nx::Int, periodic::Bool)
    ny, nz = l.ny, l.nz
    A = ny * nz
    x, r = divrem(i - 1, A)
    z, y = divrem(r, ny)
    base = A * x + 1
    ym = _wrap(y - 1, ny, periodic)
    yp = _wrap(y + 1, ny, periodic)
    zm = _wrap(z - 1, nz, periodic)
    zp = _wrap(z + 1, nz, periodic)
    return (
        x + 1 < nx ? i + A : 0,
        ym < 0 ? 0 : base + ny * z + ym,
        yp < 0 ? 0 : base + ny * z + yp,
        zm < 0 ? 0 : base + ny * zm + y,
        zp < 0 ? 0 : base + ny * zp + y,
        x > 0 ? i - A : 0,
    )
end
