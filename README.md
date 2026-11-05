# SelfOrganizedPercolation.jl

Self-organized percolation (Alencar, Andrade & Lucena, *Phys. Rev. E* 56, R2379, 1997)
with one or more immiscible fluids injected from the left wall.

```julia
using SelfOrganizedPercolation, Random

model = SOPModel(SquareLattice(2000);
    percolation = SitePercolation(), seeds = 100, p0 = 1.0, c = 0.2, f_T = 0.1, tmax = 1000)

result = simulate(model; rng = Xoshiro(1))
result.p        # p(t), one column per fluid
result.f        # f(t) = N(t)/L^(d-1), one column per fluid
result.status   # :extinct, :spanned, :tmax or :maxlength

means = ensemble(model, 64; f = r -> stationary_mean(@view r.p[:, 1]))   # start with julia -t auto
```

## Options

| Keyword | Values |
|---|---|
| lattice | `SquareLattice(W)`, `TriangularLattice(W)` |
| `percolation` | `BondPercolation()`, `SitePercolation()` |
| `rejection` | `Permanent()` (original SOP), `Persistent()` (a rejected throat can open later) |
| `right` | `FixedRight(L)`, `GrowingRight(L0; maxlength)` (the wall moves ahead of the front) |
| `periodic` | periodic transverse boundaries; `false` matches the 1997 paper |
| `species` | number of immiscible fluids: up to 254 for bond, 8 for site |
| `f_T` | a number, `PowerLawTarget(a, α)`, or any function of `t` |

Active sites are measured in units of the wall surface area, so the feedback rule is
`p(t+1) = p(t) + c[f_T - f(t)]` with `f(t) = N(t)/L^(d-1)`. The parameters of the 1997
paper convert as `c = k L^(d-1)` and `f_T = N_T / L^(d-1)`.

## Contested sites

Fluids never share a site. When several successful bonds reach an empty site in one
step, one bond is chosen uniformly at random and its fluid takes the site.

## Adding a lattice

A new lattice needs one type and one `neighbors` method. The 3D cubic lattice goes
the same way, with `transverse_sites` returning `Wy * Wz`; `x` stays the last index.

## Reproducing the paper

```sh
julia -t auto --project=reproduce -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia -t auto --project=reproduce reproduce/paper1997.jl
julia -t auto --project=reproduce reproduce/validate_extensions.jl
julia -t auto --project=reproduce reproduce/benchmark.jl
```

Figures and `stationary_p.csv` are written to `reproduce/output/`.
