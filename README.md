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
| lattice | `SquareLattice(W)`, `TriangularLattice(W)`, `CubicLattice(ny, nz)` |
| `percolation` | `BondPercolation()`, `SitePercolation()` |
| `rejection` | `Permanent()` (original SOP), `Persistent()` (a rejected throat can open later) |
| `right` | `FixedRight(L)`, `GrowingRight(L0; maxlength)` (the wall moves ahead of the front) |
| `periodic` | periodic transverse boundaries; `false` matches the 1997 paper |
| `species` | number of immiscible fluids: up to 254 for bond, 8 for site |
| `trapping` | `false` (default) or `true`: defender cut off from the outlet is trapped |
| `f_T` | a number, `PowerLawTarget(a, α)`, or any function of `t` |

Active sites are measured in units of the wall surface area, so the feedback rule is
`p(t+1) = p(t) + c[f_T - f(t)]` with `f(t) = N(t)/L^(d-1)`. The parameters of the 1997
paper convert as `c = k L^(d-1)` and `f_T = N_T / L^(d-1)`.

## Contested sites

Fluids never share a site. When several successful bonds reach an empty site in one
step, one bond is chosen uniformly at random and its fluid takes the site.

## Quenched medium

Under `Persistent`, every site and bond has a fixed threshold `x` uniform in `[0, 1)`,
and it opens once `x < p(t)`. Thresholds come from a `Medium`, a hash of a seed and the
element's index, so they do not depend on the order of invasion or on how far the
lattice grows. All fluids share one medium, and a bond has one threshold.

```julia
medium = Medium(42)
a = simulate(model_a; rng = Xoshiro(1), medium)   # same porous medium,
b = simulate(model_b; rng = Xoshiro(2), medium)   # different seeds and competition
```

By default the medium is drawn from `rng`, so runs with the same seed (including
realization `i` of `ensemble`) share it. `Permanent` gives each element one test at
first contact and still draws it from `rng`.

## Invasion percolation

`IPModel` takes the same lattice, percolation type, boundaries and seeds as `SOPModel`.
Each step it picks a fluid uniformly at random and invades the smallest threshold on
that fluid's perimeter. There is no trapping unless `trapping = true`. The perimeter is a binary heap
(Masson & Pride, *Transp. Porous Med.* 102, 301, 2014), so `M` invasions take
`O(M log M)` time. Thresholds are hashed rather than stored, so memory is about one
byte per lattice site plus 16 bytes per perimeter element.

```julia
ip = IPModel(SquareLattice(4096); percolation = SitePercolation(), periodic = true, seeds = 4096)
r = simulate(ip; rng = Xoshiro(1), medium = Medium(42))
r.x         # threshold of each invaded element, in order
r.status    # :spanned when a fluid reaches the right wall
```

Run with the same `medium` (or the same `rng` seed), `IPModel` and a `Persistent`
`SOPModel` see the same thresholds. Use `keep_history = false` on very large runs.
`reproduce/invasion.jl` checks the cost per invasion and the fractal dimension.

## Trapping

With `trapping = true` (in `SOPModel`, `IPModel` and `CVModel`), a region of
defender cut off from the outlet, the last column, is trapped: it can no longer be
invaded, and its sites read `TRAPPED` in `sites`. With `GrowingRight` the outlet is
the current last column, which stays ahead of the front. In SOP the sites of one
step are placed one at a time, so a site trapped earlier in the step is skipped.

After each invasion, one search starts from each empty neighbor of the invaded
site, and the searches advance in turn until all but one are resolved: joined to
another, at the outlet, or closed. Closed regions are trapped. The cost per
invasion is about twice that of the run without trapping on the square lattice.
In `CVModel` trapping keeps the front, and so the compliance of the menisci, from
growing without bound, as in Furuberg, Måløy & Feder (*Phys. Rev. E* 53, 966, 1996).

## Capacitive-volume model

`CVModel` is the model of slow drainage of Måløy et al. (*Phys. Rev. Lett.* 68, 2161,
1992): invasion percolation with a pressure and a time. One fluid is injected at a
constant `rate`, and the volume that does not fill pores is stored by a compliance
`C = compliance + meniscus * n_f`, the injection system plus the `n_f` menisci on the
front. The weakest element on the perimeter is invaded once the pressure reaches its
threshold; if the pressure is below it, it first builds up, which takes time and
starts a new burst. Each invaded pore lowers the pressure by `1/C`.

The invasion order is exactly that of `IPModel` on the same medium, so the model adds
only the pressure signal. Pressure is in threshold units, the `p` of `SOPModel`, and
SOP's parameters convert as `compliance = A/c` and `rate = f_T A`, with `A` the wall
area. Måløy et al. have `meniscus > 0` and `compliance = 0`; SOP has the opposite.

```julia
cv = CVModel(SquareLattice(512); compliance = 512 / 0.05, meniscus = 0.0, rate = 1.0,
    percolation = BondPercolation(), seeds = 512, periodic = true)
r = simulate(cv; rng = Xoshiro(1), medium = Medium(42))
r.p, r.t, r.nf   # pressure after, time of, and menisci before each invasion
b = bursts(r)    # size, drop, peak, start and time of each burst
```

## Three dimensions

`CubicLattice(ny, nz)` injects through an `ny × nz` wall; `CubicLattice(W)` is a
`W × W` wall. `f(t)` is measured per wall site, `N(t)/(ny nz)`, and `periodic` wraps
both transverse directions. `sites` and `times` come back as `ny × nz × L` arrays.
`BlockSeeds` is two-dimensional only.

## Adding a lattice

A new lattice needs one type and one `neighbors` method. If its wall is not a line it
also defines `transverse_sites`, `default_length` and `lattice_shape`, as
`CubicLattice` does; `x` stays the last index.

## Reproducing the paper

```sh
julia -t auto --project=reproduce -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia -t auto --project=reproduce reproduce/paper1997.jl
julia -t auto --project=reproduce reproduce/validate_extensions.jl
julia -t auto --project=reproduce reproduce/benchmark.jl
julia -t auto --project=reproduce reproduce/invasion.jl
```

Figures and `stationary_p.csv` are written to `reproduce/output/`.
