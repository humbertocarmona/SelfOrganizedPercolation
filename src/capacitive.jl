"""
    CVModel(lattice; compliance, meniscus=0.0, rate=1.0, p0=0.0, kwargs...)

Capacitive-volume model of slow drainage (Måløy, Furuberg, Feder & Jøssang, PRL
68, 2161, 1992) on the medium of [`IPModel`](@ref): invasion percolation with a
pressure and a time. One fluid is injected at a constant `rate`; the volume that
does not fill pores is stored by the compliance of the system and of the menisci,

    C(t) dp/dt = rate - (volume of the pores invaded),
    C = compliance + meniscus * n_f,

where `n_f` is the number of menisci on the front: perimeter bonds for
`BondPercolation`, perimeter sites for `SitePercolation`. The pressure `p` is in
the units of the thresholds, so it is the `p` of [`SOPModel`](@ref).

Events are invasions. The weakest element on the perimeter, of threshold `m`, is
invaded as soon as `p >= m`. If `p < m`, the pressure first builds up to `m`,
which takes a time `C (m - p) / rate`, and a new burst starts; otherwise the
element is invaded within the current burst, at no cost in time. Each invasion
fills a pore of unit volume and lowers `p` by `1/C`. The invasion order is
therefore that of invasion percolation on the same medium; the model adds the
pressure, the time and the bursts.

Units: volumes in pores, pressures in threshold units. Filling one pore lowers
`p` by `c/A` in SOP, with `A` the wall area, so a constant system compliance
corresponds to `compliance = A/c`, and SOP injects `rate = f_T A` pores per step.
Måløy et al. have `meniscus > 0` and `compliance = 0`; SOP has the opposite.

# Keywords
- `compliance`: compliance of the injection system, pores per unit of `p`.
- `meniscus`: compliance of one meniscus, pores per unit of `p` (default 0).
- `rate`: injected pores per unit time (default 1).
- `p0`: initial pressure (default 0).
- `percolation`, `right`, `periodic`, `seeds`, `layout`, `active_fraction`,
  `nmax`, `trapping`: as in `IPModel`. There is a single fluid. Furuberg, Måløy
  & Feder (PRE 53, 966, 1996) trap the defender, which keeps the front, and so
  the compliance of the menisci, from growing without bound.
"""
struct CVModel{M<:IPModel}
    ip::M
    compliance::Float64
    meniscus::Float64
    rate::Float64
    p0::Float64
end

function CVModel(lattice::Lattice; compliance::Real, meniscus::Real=0.0,
        rate::Real=1.0, p0::Real=0.0, kwargs...)
    haskey(kwargs, :species) && kwargs[:species] != 1 &&
        throw(ArgumentError("CVModel has a single fluid"))
    compliance >= 0 || throw(ArgumentError("compliance must be non-negative"))
    meniscus >= 0 || throw(ArgumentError("meniscus must be non-negative"))
    compliance + meniscus > 0 ||
        throw(ArgumentError("compliance and meniscus cannot both be zero"))
    rate > 0 || throw(ArgumentError("rate must be positive"))
    ip = IPModel(lattice; species=1, kwargs...)
    return CVModel(ip, Float64(compliance), Float64(meniscus), Float64(rate),
        Float64(p0))
end

"""
Outcome of one capacitive-volume run. For the `n`-th invasion, `x[n]` is the
threshold of the invaded element, `t[n]` the time at which it is invaded, `p[n]`
the pressure just after it, and `nf[n]` the number of menisci on the front just
before it. The pressure just before invasion `n` is `max(p[n-1], x[n])`, and a
burst starts at invasion `n` when `x[n] > p[n-1]`; see [`bursts`](@ref). The
histories are empty unless `keep_history=true`.

`status`, `mass`, `length`, `sites` and `times` are as in [`IPResult`](@ref);
`times` gives the invasion order, not the time `t`.
"""
struct CVResult
    x::Vector{Float64}
    t::Vector{Float64}
    p::Vector{Float64}
    nf::Vector{Int32}
    status::Symbol
    mass::Int
    length::Int
    p0::Float64
    sites::Union{Nothing,Array{UInt8}}
    times::Union{Nothing,Array{Int32}}
end

# Menisci that disappear when site `j` is invaded: the perimeter entries that
# point to it. In site percolation `j` is on the perimeter once. In bond
# percolation every occupied neighbor that has opened its perimeter holds a bond
# to `j`; inactive seeds never did.
menisci_at(st, m, ::SitePercolation, j, inactive) = 1
function menisci_at(st, m, ::BondPercolation, j, inactive)
    n = 0
    @inbounds for i in neighbors(m.lattice, j, st.nx, m.periodic)
        i == 0 && continue
        st.sites[i] == 0x01 && !(i in inactive) && (n += 1)
    end
    return n
end
# A trapped site drops out of the perimeter: in site percolation if the fluid had
# reached it, in bond percolation with every bond that points to it. The site `j`
# just invaded has not opened its perimeter yet, so it holds no bond to `u`.
menisci_trapped(st, m, ::SitePercolation, u, j, inactive) =
    Int(@inbounds st.tested[u] & 0x01)
menisci_trapped(st, m, ::BondPercolation, u, j, inactive) =
    menisci_at(st, m, BondPercolation(), u, inactive) -
    count(==(j), neighbors(m.lattice, u, st.nx, m.periodic))

function invade_cv!(st::IPState, cv::CVModel, inactive, hist)
    m = cv.ip
    W = transverse_sites(m.lattice)
    heap = st.heaps[1]
    nf = length(heap)            # no entry is stale before the first invasion
    p, t = cv.p0, 0.0
    n = 0
    keep = !isempty(st.times)
    while n < m.nmax
        next = next_element!(heap, st.sites)
        next === nothing && return :extinct, p
        x, j = next
        col = (j - 1) ÷ W
        if !room_right!(st, m.right, W, col)
            heap_push!(heap, next)
            return :maxlength, p
        end
        C = cv.compliance + cv.meniscus * nf
        if x > p                 # build up to the threshold: a new burst
            t += C * (x - p) / cv.rate
            p = x
        end
        p -= 1 / C
        n += 1
        if hist !== nothing
            push!(hist.x, x); push!(hist.t, t); push!(hist.p, p)
            push!(hist.nf, Int32(nf))
        end
        nf -= menisci_at(st, m, m.percolation, j, inactive)
        @inbounds st.sites[j] = 0x01
        if keep
            n <= typemax(Int32) ||
                throw(OverflowError("invasion order exceeds Int32"))
            @inbounds st.times[j] = Int32(n)
        end
        st.mass[1] += 1
        st.xmax = max(st.xmax, col)
        spans(m.right, st, col) && return :spanned, p
        if m.trapping
            for u in trap!(st.trap, st.sites, m.lattice, m.periodic, st.nx, j)
                nf -= menisci_trapped(st, m, m.percolation, u, j, inactive)
            end
        end
        before = length(heap)
        add_perimeter!(st, m, 1, j)
        nf += length(heap) - before
    end
    return :nmax, p
end

"""
    simulate(model::CVModel; rng=Xoshiro(), medium=Medium(rand(rng, UInt64)),
             keep_history=true, keep_sites=false, keep_times=false)

Runs the capacitive-volume model and returns a [`CVResult`](@ref). The medium
and the seed positions are drawn from `rng` as for `IPModel` and `SOPModel`, so
the three models run with the same seed share them.
"""
function simulate(cv::CVModel; rng=Xoshiro(),
        medium::Medium=Medium(rand(rng, UInt64)), keep_history::Bool=true,
        keep_sites::Bool=false, keep_times::Bool=false)
    m = cv.ip
    st = initial_state(m, rng, medium, keep_times)
    # Seeds that never opened their perimeter: no meniscus points away from them.
    inactive = setdiff!(Set(findall(!=(EMPTY), st.sites)), st.active)
    hist = keep_history ?
        (x=Float64[], t=Float64[], p=Float64[], nf=Int32[]) : nothing
    status, _ = invade_cv!(st, cv, inactive, hist)
    shape = lattice_shape(m.lattice, st.nx)
    h = something(hist, (x=Float64[], t=Float64[], p=Float64[], nf=Int32[]))
    return CVResult(h.x, h.t, h.p, h.nf, status, st.mass[1], st.xmax + 1,
        cv.p0, keep_sites ? reshape(st.sites, shape) : nothing,
        keep_times ? reshape(st.times, shape) : nothing)
end

"""
    bursts(result::CVResult)

Splits the invasions into bursts, as in Måløy et al.: a burst starts when the
pressure has to build up to the weakest threshold on the front, and it goes on,
at no cost in time, while the next threshold is below the pressure. Returns a
named tuple of vectors, one entry per burst: `size` (invasions), `drop` (fall
of the pressure from its peak, the first threshold, to its value at the end),
`peak`, `start` (index of the first invasion) and `time`.
"""
function bursts(r::CVResult)
    n = length(r.x)
    size, drop, peak, start, time = Int[], Float64[], Float64[], Int[], Float64[]
    prev = r.p0
    for k in 1:n
        if k == 1 || r.x[k] > prev
            push!(size, 0); push!(peak, max(prev, r.x[k])); push!(start, k)
            push!(time, r.t[k]); push!(drop, 0.0)
        end
        size[end] += 1
        drop[end] = peak[end] - r.p[k]
        prev = r.p[k]
    end
    return (; size, drop, peak, start, time)
end
