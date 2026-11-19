"Value of `sites` for a defender site cut off from the outlet: trapped."
const TRAPPED = 0xff

# Scratch space of the trapping rule. `mark[u]` is a stamp: the search that
# reached site `u` in the current call is `mark[u] - base`, and stamps of earlier
# calls are at most `base`. With trapping off, `mark` is empty.
mutable struct Trapper
    mark::Vector{UInt32}
    base::UInt32
    queues::Vector{Vector{Int}}
    heads::Vector{Int}
    parent::Vector{Int}
    open::Vector{Bool}
    out::Vector{Int}               # sites trapped by the last call
end
const MAXSEARCH = 8                # more than the largest coordination number
Trapper(n::Integer) = Trapper(zeros(UInt32, n), UInt32(0),
    [Int[] for _ in 1:MAXSEARCH], zeros(Int, MAXSEARCH), zeros(Int, MAXSEARCH),
    fill(false, MAXSEARCH), Int[])

enabled(tr::Trapper) = !isempty(tr.mark)

@inline function root(parent, a)
    @inbounds while parent[a] != a
        a = parent[a]
    end
    return a
end

@inline function live(tr::Trapper, k, r)
    @inbounds for b in 1:k
        root(tr.parent, b) == r && tr.heads[b] <= length(tr.queues[b]) &&
            return true
    end
    return false
end

"""
    trap!(tr, sites, lattice, periodic, nx, j)

Traps the defender regions that the invasion of site `j` cut off from the outlet,
the last column of the lattice, and returns them in `tr.out`. Every empty site is
assumed to reach the outlet before `j` is invaded, so every region cut off now
touches `j`. One search starts from each empty neighbor of `j`, and the searches
advance in turn, one site each. A search joins another when it reaches a site
that one has visited, and its group is open when it reaches the outlet. A group
whose searches run out of sites is closed, and trapped.

Unless `j` is on the outlet, at most one region reaches the outlet: two would need
an invaded path to the outlet between them. So once one group is left unresolved
and none is open, that group is the open one and the searches stop. The cost of a
call is then a few times the size of the regions it traps, plus the cost of
finding that the others meet.

Defender sites connect through the lattice's own neighbors, which is the external
perimeter of Furuberg, Måløy & Feder (PRE 53, 966, 1996) on the square lattice.
"""
function trap!(tr::Trapper, sites, lattice, periodic, nx, j)
    empty!(tr.out)
    W = length(sites) ÷ nx
    sink = nx - 1
    if tr.base > typemax(UInt32) - MAXSEARCH
        fill!(tr.mark, 0x00000000)
        tr.base = UInt32(0)
    end
    base = tr.base
    tr.base += UInt32(MAXSEARCH)
    k = 0
    @inbounds for u in neighbors(lattice, j, nx, periodic)
        (u == 0 || sites[u] != EMPTY || tr.mark[u] > base) && continue
        k += 1
        q = tr.queues[k]
        empty!(q)
        push!(q, u)
        tr.heads[k] = 1
        tr.parent[k] = k
        tr.open[k] = (u - 1) ÷ W == sink
        tr.mark[u] = base + UInt32(k)
    end
    k <= 1 && return tr.out
    onsink = (j - 1) ÷ W == sink
    @inbounds while true
        unresolved, anyopen = 0, false
        for a in 1:k
            root(tr.parent, a) == a || continue
            if tr.open[a]
                anyopen = true
            elseif live(tr, k, a)
                unresolved += 1
            end
        end
        unresolved == 0 && break
        unresolved == 1 && !anyopen && !onsink && break
        for a in 1:k
            tr.open[root(tr.parent, a)] && continue
            q = tr.queues[a]
            h = tr.heads[a]
            h > length(q) && continue
            v = q[h]
            tr.heads[a] = h + 1
            for u in neighbors(lattice, v, nx, periodic)
                (u == 0 || sites[u] != EMPTY) && continue
                mu = tr.mark[u]
                if mu > base
                    ra, rb = root(tr.parent, a), root(tr.parent, Int(mu - base))
                    if ra != rb
                        tr.parent[rb] = ra
                        tr.open[ra] |= tr.open[rb]
                    end
                    continue
                end
                tr.mark[u] = base + UInt32(a)
                push!(q, u)
                (u - 1) ÷ W == sink && (tr.open[root(tr.parent, a)] = true)
            end
        end
    end
    @inbounds for a in 1:k
        r = root(tr.parent, a)
        (tr.open[r] || live(tr, k, r)) && continue
        for u in tr.queues[a]
            sites[u] = TRAPPED
            push!(tr.out, u)
        end
    end
    return tr.out
end
