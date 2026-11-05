# Binary min-heap of (threshold, site) used by the Persistent rule. It is written
# here so that a perimeter of millions of bonds stays allocation-free between resizes.
const Threshold = Tuple{Float64,Int}

function heap_push!(h::Vector{Threshold}, item::Threshold)
    push!(h, item)
    i = length(h)
    @inbounds while i > 1
        parent = i >> 1
        h[parent][1] <= item[1] && break
        h[i] = h[parent]
        i = parent
    end
    @inbounds h[i] = item
    return h
end

function heap_pop!(h::Vector{Threshold})
    @inbounds top = h[1]
    last = pop!(h)
    n = length(h)
    n == 0 && return top
    i = 1
    @inbounds while true
        child = 2i
        child > n && break
        child < n && h[child+1][1] < h[child][1] && (child += 1)
        h[child][1] >= last[1] && break
        h[i] = h[child]
        i = child
    end
    @inbounds h[i] = last
    return top
end
