"""
    ensemble(model, n; seed=1, ntasks=Threads.nthreads(), f=identity)

Runs `n` independent realizations on `ntasks` threads and returns `f(result)` for
each. Realization `i` always uses the same random stream, so the output does not
depend on `ntasks`. Use `f` to keep a summary instead of the full time series.
"""
function ensemble(m::SOPModel, n::Integer; seed::Integer=1, ntasks::Integer=Threads.nthreads(), f=identity)
    n >= 1 || throw(ArgumentError("n must be positive"))
    run(i) = f(simulate(m; rng=Xoshiro(seed + i - 1)))
    results = Vector{Any}(undef, n)
    next = Threads.Atomic{Int}(1)
    @sync for _ in 1:min(ntasks, n)
        Threads.@spawn while true
            i = Threads.atomic_add!(next, 1)
            i > n && break
            results[i] = run(i)
        end
    end
    # Narrow the element type, which may be a union if `f` returns different types.
    return identity.(results)
end

"""
    stationary_mean(series; burn=length(series) ÷ 2)

Mean of the time series after dropping the first `burn` points.
"""
function stationary_mean(series::AbstractVector; burn::Integer=length(series) ÷ 2)
    tail = @view series[burn+1:end]
    isempty(tail) && throw(ArgumentError("burn removes the whole series"))
    return sum(tail) / length(tail)
end
