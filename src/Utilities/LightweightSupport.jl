module LightweightSupport

using Logging
using Dates

export Progress, next!, @showprogress, find_zero_bisection, knn_bruteforce,
    StreamLogger, TeeLogger

mutable struct Progress
    total::Int
    message::String
    current::Int
    interval::Int
end

function Progress(total::Integer, message::AbstractString = "")
    total_int = max(0, Int(total))
    interval = max(1, cld(max(total_int, 1), 20))
    return Progress(total_int, String(message), 0, interval)
end

"""
    @showprogress ex

进度显示宏。当前实现为透传（不改变原表达式语义），用于保持调用方 API 兼容；
与 `Progress` / `next!` 配合可在循环中手动输出进度。
"""
macro showprogress(ex)
    return esc(ex)
end

function next!(progress::Progress, increment::Integer = 1)
    progress.current = min(progress.total, progress.current + Int(increment))
    if !isempty(progress.message) && (progress.current == progress.total || progress.current % progress.interval == 0)
        @info progress.message progress = "$(progress.current)/$(progress.total)"
    end
    return progress.current
end

function _bisect_zero(f, left::Float64, right::Float64, fleft::Float64, fright::Float64; tol::Float64, max_iter::Int)
    for _ = 1:max_iter
        mid = (left + right) / 2
        fmid = Float64(f(mid))
        if abs(fmid) <= tol || abs(right - left) <= tol
            return mid
        end
        if signbit(fleft) == signbit(fmid)
            left = mid
            fleft = fmid
        else
            right = mid
            fright = fmid
        end
    end
    return (left + right) / 2
end

function find_zero_bisection(f, x0::Real; initial_step::Real = 1.0, tol::Real = 1e-10, max_expand::Integer = 64, max_iter::Integer = 256)
    x_center = Float64(x0)
    f_center = Float64(f(x_center))
    abs(f_center) <= tol && return x_center

    step = max(abs(Float64(initial_step)), tol)
    left = x_center
    right = x_center
    fleft = f_center
    fright = f_center

    for _ = 1:Int(max_expand)
        right = x_center + step
        fright = Float64(f(right))
        if signbit(f_center) != signbit(fright)
            return _bisect_zero(f, x_center, right, f_center, fright; tol = Float64(tol), max_iter = Int(max_iter))
        end

        left = x_center - step
        fleft = Float64(f(left))
        if signbit(fleft) != signbit(f_center)
            return _bisect_zero(f, left, x_center, fleft, f_center; tol = Float64(tol), max_iter = Int(max_iter))
        end

        step *= 2
    end

    throw(ArgumentError("Failed to bracket root around x0=$(x0)"))
end

function knn_bruteforce(target_nodes::AbstractMatrix{T}, query_nodes::AbstractMatrix{T}, k::Integer; return_distance::Bool = true) where {T<:Real}
    n_target = size(target_nodes, 2)
    n_query = size(query_nodes, 2)
    k_eff = min(Int(k), n_target)
    idxs = Vector{Vector{Int}}(undef, n_query)
    dists = Vector{Vector{Float64}}(undef, n_query)

    for j = 1:n_query
        distances = Vector{Float64}(undef, n_target)
        query = @view query_nodes[:, j]
        for i = 1:n_target
            delta = @view target_nodes[:, i]
            acc = 0.0
            @inbounds for n = 1:length(query)
                diff = Float64(delta[n] - query[n])
                acc += diff * diff
            end
            distances[i] = acc
        end
        order = partialsortperm(distances, 1:k_eff)
        idxs[j] = collect(order)
        if return_distance
            dists[j] = sqrt.(distances[order])
        else
            dists[j] = Float64[]
        end
    end

    return idxs, dists
end

struct StreamLogger <: AbstractLogger
    io::IO
    min_level::LogLevel
end

Logging.min_enabled_level(logger::StreamLogger) = logger.min_level
Logging.shouldlog(logger::StreamLogger, level, _module, group, id) = level >= logger.min_level
Logging.catch_exceptions(::StreamLogger) = false

function Logging.handle_message(logger::StreamLogger, level, message, _module, group, id, file, line; kwargs...)
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println(logger.io, "[$timestamp] [$level] [$(_module)] $(message)")
    flush(logger.io)
    return nothing
end

struct TeeLogger{T<:Tuple} <: AbstractLogger
    loggers::T
end

TeeLogger(loggers::AbstractLogger...) = TeeLogger(tuple(loggers...))

function Logging.min_enabled_level(logger::TeeLogger)
    return reduce(min, map(Logging.min_enabled_level, logger.loggers))
end

function Logging.shouldlog(logger::TeeLogger, level, _module, group, id)
    return any(inner -> Logging.shouldlog(inner, level, _module, group, id), logger.loggers)
end

Logging.catch_exceptions(::TeeLogger) = false

function Logging.handle_message(logger::TeeLogger, level, message, _module, group, id, file, line; kwargs...)
    for inner in logger.loggers
        if Logging.shouldlog(inner, level, _module, group, id)
            Logging.handle_message(inner, level, message, _module, group, id, file, line; kwargs...)
        end
    end
    return nothing
end

end
