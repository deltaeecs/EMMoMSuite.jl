"""
    select_gap_feed_edges(basis::RWGBasis; axis=3, center=0.0, atol=nothing)

Select the non-boundary RWG edges that lie on the feed gap plane of a wire-like
antenna mesh.

When an exact plane is present (for example `z=0` for an even-layer dipole mesh),
this returns that full ring. If no basis centers lie exactly on the requested
plane, it falls back to the nearest basis-center layer.
"""
function select_gap_feed_edges(
    basis::RWGBasis;
    axis::Int = 3,
    center::Real = 0.0,
    atol::Union{Nothing,Real} = nothing,
)
    axis in 1:3 || throw(ArgumentError("axis must be 1, 2, or 3"))

    N = num_basis(basis)
    N == 0 && return Int[]

    distances = Vector{Float64}(undef, N)
    max_abs_coord = 0.0
    min_dist = Inf

    @inbounds for idx = 1:N
        coord = Float64(basis.functions[idx].center[axis])
        dist = abs(coord - center)
        distances[idx] = dist
        max_abs_coord = max(max_abs_coord, abs(coord))
        min_dist = min(min_dist, dist)
    end

    plane_tol = isnothing(atol) ? 128 * eps(Float64) * max(1.0, max_abs_coord) : Float64(atol)
    threshold = min_dist + plane_tol

    selected = Int[]
    @inbounds for idx = 1:N
        bf = basis.functions[idx]
        bf.is_boundary && continue
        distances[idx] <= threshold && push!(selected, idx)
    end

    return selected
end