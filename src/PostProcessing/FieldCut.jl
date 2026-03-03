module FieldCut

using StaticArrays
using LinearAlgebra
using ...CoreModule
using ...Geometry
using ...BasisFunctions
using ..NearField: calculate_near_field

export field_cut_line, field_cut_plane

"""
    field_cut_line(p1, p2, N, basis, I_coeffs [, perms])

Sample the electric near-field along the line segment from `p1` to `p2`
at `N` uniformly-spaced observation points (endpoints included).

# Returns
- `pts::Vector{SVector{3,FT}}` — observation point coordinates
- `E::Vector{SVector{3,Complex{FT}}}` — electric field at each point

# Example
```julia
pts, E = field_cut_line(SVector(0,0,1.0), SVector(0,0,3.0), 20, basis, I)
```
"""
function field_cut_line(
    p1::SVector{3,FT},
    p2::SVector{3,FT},
    N::Integer,
    args...;          # basis [, I_coeffs [, perms]] forwarded to calculate_near_field
    kwargs...,
) where {FT<:Real}
    N >= 2 || throw(ArgumentError("N must be ≥ 2 to include both endpoints"))

    ts = range(FT(0), FT(1); length=N)
    pts = [p1 + t * (p2 - p1) for t in ts]

    E = calculate_near_field(pts, args...; kwargs...)
    return pts, E
end

"""
    field_cut_plane(origin, u, v, Nu, Nv, basis, I_coeffs [, perms])

Sample the electric near-field on a 2-D planar grid.

The grid spans:
  `P(i,j) = origin + (i-1)/(Nu-1) * u + (j-1)/(Nv-1) * v`

for `i = 1..Nu`, `j = 1..Nv`.

# Arguments
- `origin`: Corner point of the plane
- `u`: Edge vector along the first axis (full extent)
- `v`: Edge vector along the second axis (full extent)
- `Nu`, `Nv`: Number of sample points along each axis (≥ 2)

# Returns
- `pts::Matrix{SVector{3,FT}}` — shape `(Nu, Nv)`
- `E::Matrix{SVector{3,Complex{FT}}}` — shape `(Nu, Nv)`

# Example
```julia
origin = SVector(0.0, -1.0, -1.0)
u_vec  = SVector(0.0,  2.0,  0.0)
v_vec  = SVector(0.0,  0.0,  2.0)
pts, E = field_cut_plane(origin, u_vec, v_vec, 21, 21, basis, I)
```
"""
function field_cut_plane(
    origin::SVector{3,FT},
    u::SVector{3,FT},
    v::SVector{3,FT},
    Nu::Integer,
    Nv::Integer,
    args...;
    kwargs...,
) where {FT<:Real}
    Nu >= 2 || throw(ArgumentError("Nu must be ≥ 2"))
    Nv >= 2 || throw(ArgumentError("Nv must be ≥ 2"))

    # Build flat list of points (row-major over u, then v)
    pts_flat = Vector{SVector{3,FT}}(undef, Nu * Nv)
    idx = 0
    for j in 1:Nv, i in 1:Nu
        idx += 1
        su = FT(i - 1) / FT(Nu - 1)
        sv = FT(j - 1) / FT(Nv - 1)
        pts_flat[idx] = origin + su * u + sv * v
    end

    E_flat = calculate_near_field(pts_flat, args...; kwargs...)

    # Reshape to (Nu, Nv)
    pts = reshape(pts_flat, Nu, Nv)
    E   = reshape(E_flat,   Nu, Nv)
    return pts, E
end

end # module FieldCut
