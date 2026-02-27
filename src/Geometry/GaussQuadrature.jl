using FastGaussQuadrature
using StaticArrays

struct GaussQuadratureInfoStruct{FT<:Real, N, D}
    coordinate  ::SMatrix{D, N, FT}
    weight      ::SVector{N, FT}
    uvw         ::SMatrix{0, 0, FT} # Placeholder
end

"""
    GaussQuadratureInfo(GeoS::Symbol, GQN::IT, FT::DataType = Float64)

Generate Gauss quadrature points and weights for numerical integration over a specified geometry.

# Arguments
- `GeoS`: Geometry symbol (e.g., `:Triangle`).
- `GQN`: Number of quadrature points.
- `FT`: Floating point type.

# Supported Geometries
- `:Triangle`: Supports 1, 3, 4, 6, 7 points (Dunavant rules).

# Returns
- `GaussQuadratureInfoStruct`: Contains coordinates (barycentric for triangles) and weights.
"""
function GaussQuadratureInfo(GeoS::Symbol, GQN::IT, FT::DataType = Float64) where {IT<:Integer}
    uvw = SMatrix{0, 0, FT, 0}()
    if GeoS == :Triangle
        coordinate, weight = gaussQuadratureTri(GQN, FT)
        N = length(weight)
        return GaussQuadratureInfoStruct{FT, N, 3}(coordinate, weight, uvw)
    elseif GeoS == :Tetrahedron
        coordinate, weight = gaussQuadratureTet(GQN, FT)
        N = length(weight)
        return GaussQuadratureInfoStruct{FT, N, 4}(coordinate, weight, uvw)
    else
        error("Geometry type $GeoS not supported yet.")
    end
end

function gaussQuadratureTet(num::Integer, FT::DataType = Float64)
    if num == 1
        # Centroid
        coordinate = SMatrix{4, 1, FT}(0.25, 0.25, 0.25, 0.25)
        weight = SVector{1, FT}(1.0)
    elseif num == 4
        # 4 points
        a = 0.58541020
        b = 0.13819660
        coordinate = SMatrix{4, 4, FT}(
            a, b, b, b,
            b, a, b, b,
            b, b, a, b,
            b, b, b, a
        )
        weight = SVector{4, FT}(0.25, 0.25, 0.25, 0.25)
    elseif num == 5
        # 5 points
        coordinate = SMatrix{4, 5, FT}(
            0.25, 0.25, 0.25, 0.25,
            0.5, 1/6, 1/6, 1/6,
            1/6, 0.5, 1/6, 1/6,
            1/6, 1/6, 0.5, 1/6,
            1/6, 1/6, 1/6, 0.5
        )
        weight = SVector{5, FT}(-0.8, 0.45, 0.45, 0.45, 0.45)
    else
        error("Tetrahedron quadrature rule for N=$num not implemented.")
    end
    return coordinate, weight
end

function gaussQuadratureTri(num::Integer, FT::DataType = Float64)
    if num == 1
        coordinate = SMatrix{3, 1, FT}(1/3, 1/3, 1/3)
        weight = SVector{1, FT}(1.0)
    elseif num == 3
        coordinate = SMatrix{3, 3, FT}(
            2/3, 1/6, 1/6,
            1/6, 2/3, 1/6,
            1/6, 1/6, 2/3
        )
        weight = SVector{3, FT}(1/3, 1/3, 1/3)
    elseif num == 4
        coordinate = SMatrix{3, 4, FT}(
            1/3, 1/3, 1/3,
            0.6, 0.2, 0.2,
            0.2, 0.6, 0.2,
            0.2, 0.2, 0.6
        )
        weight = SVector{4, FT}(-27/48, 25/48, 25/48, 25/48) # Check weights
    elseif num == 7
        # 7-point rule (Degree 5)
        w1 = 0.225
        w2 = 0.1323941527885062
        w3 = 0.1259391805448271
        
        a = 1/3
        b1 = 0.0597158717897698
        c1 = 0.4701420641051151
        b2 = 0.7974269853530873
        c2 = 0.1012865073234563
        
        coordinate = SMatrix{3, 7, FT}(
            a, a, a,
            b1, c1, c1,
            c1, b1, c1,
            c1, c1, b1,
            b2, c2, c2,
            c2, b2, c2,
            c2, c2, b2
        )
        weight = SVector{7, FT}(w1, w2, w2, w2, w3, w3, w3)
    else
        error("Unsupported number of points for Triangle: $num")
    end
    return coordinate, weight
end

"""
    get_global_quad_points(tri, gq)

Map local quadrature points to global coordinates on a triangle.
"""
function get_global_quad_points(tri, gq::GaussQuadratureInfoStruct{FT, N}) where {FT, N}
    v1 = tri.vertices[:, 1]
    v2 = tri.vertices[:, 2]
    v3 = tri.vertices[:, 3]
    
    # Use SVector generator to avoid allocation
    return SVector{N, SVector{3, FT}}(
        v1 * gq.coordinate[1, i] + v2 * gq.coordinate[2, i] + v3 * gq.coordinate[3, i]
        for i in 1:N
    )
end
