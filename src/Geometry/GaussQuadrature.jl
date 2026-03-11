using FastGaussQuadrature
using StaticArrays

struct GaussQuadratureInfoStruct{FT<:Real,N,D}
    coordinate::SMatrix{D,N,FT}
    weight::SVector{N,FT}
    uvw::SMatrix{0,0,FT} # Placeholder
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
    uvw = SMatrix{0,0,FT,0}()
    if GeoS == :Triangle
        coordinate, weight = gaussQuadratureTri(GQN, FT)
        N = length(weight)
        return GaussQuadratureInfoStruct{FT,N,3}(coordinate, weight, uvw)
    elseif GeoS == :Tetrahedron
        coordinate, weight = gaussQuadratureTet(GQN, FT)
        N = length(weight)
        return GaussQuadratureInfoStruct{FT,N,4}(coordinate, weight, uvw)
    elseif GeoS == :Hexahedron
        coordinate, weight = gaussQuadratureHexa(GQN, FT)
        N = length(weight)
        return GaussQuadratureInfoStruct{FT,N,8}(coordinate, weight, uvw)
    elseif GeoS == :Quadrangle
        coordinate, weight = gaussQuadratureQuad(GQN, FT)
        N = length(weight)
        return GaussQuadratureInfoStruct{FT,N,4}(coordinate, weight, uvw)
    else
        error("Geometry type $GeoS not supported yet.")
    end
end

"""
    gaussQuadratureHexa(num::Integer, FT::DataType = Float64)

Generate Gauss quadrature for hexahedra in [0,1]³ using tensor-product Gauss-Legendre rules.
`num` must be a perfect cube (n³): 1, 8, 27, 64.
Returns coordinate (8×N shape function values) and weight (N-vector, sum=1).
"""
function gaussQuadratureHexa(num::Integer, FT::DataType = Float64)
    # Determine 1D order
    n1d = round(Int, cbrt(num))
    @assert n1d^3 == num "Hexahedron GQ requires N = n³ (got $num)"

    # 1D Gauss-Legendre on [-1,1], then map to [0,1]
    x_gl, w_gl = gausslegendre(n1d)
    x_01 = FT.(0.5 .+ 0.5 .* x_gl)  # Map to [0,1]
    w_01 = FT.(0.5 .* w_gl)           # Scale weights

    N = num
    # Build coordinate matrix (8×N): shape function values at each GQ point
    # Shape functions for 8-node hex on [0,1]³:
    # N₁ = (1-u)(1-v)(1-w), N₂ = u(1-v)(1-w), N₃ = uv(1-w), N₄ = (1-u)v(1-w)
    # N₅ = (1-u)(1-v)w, N₆ = u(1-v)w, N₇ = uvw, N₈ = (1-u)vw
    coords = zeros(FT, 8, N)
    weights = zeros(FT, N)

    idx = 1
    for k = 1:n1d, j = 1:n1d, i = 1:n1d
        u, v, w = x_01[i], x_01[j], x_01[k]
        coords[1, idx] = (1 - u) * (1 - v) * (1 - w)
        coords[2, idx] = u * (1 - v) * (1 - w)
        coords[3, idx] = u * v * (1 - w)
        coords[4, idx] = (1 - u) * v * (1 - w)
        coords[5, idx] = (1 - u) * (1 - v) * w
        coords[6, idx] = u * (1 - v) * w
        coords[7, idx] = u * v * w
        coords[8, idx] = (1 - u) * v * w
        weights[idx] = w_01[i] * w_01[j] * w_01[k]
        idx += 1
    end

    coordinate = SMatrix{8,N,FT}(coords)
    weight = SVector{N,FT}(weights)
    return coordinate, weight
end

"""
    gaussQuadratureQuad(num::Integer, FT::DataType = Float64)

Generate Gauss quadrature for quadrangles in [0,1]² using tensor-product Gauss-Legendre rules.
`num` must be a perfect square (n²): 1, 4, 9, 16.
Returns coordinate (4×N shape function values) and weight (N-vector, sum=1).
"""
function gaussQuadratureQuad(num::Integer, FT::DataType = Float64)
    n1d = round(Int, sqrt(num))
    @assert n1d^2 == num "Quadrangle GQ requires N = n² (got $num)"

    x_gl, w_gl = gausslegendre(n1d)
    x_01 = FT.(0.5 .+ 0.5 .* x_gl)
    w_01 = FT.(0.5 .* w_gl)

    N = num
    coords = zeros(FT, 4, N)
    weights = zeros(FT, N)

    idx = 1
    for j = 1:n1d, i = 1:n1d
        u, v = x_01[i], x_01[j]
        coords[1, idx] = (1 - u) * (1 - v)
        coords[2, idx] = u * (1 - v)
        coords[3, idx] = u * v
        coords[4, idx] = (1 - u) * v
        weights[idx] = w_01[i] * w_01[j]
        idx += 1
    end

    coordinate = SMatrix{4,N,FT}(coords)
    weight = SVector{N,FT}(weights)
    return coordinate, weight
end

"""
    gaussQuadratureHexa1D(n1d::Integer, FT::DataType = Float64)

Return 1D Gauss-Legendre points and weights on [0,1] for hexahedral quadrature.
Useful for RBF free-end mappings.
"""
function gaussQuadratureHexa1D(n1d::Integer, FT::DataType = Float64)
    x_gl, w_gl = gausslegendre(n1d)
    x_01 = FT.(0.5 .+ 0.5 .* x_gl)
    w_01 = FT.(0.5 .* w_gl)
    return x_01, w_01
end

function gaussQuadratureTet(num::Integer, FT::DataType = Float64)
    if num == 1
        # Centroid
        coordinate = SMatrix{4,1,FT}(0.25, 0.25, 0.25, 0.25)
        weight = SVector{1,FT}(1.0)
    elseif num == 4
        # 4 points
        a = 0.58541020
        b = 0.13819660
        coordinate = SMatrix{4,4,FT}(a, b, b, b, b, a, b, b, b, b, a, b, b, b, b, a)
        weight = SVector{4,FT}(0.25, 0.25, 0.25, 0.25)
    elseif num == 5
        # 5 points
        coordinate = SMatrix{4,5,FT}(
            0.25,
            0.25,
            0.25,
            0.25,
            0.5,
            1 / 6,
            1 / 6,
            1 / 6,
            1 / 6,
            0.5,
            1 / 6,
            1 / 6,
            1 / 6,
            1 / 6,
            0.5,
            1 / 6,
            1 / 6,
            1 / 6,
            1 / 6,
            0.5,
        )
        weight = SVector{5,FT}(-0.8, 0.45, 0.45, 0.45, 0.45)
    elseif num == 11
        a1 = FT(0.714285714286)
        b1 = FT(0.095238095238)
        a2 = FT(0.399403576167)
        b2 = FT(0.100596423833)
        coordinate = SMatrix{4,11,FT}(
            FT(1 / 4), FT(1 / 4), FT(1 / 4), FT(1 / 4),
            a1, b1, b1, b1,
            b1, a1, b1, b1,
            b1, b1, a1, b1,
            b1, b1, b1, a1,
            a2, a2, b2, b2,
            b2, a2, a2, b2,
            b2, b2, a2, a2,
            a2, b2, b2, a2,
            b2, a2, b2, a2,
            a2, b2, a2, b2,
        )
        weight = SVector{11,FT}(
            FT(-0.078933333333),
            FT(0.0457333333333),
            FT(0.0457333333333),
            FT(0.0457333333333),
            FT(0.0457333333333),
            FT(0.1493333333333),
            FT(0.1493333333333),
            FT(0.1493333333333),
            FT(0.1493333333333),
            FT(0.1493333333333),
            FT(0.1493333333333),
        )
    else
        error("Tetrahedron quadrature rule for N=$num not implemented.")
    end
    return coordinate, weight
end

function gaussQuadratureTri(num::Integer, FT::DataType = Float64)
    if num == 1
        coordinate = SMatrix{3,1,FT}(1 / 3, 1 / 3, 1 / 3)
        weight = SVector{1,FT}(1.0)
    elseif num == 3
        coordinate = SMatrix{3,3,FT}(2 / 3, 1 / 6, 1 / 6, 1 / 6, 2 / 3, 1 / 6, 1 / 6, 1 / 6, 2 / 3)
        weight = SVector{3,FT}(1 / 3, 1 / 3, 1 / 3)
    elseif num == 4
        coordinate =
            SMatrix{3,4,FT}(1 / 3, 1 / 3, 1 / 3, 0.6, 0.2, 0.2, 0.2, 0.6, 0.2, 0.2, 0.2, 0.6)
        weight = SVector{4,FT}(-27 / 48, 25 / 48, 25 / 48, 25 / 48) # Check weights
    elseif num == 7
        # 7-point rule (Degree 5)
        w1 = 0.225
        w2 = 0.1323941527885062
        w3 = 0.1259391805448271

        a = 1 / 3
        b1 = 0.0597158717897698
        c1 = 0.4701420641051151
        b2 = 0.7974269853530873
        c2 = 0.1012865073234563

        coordinate = SMatrix{3,7,FT}(
            a,
            a,
            a,
            b1,
            c1,
            c1,
            c1,
            b1,
            c1,
            c1,
            c1,
            b1,
            b2,
            c2,
            c2,
            c2,
            b2,
            c2,
            c2,
            c2,
            b2,
        )
        weight = SVector{7,FT}(w1, w2, w2, w2, w3, w3, w3)
    else
        error("Unsupported number of points for Triangle: $num")
    end
    return coordinate, weight
end

"""
    get_global_quad_points(tri, gq)

Map local quadrature points to global coordinates on a triangle.
"""
function get_global_quad_points(tri, gq::GaussQuadratureInfoStruct{FT,N}) where {FT,N}
    v1 = tri.vertices[:, 1]
    v2 = tri.vertices[:, 2]
    v3 = tri.vertices[:, 3]

    # Use SVector generator to avoid allocation
    return SVector{N,SVector{3,FT}}(
        v1 * gq.coordinate[1, i] + v2 * gq.coordinate[2, i] + v3 * gq.coordinate[3, i] for i = 1:N
    )
end
