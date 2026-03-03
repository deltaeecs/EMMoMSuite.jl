"""
    GeomKernel.jl

Phase 18.1 — CSG / B-Rep geometric kernel: data structures and basic solids.

Exported types:
  BRepFace     — a planar polygon face of a B-Rep solid
  BRepSolid    — a closed polyhedral solid (vertices + faces + edges + labels)
  CSGNode      — node in a Constructive Solid Geometry tree

Exported constructors / utilities:
  box_solid              — axis-aligned box
  solid_volume           — divergence-theorem volume (signed → abs)
  solid_surface_area     — total surface area
  check_manifold         — verify each edge is shared by exactly 2 faces
  convert_to_triangle_mesh — triangulate all faces → TriangleMesh
"""

using LinearAlgebra
using StaticArrays

export BRepFace, BRepSolid, CSGNode
export box_solid, solid_volume, solid_surface_area
export check_manifold, convert_to_triangle_mesh

# ─────────────────────────────────────────────────────────────────────────────
# Data types
# ─────────────────────────────────────────────────────────────────────────────

"""
    BRepFace

One planar face of a B-Rep solid.  Vertices are listed in CCW order when the
face is viewed from *outside* (outward normal direction), consistent with the
right-hand rule.
"""
struct BRepFace
    vertex_indices::Vector{Int}   # polygon winding (CCW from outside)
end

"""
    BRepSolid{FT}

A closed polyhedral solid represented as a Boundary Representation (B-Rep).

- `vertices`         : 3-D coordinates of each vertex.
- `edges`            : undirected edges `(v_lo, v_hi)` (no duplicates).
- `faces`            : polygon faces with CCW outward winding.
- `material_label`   : material name (empty = default).
- `boundary_labels`  : face-index → boundary condition name.
"""
struct BRepSolid{FT<:AbstractFloat}
    vertices        :: Vector{SVector{3,FT}}
    edges           :: Vector{Tuple{Int,Int}}
    faces           :: Vector{BRepFace}
    material_label  :: String
    boundary_labels :: Dict{Int,String}
end

"""
    CSGNode

A node in a Constructive Solid Geometry (CSG) tree.

- `op`    : `:union`, `:intersect`, `:subtract`, or `:leaf`.
- `left`  : left sub-node or leaf solid.
- `right` : right sub-node, leaf solid, or `nothing` for leaf nodes.
"""
struct CSGNode
    op     :: Symbol
    left   :: Union{CSGNode, BRepSolid{<:AbstractFloat}}
    right  :: Union{CSGNode, BRepSolid{<:AbstractFloat}, Nothing}
end

# ─────────────────────────────────────────────────────────────────────────────
# box_solid
# ─────────────────────────────────────────────────────────────────────────────

"""
    box_solid(Lx, Ly, Lz; origin=SVector(0.,0.,0.), material_label="", boundary_labels=Dict()) → BRepSolid

Create an axis-aligned box with dimensions Lx × Ly × Lz.

`origin` is the position of the minimum-coordinate corner vertex
(i.e. the box spans `[origin_x, origin_x+Lx] × ... × [origin_z, origin_z+Lz]`).

The 6 quadrilateral faces are ordered with CCW winding when viewed from outside:
  1. Bottom  (z=0,   outward n = −z)
  2. Top     (z=Lz,  outward n = +z)
  3. Front   (y=0,   outward n = −y)
  4. Back    (y=Ly,  outward n = +y)
  5. Left    (x=0,   outward n = −x)
  6. Right   (x=Lx,  outward n = +x)
"""
function box_solid(
    Lx::FT, Ly::FT, Lz::FT;
    origin::SVector{3,FT} = SVector(zero(FT), zero(FT), zero(FT)),
    material_label::String = "",
    boundary_labels::Dict{Int,String} = Dict{Int,String}(),
) where {FT<:AbstractFloat}
    o = origin
    # 8 corners (1-indexed)
    verts = SVector{3,FT}[
        o + SVector(zero(FT), zero(FT), zero(FT)),   # 1: (0,0,0)
        o + SVector(Lx,       zero(FT), zero(FT)),   # 2: (Lx,0,0)
        o + SVector(zero(FT), Ly,       zero(FT)),   # 3: (0,Ly,0)
        o + SVector(Lx,       Ly,       zero(FT)),   # 4: (Lx,Ly,0)
        o + SVector(zero(FT), zero(FT), Lz      ),   # 5: (0,0,Lz)
        o + SVector(Lx,       zero(FT), Lz      ),   # 6: (Lx,0,Lz)
        o + SVector(zero(FT), Ly,       Lz      ),   # 7: (0,Ly,Lz)
        o + SVector(Lx,       Ly,       Lz      ),   # 8: (Lx,Ly,Lz)
    ]

    # 6 quad faces with outward CCW winding
    # (verified by divergence theorem: Σ V1·(V2×V3)/6 = Lx*Ly*Lz)
    faces = BRepFace[
        BRepFace([1, 3, 4, 2]),   # 1 bottom (z=0,  n=-z)
        BRepFace([5, 6, 8, 7]),   # 2 top    (z=Lz, n=+z)
        BRepFace([1, 2, 6, 5]),   # 3 front  (y=0,  n=-y)
        BRepFace([3, 7, 8, 4]),   # 4 back   (y=Ly, n=+y)
        BRepFace([1, 5, 7, 3]),   # 5 left   (x=0,  n=-x)
        BRepFace([2, 4, 8, 6]),   # 6 right  (x=Lx, n=+x)
    ]

    edges = _extract_edges(faces)

    return BRepSolid{FT}(verts, edges, faces, material_label, boundary_labels)
end

# Allow Integer arguments to be auto-promoted to Float64
function box_solid(Lx::Real, Ly::Real, Lz::Real; kwargs...)
    return box_solid(Float64(Lx), Float64(Ly), Float64(Lz); kwargs...)
end

# ─────────────────────────────────────────────────────────────────────────────
# solid_volume
# ─────────────────────────────────────────────────────────────────────────────

"""
    solid_volume(solid::BRepSolid) → Float64

Compute the (positive) volume of a closed solid via the divergence theorem:

    V = (1/6) |Σ_{triangles} v₁ · (v₂ × v₃)|

Each face is triangulated by fan triangulation from its first vertex.
The face winding must be outward-CCW for the formula to give positive volume.
"""
function solid_volume(solid::BRepSolid{FT}) where {FT}
    vs = solid.vertices
    vol = zero(FT)
    for face in solid.faces
        idx = face.vertex_indices
        n   = length(idx)
        # Fan triangulation: (idx[1], idx[k], idx[k+1]) for k = 2..n-1
        v0 = vs[idx[1]]
        for k in 2:n-1
            va = vs[idx[k]]
            vb = vs[idx[k+1]]
            vol += dot(v0, cross(va, vb))
        end
    end
    return abs(vol) / 6
end

# ─────────────────────────────────────────────────────────────────────────────
# solid_surface_area
# ─────────────────────────────────────────────────────────────────────────────

"""
    solid_surface_area(solid::BRepSolid) → Float64

Total surface area of the solid (sum over all faces, triangulated).
"""
function solid_surface_area(solid::BRepSolid{FT}) where {FT}
    vs   = solid.vertices
    area = zero(FT)
    for face in solid.faces
        idx = face.vertex_indices
        n   = length(idx)
        v0  = vs[idx[1]]
        for k in 2:n-1
            va = vs[idx[k]]
            vb = vs[idx[k+1]]
            area += norm(cross(va - v0, vb - v0)) / 2
        end
    end
    return area
end

# ─────────────────────────────────────────────────────────────────────────────
# check_manifold
# ─────────────────────────────────────────────────────────────────────────────

"""
    check_manifold(solid::BRepSolid; warn=true) → Bool

Return `true` if every edge is shared by exactly 2 faces (2-manifold condition).
Edges shared by 0, 1, or ≥ 3 faces trigger a `@warn` (unless `warn=false`) and
cause the function to return `false`.
"""
function check_manifold(solid::BRepSolid; warn::Bool = true) :: Bool
    edge_count = Dict{Tuple{Int,Int}, Int}()
    for face in solid.faces
        idx = face.vertex_indices
        n   = length(idx)
        for k in 1:n
            a = idx[k]
            b = idx[mod1(k+1, n)]
            key = a < b ? (a, b) : (b, a)
            edge_count[key] = get(edge_count, key, 0) + 1
        end
    end

    ok = true
    for (e, cnt) in edge_count
        if cnt != 2
            if warn
                @warn "Non-manifold edge $e shared by $cnt faces (expected 2)"
            end
            ok = false
        end
    end
    return ok
end

# ─────────────────────────────────────────────────────────────────────────────
# convert_to_triangle_mesh
# ─────────────────────────────────────────────────────────────────────────────

"""
    convert_to_triangle_mesh(solid::BRepSolid) → TriangleMesh

Triangulate all faces of `solid` (fan triangulation from first vertex) and
return a `TriangleMesh`.

- Node matrix: columns are vertex coordinates from `solid.vertices`.
- Triangle connectivity: preserves the original vertex indices (no node merge).
- Tags: the 1-based face index is stored as the triangle tag.
"""
function convert_to_triangle_mesh(solid::BRepSolid{FT}) where {FT}
    vs = solid.vertices
    Nv = length(vs)

    # Node matrix: 3 × Nv
    node = Matrix{FT}(undef, 3, Nv)
    for i in 1:Nv
        node[:, i] .= vs[i]
    end

    # Count total number of triangles
    ntri = sum(max(length(f.vertex_indices) - 2, 0) for f in solid.faces)

    triangles = Matrix{Int}(undef, 3, ntri)
    tags      = Vector{Int}(undef, ntri)

    ti = 1
    for (fi, face) in enumerate(solid.faces)
        idx = face.vertex_indices
        n   = length(idx)
        for k in 2:n-1
            triangles[:, ti] = [idx[1], idx[k], idx[k+1]]
            tags[ti]          = fi
            ti += 1
        end
    end

    return TriangleMesh(ntri, node, triangles, tags)
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

# Extract all undirected edges from a face list (no duplicates, sorted).
function _extract_edges(faces::Vector{BRepFace}) :: Vector{Tuple{Int,Int}}
    edge_set = Set{Tuple{Int,Int}}()
    for face in faces
        idx = face.vertex_indices
        n   = length(idx)
        for k in 1:n
            a = idx[k]
            b = idx[mod1(k+1, n)]
            push!(edge_set, a < b ? (a, b) : (b, a))
        end
    end
    return sort!(collect(edge_set))
end
