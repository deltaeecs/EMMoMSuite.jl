"""
    BooleanOps.jl

Phase 18.2a — Boolean operations on convex polyhedral B-Rep solids.

Only **convex** polyhedra are fully supported for `intersect_solids`.
`union_solids` and `subtract_solid` return a `CSGNode` (deferred evaluation);
use `csg_volume` for inclusion-exclusion volume computation.

Exported functions:
  intersect_solids  — Sutherland-Hodgman 3-D convex polyhedron clipping
  union_solids      — return CSGNode(:union, a, b)
  subtract_solid    — return CSGNode(:subtract, a, b)
  csg_volume        — volume evaluation via inclusion-exclusion (convex CSG trees)
"""

using LinearAlgebra
using StaticArrays

export intersect_solids, union_solids, subtract_solid, csg_volume

# ─────────────────────────────────────────────────────────────────────────────
# Public: union / subtract (deferred as CSGNode)
# ─────────────────────────────────────────────────────────────────────────────

"""
    union_solids(a::BRepSolid, b::BRepSolid) → CSGNode

Return a CSG tree node representing the set-union A ∪ B.
Use `csg_volume` to evaluate its volume.
"""
union_solids(a::BRepSolid, b::BRepSolid) = CSGNode(:union, a, b)

"""
    subtract_solid(a::BRepSolid, b::BRepSolid) → CSGNode

Return a CSG tree node representing the set-difference A \\ B.
Use `csg_volume` to evaluate its volume.
"""
subtract_solid(a::BRepSolid, b::BRepSolid) = CSGNode(:subtract, a, b)

# ─────────────────────────────────────────────────────────────────────────────
# Public: intersect_solids
# ─────────────────────────────────────────────────────────────────────────────

"""
    intersect_solids(a::BRepSolid, b::BRepSolid) → BRepSolid

Compute the geometric intersection A ∩ B by clipping `a` against each face
half-space of `b` using the 3-D Sutherland-Hodgman algorithm.

Float types are automatically promoted to `promote_type(FA, FB)`.

!!! note
    Both `a` and `b` must be **convex** polyhedra.  The result is also convex.
    If the intersection is empty, returns a degenerate solid with no faces.
"""
function intersect_solids(a::BRepSolid{FA}, b::BRepSolid{FB}) where {FA, FB}
    FT = promote_type(FA, FB)
    a2 = _convert_solid(FT, a)
    b2 = _convert_solid(FT, b)
    return _intersect_impl(a2, b2)
end

function _intersect_impl(a::BRepSolid{FT}, b::BRepSolid{FT}) where {FT}
    result = a
    for face in b.faces
        n, d = _face_halfspace(b.vertices, face)
        clipped = _clip_by_halfspace(result, n, d)
        if isnothing(clipped)
            return _empty_solid(a)
        end
        result = clipped
    end
    return result
end

# ─────────────────────────────────────────────────────────────────────────────
# Public: csg_volume
# ─────────────────────────────────────────────────────────────────────────────

"""
    csg_volume(node::CSGNode) → Float64

Evaluate the volume of a CSG tree using set-theoretic inclusion-exclusion:
- `:leaf`      → volume of the leaf solid
- `:union`     → V(L) + V(R) - V(L ∩ R)
- `:subtract`  → V(L) - V(L ∩ R)
- `:intersect` → V(L ∩ R)

!!! note
    Only convex `BRepSolid` leaf nodes are supported in Phase 18.2a.
    Nested `CSGNode` children must themselves be `:leaf` nodes.
"""
function csg_volume(node::CSGNode)
    node.op == :leaf && return _csg_vol(node.left)

    v_left  = _csg_vol(node.left)
    iv      = _csg_intersection_volume(node.left, node.right)

    node.op == :intersect && return iv
    node.op == :union     && return v_left + _csg_vol(node.right) - iv
    node.op == :subtract  && return v_left - iv

    error("csg_volume: unknown operation :$(node.op)")
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

# Volume dispatch
_csg_vol(s::BRepSolid) = solid_volume(s)
_csg_vol(n::CSGNode)   = csg_volume(n)

# Intersection volume (only BRepSolid leaves supported in 18.2a)
function _csg_intersection_volume(a, b)
    (a isa BRepSolid && b isa BRepSolid) || return 0.0
    return solid_volume(intersect_solids(a, b))
end

# Compute the outward unit normal and plane constant d for a face.
# Assumes CCW winding when viewed from outside (right-hand rule).
function _face_halfspace(vertices::Vector{SVector{3,FT}}, face::BRepFace) where {FT}
    idx = face.vertex_indices
    v0  = vertices[idx[1]]
    v1  = vertices[idx[2]]
    v2  = vertices[idx[3]]
    n   = normalize(cross(v1 - v0, v2 - v0))
    d   = dot(n, v0)
    return n, d
end

# Return an empty BRepSolid with the same labels as template.
function _empty_solid(template::BRepSolid{FT}) where {FT}
    return BRepSolid{FT}(
        SVector{3,FT}[],
        Tuple{Int,Int}[],
        BRepFace[],
        template.material_label,
        Dict{Int,String}(),
    )
end

# Convert a BRepSolid to a different float type.
function _convert_solid(::Type{FT2}, s::BRepSolid{FT1}) where {FT2, FT1}
    return BRepSolid{FT2}(
        [SVector{3,FT2}(v) for v in s.vertices],
        copy(s.edges),
        copy(s.faces),
        s.material_label,
        copy(s.boundary_labels),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Sutherland-Hodgman 3-D convex clipping
# ─────────────────────────────────────────────────────────────────────────────

"""
    _clip_by_halfspace(solid, n, d) → BRepSolid or nothing

Clip `solid` against the halfspace  {x : n·x ≤ d}.
Returns `nothing` if the result is empty (entire solid is outside).

Algorithm:
1. Classify vertices as inside (n·v ≤ d) or outside.
2. Clip each face polygon via Sutherland-Hodgman in 2D (on the face plane).
3. Collect all edge-plane intersection vertices and build the cap face.
"""
function _clip_by_halfspace(solid::BRepSolid{FT}, n::SVector{3,FT}, d::FT) where {FT}
    vs  = solid.vertices
    N   = length(vs)

    # Signed distance from clipping plane (positive = outside)
    tol = FT(1e-9)
    sd  = [dot(n, v) - d for v in vs]
    inside = [s <= tol for s in sd]

    all(inside)  && return solid   # nothing to clip
    !any(inside) && return nothing # fully outside

    # ── Build new vertex list ──────────────────────────────────────────────
    new_verts = SVector{3,FT}[]
    v_map     = fill(0, N)          # old index → new index (0 = not inside)

    for i in 1:N
        if inside[i]
            push!(new_verts, vs[i])
            v_map[i] = length(new_verts)
        end
    end

    # Storage for cap face vertices (intersection of edges with clipping plane)
    cap_pos     = SVector{3,FT}[]  # positions (for dedup)
    cap_new_idx = Int[]            # corresponding new_verts indices

    # Helper: find or create intersection vertex on the clipping plane
    function get_cap_vertex(i::Int, j::Int)::Int
        # t ∈ (0,1) s.t. dot(n, vs[i]+t*(vs[j]-vs[i])) = d
        t     = -sd[i] / (sd[j] - sd[i])
        v_int = vs[i] + t * (vs[j] - vs[i])
        # Dedup by position
        for (ci, cv) in enumerate(cap_pos)
            norm(cv - v_int) < FT(1e-8) && return cap_new_idx[ci]
        end
        push!(new_verts, v_int)
        idx = length(new_verts)
        push!(cap_pos, v_int)
        push!(cap_new_idx, idx)
        return idx
    end

    # ── Sutherland-Hodgman: process each face ─────────────────────────────
    new_faces = BRepFace[]

    for face in solid.faces
        idx = face.vertex_indices
        nf  = length(idx)
        clipped = Int[]

        for k in 1:nf
            i = idx[k]
            j = idx[mod1(k + 1, nf)]

            if inside[i]
                push!(clipped, v_map[i])
                !inside[j] && push!(clipped, get_cap_vertex(i, j))
            else
                inside[j] && push!(clipped, get_cap_vertex(i, j))
            end
        end

        length(clipped) >= 3 && push!(new_faces, BRepFace(clipped))
    end

    # ── Build cap face (CCW orientation when viewed from +n) ──────────────
    if length(cap_new_idx) >= 3
        ordered = _order_polygon_ccw(new_verts, cap_new_idx, n)
        push!(new_faces, BRepFace(ordered))
    end

    isempty(new_faces) && return nothing

    edges = _extract_edges(new_faces)
    return BRepSolid{FT}(
        new_verts, edges, new_faces,
        solid.material_label, copy(solid.boundary_labels),
    )
end

"""
    _order_polygon_ccw(verts, indices, n) → Vector{Int}

Order the vertices at `verts[indices]` in counter-clockwise order when viewed
from the direction of `n` (outward normal of the polygon).

Uses a planar projection onto the tangent plane perpendicular to `n` and sorts
by azimuthal angle.  Only valid for **convex** polygon vertex sets.
"""
function _order_polygon_ccw(
    verts::Vector{SVector{3,FT}},
    indices::Vector{Int},
    n::SVector{3,FT},
) where {FT}
    length(indices) <= 2 && return copy(indices)

    centroid = sum(verts[i] for i in indices) / length(indices)

    # Tangent basis: u⊥n, v = n̂×u  →  (u, v, n̂) is right-handed
    n̂ = normalize(n)
    # Pick u as a cardinal axis not parallel to n̂
    u_tmp = abs(n̂[1]) < FT(0.9) ? SVector{3,FT}(1, 0, 0) :
                                    SVector{3,FT}(0, 1, 0)
    u = normalize(u_tmp - dot(u_tmp, n̂) * n̂)
    v = cross(n̂, u)   # (u, v, n̂) right-handed: u×v = n̂

    # Sort indices by angle in the (u,v) plane
    angles = [atan(dot(v, verts[i] - centroid),
                   dot(u, verts[i] - centroid))  for i in indices]
    return indices[sortperm(angles)]
end
